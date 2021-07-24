#include "leastsquaressolver.h"
#include "globals.h"

double M_[MAX_EQNS][MAX_EQNS],
M_0[MAX_EQNS][MAX_EQNS],
MWM[MAX_EQNS][MAX_EQNS],
x_m[MAX_EQNS],
b_m[MAX_EQNS],
w[MAX_EQNS],
mwb[MAX_EQNS],
LU[MAX_EQNS][MAX_EQNS],
Inv[MAX_EQNS][MAX_EQNS];
int ps[MAX_EQNS];


double get_w(double r2,double delta)
{
    // return exp(-r2/(delta*delta));
    //delta*delta/(r2+delta*delta);
    // return delta/(sqrt(r2)+delta);
    return delta*delta/(r2+delta*delta);

}
leastSquaresSolver::leastSquaresSolver()
{
    pww=0.49;
    rad=0.0;//0.125;
}


void leastSquaresSolver::LU_decompose(int size)
{
    int i,j,k,pivotindex;
    double scales[MAX_EQNS];
    double normrow,pivot,size1,biggest,mult;

    for (i=0;i<size;i++) //заполнение начальными данными
    {
        ps[i]=i;//маппинг изначального порядка на переставленный.
        normrow=0;//максимум в итой строке

        for (j=0;j<size;j++)
        {
            LU[i][j]=M_[i][j];
            if (normrow<fabs(LU[i][j]))
                normrow=fabs(LU[i][j]);
        }
        if (normrow!=0)
            scales[i]=1.0/normrow;//для общих множителей
        else
        {
            scales[i]=0.0;
            //     err_code(DIV_ZERO);
        }
    }
    //метод гаусса с частичным упорядочиванием

    for (k=0;k<size-1;k++)
    {
        biggest=0;
        for (i=k; i<size;i++)
        {
            size1=fabs(LU[ps[i]][k])*scales[ps[i]];
            if (biggest<size1)
            {
                biggest=size1;
                pivotindex=i;
            }
        }

        if (biggest==0)
        {
            //	err_code(1);
            pivotindex=0;
        }

        if (pivotindex!=k)
        {
            j=ps[k];
            ps[k]=ps[pivotindex];
            ps[pivotindex]=j;
        }

        pivot=LU[ps[k]][k];

        for (i=k+1;i<size;i++)
        {
            mult=LU[ps[i]][k]/pivot;
            LU[ps[i]][k]=mult;

            if (mult!=0.0)
            {
                for (j=k+1; j<size;j++)
                    LU[ps[i]][j]-=mult*LU[ps[k]][j];
            }
        }
    }
    //      if (LU[ps[VAR_NUM-1]][VAR_NUM-1]==0.0) err_code(1);
}


void leastSquaresSolver::m_solve(int num)
{
    int i,j;
    double dot;

    for (i=0;i<num;i++)
    {
        dot=0;
        for (j=0;j<i;j++)
            dot+=LU[ps[i]][j]*x_m[j];

        x_m[i]=b_m[ps[i]]-dot;
    }

    for (i=num-1; i>=0;i--)
    {
        dot=0.0;

        for (j=i+1;j<num;j++)
            dot+=LU[ps[i]][j]*x_m[j];

        x_m[i]=(x_m[i]-dot)/LU[ps[i]][i];
    }
}

void leastSquaresSolver::m_invert(int size)
{

    int i,j;
    //err_code(1);
    LU_decompose(size);

    for (j=0;j<size;j++)
    {

        for (i=0;i<size;i++)
        {
            if (i==j)
                b_m[i]=1;
            else
                b_m[i]=0;
        }

        m_solve(size);

        for (i=0;i<size;i++)
            Inv[i][j]=x_m[i];
    }
}

void leastSquaresSolver::nullify_m(double m[MAX_EQNS][MAX_EQNS], int in, int im)
{
    for (int i=0;i<in;i++)
        for(int j=0;j<im;j++)
            m[i][j]=0;
}

void leastSquaresSolver::nullify_v(double v[MAX_EQNS], int in)
{
    for (int i=0;i<in;i++)
        v[i]=0;
}

double leastSquaresSolver::getDivA(posVelAccelVort &p, double delta)
{
    int var_num=VAR_NUM;
    int eq_num=p.neighbours.size();
    // printf("eq_num=%d \n",eq_num);
    //first index is a row number
    //second index is a column number  M[eq_num][var_num]
    // get Mt*W*M
    //nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=Layer.currValues[p.neighbours.at(i)].x-p.x;
        dyl=Layer.currValues[p.neighbours.at(i)].y-p.y;
        dzl=Layer.currValues[p.neighbours.at(i)].z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=get_w(r2,delta);//exp(-r2/(delta*delta));//1.0;

        b_m[i]=0.0;
    }

    nullify_m(M_, var_num, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }

    //ax:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].ax;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double dax = x_m[FX];

    //ay:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].ay;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double day = x_m[FY];
    //az:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].az;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double daz = x_m[FZ];
    return dax+day+daz;
}

v3 leastSquaresSolver::getVort(posVelAccelVort &p, double delt)
{
    int var_num=VAR_NUM;
    int eq_num=p.neighbours.size();
    // printf("eq_num=%d \n",eq_num);
    //first index is a row number
    //second index is a column number  M[eq_num][var_num]
    // get Mt*W*M
    //nullify_m(M_0);

    double x0 = -50;
    double x1 = 50;
    double y0 = -25;
    double y1 = 25;
    double z0 = 0.01;
    double z1 = 30.01;
    double x,y,z;
    x=p.x; y=p.y; z=p.z;

    double d[5];
    d[0]=fabs(x-x0);
    d[1]=fabs(x-x1);
    d[2]=fabs(y-y0);
    d[3]=fabs(y-y1);
    d[4]=fabs(z-z1);

    double dist=d[0];

    for (int i=1;i<5;i++)
    {
        if (dist>d[i]) dist=d[i];
    }
    double ggg=1.86;
    double dd=ggg; //
    dist=fmax(dd-dist,0.0)/dd+1;
    double delta=delt*(dist);

    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=Layer.currValues[p.neighbours.at(i)].x-p.x;
        dyl=Layer.currValues[p.neighbours.at(i)].y-p.y;
        dzl=Layer.currValues[p.neighbours.at(i)].z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=get_w(r2,delta);//exp(-r2/(delta*delta));//1.0;

        b_m[i]=0.0;
    }

    nullify_m(M_, var_num, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }

    //u:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].u;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double dUDy = x_m[FY];
    double dUDz = x_m[FZ];

    //v:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].v;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double dVDx = x_m[FX];
    double dVDz = x_m[FZ];

    //w:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].w;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double dWDx = x_m[FX];
    double dWDy = x_m[FY];

    return v3(dWDy - dVDz ,dUDz - dWDx,dVDx - dUDy);
}


v3 leastSquaresSolver::getRotVort(posVelAccelVort &p, double delta)
{
    int var_num=VAR_NUM;
    int eq_num=p.neighbours.size();
    // printf("eq_num=%d \n",eq_num);
    //first index is a row number
    //second index is a column number  M[eq_num][var_num]
    // get Mt*W*M
    //nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=Layer.currValues[p.neighbours.at(i)].x-p.x;
        dyl=Layer.currValues[p.neighbours.at(i)].y-p.y;
        dzl=Layer.currValues[p.neighbours.at(i)].z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=get_w(r2,delta);//exp(-r2/(delta*delta));//1.0;

        b_m[i]=0.0;
    }

    nullify_m(M_, var_num, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }

    //vortx:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].vortx;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double dVotrxDy = x_m[FY];
    double dVotrxDz = x_m[FZ];

    //vorty:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].vorty;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double dVotryDx = x_m[FX];
    double dVotryDz = x_m[FZ];

    //vortz:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].vortz;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double dVotrzDx = x_m[FX];
    double dVotrzDy = x_m[FY];

    return v3(dVotrzDy - dVotryDz ,dVotrxDz - dVotrzDx,dVotryDx - dVotrxDy);
}

posVelAccelVort leastSquaresSolver::interpFull(posVelAccelVort &p, double delt)
{
    int var_num=VAR_NUM;
    int eq_num=p.neighbours.size();
    // printf("eq_num=%d \n",eq_num);
    //first index is a row number
    //second index is a column number  M[eq_num][var_num]
    // get Mt*W*M
    //nullify_m(M_0);
    posVelAccelVort res;

    double x0 = -50;
    double x1 = 50;
    double y0 = -25;
    double y1 = 25;
    double z0 = 0.01;
    double z1 = 30.01;
    double x,y,z;
    x=p.x; y=p.y; z=p.z;

    double d[5];
    d[0]=fabs(x-x0);
    d[1]=fabs(x-x1);
    d[2]=fabs(y-y0);
    d[3]=fabs(y-y1);
    d[4]=fabs(z-z1);






    double dist=d[0];

    for (int i=1;i<5;i++)
    {
        if (dist>d[i]) dist=d[i];
    }
    double ggg=1.86;
    double dd=ggg; //
    dist=fmax(dd-dist,0.0)/dd+1;
    double delta=delt*(dist);


    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=Layer.currValues[p.neighbours.at(i)].x-p.x;
        dyl=Layer.currValues[p.neighbours.at(i)].y-p.y;
        dzl=Layer.currValues[p.neighbours.at(i)].z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=get_w(r2,delta);;//exp(-r2/(delta*delta));//1.0;

        b_m[i]=0.0;
    }

    nullify_m(M_, var_num, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }
    LU_decompose(var_num);

    ////U
    for (int i=0;i<eq_num;i++) {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].u;
    }
    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++) {
        for (int n=0;n<eq_num;n++) {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }
    for (int i=0;i<var_num;i++) {
        b_m[i]=mwb[i];  //its mvb
    }
    m_solve(var_num);
    res.u = x_m[F];
    res.dux = x_m[FX];
    res.duy = x_m[FY];
    res.duz = x_m[FZ];

    ////V
    for (int i=0;i<eq_num;i++) {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].v;
    }
    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++) {
        for (int n=0;n<eq_num;n++) {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }
    for (int i=0;i<var_num;i++) {
        b_m[i]=mwb[i];  //its mvb
    }
    m_solve(var_num);
    res.v = x_m[F];
    res.dvx = x_m[FX];
    res.dvy = x_m[FY];
    res.dvz = x_m[FZ];

    ////W
    for (int i=0;i<eq_num;i++) {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].w;
    }
    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++) {
        for (int n=0;n<eq_num;n++) {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }
    for (int i=0;i<var_num;i++) {
        b_m[i]=mwb[i];  //its mvb
    }
    m_solve(var_num);
    res.w = x_m[F];
    res.dwx = x_m[FX];
    res.dwy = x_m[FY];
    res.dwz = x_m[FZ];

    ////P
    for (int i=0;i<eq_num;i++) {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].f;
    }
    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++) {
        for (int n=0;n<eq_num;n++) {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }
    for (int i=0;i<var_num;i++) {
        b_m[i]=mwb[i];  //its mvb
    }
    m_solve(var_num);

    ////ax
    for (int i=0;i<eq_num;i++) {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].ax;
    }
    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++) {
        for (int n=0;n<eq_num;n++) {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }
    for (int i=0;i<var_num;i++) {
        b_m[i]=mwb[i];  //its mvb
    }
    m_solve(var_num);
    res.ax = x_m[F];

    ////ay
    for (int i=0;i<eq_num;i++) {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].ay;
    }
    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++) {
        for (int n=0;n<eq_num;n++) {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }
    for (int i=0;i<var_num;i++) {
        b_m[i]=mwb[i];  //its mvb
    }
    m_solve(var_num);
    res.ay = x_m[F];

    ////az
    for (int i=0;i<eq_num;i++) {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].az;
    }
    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++) {
        for (int n=0;n<eq_num;n++) {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }
    for (int i=0;i<var_num;i++) {
        b_m[i]=mwb[i];  //its mvb
    }
    m_solve(var_num);
    res.az = x_m[F];

    ////vortx
    for (int i=0;i<eq_num;i++) {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].vortx;
    }
    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++) {
        for (int n=0;n<eq_num;n++) {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }
    for (int i=0;i<var_num;i++) {
        b_m[i]=mwb[i];  //its mvb
    }
    m_solve(var_num);
    res.vortx =  x_m[F];//res.dwy - res.dvz;//

    ////vorty
    for (int i=0;i<eq_num;i++) {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].vorty;
    }
    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++) {
        for (int n=0;n<eq_num;n++) {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }
    for (int i=0;i<var_num;i++) {
        b_m[i]=mwb[i];  //its mvb
    }
    m_solve(var_num);
    res.vorty = x_m[F];//res.duz - res.dwx;//
    ////vortz
    for (int i=0;i<eq_num;i++) {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].vortz;
    }
    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++) {
        for (int n=0;n<eq_num;n++) {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }
    for (int i=0;i<var_num;i++) {
        b_m[i]=mwb[i];  //its mvb
    }
    m_solve(var_num);
    res.vortz = x_m[F];//res.dvx - res.duy;//

    ////rotvortx
    for (int i=0;i<eq_num;i++) {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].rotvortx;
    }
    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++) {
        for (int n=0;n<eq_num;n++) {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }
    for (int i=0;i<var_num;i++) {
        b_m[i]=mwb[i];  //its mvb
    }
    m_solve(var_num);
    res.rotvortx =  x_m[F];//res.dwy - res.dvz;//

    ////rotvorty
    for (int i=0;i<eq_num;i++) {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].rotvorty;
    }
    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++) {
        for (int n=0;n<eq_num;n++) {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }
    for (int i=0;i<var_num;i++) {
        b_m[i]=mwb[i];  //its mvb
    }
    m_solve(var_num);
    res.rotvorty = x_m[F];//res.duz - res.dwx;//

    ////rotvortz
    for (int i=0;i<eq_num;i++) {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].rotvortz;
    }
    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++) {
        for (int n=0;n<eq_num;n++) {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }
    for (int i=0;i<var_num;i++) {
        b_m[i]=mwb[i];  //its mvb
    }
    m_solve(var_num);
    res.rotvortz = x_m[F];//res.dvx - res.duy;//


    return res;
}

void leastSquaresSolver::interp(posVelAccelVort &p, double delta)
{
    int var_num=VAR_NUM;
    int eq_num=p.neighbours.size();
    // printf("eq_num=%d \n",eq_num);
    //first index is a row number
    //second index is a column number  M[eq_num][var_num]
    // get Mt*W*M
    //nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=Layer.currValues[p.neighbours.at(i)].x-p.x;
        dyl=Layer.currValues[p.neighbours.at(i)].y-p.y;
        dzl=Layer.currValues[p.neighbours.at(i)].z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=get_w(r2,delta);//exp(-r2/(delta*delta));//1.0;

        b_m[i]=0.0;
    }

    nullify_m(M_, var_num, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }

    //ax:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].ax;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double dp = x_m[F];


    //ay:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].ay;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double dpref = x_m[F];
    //az:


    //vortx:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].vortx;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double vortx = x_m[F];

    //vorty:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].vorty;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double vorty = x_m[F];

    //vortz:
    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].vortz;
    }

    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    double vortz = x_m[F];

    p.ax=dp;
    p.ay=dpref;
    p.vortx = vortx;
    p.vorty = vorty;
    p.vortz = vortz;

}

void leastSquaresSolver::get_vortEq_combo(posVelAccelVort &p, deriv3D &res, double delta)
{
    //we do not chack if p is a boundary point this should be done before calling this function

    int var_num=VAR_NUM;
    int eq_num=p.neighbours.size();

    //nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=Layer.currValues[p.neighbours.at(i)].x-p.x;
        dyl=Layer.currValues[p.neighbours.at(i)].y-p.y;
        dzl=Layer.currValues[p.neighbours.at(i)].z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=get_w(r2,delta);//exp(-r2/(delta*delta));//1.0;  todo: make w anisotrpopic in the direction normal to the wall

        b_m[i]=0.0;
    }
    int i=eq_num; //constraint for laplacian
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0;             M_0[i][2]=0;       M_0[i][3]=0;
    M_0[i][4]=1;             M_0[i][5]=0;       M_0[i][6]=0;
    M_0[i][7]=1;             M_0[i][8]=0;       M_0[i][9]=1;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;
    /////////////////////


    nullify_m(M_, var_num, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }

    for (int i=0;i<eq_num-1;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].f;//from prev interation
    }

    b_m[eq_num-1]=p.rhs;//from prev interation


    nullify_v(mwb, var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose(var_num);
    m_solve(var_num);

    p.f=x_m[0];//res.d[F];
}

void leastSquaresSolver::get_vortEq_matr(posVelAccelVort &p,  double delta)
{
    int var_num=VAR_NUM;
    int eq_num=p.neighbours.size();

    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=Layer.currValues[p.neighbours.at(i)].x-p.x;
        dyl=Layer.currValues[p.neighbours.at(i)].y-p.y;
        dzl=Layer.currValues[p.neighbours.at(i)].z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=get_w(r2,delta);//exp(-r2/(delta*delta));//1.0;  todo: make w anisotrpopic in the direction normal to the wall

        b_m[i]=0.0;
    }
    int i=eq_num; //constraint for laplacian
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0;             M_0[i][2]=0;       M_0[i][3]=0;
    M_0[i][4]=1;             M_0[i][5]=0;       M_0[i][6]=0;
    M_0[i][7]=1;             M_0[i][8]=0;       M_0[i][9]=1;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;
    /////////////////////

    nullify_m(M_,var_num,var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }

    m_invert(var_num);

    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {

            p.INV[i][j]=Inv[i][j];

        }
    }
}

void leastSquaresSolver::get_vortEq_nomatr(posVelAccelVort &p, double delta)
{
    //we do not chack if p is a boundary point this should be done before calling this function
    int var_num=VAR_NUM;
    int eq_num=p.neighbours.size();

    //nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=Layer.currValues[p.neighbours.at(i)].x-p.x;
        dyl=Layer.currValues[p.neighbours.at(i)].y-p.y;
        dzl=Layer.currValues[p.neighbours.at(i)].z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=get_w(r2,delta);//exp(-r2/(delta*delta));//1.0;  todo: make w anisotrpopic in the direction normal to the wall

        b_m[i]=0.0;
    }
    int i=eq_num; //constraint for laplacian
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0;             M_0[i][2]=0;       M_0[i][3]=0;
    M_0[i][4]=1;             M_0[i][5]=0;       M_0[i][6]=0;
    M_0[i][7]=1;             M_0[i][8]=0;       M_0[i][9]=1;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;
    /////////////////////


    for (int i=0;i<eq_num-1;i++)
    {
        b_m[i]=Layer.currValues[p.neighbours.at(i)].f;//from prev interation
    }

    b_m[eq_num-1]=p.rhs;//from prev interation
    //b_m[eq_num-1]=p.u;//from prev interation


    nullify_v(mwb,var_num);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    p.f=0.0;
    for (int j=0;j<var_num;j++)
    {
        p.f+=p.INV[0][j]*mwb[j];
    }

}

///debug algorithms:
void leastSquaresSolver::fill_u()
{

    for (int i=0;i<Layer.currValues.size();i++)
    {

        double dx=Layer.currValues[i].x;
        double dy=Layer.currValues[i].y;
        double dz=Layer.currValues[i].z;

        double r2=dx*dx + dy*dy + dz*dz;

        Layer.currValues[i].u0=-dz/(r2+0.01);
        Layer.currValues[i].v0=0.0;
        Layer.currValues[i].w0=dx/(r2+0.01);

    }
    getVort_FromVel_LS();
    getVel_FromVort();


}
void leastSquaresSolver::fill_omega()
{
    Layer.currValues.clear();
    for (int i=0;i<30;i++)
    {
        posVelAccelVort n;
        n.x=1.0*(rand()*1.0/RAND_MAX-0.5);
        n.y=0.01*(rand()*1.0/RAND_MAX-0.5);
        n.z=1.0*(rand()*1.0/RAND_MAX-0.5);;

        n.vortx=0.001*(rand()*1.0/RAND_MAX-0.5);
        n.vorty=rand()*1.0/RAND_MAX;
        n.vortz=0.001*(rand()*1.0/RAND_MAX-0.5);
        Layer.currValues.push_back(n);
    }
    getVel_FromVort();
}

void leastSquaresSolver::getVel_FromVort()
{
    //Ax=b where x is omega and b is vel
    //three vectro components are going in sequence for each i
    //velocity at i'th points from j'th points
    int num=Layer.currValues.size();

    for (int i=0;i<num;i++)
    {
        //the
        x_m[i*3]=Layer.currValues[i].vortx;
        x_m[i*3+1]=Layer.currValues[i].vorty;
        x_m[i*3+2]=Layer.currValues[i].vortz;

        //b_m[i]+=M_[n][i]*x_m[n];
        //ux[i]=M_[n][i]*wx[n] + M_[n+1][i]*wy[n] +M_[n+2][i]*wz[n]
        // ux=zr3*(dy*wz - dz*wy)
        // uy=zr3*(dz*wx - dx*wz)
        // uz=zr3*(dx*wy - dy*wx)
        for (int j=0;j<num;j++)
        {
            double dx,dy,dz;
            dx=Layer.currValues[i].x-Layer.currValues[j].x;
            dy=Layer.currValues[i].y-Layer.currValues[j].y;
            dz=Layer.currValues[i].z-Layer.currValues[j].z;
            double r2=dx*dx + dy*dy + dz*dz;
            double zr3=1.0/(r2*sqrt(r2)+0.001);
            //ux from wx                 ux from wy                  ux from wz
            M_[i*3][j*3]=0.0;            M_[i*3 + 1][j*3]=-zr3*dz;   M_[i*3 + 2][j*3]=zr3*dy;

            //uy from wx                 uy from wy                      uy from wz
            M_[i*3][j*3 + 1]=zr3*dz;     M_[i*3 + 1][j*3 + 1]=0.0;   M_[i*3 + 2][j*3 + 1]=-zr3*dx;

            //uz from wx                 uz from wy                      uz from wz
            M_[i*3][j*3 + 2]=-zr3*dy;     M_[i*3 + 1][j*3 + 2]=zr3*dx;   M_[i*3 + 2][j*3 + 2]=0.0;

        }
    }

    for (int i=0;i<num*3;i++)
    {
        b_m[i]=0;
        for (int n=0;n<num*3;n++)
        {
            b_m[i]+=M_[n][i]*x_m[n];  //Ax=b
        }
    }

    for (int i=0;i<num;i++)
    {
        Layer.currValues[i].u=b_m[i*3];
        Layer.currValues[i].v=b_m[i*3+1];
        Layer.currValues[i].w=b_m[i*3+2];

        //printf("u=%f v=%f w=%f \n",m_p[i].u,m_p[i].v,m_p[i].w);
    }
}

void leastSquaresSolver::getVort_FromVel()
{
    //Ax=b where x is omega and b is vel
    //three vectro components are going in sequence for each i
    //velocity at i'th points from j'th points
    int num=Layer.currValues.size();
    for (int i=0;i<num;i++)
    {
        //the direct inversion method
        b_m[i*3]=Layer.currValues[i].u0;
        b_m[i*3+1]=Layer.currValues[i].v0;
        b_m[i*3+2]=Layer.currValues[i].w0;

        //b_m[i]+=M_[n][i]*x_m[n];
        //ux[i]=M_[n][i]*wx[n] + M_[n+1][i]*wy[n] +M_[n+2][i]*wz[n]
        // ux=zr3*(dy*wz - dz*wy)
        // uy=zr3*(dz*wx - dx*wz)
        // uz=zr3*(dx*wy - dy*wx)
        for (int j=0;j<num;j++)
        {
            double dx,dy,dz;
            dx=Layer.currValues[i].x-Layer.currValues[j].x;
            dy=Layer.currValues[i].y-Layer.currValues[j].y;
            dz=Layer.currValues[i].z-Layer.currValues[j].z;
            double r2=dx*dx + dy*dy + dz*dz;
            double zr3=1.0/(r2*sqrt(r2)+0.001);
            //ux from wx                 ux from wy                  ux from wz
            M_[i*3][j*3]=0.0;            M_[i*3 + 1][j*3]=-zr3*dz;   M_[i*3 + 2][j*3]=zr3*dy;

            //uy from wx                 uy from wy                      uy from wz
            M_[i*3][j*3 + 1]=zr3*dz;     M_[i*3 + 1][j*3 + 1]=0.0;   M_[i*3 + 2][j*3 + 1]=-zr3*dx;

            //uz from wx                 uz from wy                      uz from wz
            M_[i*3][j*3 + 2]=-zr3*dy;     M_[i*3 + 1][j*3 + 2]=zr3*dx;   M_[i*3 + 2][j*3 + 2]=0.0;

        }
    }

    LU_decompose(num*3);
    m_solve(num*3);


    for (int i=0;i<num;i++)
    {
        Layer.currValues[i].vortx=x_m[i*3];
        Layer.currValues[i].vorty=x_m[i*3+1];
        Layer.currValues[i].vortz=x_m[i*3+2];
    }
}

void leastSquaresSolver::getVort_FromVel_LS()
{
    //solving in least squares sense
    //Ax=b where x is omega and b is vel
    //three vectro components are going in sequence for each i
    //velocity at i'th points from j'th points

    int num=Layer.currValues.size();

    for (int i=0;i<num;i++)
    {
        //the Least squares inversion method
        b_m[i*3]=  Layer.currValues[i].u0;
        b_m[i*3+1]=Layer.currValues[i].v0;
        b_m[i*3+2]=Layer.currValues[i].w0;

        //b_m[i]+=M_[n][i]*x_m[n];
        //ux[i]=M_[n][i]*wx[n] + M_[n+1][i]*wy[n] +M_[n+2][i]*wz[n]
        // ux=zr3*(dy*wz - dz*wy)
        // uy=zr3*(dz*wx - dx*wz)
        // uz=zr3*(dx*wy - dy*wx)
        for (int j=0;j<num;j++)
        {
            double dx,dy,dz;
            dx=Layer.currValues[i].x - Layer.currValues[j].x;
            dy=Layer.currValues[i].y - Layer.currValues[j].y;
            dz=Layer.currValues[i].z - Layer.currValues[j].z;
            double r2=dx*dx + dy*dy + dz*dz;
            double zr3=1.0/(r2*sqrt(r2)+0.001);
            //ux from wx                 ux from wy                  ux from wz
            M_0[i*3][j*3]=0.0;            M_0[i*3 + 1][j*3]=-zr3*dz;   M_0[i*3 + 2][j*3]=zr3*dy;

            //uy from wx                 uy from wy                      uy from wz
            M_0[i*3][j*3 + 1]=zr3*dz;     M_0[i*3 + 1][j*3 + 1]=0.0;   M_0[i*3 + 2][j*3 + 1]=-zr3*dx;

            //uz from wx                 uz from wy                      uz from wz
            M_0[i*3][j*3 + 2]=-zr3*dy;     M_0[i*3 + 1][j*3 + 2]=zr3*dx;   M_0[i*3 + 2][j*3 + 2]=0.0;

        }
    }

    for (int i=0;i<num*3;i++)
    {
        for (int j=0;j<num*3;j++)
        {
            M_[i][j]=0.0;
            for (int n=0;n<num*3;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j];  //its mvm
            }
        }
    }

    for (int i=0;i<num*3;i++)
    {
        mwb[i]=0.0;
        for (int n=0;n<num*3;n++)
        {
            mwb[i]+=M_0[n][i]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<num*3;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }

    LU_decompose(num*3);
    m_solve(num*3);


    for (int i=0;i<num;i++)
    {
        Layer.currValues[i].vortx=x_m[i*3];
        Layer.currValues[i].vorty=x_m[i*3+1];
        Layer.currValues[i].vortz=x_m[i*3+2];
    }
}
