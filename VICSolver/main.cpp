#include "globals.h"
#include "leastsquaressolver.h"

void threadOptimize(int threadIdx, int startIdx, int endIdx);
void display(void);
void init();

int iCoarse=0;
int jCoarse=0;
int kCoarse=0;

int pnum_in=0;

void findNeighbors_coarse_smart();
void solvePoisson();
void solveVortEq();
void solvePoisson_nomatr();
void solveVortEq_nomatr();
void calcVort();
void calcRotVort();

posVelAccelVort interpolateFull(double x, double y, double z);
posVelAccelVort interpolate_smart(double x,double y, double z);
int removeWorst_free(mLayer& m, int* arr, int num,double x, double y,double z,int desiredNum);
int desiredStencilNum=24;//19

leastSquaresSolver lssA;

double reference[251][126][16];
enum ref
{
    X_=0,
    Y_=1,
    Z_=2,
    U_=3,
    V_=4,
    W_=5,
    U1_=6,
    V1_=7,
    W1_=8,
    AX_=9,
    AY_=10,
    AZ_=11,
    P_=12,
    VortX_=13,
    VortY_=14,
    VortZ_=15
};

void loadReference()
{
    FILE* f=fopen("../MSTfield00801_cross.dat","r");
    char str[2048];
    for (int i=0;i<17;i++) fgets(str,2048,f);

    int mass=0;
    for (int k=0;k<126;k++)
        for (int i=0;i<251;i++)
        {
            double x_,y_,z_,u_,v_,w_,u1_,v1_,w1_,ax_,ay_,az_,p_;
            fgets(str,2048,f);
            sscanf(str,"%lf  %lf %lf %lf %lf %lf %lf  %lf  %lf  %lf  %lf  %lf %lf",
                   &x_,&y_, &z_,&u_,&v_,&w_,&u1_,&v1_,&w1_,&ax_,&ay_,&az_,&p_);
            reference[i][k][X_]=x_;
            reference[i][k][Y_]=y_;
            reference[i][k][Z_]=z_;
            reference[i][k][U_]=u_;
            reference[i][k][V_]=v_;
            reference[i][k][W_]=w_;
            reference[i][k][U1_]=u1_;
            reference[i][k][V1_]=v1_;
            reference[i][k][W1_]=w1_;
            reference[i][k][AX_]=ax_;
            reference[i][k][AY_]=ay_;
            reference[i][k][AZ_]=az_;
            reference[i][k][P_]=p_;
            mass+=1;
        }

    for (int k=0;k<126-1;k++)
        for (int i=0;i<251-1;i++)
        {
            reference[i][k][VortX_]=0.0;
            reference[i][k][VortY_]=(reference[i][k+1][U_] - reference[i][k][U_]) / (reference[i][k+1][Z_] - reference[i][k][Z_])
                    - (reference[i+1][k][W_] - reference[i][k][W_]) / (reference[i+1][k][X_] - reference[i][k][X_]);
            reference[i][k][VortZ_]=0.0;
        }
}

void display(void)
{
    double orient_x=0.0;
    double orient_y=0.0;
    double orient_z=5.0;

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    o_x=cos(-ry/180*M_PI)*cos(-rx/180*M_PI);
    o_y=cos(-ry/180*M_PI)*sin(-rx/180*M_PI);
    o_z=sin(-ry/180*M_PI);
    orient_x=view_x+o_x;
    orient_y=view_y+o_y;
    orient_z=view_z+o_z;
    gluLookAt(view_x,view_y,view_z,orient_x,orient_y,orient_z,0,0,1);

    double mean_curr=0.0;
    int mass=0;
    for( size_t i = 0; i < 250; i++)
    {
        for( size_t k = 0; k < 126; k++)
        {
            int j=120;
            mean_curr+=reference[i][k][U_] ;
            mass++;
        }
    }
    mean_curr/=mass;

    {
        glPointSize(5);
        glBegin(GL_POINTS);
        for( size_t i = 0; i < Layer.currValues.size(); i+=1 )
        {
            if( Layer.currValues[i].y<Layer.y0_+(54*1.0/100)*(Layer.y1_-Layer.y0_)
                    && Layer.currValues[i].y>Layer.y0_+(46*1.0/100)*(Layer.y1_-Layer.y0_)) {
                double val = (Layer.currValues[i].u-mean_curr)*scale;
                glColor3f(-val,-val,val);
                if (fabs(Layer.currValues[i].ax) < 1e-5)
                    glColor3f(1,0,0);
                glVertex3f(Layer.currValues[i].x, Layer.currValues[i].y, Layer.currValues[i].z - Layer.z1_);
            }
        }
        glEnd();
    }

    if(drawVelocity)
    {
        glLineWidth(1);
        glBegin(GL_LINES);
        for( size_t i = 0; i < Layer.currValues.size(); i+=1 )
        {
            glColor3f(1,1,1);
            glVertex3f(Layer.currValues[i].x, Layer.currValues[i].y, Layer.currValues[i].z);
            glColor3f(0,0,0);
            glVertex3f(Layer.currValues[i].x+scale*5.0*Layer.currValues[i].u, Layer.currValues[i].y+scale*5.0*Layer.currValues[i].v, Layer.currValues[i].z+scale*5.0*Layer.currValues[i].w);
        }
        glEnd();
    }

    if(drawAcceleration)
    {
        glLineWidth(2);
        glBegin(GL_LINES);
        for( size_t i = 0; i < Layer.currValues.size(); i+=1 )
        {
            glColor3f(1,0,1);
            glVertex3f(Layer.currValues[i].x, Layer.currValues[i].y, Layer.currValues[i].z);
            glColor3f(0,0,0);
            glVertex3f(Layer.currValues[i].x+scale*1.0*Layer.currValues[i].ax, Layer.currValues[i].y+scale*1.0*Layer.currValues[i].ay, Layer.currValues[i].z+scale*1.0*Layer.currValues[i].az);
        }
        glEnd();
    }

    for( size_t i = 0; i < 250; i++)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for( size_t k = 0; k < 126; k++)
        {
            int j=120;
            double x,y,z;
            x=reference[i][k][X_];
            y=reference[i][k][Y_];
            z=Layer.z1_*2+0.002+reference[i][k][Z_];

            double U = (reference[i][k][U_] - mean_curr)*scale;
            glColor3f(-U ,-U ,U );
            glVertex3f(x,y,z);

            x=reference[i+1][k][X_];
            y=reference[i+1][k][Y_];
            z=Layer.z1_*2+0.002+reference[i+1][k][Z_];

            U  = (reference[i+1][k][U_] - mean_curr)*scale;
            glColor3f(-U ,-U ,U);
            glVertex3f(x,y,z);
        }
        glEnd();
    }
    if(drawArr)
    {
        posVelAccelVort cell;
        static double arr[240][240];

        for( size_t i = 0; i < 240; i++)
        {
            for( size_t k = 0; k < 80; k++)
            {
                int j=120;
                double x,y,z;
                x=Layer.x0_+(i*1.0/240)*(Layer.x1_-Layer.x0_);
                y=Layer.y0_+(j*1.0/240)*(Layer.y1_-Layer.y0_);
                z=Layer.z0_+(k*1.0/80)*(Layer.z1_-Layer.z0_);

                cell = interpolateFull(x,y,z);
                arr[i][k]=-cell.rotvorty;
            }
        }


        glPointSize(7);
        glLineWidth(3);
        for( size_t i = 1; i < 239; i++) {
            glBegin(GL_TRIANGLE_STRIP);
            for( size_t k = 0; k < 80; k++) {
                int j=120;
                double x,y,z;
                x=Layer.x0_+(i*1.0/240)*(Layer.x1_-Layer.x0_);
                y=Layer.y0_+(j*1.0/240)*(Layer.y1_-Layer.y0_);
                z=Layer.z1_+(k*1.0/80)*(Layer.z1_-Layer.z0_)+0.001;

                glColor3f(arr[i][k]*scale,arr[i][k]*scale,-arr[i][k]*scale);
                glVertex3f(x,y,z);

                x=Layer.x0_+((i+1)*1.0/240)*(Layer.x1_-Layer.x0_);
                y=Layer.y0_+(j*1.0/240)*(Layer.y1_-Layer.y0_);
                z=Layer.z1_+(k*1.0/80)*(Layer.z1_-Layer.z0_)+0.001;

                glColor3f(arr[i+1][k]*scale,arr[i+1][k]*scale,-arr[i+1][k]*scale);
                glVertex3f(x,y,z);
            }
            glEnd();
        }
    }

    if( redr) {
        solveVortEq();
        glutPostRedisplay();
    }
    glutSwapBuffers();
}

void m_m(int x,int y)
{
    if (rotate==1) {
        rx=rx0+0.1*(x-mx0);
        ry=ry0+0.1*(y-my0);
    }
    glutPostRedisplay();
}

void m_d(int button, int state,int x, int y)
{
    if (state==GLUT_UP) {
        rotate=0;
        rx0=rx;
        ry0=ry;
    }
    if (state==GLUT_DOWN) {
        rotate=1;
        mx0=x;
        my0=y;
    }
    mouse_x=(1.0*x)/W_WIDTH;
    mouse_y=(W_HEIGHT-(1.0*y))/W_HEIGHT;
    glutPostRedisplay();
}

void setNiNjNk(int opticell, double j_in_i,double k_in_i) //opricell is the desired number of points in the cell
{
    int num_cells=Layer.currValues.size()/opticell;

    double ni= pow(num_cells/(j_in_i*k_in_i),0.3333333);
    double nj= ni*j_in_i;
    double nk= ni*k_in_i;

    Layer.ni=int (ni);
    Layer.nj=int (nj);
    Layer.nk=int (nk);

    printf("ni=%d nj=%d nk=%d \n",Layer.ni,Layer.nj,Layer.nk);
}

void loadCurrValuesFromDNS(int i)
{
    xmin = 1e10;
    xmax = -1e10;
    ymin = 1e10;
    ymax = -1e10;
    zmin = 1e10;
    zmax = -1e10;
    int n=0; int nn=0;
    double xx, yy, zz;
    double ux, uy, uz;
    double ax, ay, az;
    double pp;
    int num;

    Layer.currValues.clear();
    {
        FILE *file_data = fopen("../before/0_005_41_.dat", "r");
        //FILE *file_data = fopen("../before/0_05_41_.dat", "r");
        //FILE *file_data = fopen("../before/0_12_41_.dat", "r");
        //FILE *file_data = fopen("../before/0_2_41_.dat", "r");

        printf ("HERe! \n");

        int mass=0;

        while((!feof (file_data)))
        {
            fscanf(file_data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", &xx, &yy, &zz, &ux, &uy, &uz, &ax, &ay, &az, &pp, &num);
            if(zz<20.0) {
                mass++;
                posVelAccelVort pVec(xx,yy,zz,
                                     ux, uy, uz,
                                     ax, ay, az, 0,0,0);

                Layer.currValues.push_back(pVec);

                xmin = xx < xmin ? xx : xmin;
                xmax = xx > xmax ? xx : xmax;
                ymin = yy < ymin ? yy : ymin;
                ymax = yy > ymax ? yy : ymax;
                zmin = zz < zmin ? zz : zmin;
                zmax = zz > zmax ? zz : zmax;
                nn++;
            }
            n++;
        }
        fclose(file_data);
        pnum_in=Layer.currValues.size();
    }

    int m=0;
    int dn=nn/16;
    while(  m<dn)
    {
        xx=xmin+(rand()*1.0/RAND_MAX)*(xmax-xmin);
        yy=ymin+(rand()*1.0/RAND_MAX)*(ymax-ymin);
        zz=0.0;
        ux=0.0; uy=0.0; uz=0.0; ax=0.0; ay=0.0;az=0.0;
        posVelAccelVort pVec(xx, yy, zz,
                             ux, uy, uz,
                             ax, ay, az,0,0,0);

        pVec.isBound = true;

        Layer.currValues.push_back(pVec);

        xmin = xx < xmin ? xx : xmin;
        xmax = xx > xmax ? xx : xmax;
        ymin = yy < ymin ? yy : ymin;
        ymax = yy > ymax ? yy : ymax;
        zmin = zz < zmin ? zz : zmin;
        zmax = zz > zmax ? zz : zmax;
        m++;
        nn++;
    }

    printf(" xMin=%f  xMax=%f yMin=%f yMax=%f zMin=%f zMax=%f n=%d  \n", xmin,xmax,ymin,ymax,zmin,zmax,n);

    Layer.x0_=xmin-0.001*(xmax-xmin);
    Layer.x1_=xmax+0.001*(xmax-xmin);

    Layer.y0_=ymin-0.001*(ymax-ymin);
    Layer.y1_=ymax+0.001*(ymax-ymin);

    Layer.z0_=zmin-0.001*(zmax-zmin);
    Layer.z1_=zmax+0.001*(zmax-zmin);

    setNiNjNk(24, 16.0/30.0 , 6.0/30.0);
    findNeighbors_coarse_smart();
    printf("Layer.rad = %e\n",Layer.rad);

    calcVort();
    findNeighbors_coarse_smart();
    calcRotVort();
}

int removeWorst(int* arr, int num,int base,int desiredNum)
{
    if (num<=desiredNum)
    {
        printf("num<desiredNum = %d \n",num);
        return num;
    }
    static double f1[600]; //field function from each other
    static double f2[600]; //field function from base

    for (int i=0;i<num;i++)
    {
        int in1=arr[i];
        double dx=Layer.currValues[base].x-Layer.currValues[in1].x;
        double dy=Layer.currValues[base].y-Layer.currValues[in1].y;
        double dz=Layer.currValues[base].z-Layer.currValues[in1].z;
        f2[i] = (dx*dx + dy*dy +dz*dz);
    }

    double fmax=0.0;
    int i_max;

    for (int i=0;i<num;i++) {
        f1[i]=0.0;
        for (int j=0;j<num;j++){
            if (i!=j) {
                int in1=arr[i];
                int in2=arr[j];
                double dx=Layer.currValues[in2].x-Layer.currValues[in1].x;
                double dy=Layer.currValues[in2].y-Layer.currValues[in1].y;
                double dz=Layer.currValues[in2].z-Layer.currValues[in1].z;
                f1[i]+= 1.0/(dx*dx + dy*dy +dz*dz);
            }
        }
        if (f1[i]*f2[i]>fmax) {
            fmax=f1[i]*f2[i];
            i_max=i;
        }
    }
    int curr_num=num;
    while (curr_num>desiredNum)
    {
        int in_max=arr[i_max];
        arr[i_max]=arr[curr_num-1];
        f1[i_max]=f1[curr_num-1];
        f2[i_max]=f2[curr_num-1];
        curr_num--;
        fmax=0.0;
        for (int i=0;i<curr_num;i++)
        {
            int in1=arr[i];
            double dx=Layer.currValues[in_max].x-Layer.currValues[in1].x;
            double dy=Layer.currValues[in_max].y-Layer.currValues[in1].y;
            double dz=Layer.currValues[in_max].z-Layer.currValues[in1].z;
            f1[i]-= 1.0/(dx*dx + dy*dy +dz*dz);
            if (f1[i]*f2[i]>fmax) {
                fmax=f1[i]*f2[i];
                i_max=i;
            }
        }
    }
    return curr_num;
}


int removeWorst_free(int* arr, int num,double x, double y,double z,int desiredNum)
{
    if (num<=desiredNum){
        printf("num<desiredNum = %d \n",num);
        return num;
    }

    static double f1[600]; //field function from each other
    static double f2[600]; //field function from base

    for (int i=0;i<num;i++) {
        int in1=arr[i];
        double dx=x-Layer.currValues[in1].x;
        double dy=y-Layer.currValues[in1].y;
        double dz=z-Layer.currValues[in1].z;
        f2[i] = (dx*dx + dy*dy +dz*dz);
    }

    double fmax=0.0;
    int i_max;

    for (int i=0;i<num;i++) {
        f1[i]=0.0;
        for (int j=0;j<num;j++) {
            if (i!=j) {
                int in1=arr[i];
                int in2=arr[j];
                double dx=Layer.currValues[in2].x-Layer.currValues[in1].x;
                double dy=Layer.currValues[in2].y-Layer.currValues[in1].y;
                double dz=Layer.currValues[in2].z-Layer.currValues[in1].z;
                f1[i]+= 1.0/(dx*dx + dy*dy +dz*dz);
            }
        }
        if (f1[i]*f2[i]>fmax) {
            fmax=f1[i]*f2[i];
            i_max=i;
        }
    }

    int curr_num=num;
    while (curr_num>desiredNum)
    {
        int in_max=arr[i_max];
        arr[i_max]=arr[curr_num-1];
        f1[i_max]=f1[curr_num-1];
        f2[i_max]=f2[curr_num-1];
        curr_num--;
        fmax=0.0;
        for (int i=0;i<curr_num;i++)
        {
            int in1=arr[i];
            double dx=Layer.currValues[in_max].x-Layer.currValues[in1].x;
            double dy=Layer.currValues[in_max].y-Layer.currValues[in1].y;
            double dz=Layer.currValues[in_max].z-Layer.currValues[in1].z;
            f1[i]-= 1.0/(dx*dx + dy*dy +dz*dz);
            if (f1[i]*f2[i]>fmax) {
                fmax=f1[i]*f2[i];
                i_max=i;
            }
        }
    }
    return curr_num;
}

void findNeighbors_coarse_smart() //11 template
{
    double grid_dx = (Layer.x1_ - Layer.x0_) * 1.0 / (Layer.ni);
    double grid_dy = (Layer.y1_ - Layer.y0_) * 1.0 / (Layer.nj);
    double grid_dz = (Layer.z1_ - Layer.z0_) * 1.0 / (Layer.nk);

    Layer.rad=0.25*pow(grid_dx*grid_dy*grid_dz,0.33333);//0.45

    printf("rad=%f \n",Layer.rad);

    for (size_t i = 0; i < Layer.ni; i++)
        for (size_t j = 0; j < Layer.nj; j++)
            for (size_t k = 0; k < Layer.nk; k++)
                Layer.particlesGrid[i][j][k].clear();

    for (int i = 0; i < Layer.currValues.size(); i++)
    {
        int xIdx = int((Layer.currValues[i].x - Layer.x0_)/grid_dx);
        int yIdx = int((Layer.currValues[i].y - Layer.y0_)/grid_dy);
        int zIdx = int((Layer.currValues[i].z - Layer.z0_)/grid_dz);

        if(xIdx> Layer.ni-1)
            xIdx =Layer.ni-1;
        if(yIdx> Layer.nj-1)
            yIdx =Layer.nj-1;
        if(zIdx> Layer.nk-1)
            zIdx =Layer.nk-1;

        Layer.particlesGrid[xIdx][yIdx][zIdx].push_back(i);
    }

    int size_max=0;
    for (int s = 0; s < Layer.currValues.size(); s++) {
        Layer.currValues[s].neighbours.clear();
        int xIdx = int((Layer.currValues[s].x - Layer.x0_)/grid_dx);
        int yIdx = int((Layer.currValues[s].y - Layer.y0_)/grid_dy);
        int zIdx = int((Layer.currValues[s].z - Layer.z0_)/grid_dz);
        int im = fmax(0, xIdx-1);
        int ip = fmin(Layer.ni - 1, xIdx+1);
        int jm = fmax(0, yIdx-1);
        int jp = fmin(Layer.nj - 1, yIdx+1);
        int km = fmax(0, zIdx-1);
        int kp = fmin(Layer.nk - 1, zIdx+1);

        int indices[500]; //indices for sorting;
        int indNum=0;

        for (int i = im; i <= ip ; i++)
            for (int j = jm; j <=jp; j++)
                for (int k = km; k <=kp; k++)
                    for(int n = 0; n < Layer.particlesGrid[i][j][k].size(); n++)
                    {
                        int i_nb=Layer.particlesGrid[i][j][k].at(n);
                        double r2 = (Layer.currValues[s].x - Layer.currValues[i_nb].x) * (Layer.currValues[s].x - Layer.currValues[i_nb].x)
                                +   (Layer.currValues[s].y - Layer.currValues[i_nb].y) * (Layer.currValues[s].y - Layer.currValues[i_nb].y)
                                +   (Layer.currValues[s].z - Layer.currValues[i_nb].z) * (Layer.currValues[s].z - Layer.currValues[i_nb].z);
                        if(r2 < fmax((16.0*0.05*0.05),16*Layer.rad * Layer.rad))//if(r2 <((4.0*Layer.rad) * (4.0*Layer.rad)))
                        {
                            indices[indNum]=i_nb;
                            indNum++;
                        }
                    }

        indNum=removeWorst(indices,indNum,s,desiredStencilNum);
        if (s%1000==1)
            printf("%d \% \n",int(s*100.0/Layer.currValues.size()));

        for (int i=0;i<indNum;i++)
        {
            int i_nb=indices[i];
            Layer.currValues[s].neighbours.push_back(i_nb);
        }

        lssA.get_vortEq_matr(Layer.currValues[s],Layer.rad);

        if(Layer.currValues[s].neighbours.size()>size_max) size_max=Layer.currValues[s].neighbours.size();
        if(Layer.currValues[s].neighbours.size()<11) printf("SIZE<11 ----------------------------------------- \n");
        if(Layer.currValues[s].neighbours.size()>=MAX_EQNS) printf("SIZE>=MAX_EQNS ----------------------------------------- \n");
    }
    printf("done neigh max size=%d \n",size_max);
}

void calcVort()
{
    for (int s = 0; s < Layer.currValues.size(); s++)
    {
        posVelAccelVort vort = lssA.interpFull(Layer.currValues[s],Layer.rad);
        Layer.currValues[s].vortx = vort.dwy-vort.dvz;
        Layer.currValues[s].vorty = vort.duz-vort.dwx;
        Layer.currValues[s].vortz = vort.dvx-vort.duy;

        if (Layer.currValues[s].neighbours.size()<11){
            printf("s=%d x=%f y=%f z=%f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  < 11 \n",s,Layer.currValues[s].x,Layer.currValues[s].y,Layer.currValues[s].z);
            Layer.currValues[s].vortx=0.0;
            Layer.currValues[s].vorty=0.0;
            Layer.currValues[s].vortz=0.0;
        }
        if (Layer.currValues[s].neighbours.size()>5000){
            printf("s=%d x=%f y=%f z=%f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  > 5000 \n",s,Layer.currValues[s].x,Layer.currValues[s].y,Layer.currValues[s].z);
            Layer.currValues[s].vortx=0.0;
            Layer.currValues[s].vorty=0.0;
            Layer.currValues[s].vortz=0.0;
        }
        if (s%5000==0)
            printf("S=%d \n",s);
    }

    printf("vort calculating finished\n");
}

void calcRotVort()
{
    for (int s = 0; s < Layer.currValues.size(); s++)
    {
        v3 rotvort = lssA.getRotVort(Layer.currValues[s],Layer.rad);
        Layer.currValues[s].rotvortx = rotvort.x;
        Layer.currValues[s].rotvorty = rotvort.y;
        Layer.currValues[s].rotvortz = rotvort.z;

        if (Layer.currValues[s].neighbours.size()<11){
            printf("s=%d x=%f y=%f z=%f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  < 11 \n",s,Layer.currValues[s].x,Layer.currValues[s].y,Layer.currValues[s].z);
            Layer.currValues[s].rotvortx =0.0;
            Layer.currValues[s].rotvorty =0.0;
            Layer.currValues[s].rotvortz =0.0;
        }
        if (Layer.currValues[s].neighbours.size()>5000){
            printf("s=%d x=%f y=%f z=%f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  > 5000 \n",s,Layer.currValues[s].x,Layer.currValues[s].y,Layer.currValues[s].z);
            Layer.currValues[s].rotvortx =0.0;
            Layer.currValues[s].rotvorty =0.0;
            Layer.currValues[s].rotvortz =0.0;
        }
        if (s%5000==0)
            printf("S=%d \n",s);
    }
    printf("rotVort calculating finished\n");
}

posVelAccelVort interpolateFull(double x,double y, double z)
{
    double grid_dx = (Layer.x1_ - Layer.x0_) * 1.0 / (Layer.ni);
    double grid_dy = (Layer.y1_ - Layer.y0_) * 1.0 / (Layer.nj);
    double grid_dz = (Layer.z1_ - Layer.z0_) * 1.0 / (Layer.nk);
    int xIdx = int((x - Layer.x0_)/grid_dx);
    int yIdx = int((y - Layer.y0_)/grid_dy);
    int zIdx = int((z - Layer.z0_)/grid_dz);
    int im = fmax(0, xIdx-1);
    int ip = fmin(Layer.ni - 1, xIdx+1);
    int jm = fmax(0, yIdx-1);
    int jp = fmin(Layer.nj - 1, yIdx+1);
    int km = fmax(0, zIdx-1);
    int kp = fmin(Layer.nk - 1, zIdx+1);

    double xr=1.0;
    double yr=1.0;
    double zr=1.0;
    int indices[500]; //indices for sorting;
    int indNum=0;

    for (int i = im; i <= ip ; i++)
        for (int j = jm; j <=jp; j++)
            for (int k = km; k <=kp; k++)
                for(int n = 0; n < Layer.particlesGrid[i][j][k].size(); n++)
                {
                    int i_nb=Layer.particlesGrid[i][j][k].at(n);
                    double r2 = (x - Layer.currValues[i_nb].x)/xr * (x - Layer.currValues[i_nb].x)/xr
                            +   (y - Layer.currValues[i_nb].y)/yr * (y - Layer.currValues[i_nb].y)/yr
                            +   (z - Layer.currValues[i_nb].z)/zr * (z - Layer.currValues[i_nb].z)/zr;
                    if(r2 < fmax((16.0*0.05*0.05),16*Layer.rad * Layer.rad))
                    {
                        indices[indNum]=i_nb;
                        indNum++;
                    }
                }
    indNum=removeWorst_free(indices,indNum,x,y,z,desiredStencilNum);
    posVelAccelVort n;
    for (int i=0;i<indNum;i++)
    {
        int i_nb=indices[i];
        n.neighbours.push_back(i_nb);
    }

    n.x =x;
    n.y = y;
    n.z = z;
    return lssA.interpFull(n,Layer.rad);
}

void saveInTecplot(char* fname)
{
    int XN = 251;
    int YN = 51;
    int ZN = 126;

    double xmin_ = 1.75;
    double xmax_ = 6.75;
    double ymin_ = -1.25;
    double ymax_ = 1.25;
    double zmin_ = 0.001;
    double zmax_ = 1.0;

    FILE *file = fopen(fname, "w");

    fprintf(file, "TITLE =\"  \" \n VARIABLES=\"x\" \n \"y\"\n \"z\" \n \"Ux\" \n \"Uy\" \n \"Uz\" \n \"Ux_x\" \n \"Ux_y\" \n \"Ux_z\" \n \"Uy_x\" \n \"Uy_y\" \n \"Uy_z\" \n \"Uz_x\" \n \"Uz_y\" \n \"Uz_z\" \n \"ax\" \n \"ay\" \n \"az\" \n\"P\" \n ZONE T=\"  \" \n");
    fprintf(file," I=%d J=%d K=%d F=POINT \n", XN,YN, ZN);

    //calc_mean_pressure
    double mean_p_save=0.0;
    double mass=0.0;

    for( int k = 0; k < ZN; k++ ){
        int j=26;
        for( int i = 0; i < XN; i++ ){

            double x = xmin_ + i * (xmax_ - xmin_) / (XN-1);
            double y = ymin_ + j * (ymax_ - ymin_) / (YN-1);
            double z = zmin_ + k * (zmax_ - zmin_) / (ZN-1);
            posVelAccelVort res = interpolateFull(x, y, z);
            mass++;
        }
    }
    mean_p_save/=mass;

    for( int k = 0; k < ZN; k++ ){
        printf("k=%d\n", k);
        for( int j = 0; j < YN; j++ ){
            for( int i = 0; i < XN; i++ ){
                double x = xmin_ + i * (xmax_ - xmin_) / (XN-1);
                double y = ymin_ + j * (ymax_ - ymin_) / (YN-1);
                double z = zmin_ + k * (zmax_ - zmin_) / (ZN-1);
                posVelAccelVort res = interpolateFull( x, y, z);

                fprintf(file,"%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",
                        x, y, z, res.u, res.v, res.w,0,0,0,0,0,0,0,0,0, res.ax, res.ay, res.az);
            }
        }
    }
    fclose(file);

}

posVelAccelVort interpolate_smart(double x,double y, double z)
{
    double grid_dx = (Layer.x1_ - Layer.x0_) * 1.0 / (Layer.ni);
    double grid_dy = (Layer.y1_ - Layer.y0_) * 1.0 / (Layer.nj);
    double grid_dz = (Layer.z1_ - Layer.z0_) * 1.0 / (Layer.nk);
    int xIdx = int((x - Layer.x0_)/grid_dx);
    int yIdx = int((y - Layer.y0_)/grid_dy);
    int zIdx = int((z - Layer.z0_)/grid_dz);
    int im = fmax(0, xIdx-1);
    int ip = fmin(Layer.ni - 1, xIdx+1);
    int jm = fmax(0, yIdx-1);
    int jp = fmin(Layer.nj - 1, yIdx+1);
    int km = fmax(0, zIdx-1);
    int kp = fmin(Layer.nk - 1, zIdx+1);
    vector<posVelAccelVort> m_p;

    double xr=1.0;
    double yr=1.0;
    double zr=1.0;
    int indices[500]; //indices for sorting;
    int indNum=0;

    for (int i = im; i <= ip ; i++)
        for (int j = jm; j <=jp; j++)
            for (int k = km; k <=kp; k++)
                for(int n = 0; n < Layer.particlesGrid[i][j][k].size(); n++)
                {
                    int i_nb=Layer.particlesGrid[i][j][k].at(n);
                    double r2 = (x - Layer.currValues[i_nb].x)/xr * (x - Layer.currValues[i_nb].x)/xr
                            +   (y - Layer.currValues[i_nb].y)/yr * (y - Layer.currValues[i_nb].y)/yr
                            +   (z - Layer.currValues[i_nb].z)/zr * (z - Layer.currValues[i_nb].z)/zr;
                    if(r2 < fmax((16.0*0.05*0.05),16*Layer.rad * Layer.rad)) {
                        indices[indNum]=i_nb;
                        indNum++;
                    }
                }

    indNum=removeWorst_free(indices,indNum,x,y,z,desiredStencilNum);
    posVelAccelVort n;
    for (int i=0;i<indNum;i++)
    {
        int i_nb=indices[i];
        n.neighbours.push_back(i_nb);
    }


    n.x =x;
    n.y = y;
    n.z = z;

    lssA.interp(n,Layer.rad);
    posVelAccelVort res(x,y,z,0,0,0,0,0,0,0,0,0);
    res.vortx = n.vortx;
    res.vorty = n.vorty;
    res.vortz = n.vortz;
    return res;
}

void solveVortEq()
{
    posVelAccelVort n;
    deriv3D dd;
    for (int s = 0; s < Layer.currValues.size(); s++){
        Layer.updateNeighbourFByVar(Layer.currValues[s], mLayer::kU);

        n.x = Layer.currValues[s].x;
        n.y = Layer.currValues[s].y;
        n.z = Layer.currValues[s].z;
        n.f = Layer.currValues[s].u;
        n.rhs=-Layer.currValues[s].rotvortx;
        n.neighbours = Layer.currValues[s].neighbours;

        lssA.get_vortEq_combo(n,dd,Layer.rad);
        Layer.currValues[s].u =  n.f;
    }
    //printf("done\n");
}

int iii;
void solveVortEq_nomatr()
{
    iii++;
    posVelAccelVort n;
    double res=0.0;
    double res2=0.0;//all below is for residuals
    for (int s = 0; s < Layer.currValues.size(); s++)
    {
        if(Layer.currValues[s].isBound)
            continue;
        Layer.updateNeighbourFByVar(Layer.currValues[s], mLayer::kU);

        n.x = Layer.currValues[s].x;
        n.y = Layer.currValues[s].y;
        n.z = Layer.currValues[s].z;
        n.f = Layer.currValues[s].u;
        n.rhs=-Layer.currValues[s].rotvortx;
        n.neighbours = Layer.currValues[s].neighbours;
    }

    for (int s = 0; s < Layer.currValues.size(); s++)
    {
        if(Layer.currValues[s].isBound)
            continue;
        Layer.updateNeighbourFByVar(Layer.currValues[s], mLayer::kU);

        n.x = Layer.currValues[s].x;
        n.y = Layer.currValues[s].y;
        n.z = Layer.currValues[s].z;
        n.f = Layer.currValues[s].u;
        n.rhs=-Layer.currValues[s].rotvortx;
        n.neighbours = Layer.currValues[s].neighbours;

        for (int i=0;i<VAR_NUM;i++)
            for (int j=0;j<VAR_NUM;j++)
                n.INV[i][j]=Layer.currValues[s].INV[i][j];

        lssA.get_vortEq_nomatr(n,Layer.rad);
        res+=abs(Layer.currValues[s].u - n.f);
        Layer.currValues[s].u = n.f;
    }
    res/=Layer.currValues.size();
    res2/=Layer.currValues.size();
    printf("i=%d res=%e res2=%e \n",iii, res, res2);
    /*FILE* ff = fopen("res_0.2.txt","a");
    fprintf(ff,"%d %e %e \n",iii, res, res2);
    fclose(ff);*/
    /*double p0=Layer.currValues[0].p;
    for (int s = 0; s < Layer.currValues.size(); s++)
    {
        Layer.currValues[s].p-=p0;
    }*/
}

void kb(unsigned char key, int x, int y){
    if (key=='1') {
        drawArr = !drawArr;
    }
    if (key=='2') {

    }
    if (key=='3') {

    }
    if (key=='4') {
        drawVelocity = !drawVelocity;
    }
    if (key=='5') {

    }
    if (key==' ') {
        redr=!redr;
    }
    if  (key=='x'){
        saveInTecplot("outTecFile.dat");
    }
    if (key=='['){
        scale/=1.2;
    }
    if (key==']') {
        scale*=1.2;
    }
    if (key=='w') {
        view_x+=(o_x)*0.05;
        view_y+=(o_y)*0.05;
        view_z+=(o_z)*0.05;
    }
    if (key=='s') {
        view_x-=(o_x)*0.05;
        view_y-=(o_y)*0.05;
        view_z-=(o_z)*0.05;
    }
    if (key=='q') {
        view_z+=0.05;
    }
    if (key=='e') {
        view_z-=0.05;
    }
    if (key=='a') {
        double l2=sqrt(o_y*o_y+o_x*o_x);
        view_y+=(o_x)*0.05/l2;
        view_x+=-(o_y)*0.05/l2;
    }
    if (key=='d') {
        double    l2=sqrt(o_y*o_y+o_x*o_x);
        view_x+=(o_y)*0.05/l2;
        view_y+=-(o_x)*0.05/l2;
    }
    glutPostRedisplay();
}

void init()
{
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluPerspective(45.0f, W_WIDTH*1.0/W_HEIGHT, 0.1f, 2500.0f);
    glMatrixMode (GL_MODELVIEW);

    loadReference();

    loadCurrValuesFromDNS(currTime);

    view_x=0.5*(Layer.x1_+Layer.x0_);
    view_y=Layer.y0_ -2.0*(Layer.y1_-Layer.y0_);
    view_z=0.5*(Layer.z1_+Layer.z0_);
    printf("Viewx=%f viewy=%f viewz=%f \n",view_x,view_y,view_z);
}

int main(int argc, char** argv)
{
    glutInit(&argc,argv);
    //glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);

    glutInitWindowSize(W_HEIGHT*(1-0)/(1-0),W_HEIGHT);
    glutInitWindowPosition(0,0);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
    glutMotionFunc(m_m);
    glutMouseFunc(m_d);
    glutKeyboardFunc(kb);
    init();
    glutMainLoop();
}
