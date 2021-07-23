#ifndef LEASTSQUARESSOLVER_H
#define LEASTSQUARESSOLVER_H
#include <vector>
#include <globals.h>

typedef struct
{double d[VAR_NUM];} deriv3D;

typedef struct
{double nx,ny,nz;
    double d;//normal distance from  coordinate origin  i.e.  nx*x+ny*y+nz*z=d
} boundary_plane;

class leastSquaresSolver
{
private:


public:
    enum deriv_order
    {
        F,
        FX,
        FY,
        FZ,
        FXX,
        FXY,
        FXZ,
        FYY,
        FYZ,
        FZZ
    };

    leastSquaresSolver();

    //std::vector<node3d> *m_p; //points cloud

    std::vector<boundary_plane> walls; //all boundaries are here for now
    double pww;
    double rad;
    void LU_decompose(int size);
    void LU_decompose2(int size);
    void m_solve(int num);
    void m_invert(int size);

    void nullify_m(double m[MAX_EQNS][MAX_EQNS], int in, int im);
    void nullify_v(double v[MAX_EQNS], int in);
    double getDivA(posVelAccelVort &p, double delta);
    v3 getVort(posVelAccelVort &p, double delta);
    v3 getRotVort(posVelAccelVort &p, double delta);
    void interp(posVelAccelVort &p, double delta);
    posVelAccelVort interpFull(posVelAccelVort &p, double delta);
    double dist(posVelAccelVort &p1,posVelAccelVort &p2);
    double distDx(posVelAccelVort &p1,posVelAccelVort &p2);

    void get_vortEq_combo(posVelAccelVort &p, deriv3D &res, double delta);

    void get_vortEq_nomatr(posVelAccelVort &p,  double delta);
    void get_vortEq_matr(posVelAccelVort &p, double delta);
};


#endif // LEASTSQUARESSOLVER_H
