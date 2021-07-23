#ifndef GLOBALS_H
#define GLOBALS_H

#ifdef __linux__
#include  <GL/gl.h>
#include  <GL/glu.h>
#include  <GL/glut.h>
#elif _WIN32
#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <string.h>
#include <fstream>
#include "dirent.h"
#include <thread>
#include <sys/time.h>


#define W_WIDTH 1600
#define W_HEIGHT 1300
#define NX 251
#define NY 126
#define NZ 76


//Lamb:
//#define NUMCELLS_X 8
//#define NUMCELLS_Y 8
//#define NUMCELLS_Z 1


//vlad_high:
//#define NUMCELLS_X 45
//#define NUMCELLS_Y 30
//#define NUMCELLS_Z 6

//vlad_med: // mod5
/*#define NUMCELLS_X 30
#define NUMCELLS_Y 16
#define NUMCELLS_Z 6*/


#define NUMCELLS_X 150
#define NUMCELLS_Y 80
#define NUMCELLS_Z 30



//vald_low //mod 16
/*#define NUMCELLS_X 15
#define NUMCELLS_Y 8
#define NUMCELLS_Z 3*/

/*#define NUMCELLS_X 10
#define NUMCELLS_Y 5
#define NUMCELLS_Z 2*/

//piv_high:
//#define NUMCELLS_X 20
//#define NUMCELLS_Y 15
//#define NUMCELLS_Z 10


//piv_med:
/*#define NUMCELLS_X 20
#define NUMCELLS_Y 10
#define NUMCELLS_Z 5*/


//piv_low:
/*#define NUMCELLS_X 12
#define NUMCELLS_Y 6
#define NUMCELLS_Z 4*/


/*#define x0_ 1.75//-2.5//-50
#define x1_ 6.75//2.5//50
#define y0_ -1.25//-25
#define y1_ 1.25//25
#define z0_ 0.04//0.01
#define z1_ 1.0//30.01

#define dx ((x1_-x0_) * 1.0 / NX)
#define dy ((y1_-y0_) * 1.0 / NY)
#define dz ((z1_-z0_) * 1.0 / NZ)
*/
#define THREADNUM 16
#define PATH "/home/user/PIV_CHALLENGE_2020/paper_PIV/data/before/"//"DA_ppp_0_025/"//"E:/projects/splineViewer/0.080ppp.txt"//"DA_ppp_0_160/"//"NoisyTracks/"//"ExactTracks/"//"NoisyTracks/"//"ExactTracks/"//"NoisyTracks/"//"ExactTracks/"//"NoisyTracks/"//"E:/projects/splineViewer/ActiveLongTracks50_0.080ppp_best.txt"
//"ExactTracks/"//"E:/projects/splineViewer/ActiveLongTracks50_0.050ppp_best.txt"
//"DA_ppp_0_025/""DA_ppp_0_160/"//"DA_ppp_0_005/""test/""testNoise10x/""testNoise100x/"

//number of variables
#define VAR_NUM 10

//maximum number of points
#define MAX_EQNS 1000

using namespace std;

struct v3
{
    double x, y, z;
    v3(double ix, double iy, double iz)
        : x(ix), y(iy), z(iz) {}
};

struct posVelAccelVort
{
    posVelAccelVort(){}
    posVelAccelVort(double ix, double iy, double iz, double iu, double iv, double iw, double iax, double iay, double iaz, double ivortx, double ivorty, double ivortz)
        : x(ix), y(iy), z(iz), u(iu), v(iv), w(iw), ax(iax), ay(iay), az(iaz), vortx(ivortx), vorty(ivorty), vortz(ivortz) {}
    double x, y, z;
    double u, v, w;
    double ax, ay, az;
    double vortx, vorty, vortz;
    double rotvortx, rotvorty, rotvortz;
    double dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz;
    vector<int> neighbours;
    double f;
    double rhs;
    int my_ind;//this point location in m_p;
    double INV[VAR_NUM][VAR_NUM];
    bool isBound = false;
};

struct mLayer
{
    int ni,nj,nk;
    double rad;
    vector<int> particlesGrid[NUMCELLS_X][NUMCELLS_Y][NUMCELLS_Z];
    vector<posVelAccelVort> currValues;
    double x0_,y0_,z0_,x1_,y1_,z1_;
    enum varName
    {

        kU,
        kV,
        kW,
        kAX,
        kAY,
        kAZ,
        kVortX,
        kVortY,
        kVortZ
    };
    void updateNeighbourFByVar(posVelAccelVort iN, varName iName);
};

double get_time(void);
void get_color(double gval, double min, double max);

extern mLayer Layer;
extern double xmin;
extern double xmax;
extern double ymin;
extern double ymax;
extern double zmin;
extern double zmax;

extern int clear_w;
extern int mx0,my0;
extern int rotate;
extern float rx0;
extern float ry0;
extern float rx;
extern float ry;
extern double mouse_x,mouse_y;
extern int redr;
extern double ck;
extern double scale;
extern double view_x;
extern double view_y;
extern double view_z;
extern double o_x;
extern double o_y;
extern double o_z;
extern double power;

extern int iNum;
extern int jNum;
extern int kNum;
extern int currTime;
extern int itn;
extern int kCur;
extern double sc;
extern double cv;


extern bool drawVelocity;
extern bool drawAcceleration;
extern bool drawArr;

extern double minVel;
extern double maxVel;

extern double minAccel;
extern double maxAccel;

extern double dx,dy,dz;

extern double minPressure;
extern double maxPressure;

extern bool rangeCalculated;

extern int isNoOptimized[THREADNUM];


#endif // GLOBALS_H
