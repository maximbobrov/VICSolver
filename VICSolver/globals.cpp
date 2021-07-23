#include "globals.h"



mLayer Layer;
double xmin = 1e10;
double xmax = -1e10;
double ymin = 1e10;
double ymax = -1e10;
double zmin = 1e10;
double zmax = -1e10;

int clear_w = 1.0;

int mx0,my0;
int rotate = 0;
float rx0 = -90.0;
float ry0 = 0.0;
float rx = rx0;
float ry = ry0;
double mouse_x,mouse_y;
int redr=0;
double ck=0.1;
double scale  = 8.5;
double view_x=14.507903;
double view_y=8.300000;
double view_z=24.2;
double o_x=0.0;
double o_y=0.0;
double o_z=0.0;
double power = 0.5;

int iNum = NX-1;
int jNum = NY-1;
int kNum = NZ-1;
int currTime = 7;
int itn=0;
int kCur= 0;
double sc=1;
double cv=0.001;

double dx,dy,dz;

bool drawVelocity = false;
bool drawAcceleration = false;
bool drawArr = false;

double minVel = 1e100;
double maxVel = -1e100;

double minAccel = 1e100;
double maxAccel = -1e100;

double minPressure = 1e100;
double maxPressure = -1e100;

bool rangeCalculated = false;

int isNoOptimized[THREADNUM];


double get_time(void) {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    return ((double)(tv.tv_sec+tv.tv_usec*1.0e-6));
}

void get_color(double gval, double min, double max)
{
    const int nn=4;
    int i;
    double val;
    val=gval;
    if (val>max) val=max;
    if (val<min) val=min;

    typedef struct {
        double x,y,z;
    } XYZ;

    XYZ col_table[5];

    col_table[0].x = 0.0; col_table[0].y = 0.0; col_table[0].z = 1.0;
    col_table[1].x = 0.0; col_table[1].y = 1.0; col_table[1].z = 1.0;
    col_table[2].x = 0.0; col_table[2].y = 1.0; col_table[2].z = 0.0;
    col_table[3].x = 1.0; col_table[3].y = 1.0; col_table[3].z = 0.0;
    col_table[4].x = 1.0; col_table[4].y = 0.0; col_table[4].z = 0.0;

    double alpha;
    if ((max-min) > 1e-35)
    {
        alpha=(val-min)/(max-min)*nn;
        i=(int)(alpha);
        alpha=alpha-i;
    }
    else
    {
        alpha=0.0;
        i=2;
    }
    glColor3f(col_table[i].x * (1 - alpha) + col_table[i+1].x * alpha, col_table[i].y * (1 - alpha) + col_table[i+1].y * alpha, col_table[i].z * (1 - alpha) + col_table[i+1].z * alpha);
}

void mLayer::updateNeighbourFByVar(posVelAccelVort iN, mLayer::varName iName)
{
    for (int i=0;i<iN.neighbours.size();i++)
    {
        if(iName == kU)
            Layer.currValues[iN.neighbours[i]].f=Layer.currValues[iN.neighbours[i]].u;
        if(iName == kV)
            Layer.currValues[iN.neighbours[i]].f=Layer.currValues[iN.neighbours[i]].v;
        if(iName == kW)
            Layer.currValues[iN.neighbours[i]].f=Layer.currValues[iN.neighbours[i]].w;
        if(iName == kAX)
            Layer.currValues[iN.neighbours[i]].f=Layer.currValues[iN.neighbours[i]].ax;
        if(iName == kAY)
            Layer.currValues[iN.neighbours[i]].f=Layer.currValues[iN.neighbours[i]].ay;
        if(iName == kAZ)
            Layer.currValues[iN.neighbours[i]].f=Layer.currValues[iN.neighbours[i]].az;
        if(iName == kVortX)
            Layer.currValues[iN.neighbours[i]].f=Layer.currValues[iN.neighbours[i]].vortx;
        if(iName == kVortY)
            Layer.currValues[iN.neighbours[i]].f=Layer.currValues[iN.neighbours[i]].vorty;
        if(iName == kVortZ)
            Layer.currValues[iN.neighbours[i]].f=Layer.currValues[iN.neighbours[i]].vortz;
    }
}






















