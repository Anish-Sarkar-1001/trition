#include "main.h"
#include "globals.h"

typedef struct Data_{
    /*
    *Main array containing variable values 
    *in the order var;k;j;i*/
    double ****vCentroid;
    double ****vStaggered;
    double ****advecFlux;
}Data;

typedef struct timeData_ {
    double ****vCentroid;
    double ****vStaggered;
    int iPos;
    int jPos;
    int kPos;
}tData;

typedef struct midT_ {
    double ****vS;
    double ****vC;
    int iPos;
    int jPos;
    int kPos;
}midT;

typedef struct Grid_{
    int npTotal[3]; /*Array for number of total points in each direction*/
    double xBeg[3]; /*Array for Beginning points in each direction*/
    double xEnd[3]; /*Array for Ending points in each direciton*/
    int xSpacing[3]; /*Type of grid spacing*/
    double **xPoints;
    double *xCentroids[3];
    double *xArea[3];
    double *dCentroid[3];
    double *dPoints[3];
    double ****Area;
    #if STAGGERED_GRID == YES
    double ****v1StagArea;
    double ****v2StagArea;
    double ****v3StagArea;
    #else
    #endif
    int iPos;
    int jPos;
    int kPos;
}Grid;

#if STAGGERED_GRID == YES
typedef struct Staggered_Grid_{
    double **xStaggeredPoints;
    double *dStaggeredPoints[3];
}StagGrid;
#endif

typedef struct Boundary_{
    double ****bCentroids; /*FIrst index points to variable,
                            2nd index points to face number,
                            3rd index points to perpendicular 1 value
                            4th index points to perpendicular 2 value*/
    #if SATGGERED_GRID == YES
        double ****b0Stag;
        double ****b1_2Stag;
        double ****b1Stag;
    #else
    double ****b0Stag;
    double ****b1_2Stag;
    double ****b1Stag;
    #endif
}Boundary;
