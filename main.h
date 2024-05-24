#define  PHYSICS            UNSTEADY_FLUID_FLOW
#define  COMPRESSIBLE       NO
#define  VISCOSITY          YES
#define  STAGGERED_GRID     NO

/*X1 Grid information*/
#define  X1BEG          1.0
#define  X1END          2.0
#define  X1POINTS       80
#define  X1SCALING      UNI

/*X2 Grid information*/
#define  X2BEG          1.0
#define  X2END          2.0
#define  X2POINTS       80
#define  X2SCALING      UNI

/*X3 Grid Information*/
#define  X3BEG          1.0
#define  X3END          2.0
#define  X3POINTS       1
#define  X3SCALING      UNI

/*Dimensions*/
#define  DIMENSIONS     2

/*Solver Information*/
#define  SOLVER         FINITE_VOL

/*TIme Information*/
#define  TIMESTEPPING   EXPLICIT
#define  DT             0.01
#define  TSTOP          10.0

/*Boundray Information*/
//Pressure
#define X1PBEGBOUNDARY   ZERO_GRAD
#define X1PENDBOUNDARY   ZERO_GRAD
#define X2PBEGBOUNDARY   ZERO_GRAD
#define X2PENDBOUNDARY   CONST
#define X3PBEGBOUNDARY   ZERO_GRAD
#define X3PENDBOUNDARY   ZERO_GRAD

//Velocity along x
#define X1UBEGBOUNDARY   CONST
#define X1UENDBOUNDARY   CONST
#define X2UBEGBOUNDARY   CONST
#define X2UENDBOUNDARY   CONST
#define X3UBEGBOUNDARY   ZERO_GRAD
#define X3UENDBOUNDARY   ZER0_GRAD

//Velocity along y
#define X1VBEGBOUNDARY   CONST
#define X1VENDBOUNDARY   CONST
#define X2VBEGBOUNDARY   CONST
#define X2VENDBOUNDARY   CONST
#define X3VBEGBOUNDARY   ZERO_GRAD
#define X3VENDBOUNDARY   ZERO_GRAD

//Velocity along z
#define X1WBEGBOUNDARY   CONST
#define X1WENDBOUNDARY   CONST
#define X2WBEGBOUNDARY   CONST
#define X2WENDBOUNDARY   CONST
#define X3WBEGBOUNDARY   CONST
#define X3WENDBOUNDARY   CONST

/*Boundary values*/
//Pressure
#define BXP1BEG          0.0
#define BXP1END          0.0
#define BXP2BEG          0.0
#define BXP2END          0.0
#define BXP3BEG          0.0
#define BXP3END          0.0

//Velocity along x
#define BXU1BEG          0.0
#define BXU1END          0.0
#define BXU2BEG          0.0
#define BXU2END          1.0
#define BXU3BEG          0.0
#define BXU3END          0.0

//Velocity along y
#define BXV1BEG          0.0
#define BXV1END          0.0
#define BXV2BEG          0.0
#define BXV2END          0.0
#define BXV3BEG          0.0
#define BXV3END          0.0

//Velocity along z
#define BXW1BEG          0.0
#define BXW1END          0.0
#define BXW2BEG          0.0
#define BXW2END          0.0
#define BXW3BEG          0.0
#define BXW3END          0.0

#define REYNOLDS_NUMBER    400.0

#define TOLLERANCE      1.e-5

#define LOG_OUTPUT      1