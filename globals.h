#include "main.h"
/*Contains all the global variables that will be used during runtime*/

/*Global variables for time stepping*/
int g_timeStepping; /*Type of time stepping algorithm used*/
int g_step;
double g_dt; /*Current integration time step*/
double g_tStart; /*Interation start time*/
double g_tStop; /*Integration stop time*/
#ifdef REYNOLDS_NUMBER
double g_Re;
#endif

/*Grid information*/

/*Solver information*/
int g_solver; /*Type of solver -> Finite difference/finite volumes*/

int IBEG;
int IEND;
int JBEG;
int JEND;
int KBEG;
int KEND;