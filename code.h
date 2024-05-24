#include "main.h"

#define NO  0
#define YES 1

#define X1  0
#define X2  1
#define X3  2

#define UNI  1
#define LOGP 2
#define LOGM 3

#define FINITE_DIFF 1
#define FINITE_VOL  2

#define EXPLICIT        1
#define IMPLICIT        2
#define SEMIIMPLICIT    3

#define STEADY_HEAT_CONDUCTION      1
#define UNSTEADY_FLUID_FLOW         2

#if PHYSICS == STEADY_HEAT_CONDUCTION
    #define TMP         0
#endif

#if PHYSICS == UNSTEADY_FLUID_FLOW
#define PRS         0
    #if STAGGERED_GRID == YES
        #define V1S         0
        #define V2S         1
        #define V3S         2
    #else
        #define V1          1
        #define V2          2
        #define V3          3
        #define V_1         0
        #define V_2         1
        #define V_3         2
    #endif
#endif

#define X1BEGBOUND      0
#define X1ENDBOUND      1
#define X2BEGBOUND      2   
#define X2ENDBOUND      3 
#define X3BEGBOUND      4 
#define X3ENDBOUND      5   

#define ZERO_GRAD       0
#define GRAD            1
#define CONST           2
#define FLUX            3
#define CONVECTIVE      4
#define REFLECTIVE      5

#define E               0
#define W               1
#define N               2
#define S               3
#define F               4
#define B               5

#define FACE            6

#ifndef BX1BEG_H
    #define BX1BEG_H    0
#endif

#ifndef BX1BEG_T
    #define BX1BEG_T    0
#endif

#ifndef BX1END_H
    #define BX1END_H    0
#endif

#ifndef BX1END_T
    #define BX1END_T    0
#endif

#ifndef BX2BEG_H
    #define BX2BEG_H    0
#endif

#ifndef BX2BEG_T
    #define BX2BEG_T    0
#endif

#ifndef BX2END_H
    #define BX2END_H    0
#endif

#ifndef BX2END_T
    #define BX2END_T    0
#endif

#ifndef BX3BEG_H
    #define BX3BEG_H    0
#endif

#ifndef BX3BEG_T
    #define BX3BEG_T    0
#endif

#ifndef BX3END_H
    #define BX3END_H    0
#endif

#ifndef BX3END_T
    #define BX3END_T    0
#endif

#ifndef LOG_OUTPUT
    #define LOG_OUTPUT = 100
#endif

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) > (y)) ? (y) : (x))

extern void initialize();
extern void initialize_boundary();
extern void update_pressure();
extern void update_fractional_velocity();
extern void update_next_velocity();
extern void allocate_boundary();
extern double **** allocate_array_uniform_4D();
extern double **** allocate_array_semi_uniform_4D();
extern double *** allocate_array_3D();
extern double ** allocate_array_2D();
extern double * allocate_array_1D();
extern void generate_uniform_array();
void generate_logplus_array();
void generate_logminus_array();

#if PHYSICS == STEADY_HEAT_CONDUCTION
    extern void calculate_temperature();
    extern void converge_solution();
#elif PHYSICS == UNSTEADY_FLUID_FLOW
    extern void calculate_flux();
    extern void find_solution();
#endif
extern void write_output();
