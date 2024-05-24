#include <stdio.h>
#include <time.h>

#include "main.h"
#include "code.h"
#include "struct.h"
#include "globals.h"

void main(int argc, char *argv[])
{   
    clock_t begin = clock();

    Grid grid[4];
    Data data[4];
    Boundary bound[3];
    tData t_data[4];
    midT mid_time[4];

    IBEG = 1;
    IEND = X1POINTS;
    JBEG = 1;
    JEND = X2POINTS;
    KBEG = 1;
    KEND = X3POINTS;
    g_tStart = 0.0;
    g_Re = REYNOLDS_NUMBER;
    g_step = 0;


    initialize(1, "0", grid, data, t_data, mid_time);
    allocate_boundary(bound, data, grid, t_data, mid_time);
    initialize_boundary(bound, data, grid, t_data, mid_time);
    #if PHYSICS == STEADY_HEAT_CONDUCTION
    converge_solution(data, grid, bound);
    #elif PHYSICS == UNSTEADY_FLUID_FLOW
    find_solution(data, grid, bound, t_data, mid_time);
    #endif

    write_output(data, grid, bound, t_data);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time: %g s", time_spent);

    printf("\nDONE!");
    
}