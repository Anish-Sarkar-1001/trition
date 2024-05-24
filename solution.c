#include <stdio.h>
#include <math.h>

#include "code.h"
#include "struct.h"
#include "globals.h"

#if PHYSICS == STEADY_HEAT_CONDUCTION
void converge_solution(Data* , Grid*, Boundary*);

void converge_solution(Data* data, Grid* grid, Boundary* bound) {

    double initial_tmp, final_tmp;
    initial_tmp = data->vCentroid[TMP][(int)(KEND+KBEG)/2][(int)JBEG][(int)IBEG];
    final_tmp = initial_tmp+1.;

    while ((final_tmp-initial_tmp)>TOLLERANCE) { /**/
        initialize_boundary(bound, data, grid);
        initial_tmp = data->vCentroid[TMP][(int)(KEND+KBEG)/2][(int)JBEG][(int)IEND];

        calculate_temperature(data, grid);
        final_tmp = data->vCentroid[TMP][(int)(KEND+KBEG)/2][(int)JBEG][(int)IEND];
    }
}

#elif PHYSICS == UNSTEADY_FLUID_FLOW
void find_solution(Data*, Grid*, Boundary*, tData*, midT*);

void find_solution(Data* data, Grid* grid, Boundary* bound, tData* tdata, midT* midt) {
    double time = 0.0;
    double t_x, t_y, t_z;
    double minDt;
    g_dt = DT;

    printf("::::::::::::::::::::::::::::::::::::::::::\n");
    printf("::::::::::                      ::::::::::\n");
    printf(":::::::::: INTEGRATION STARTING ::::::::::\n");
    printf("::::::::::                      ::::::::::\n");
    printf("::::::::::::::::::::::::::::::::::::::::::\n\n");

    initialize_boundary(bound, data, grid, tdata, midt);
    while (time <= g_tStop) {
        g_step+=1;
        if (g_step % LOG_OUTPUT == 0) {
            printf("STEP: %d\n", g_step);
            printf("Time: %0.6f\nIntegration: %0.2f %%\n\n", time, time*100/TSTOP);
        }
        
        calculate_flux(data, grid, tdata, midt, bound);
        minDt = 100.0;
        time+=g_dt;
        for (int k=KBEG; k<=KEND; k++) {
            for (int j=JBEG; j<=JEND; j++) {
                for (int i=IBEG; i<=IEND; i++) {
                    #if STAGGERED_GRID == YES
                    t_x = fabs(tdata->vStaggered[V1S][k][j][i] == 0 ? g_dt : grid->dPoints[X1][i]/tdata->vStaggered[V1S][k][j][i]);
                    t_y = fabs(tdata->vStaggered[V2S][k][j][i] == 0 ? g_dt : grid->dPoints[X2][j]/tdata->vStaggered[V2S][k][j][i]);
                    t_z = fabs(tdata->vStaggered[V3S][k][j][i] == 0 ? g_dt : grid->dPoints[X3][k]/tdata->vStaggered[V3S][k][j][i]);
                    #else
                    t_x = fabs(tdata->vCentroid[V1][k][j][i] == 0 ? g_dt : grid->dPoints[X1][i]/tdata->vCentroid[V1][k][j][i]);
                    t_y = fabs(tdata->vCentroid[V2][k][j][i] == 0 ? g_dt : grid->dPoints[X2][j]/tdata->vCentroid[V2][k][j][i]);
                    t_z = fabs(tdata->vCentroid[V3][k][j][i] == 0 ? g_dt : grid->dPoints[X3][k]/tdata->vCentroid[V3][k][j][i]);
                    #endif
                    minDt = MIN(MIN(MIN(t_x, t_y), t_z), minDt);
                }
            }
        }
        //printf("dt: %g\n", minDt);
        if (minDt != g_dt) {
            minDt*=0.1;
        }
        
        /*for (int k=KBEG; k<=KEND; k++) {
            for (int j=JBEG-1; j<=JEND+1; j++) {
                for (int i=IBEG-1; i<=IEND+1; i++) {
                    printf("%g,", midt->vC[V_2][k][j][i]);
                } printf("\n");
            } printf("\n\n");
        }*/
        g_dt = minDt;
        for (int k=KBEG-1; k<=KEND+1; k++) {
            for (int j=JBEG-1; j<=JEND+1; j++) {
                for (int i=IBEG-1; i<=IEND+1; i++) {
                    #if STAGGERED_GRID == YES
                    data->vStaggered[V1S][k][j][i] = tdata->vStaggered[V1S][k][j][i];
                    data->vStaggered[V2S][k][j][i] = tdata->vStaggered[V2S][k][j][i];
                    data->vStaggered[V3S][k][j][i] = tdata->vStaggered[V3S][k][j][i];
                    #else
                    data->vCentroid[V1][k][j][i] = tdata->vCentroid[V1][k][j][i];
                    data->vCentroid[V2][k][j][i] = tdata->vCentroid[V2][k][j][i];
                    data->vCentroid[V3][k][j][i] = tdata->vCentroid[V3][k][j][i];
                    #endif
                }
            }
        }
    }
}
#endif
