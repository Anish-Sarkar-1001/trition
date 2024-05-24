#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "globals.h"
#include "struct.h"
#include "code.h"
#include "main.h"


void configure_definitions(Grid *);
void create_grid(Grid *);
void cell_centroid(Grid *);
void cell_length(Grid *);
void centroid_distance(Grid *);
void cell_interface(Grid *);
void staggered_area(Grid *);
void compute_area(Grid *);
void initialize_computational_domain(Data *, tData *, midT *);


void initialize(int argc, char *argv[], Grid *grid, Data *data, tData *t_data, midT *mt) {
    configure_definitions(grid);
    create_grid(grid);

    #if SOLVER == FINITE_VOL
    cell_centroid(grid);
    cell_length(grid);
    centroid_distance(grid);
    cell_interface(grid);
    #if STAGGERED_GRID == YES
    staggered_area(grid);
    #else
    compute_area(grid);
    #endif
    initialize_computational_domain(data, t_data, mt);
    #endif
}

void configure_definitions(Grid *g) {
    printf("Initializing:\n");

    g->xBeg[X1] = X1BEG;
    g->xEnd[X1] = X1END;
    g->npTotal[X1] = X1POINTS;
    g->xSpacing[X1] = X1SCALING;
    g->xBeg[X2] = X2BEG;
    g->xEnd[X2] = X2END;
    g->npTotal[X2] = X2POINTS;
    g->xSpacing[X2] = X2SCALING;
    g->xBeg[X3] = X3BEG;
    g->xEnd[X3] = X3END;
    g->npTotal[X3] = X3POINTS;
    g->xSpacing[X3] = X3SCALING;
    g_solver = SOLVER;
    g_timeStepping = TIMESTEPPING;
    g_dt = DT;
    g_tStop = TSTOP;


    printf("\n*** X1 Grid Information ***\n\n");
    printf("X1-Beginning:       %g\n", g->xBeg[X1]);
    printf("X1-End:             %g\n", g->xEnd[X1]);
    printf("X1-Number-of-cells: %d\n", g->npTotal[X1]);
    printf("X1-Grid-spacing:    %d\n", g->xSpacing[X1]);

    printf("\n*** X2 Grid Information ***\n\n");
    printf("X2-Beginning:       %g\n", g->xBeg[X2]);
    printf("X2-End:             %g\n", g->xEnd[X2]);
    printf("X2-Number-of-cells: %d\n", g->npTotal[X2]);
    printf("X2-Grid-spacing:    %d\n", g->xSpacing[X2]);

    printf("\n*** X3 Grid Information ***\n\n");
    printf("X3-Beginning:       %g\n", g->xBeg[X3]);
    printf("X3-End:             %g\n", g->xEnd[X3]);
    printf("X3-Number-of-cells: %d\n", g->npTotal[X3]);
    printf("X3-Grid-spacing:    %d\n", g->xSpacing[X3]);

    printf("\n*** Solver information ***\n\n");
    printf("Solver type:        %d\n", g_solver);

    printf("\n*** Time Information *** \n\n");
    printf("Time stepping type: %d\n", g_timeStepping);
    printf("Time step:          %g\n", g_dt);
    printf("Stop time:          %g\n", g_tStop);
    
}

void create_grid(Grid *g) {
    int np[3];
    for (int i=0; i<3; i++) {
        np[i] = g->npTotal[i]+1;
    }
    g->xPoints = allocate_array_2D(g->xPoints, 3, np);

    if (g->xPoints[X1]!=NULL) {
        printf("\nX1 grid memory allocation successful\n");
    }
    else {
        printf("\nX1 grid memory allocation failed\n");
    }

    if (g->xPoints[X2]!=NULL) {
        printf("\nX2 grid memory allocation successful\n");
    }
    else {
        printf("\nX2 grid memory allocation failed\n");
    }

    if (g->xPoints[X3]!=NULL) {
        printf("\nX3 grid memory allocation successful\n");
    }
    else {
        printf("\nX3 grid memory allocation failed\n");
    }

    if (g->xSpacing[X1] == 1) {
        generate_uniform_array(g, X1);
    }
    if (g->xSpacing[X1] == 2) {
        generate_logplus_array(g, X1);
    }
    if (g->xSpacing[X1] == 3) {
        generate_logminus_array(g, X1);
    }
    printf("X1:\t[");
    for (int index=0; index<=g->npTotal[X1]; index++) {
        printf("%g,", g->xPoints[X1][index]);
    }
    printf("]\n");

    if (g->xSpacing[X2] == 1) {
        generate_uniform_array(g, X2);
    }
    if (g->xSpacing[X2] == 2) {
        generate_logplus_array(g, X2);
    }
    if (g->xSpacing[X2] == 3) {
        generate_logminus_array(g, X2);
    }
    printf("X2:\t[");
    for (int index=0; index<=g->npTotal[X2]; index++) {
        printf("%g,", g->xPoints[X2][index]);
    }
    printf("]\n");

    if (g->xSpacing[X3] == 1) {
        generate_uniform_array(g, X3);
    }
    if (g->xSpacing[X3] == 2) {
        generate_logplus_array(g, X3);
    }
    if (g->xSpacing[X3] == 3) {
        generate_logminus_array(g, X3);
    }
    printf("X3:\t[");
    for (int index=0; index<=g->npTotal[X3]; index++) {
        printf("%g,", g->xPoints[X3][index]);
    }
    printf("]\n");
}

#if SOLVER == FINITE_VOL
void cell_centroid(Grid* g) {
    printf("Defining centroids\n");
    for (int pt=0; pt<3; pt++) {
        g->xCentroids[pt] = allocate_array_1D(g->xCentroids[pt], g->npTotal[pt]);
        for (int i = 0; i<g->npTotal[pt]; i++) {
            g->xCentroids[pt][i] = (g->xPoints[pt][i]+g->xPoints[pt][i+1])/2.; 
        }
        if (g->xCentroids != NULL) {
            printf("X%d Centroid array allocation successful\n", pt+1);
            printf("X%d cell centroids:\t[", pt+1);
            for (int i=0; i<g->npTotal[pt]; i++) {
                printf("%g,", g->xCentroids[pt][i]);
            } printf("]\n");
        }
    }
}

void cell_length(Grid* g) {
    for (int i=0; i<=2; i++) {
        g->dPoints[i] = allocate_array_1D(g->dPoints[i], g->npTotal[i]+2);
        for (int j=1; j<=g->npTotal[i]; j++) {
            g->dPoints[i][j] = g->xPoints[i][j] - g->xPoints[i][j-1];
        }
        g->dPoints[i][0] = g->dPoints[i][1];
        g->dPoints[i][g->npTotal[i]+1] = g->dPoints[i][g->npTotal[i]];
    }
    if (g->dPoints != NULL) {
        printf("Cell lengths calculated\n");
        for (int i=0; i<=2; i++) {
            printf("X%d cell lengths are:\t[", i+1);
            for (int j=0; j<=g->npTotal[i]+1; j++) {
                printf("%g,", g->dPoints[i][j]);
            } printf("]\n");
        }
    }
}

void centroid_distance(Grid* g) {
    for (int i=0; i<=2; i++) {
        g->dCentroid[i] = allocate_array_1D(g->dCentroid[i], g->npTotal[i]+1);
        if (g->npTotal[i] != 1) {
            for (int j=1; j<g->npTotal[i]; j++) {
                g->dCentroid[i][j] = g->xCentroids[i][j] - g->xCentroids[i][j-1];
            }
        g->dCentroid[i][0] = g->dPoints[i][1];
        g->dCentroid[i][g->npTotal[i]] = g->dPoints[i][g->npTotal[i]];
        }
        else {
            g->dCentroid[i][0] = g->dPoints[i][1];
            g->dCentroid[i][1] = g->dPoints[i][1];
        }
    }
    if (g->dCentroid != NULL) {
        printf("Centroid distances calculated\n");
        for (int i=0; i<=2; i++) {
            printf("X%d centroid distances are:\t[", i+1);
            for (int j=0; j<=g->npTotal[i]; j++) {
                printf("%g,", g->dCentroid[i][j]);
            } printf("]\n");
        }
    }
}

void cell_interface(Grid* g) {
    #if DIMENSIONS == 1
    double area;
    area = 1.0;
    for (int i=0; i<=2; i++) {
        g->xArea[i] = allocate_array_1D(g->xArea[i], g->npTotal[i]+1);
        for (int j=0; j<=g->npTotal[i]; j++) {
            g->xArea[i][j] = area;
        }
    }
    if (g->xArea != NULL) {
        printf("Cell interfaces processed\n");
        for (int i=0; i<=2; i++) {
            printf("Cell interface in X%d:\t[", i+1);
            for (int j=0; j<=g->npTotal[i]; j++) {
                printf("%g,", g->xArea[i][j]);
            } printf("]\n");
        }
    }
    #endif

    #if DIMENSIONS == 2
    double area;
    for (int i=0; i<=2; i++) {
        g->xArea[i] = allocate_array_1D(g->xArea[i], g->npTotal[i]+1);
        for (int j=0; j<=2; j++) {
            if (j != i && g->npTotal[j]!=1) {
                for (int k=0; k<=g->npTotal[i]; k++);
            }
        }
        for (int j=0; j<=g->npTotal[i]; j++) {
            g->xArea[i][j] = area;
        }
    }
    #endif

    #if DIMENSIONS == 3
    #endif
}
#endif

#if STAGGERED_GRID == YES
void initialize_staggered_grid(Grid *g, StagGrid *sg) {
    sg->xStaggeredPoints = allocate_array_2D(g->xPoints, 3, g->npTotal);
    for (int i=0; i<3; i++) {
        for (int j=0; i<g->npTotal[i]; j++) {
            sg->xStaggeredPoints[i][j] = g->xPoints[i][j+1];
        }
    }
}
void staggered_area(Grid *g) {
    g->v1StagArea = allocate_array_uniform_4D(g->v1StagArea, X3POINTS+2, X2POINTS+2, X1POINTS+2, FACE);
    g->v2StagArea = allocate_array_uniform_4D(g->v2StagArea, X3POINTS+2, X2POINTS+2, X1POINTS+2, FACE);
    g->v3StagArea = allocate_array_uniform_4D(g->v3StagArea, X3POINTS+2, X2POINTS+2, X1POINTS+2, FACE);
    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                for (int x=0; x<FACE; x++) {
                    g->v1StagArea[k][j][i][E] = g->dPoints[X2][j]*g->dPoints[X3][k];
                    g->v1StagArea[k][j][i][W] = g->dPoints[X2][j]*g->dPoints[X3][k];
                    g->v1StagArea[k][j][i][N] = g->dPoints[X3][k]*(g->dPoints[X1][i]+g->dPoints[X1][i+1])/2.;
                    g->v1StagArea[k][j][i][S] = g->dPoints[X3][k]*(g->dPoints[X1][i]+g->dPoints[X1][i+1])/2.;
                    g->v1StagArea[k][j][i][F] = g->dPoints[X2][j]*(g->dPoints[X1][i]+g->dPoints[X1][i+1])/2.;
                    g->v1StagArea[k][j][i][B] = g->dPoints[X2][j]*(g->dPoints[X1][i]+g->dPoints[X1][i+1])/2.;

                    g->v2StagArea[k][j][i][E] = (g->dPoints[X2][j] + g->dPoints[X2][j+1])*g->dPoints[X3][k]/2.;
                    g->v2StagArea[k][j][i][W] = (g->dPoints[X2][j] + g->dPoints[X2][j+1])*g->dPoints[X3][k]/2.;
                    g->v2StagArea[k][j][i][N] = g->dPoints[X3][k]*g->dPoints[X1][i];
                    g->v2StagArea[k][j][i][S] = g->dPoints[X3][k]*g->dPoints[X1][i];
                    g->v2StagArea[k][j][i][F] = (g->dPoints[X2][j] + g->dPoints[X2][j+1])*g->dPoints[X1][i]/2.;
                    g->v2StagArea[k][j][i][B] = (g->dPoints[X2][j] + g->dPoints[X2][j+1])*g->dPoints[X1][i]/2.;

                    g->v3StagArea[k][j][i][E] = (g->dPoints[X3][k] + g->dPoints[X3][k+1])*g->dPoints[X2][j]/2.;
                    g->v3StagArea[k][j][i][W] = (g->dPoints[X3][k] + g->dPoints[X3][k+1])*g->dPoints[X2][j]/2.;
                    g->v3StagArea[k][j][i][N] = (g->dPoints[X3][k] + g->dPoints[X3][k+1])*g->dPoints[X1][i]/2.;
                    g->v3StagArea[k][j][i][S] = (g->dPoints[X3][k] + g->dPoints[X3][k+1])*g->dPoints[X1][i]/2.;
                    g->v3StagArea[k][j][i][F] = g->dPoints[X1][i]*g->dPoints[X2][j];
                    g->v3StagArea[k][j][i][B] = g->dPoints[X1][i]*g->dPoints[X2][j];
                }
            }
        }
    }
}
#else
void compute_area(Grid *g) {
    g->Area = allocate_array_uniform_4D(g->Area, X3POINTS+2, X2POINTS+2, X1POINTS+2, FACE);
    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                for (int x=0; x<FACE; x++) {
                    g->Area[k][j][i][E] = g->dPoints[X2][j]*g->dPoints[X3][k];
                    g->Area[k][j][i][W] = g->dPoints[X2][j]*g->dPoints[X3][k];
                    g->Area[k][j][i][N] = g->dPoints[X3][k]*g->dPoints[X1][i];
                    g->Area[k][j][i][S] = g->dPoints[X3][k]*g->dPoints[X1][i];
                    g->Area[k][j][i][F] = g->dPoints[X2][j]*g->dPoints[X1][i];
                    g->Area[k][j][i][B] = g->dPoints[X2][j]*g->dPoints[X1][i];
                }
            }
        }
    }
}
#endif

void initialize_computational_domain(Data* d, tData* td, midT* mt) {
    #if PHYSICS == STEADY_HEAT_CONDUCTION
    d->vCentroid = allocate_array_uniform_4D(d->vCentroid, 1, X3POINTS+2, X2POINTS+2, X1POINTS+2);
    printf("Computational domain memmory allocation successful\n");
    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            for (int i=0; i<=X1POINTS+1; i++) {
                d->vCentroid[TMP][k][j][i] = 0.0;
            }
        }
    }
    #elif PHYSICS == UNSTEADY_FLUID_FLOW
    #if STAGGERED_GRID == YES
    d->vCentroid = allocate_array_uniform_4D(d->vCentroid, 1, X3POINTS+2, X2POINTS+2, X1POINTS+2);
    printf("Computational domain memmory allocation successful\n");
    /*
    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            for (int i=0; i<=X1POINTS+1; i++) {
                d->vCentroid[PRS][k][j][i] = 0.0;
            }
        }
    }*/
    d->vStaggered = allocate_array_uniform_4D(d->vStaggered, 3, X3POINTS+2, X2POINTS+2, X1POINTS+2);

    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            for (int i=0; i<=X1POINTS+1; i++) {
                d->vStaggered[V1S][k][j][i] = 0.0;
                d->vStaggered[V2S][k][j][i] = 0.0;
                d->vStaggered[V3S][k][j][i] = 0.0;
            }
        }
    }

    mt->vS = allocate_array_uniform_4D(mt->vS, 3, X3POINTS+2, X2POINTS+2, X1POINTS+2);

    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            for (int i=0; i<=X1POINTS+1; i++) {
                mt->vS[V1S][k][j][i] = 0.0;
                mt->vS[V2S][k][j][i] = 0.0;
                mt->vS[V3S][k][j][i] = 0.0;
            }
        }
    }

    #else
    d->vCentroid = allocate_array_uniform_4D(d->vCentroid, 4, X3POINTS+2, X2POINTS+2, X1POINTS+2);
    d->advecFlux = allocate_array_uniform_4D(d->advecFlux, 3, X3POINTS+2, X2POINTS+2, X1POINTS+2);
    printf("Computational domain memmory allocation successful\n");
    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            for (int i=0; i<=X1POINTS+1; i++) {
                d->vCentroid[PRS][k][j][i] = 0.0;
                d->vCentroid[V1][k][j][i] = 0.0;
                d->vCentroid[V2][k][j][i] = 0.0;
                d->vCentroid[V3][k][j][i] = 0.0;
            }
        }
    }
    mt->vC = allocate_array_uniform_4D(mt->vC, 3, X3POINTS+2, X2POINTS+2, X1POINTS+2);

    #endif

    #if STAGGERED_GRID == YES
    td->vCentroid = allocate_array_uniform_4D(td->vCentroid, 1, X3POINTS+2, X2POINTS+2, X1POINTS+2);
    printf("Computational domain memmory allocation successful\n");
    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            for (int i=0; i<=X1POINTS+1; i++) {
                td->vCentroid[PRS][k][j][i] = 0.0;
            }
        }
    }
    td->vStaggered = allocate_array_uniform_4D(td->vStaggered, 3, X3POINTS+2, X2POINTS+2, X1POINTS+2);

    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            for (int i=0; i<=X1POINTS+1; i++) {
                td->vStaggered[V1S][k][j][i] = 0.0;
                td->vStaggered[V2S][k][j][i] = 0.0;
                td->vStaggered[V3S][k][j][i] = 0.0;
            }
        }
    }

    #else
    td->vCentroid = allocate_array_uniform_4D(td->vCentroid, 4, X3POINTS+2, X2POINTS+2, X1POINTS+2);
    printf("Computational domain memmory allocation successful\n");
    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            for (int i=0; i<=X1POINTS+1; i++) {
                td->vCentroid[PRS][k][j][i] = 0.0;
                td->vCentroid[V1][k][j][i] = 0.0;
                td->vCentroid[V2][k][j][i] = 0.0;
                td->vCentroid[V3][k][j][i] = 0.0;
            }
        }
    }
    #endif

    if (d->vCentroid != NULL) {
        printf("Data array allocation successful\n\n");
    }

    /*printf("Grid:\n");
    for (int x=0; x<1; x++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            for (int j=0; j<=X2POINTS+1; j++) {
                for (int i=0; i<=X1POINTS+1; i++) {
                    printf("%g,", d->vCentroid[x][k][j][i]);
                } printf("\n");
            } printf("\n\n");
        }
    }
    for (int x=0; x<3; x++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            for (int j=0; j<=X2POINTS+1; j++) {
                for (int i=0; i<=X1POINTS+1; i++) {
                    printf("%g,", d->vStaggered[x][k][j][i]);
                } printf("\n");
            } printf("\n\n");
        }
    }
    for (int x=0; x<1; x++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            for (int j=0; j<=X2POINTS+1; j++) {
                for (int i=0; i<=X1POINTS+1; i++) {
                    printf("%g,", td->vCentroid[x][k][j][i]);
                } printf("\n");
            } printf("\n\n");
        }
    }
    for (int x=0; x<3; x++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            for (int j=0; j<=X2POINTS+1; j++) {
                for (int i=0; i<=X1POINTS+1; i++) {
                    printf("%g,", td->vStaggered[x][k][j][i]);
                } printf("\n");
            } printf("\n\n");
        }
    }*/
    #endif
}
