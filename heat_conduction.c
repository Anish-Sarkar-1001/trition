#include <math.h>

#include "main.h"
#include "code.h"
#include "struct.h"
#include "globals.h"

#if PHYSICS == STEADY_HEAT_CONDUCTION
void calculate_temperature(Data *, Grid *);

void calculate_temperature(Data* d, Grid* g) {
    double ir;
    double il;
    double jr;
    double jl;
    double kr;
    double kl;
    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                ir = g->dPoints[X2][j]*g->dPoints[X3][k]/g->dCentroid[X1][i];
                il = g->dPoints[X2][j]*g->dPoints[X3][k]/g->dCentroid[X1][i-1];
                jr = g->dPoints[X1][i]*g->dPoints[X3][k]/g->dCentroid[X2][j];
                jl = g->dPoints[X1][i]*g->dPoints[X3][k]/g->dCentroid[X2][j-1];
                kr = g->dPoints[X1][i]*g->dPoints[X2][j]/g->dCentroid[X3][k];
                kl = g->dPoints[X1][i]*g->dPoints[X2][j]/g->dCentroid[X3][k-1];
                d->vCentroid[TMP][k][j][i] = (d->vCentroid[TMP][k][j][i+1]*ir+
                                            d->vCentroid[TMP][k][j][i-1]*il+
                                            d->vCentroid[TMP][k][j+1][i]*jr+
                                            d->vCentroid[TMP][k][j-1][i]*jl+
                                            d->vCentroid[TMP][k+1][j][i]*kr+
                                            d->vCentroid[TMP][k-1][j][i]*kl)/
                                            (ir+il+jr+jl+kr+kl);
                
            }
        }
    }
}
#endif