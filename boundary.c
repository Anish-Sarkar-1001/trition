#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "globals.h"
#include "struct.h"
#include "code.h"
#include "main.h"

void initialize_boundary(Boundary *, Data *, Grid *, tData *, midT *);
void allocate_boundary(Boundary *, Data *, Grid *, tData *, midT *);
void update_pressure(Boundary *, Data *, Grid *, tData *, midT *);
void update_fractional_velocity(Boundary *, Data *, Grid *, tData *, midT *);
void update_next_velocity(Boundary *, Data *, Grid *, tData *, midT *);

void allocate_boundary(Boundary* b, Data* d, Grid* g, tData* td, midT* mt) {
    int boundPoints1[6] = {X2POINTS+2, X2POINTS+2, X3POINTS+2, X3POINTS+2, X1POINTS+2, X1POINTS+2};
    int boundPoints2[6] = {X3POINTS+2, X3POINTS+2, X1POINTS+2, X1POINTS+2, X2POINTS+2, X2POINTS+2};

    #if PHYSICS == STEADY_HEAT_CONDUCTION
    b->bCentroids = allocate_array_semi_uniform_4D(b->bCentroids, 1, 6, boundPoints2, boundPoints1);
    #elif PHYSICS == UNSTEADY_FLUID_FLOW
        #if STAGGERED_GRID == YES
        b->bCentroids = allocate_array_semi_uniform_4D(b->bCentroids, 1, 6, boundPoints2, boundPoints1);
        b->b0Stag = allocate_array_semi_uniform_4D(b->b0Stag, 3, 6, boundPoints2, boundPoints1);
        b->b1_2Stag = allocate_array_semi_uniform_4D(b->b1_2Stag, 3, 6, boundPoints2, boundPoints1);
        b->b1Stag = allocate_array_semi_uniform_4D(b->b1Stag, 3, 6, boundPoints2, boundPoints1);
        #else
        b->bCentroids = allocate_array_semi_uniform_4D(b->bCentroids, 1, 6, boundPoints2, boundPoints1);
        b->b0Stag = allocate_array_semi_uniform_4D(b->b0Stag, 3, 6, boundPoints2, boundPoints1);
        b->b1_2Stag = allocate_array_semi_uniform_4D(b->b1_2Stag, 3, 6, boundPoints2, boundPoints1);
        b->b1Stag = allocate_array_semi_uniform_4D(b->b1Stag, 3, 6, boundPoints2, boundPoints1);
        #endif
    #endif
}

void initialize_boundary(Boundary* b, Data* d, Grid* g, tData* td, midT* mt) {
    int boundPoints1[6] = {X2POINTS+2, X2POINTS+2, X3POINTS+2, X3POINTS+2, X1POINTS+2, X1POINTS+2};
    int boundPoints2[6] = {X3POINTS+2, X3POINTS+2, X1POINTS+2, X1POINTS+2, X2POINTS+2, X2POINTS+2};

    #if PHYSICS == STEADY_HEAT_CONDUCTION

    #if X1BEGBOUNDARY == ZERO_GRAD
    for (int k=0; k<X3POINTS; k++) {
        for (int j=0; j<X2POINTS; j++) {
            double T_b = d->vCentroid[TMP][k][j][IBEG];
            b->bCentroids[TMP][X1BEGBOUND][k][j] = 2.*T_b-d->vCentroid[TMP][k][j][IBEG];
        }
    }
    #endif

    #if X1BEGBOUNDARY == GRAD
    for (int k=0; k<X3POINTS; k++) {
        for (int j=0; j<X2POINTS; j++) {
            b->bCentroids[TMP][X1BEGBOUND][k][j] = d->vCentroid[TMP][k][j][IBEG]+(double)BX1BEG*g->dPoints[X1][IBEG];
        }
    }
    #endif

    #if X1BEGBOUNDARY == CONST
    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            double T_b = (double)BX1BEG;
            b->bCentroids[TMP][X1BEGBOUND][k][j] = 2.*T_b-d->vCentroid[TMP][k][j][IBEG];
            
        }
    }
    #endif
    
    #if X1BEGBOUNDARY == FLUX
    for (int k=0; k<X3POINTS; k++) {
        for (int j=0; j<X2POINTS; j++) {
            b->bCentroids[TMP][X1BEGBOUND][k][j] = d->vCentroid[TMP][k][j][IBEG]-
                                                    (double)X1BEG*g->dPoints[X1][IBEG]/CONDUCTIVITY;
        }
    }
    #endif

    #if X1BEGBOUNDARY == CONVECTIVE
    for (int k=0; k<X3POINTS; k++) {
        for (int j=0; j<X2POINTS; j++) {
            b->bCentroids[TMP][X1BEGBOUND][k][j] = ((2.*CONDUCTIVITY*d->vCentroid[TMP][k][j][IBEG]
            /g->dPoints[X1][IBEG])+BX1BEG_H*BX1BEG_T)/(BX1BEG_H+(2.*CONDUCTIVITY/g->dPoints[X1][IBEG]));
        }
    }
    #endif

    #if X1ENDBOUNDARY == ZERO_GRAD
    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            double T_b = d->vCentroid[TMP][k][j][IEND];
            b->bCentroids[TMP][X1ENDBOUND][k][j] = T_b;
        }
    }
    #endif

    #if X1ENDBOUNDARY == GRAD
    for (int k=0; k<X3POINTS; k++) {
        for (int j=0; j<X2POINTS; j++) {
            b->bCentroids[TMP][X1ENDBOUND][k][j] = d->vCentroid[TMP][k][j][IEND]+(double)BX1END;
        }
    }
    #endif

    #if X1ENDBOUNDARY == CONST
    for (int k=0; k<X3POINTS; k++) {
        for (int j=0; j<X2POINTS; j++) {
            b->bCentroids[TMP][X1ENDBOUND][k][j] = (double)BX1END;
        }
    }
    #endif
    
    #if X1ENDBOUNDARY == FLUX
    for (int k=0; k<X3POINTS; k++) {
        for (int j=0; j<X2POINTS; j++) {
            b->bCentroids[TMP][X1ENDBOUND][k][j] = d->vCentroid[TMP][k][j][IEND]-
                                                    (double)X1END*g->dPoints[X1][IEND]/CONDUCTIVITY;
        }
    }
    #endif

    #if X1ENDBOUNDARY == CONVECTIVE
    for (int k=0; k<X3POINTS; k++) {
        for (int j=0; j<X2POINTS; j++) {
            b->bCentroids[TMP][X1ENDBOUND][k][j] = ((2.*CONDUCTIVITY*d->vCentroid[TMP][k][j][IEND]
            /g->dPoints[X1][IEND])+BX1END_H*BX1END_T)/(BX1END_H+(2.*CONDUCTIVITY/g->dPoints[X1][IEND]));
        }
    }
    #endif

    #if X2BEGBOUNDARY == ZERO_GRAD
    for (int i=0; i<X1POINTS; i++) {
        for (int k=0; k<X3POINTS; k++) {
            double T_b = d->vCentroid[TMP][k][JBEG][i];
            b->bCentroids[TMP][X2BEGBOUND][i][k] = 2.*T_b-d->vCentroid[TMP][k][JBEG][i];
        }
    }
    #endif

    #if X2BEGBOUNDARY == GRAD
    for (int i=0; i<X1POINTS; i++) {
        for (int k=0; k<X3POINTS; k++) {
            b->bCentroids[TMP][X2BEGBOUND][i][k] = d->vCentroid[TMP][k][JBEG][i]+(double)BX2BEG;
        }
    }
    #endif

    #if X2BEGBOUNDARY == CONST
    for (int i=0; i<X1POINTS; i++) {
        for (int k=0; k<X3POINTS; k++) {
            b->bCentroids[TMP][X2BEGBOUND][i][k] = (double)BX2BEG;
        }
    }
    #endif
    
    #if X2BEGBOUNDARY == FLUX
    for (int i=0; i<=X1POINTS+1; i++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            double T_b = d->vCentroid[TMP][k][JBEG][i]+(double)BX2BEG*g->dPoints[X2][JBEG]/(CONDUCTIVITY);
            b->bCentroids[TMP][X2BEGBOUND][i][k] = 2.*T_b-d->vCentroid[TMP][k][JBEG][i];
        }
    }
    #endif

    #if X2BEGBOUNDARY == CONVECTIVE
    for (int i=0; i<=X1POINTS+1; i++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            double T_b = ((2.*CONDUCTIVITY*d->vCentroid[TMP][k][JBEG][i]
            /g->dPoints[X2][JBEG])+BX2BEG_H*BX2BEG_T)/(BX2BEG_H+(2.*CONDUCTIVITY/g->dPoints[X2][JBEG]));
            b->bCentroids[TMP][X2BEGBOUND][i][k] = 2.*T_b-d->vCentroid[TMP][k][JBEG][i];
        }
    }
    #endif

    #if X2ENDBOUNDARY == ZERO_GRAD
    for (int i=0; i<X1POINTS; i++) {
        for (int k=0; k<X3POINTS; k++) {
            double T_b = d->vCentroid[TMP][k][JEND][i];
            b->bCentroids[TMP][X2ENDBOUND][i][k] = 2.*T_b-d->vCentroid[TMP][k][JEND][i];
        }
    }
    #endif

    #if X2ENDBOUNDARY == GRAD
    for (int i=0; i<X1POINTS; i++) {
        for (int k=0; k<X3POINTS; k++) {
            b->bCentroids[TMP][X2ENDBOUND][i][k] = d->vCentroid[TMP][k][JEND][i]+(double)BX2END;
        }
    }
    #endif

    #if X2ENDBOUNDARY == CONST
    for (int i=0; i<X1POINTS; i++) {
        for (int k=0; k<X3POINTS; k++) {
            b->bCentroids[TMP][X2ENDBOUND][i][k] = (double)BX2END;
        }
    }
    #endif
    
    #if X2ENDBOUNDARY == FLUX
    for (int i=0; i<=X1POINTS+1; i++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            double T_b = d->vCentroid[TMP][k][JEND][i]+(double)BX2END*g->dPoints[X2][JEND]/CONDUCTIVITY;
            b->bCentroids[TMP][X2ENDBOUND][i][k] = 2.*T_b-d->vCentroid[TMP][k][JEND][i];
        }
    }
    #endif

    #if X2ENDBOUNDARY == CONVECTIVE
    for (int i=0; i<=X1POINTS+1; i++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            double T_b = ((2*CONDUCTIVITY*d->vCentroid[TMP][k][JEND][i]
            /g->dPoints[X2][JEND])+BX2END_H*BX2END_T)/(BX2END_H+(2*CONDUCTIVITY/g->dPoints[X2][JEND]));
            b->bCentroids[TMP][X2ENDBOUND][i][k] = 2.*T_b-d->vCentroid[TMP][k][JEND][i];
        }
    }
    #endif

    #if X3BEGBOUNDARY == ZERO_GRAD
    
    for (int j=0; j<=X2POINTS+1; j++) {
        for (int i=0; i<=X1POINTS+1; i++) {
            double T_b = d->vCentroid[TMP][KBEG][j][i];
            b->bCentroids[TMP][X3BEGBOUND][j][i] = 2.*T_b-d->vCentroid[TMP][KBEG][j][i];
            
        }
    }
    #endif

    #if X3BEGBOUNDARY == GRAD
    for (int j=0; j<X2POINTS; j++) {
        for (int i=0; i<X1POINTS; i++) {
            b->bCentroids[TMP][X3BEGBOUND][j][i] = d->vCentroid[TMP][KBEG][j][i]+(double)BX3BEG;
        }
    }
    #endif

    #if X3BEGBOUNDARY == CONST
    for (int j=0; j<X2POINTS; j++) {
        for (int i=0; i<X1POINTS; i++) {
            b->bCentroids[TMP][X3BEGBOUND][j][i] = (double)BX3BEG;
        }
    }
    #endif
    
    #if X3BEGBOUNDARY == FLUX
    for (int j=0; j<X2POINTS; j++) {
        for (int i=0; i<X1POINTS; i++) {
            b->bCentroids[TMP][X3BEGBOUND][j][i] = d->vCentroid[TMP][KBEG][j][i]-
                                                    (double)X3BEG*g->dPoints[X1][IBEG]/CONDUCTIVITY;
        }
    }
    #endif

    #if X3BEGBOUNDARY == CONVECTIVE
    for (int j=0; j<X2POINTS; j++) {
        for (int i=0; i<X1POINTS; i++) {
            b->bCentroids[TMP][X3BEGBOUND][j][i] = ((2.*CONDUCTIVITY*d->vCentroid[TMP][KBEG][j][i]
            /g->dPoints[X1][IBEG])+BX3BEG_H*BX3BEG_T)/(BX3BEG_H+(2.*CONDUCTIVITY/g->dPoints[X1][IBEG]));
        }
    }
    #endif

    #if X3ENDBOUNDARY == ZERO_GRAD
    
    for (int j=0; j<=X2POINTS+1; j++) {
        for (int i=0; i<=X1POINTS+1; i++) {
            double T_b = d->vCentroid[TMP][KEND][j][i];
            b->bCentroids[TMP][X3ENDBOUND][j][i] = 2.*T_b-d->vCentroid[TMP][KEND][j][i];
        }
    }
    #endif

    #if X3ENDBOUNDARY == GRAD
    for (int j=0; j<X2POINTS; j++) {
        for (int i=0; i<X1POINTS; i++) {
            b->bCentroids[TMP][X3ENDBOUND][j][i] = d->vCentroid[TMP][KEND][j][i]+(double)BX3END;
        }
    }
    #endif

    #if X3ENDBOUNDARY == CONST
    for (int j=0; j<X2POINTS; j++) {
        for (int i=0; i<X1POINTS; i++) {
            b->bCentroids[TMP][X3ENDBOUND][j][i] = (double)BX3END;
        }
    }
    #endif
    
    #if X3ENDBOUNDARY == FLUX
    for (int j=0; j<X2POINTS; j++) {
        for (int i=0; i<X1POINTS; i++) {
            b->bCentroids[TMP][X3ENDBOUND][j][i] = d->vCentroid[TMP][KEND][j][i]-
                                                    (double)X3END*g->dPoints[X1][IBEG]/CONDUCTIVITY;
        }
    }
    #endif

    #if X3ENDBOUNDARY == CONVECTIVE
    for (int j=0; j<X2POINTS; j++) {
        for (int i=0; i<X1POINTS; i++) {
            b->bCentroids[TMP][X3ENDBOUND][j][i] = ((2.*CONDUCTIVITY*d->vCentroid[TMP][KEND][j][i]
            /g->dPoints[X1][IBEG])+BX3END_H*BX3END_T)/(BX3END_H+(2.*CONDUCTIVITY/g->dPoints[X1][IBEG]));
        }
    }
    #endif
/*
    for (int j=0; j<X2POINTS; j++) {
        for (int i=0; i<X1POINTS; i++) {
            printf("%g", b->bCentroids[TMP][X3ENDBOUND][j][i]);
        }
    }
*/
    for (int j=0; j<=X2POINTS+1; j++) {
        for (int i=0; i<=X1POINTS+1; i++) {
            d->vCentroid[TMP][KBEG-1][j][i] = b->bCentroids[TMP][X3BEGBOUND][j][i];
            d->vCentroid[TMP][KEND+1][j][i] = b->bCentroids[TMP][X3ENDBOUND][j][i];
        }
    }

    for (int i=0; i<=X1POINTS+1; i++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            d->vCentroid[TMP][k][JBEG-1][i] = b->bCentroids[TMP][X2BEGBOUND][i][k];
            d->vCentroid[TMP][k][JEND+1][i] = b->bCentroids[TMP][X2ENDBOUND][i][k];
        }
    }

    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            d->vCentroid[TMP][k][j][IBEG-1] = b->bCentroids[TMP][X1BEGBOUND][k][j];
            d->vCentroid[TMP][k][j][IEND+1] = b->bCentroids[TMP][X1ENDBOUND][k][j];
        }
    }
    /*
    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            for (int i=0; i<=X1POINTS+1; i++) {
                printf("%g,", d->vCentroid[TMP][k][j][i]);
            } printf("\n");
        } printf("\n\n");
    }*/

    #elif PHYSICS == UNSTEADY_FLUID_FLOW

    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            #if X1PBEGBOUNDARY == ZERO_GRAD
            b->bCentroids[PRS][X1BEGBOUND][k][j] = td->vCentroid[PRS][k][j][IBEG];
            #elif X1PBEGBOUNDARY == GRAD
            b->bCentroids[PRS][X1BEGBOUND][k][j] = td->vCentroid[PRS][k][j][IBEG] - 
                (double)BXP1BEG*g->dPoints[X1][IBEG];
            #elif X1PBEGBOUNDARY == CONST
            b->bCentroids[PRS][X1BEGBOUND][k][j] = 2.*(double)BXP1BEG - td->vCentroid[PRS][k][j][IBEG];
            #endif

            #if X1PENDBOUNDARY == ZERO_GRAD
            b->bCentroids[PRS][X1ENDBOUND][k][j] = td->vCentroid[PRS][k][j][IEND];
            #elif X1PENDBOUNDARY == GRAD
            b->bCentroids[PRS][X1ENDBOUND][k][j] = td->vCentroid[PRS][k][j][IEND] + 
                (double)BXP1END*g->dPoints[X1][IEND];
            #elif X1PENDBOUNDARY == CONST
            b->bCentroids[PRS][X1ENDBOUND][k][j] = 2.*(double)BXP1END - td->vCentroid[PRS][k][j][IEND];
            #endif

            #if STAGGERED_GRID == YES
                #if X1UBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V1S][X1BEGBOUND][k][j] = d->vStaggered[V1S][k][j][IBEG];
                b->b1_2Stag[V1S][X1BEGBOUND][k][j] = mt->vS[V1S][k][j][IBEG];
                b->b1Stag[V1S][X1BEGBOUND][k][j] = td->vStaggered[V1S][k][j][IBEG];
                #elif X1UBEGBOUNDARY == GRAD
                b->b0Stag[V1S][X1BEGBOUND][k][j] = d->vStaggered[V1S][k][j][IBEG] - 
                    (double)BXU1BEG*g->dPoints[X1][IBEG];
                b->b1_2Stag[V1S][X1BEGBOUND][k][j] = mt->vS[V1S][k][j][IBEG] - 
                    (double)BXU1BEG*g->dPoints[X1][IBEG];
                b->b1Stag[V1S][X1BEGBOUND][k][j] = td->vStaggered[V1S][k][j][IBEG] - 
                    (double)BXU1BEG*g->dPoints[X1][IBEG];
                #elif X1UBEGBOUNDARY == CONST
                b->b0Stag[V1S][X1BEGBOUND][k][j] = (double)BXU1BEG;
                b->b1_2Stag[V1S][X1BEGBOUND][k][j] = (double)BXU1BEG;
                b->b1Stag[V1S][X1BEGBOUND][k][j] = (double)BXU1BEG;
                #endif

                #if X1VBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V2S][X1BEGBOUND][k][j] = d->vStaggered[V2S][k][j][IBEG];
                b->b1_2Stag[V2S][X1BEGBOUND][k][j] = mt->vS[V2S][k][j][IBEG];
                b->b1Stag[V2S][X1BEGBOUND][k][j] = td->vStaggered[V2S][k][j][IBEG];
                #elif X1VBEGBOUNDARY == GRAD
                b->b0Stag[V2S][X1BEGBOUND][k][j] = d->vStaggered[V2S][k][j][IBEG] - 
                    (double)BXV1BEG*g->dPoints[X1][IBEG];
                b->b1_2Stag[V2S][X1BEGBOUND][k][j] = mt->vS[V2S][k][j][IBEG] - 
                    (double)BXV1BEG*g->dPoints[X1][IBEG];
                b->b1Stag[V2S][X1BEGBOUND][k][j] = td->vStaggered[V2S][k][j][IBEG] - 
                    (double)BXV1BEG*g->dPoints[X1][IBEG];
                #elif X1VBEGBOUNDARY == CONST
                b->b0Stag[V2S][X1BEGBOUND][k][j] = 2.*(double)BXV1BEG - d->vStaggered[V2S][k][j][IBEG];
                b->b1_2Stag[V2S][X1BEGBOUND][k][j] = 2.*(double)BXV1BEG - mt->vS[V2S][k][j][IBEG];
                b->b1Stag[V2S][X1BEGBOUND][k][j] = 2.*(double)BXV1BEG - td->vStaggered[V2S][k][j][IBEG];
                #endif

                #if X1WBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V3S][X1BEGBOUND][k][j] = d->vStaggered[V3S][k][j][IBEG];
                b->b1_2Stag[V3S][X1BEGBOUND][k][j] = mt->vS[V3S][k][j][IBEG];
                b->b1Stag[V3S][X1BEGBOUND][k][j] = td->vStaggered[V3S][k][j][IBEG];
                #elif X1WBEGBOUNDARY == GRAD
                b->b0Stag[V3S][X1BEGBOUND][k][j] = d->vStaggered[V3S][k][j][IBEG] - 
                    (double)BXW1BEG*g->dPoints[X1][IBEG];
                b->b1_2Stag[V3S][X1BEGBOUND][k][j] = mt->vS[V3S][k][j][IBEG] - 
                    (double)BXW1BEG*g->dPoints[X1][IBEG];
                b->b1Stag[V3S][X1BEGBOUND][k][j] = td->vStaggered[V3S][k][j][IBEG] - 
                    (double)BXW1BEG*g->dPoints[X1][IBEG];
                #elif X1WBEGBOUNDARY == CONST
                b->b0Stag[V3S][X1BEGBOUND][k][j] = (double)BXW1BEG;
                b->b1_2Stag[V3S][X1BEGBOUND][k][j] = (double)BXW1BEG;
                b->b1Stag[V3S][X1BEGBOUND][k][j] = (double)BXW1BEG;
                #endif

                #if X1UENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V1S][X1ENDBOUND][k][j] = d->vStaggered[V1S][k][j][IEND];
                b->b1_2Stag[V1S][X1ENDBOUND][k][j] = mt->vS[V1S][k][j][IEND];
                b->b1Stag[V1S][X1ENDBOUND][k][j] = td->vStaggered[V1S][k][j][IEND];
                #elif X1UENDBOUNDARY == GRAD
                b->b0Stag[V1S][X1ENDBOUND][k][j] = d->vStaggered[V1S][k][j][IEND] + 
                    (double)BXU1END*g->dPoints[X1][IEND];
                b->b1_2Stag[V1S][X1ENDBOUND][k][j] = mt->vS[V1S][k][j][IEND] + 
                    (double)BXU1END*g->dPoints[X1][IEND];
                b->b1Stag[V1S][X1ENDBOUND][k][j] = td->vStaggered[V1S][k][j][IEND] + 
                    (double)BXU1END*g->dPoints[X1][IEND];
                #elif X1UENDBOUNDARY == CONST
                b->b0Stag[V1S][X1ENDBOUND][k][j] = 2.*(double)BXU1END - d->vStaggered[V1S][k][j][IEND-1];
                b->b1_2Stag[V1S][X1ENDBOUND][k][j] = 2.*(double)BXU1END - mt->vS[V1S][k][j][IEND-1];
                b->b1Stag[V1S][X1ENDBOUND][k][j] = 2.*(double)BXU1END - td->vStaggered[V1S][k][j][IEND-1];
                #endif

                #if X1VENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V2S][X1ENDBOUND][k][j] = d->vStaggered[V2S][k][j][IEND];
                b->b1_2Stag[V2S][X1ENDBOUND][k][j] = mt->vS[V2S][k][j][IEND];
                b->b1Stag[V2S][X1ENDBOUND][k][j] = td->vStaggered[V2S][k][j][IEND];
                #elif X1VENDBOUNDARY == GRAD
                b->b0Stag[V2S][X1ENDBOUND][k][j] = d->vStaggered[V2S][k][j][IEND] + 
                    (double)BXV1END*g->dPoints[X1][IEND];
                b->b1_2Stag[V2S][X1ENDBOUND][k][j] = mt->vS[V2S][k][j][IEND] + 
                    (double)BXV1END*g->dPoints[X1][IEND];
                b->b1Stag[V2S][X1ENDBOUND][k][j] = td->vStaggered[V2S][k][j][IEND] + 
                    (double)BXV1END*g->dPoints[X1][IEND];
                #elif X1VENDBOUNDARY == CONST
                b->b0Stag[V2S][X1ENDBOUND][k][j] = 2.*(double)BXV1END - d->vStaggered[V2S][k][j][IEND];
                b->b1_2Stag[V2S][X1ENDBOUND][k][j] = 2.*(double)BXV1END - mt->vS[V2S][k][j][IEND];
                b->b1Stag[V2S][X1ENDBOUND][k][j] = 2.*(double)BXV1END - td->vStaggered[V2S][k][j][IEND];
                #endif

                #if X1WENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V3S][X1ENDBOUND][k][j] = d->vStaggered[V3S][k][j][IEND];
                b->b1_2Stag[V3S][X1ENDBOUND][k][j] = mt->vS[V3S][k][j][IEND];
                b->b1Stag[V3S][X1ENDBOUND][k][j] = td->vStaggered[V3S][k][j][IEND];
                #elif X1WENDBOUNDARY == GRAD
                b->b0Stag[V3S][X1ENDBOUND][k][j] = d->vStaggered[V3S][k][j][IEND] + 
                    (double)BXW1END*g->dPoints[X1][IEND];
                b->b1_2Stag[V3S][X1ENDBOUND][k][j] = mt->vS[V3S][k][j][IEND] + 
                    (double)BXW1END*g->dPoints[X1][IEND];
                b->b1Stag[V3S][X1ENDBOUND][k][j] = td->vStaggered[V3S][k][j][IEND] + 
                    (double)BXW1END*g->dPoints[X1][IEND];
                #elif X1WENDBOUNDARY == CONST
                b->b0Stag[V3S][X1ENDBOUND][k][j] = 2.*(double)BXW1END - d->vStaggered[V3S][k][j][IEND-1];
                b->b1_2Stag[V3S][X1ENDBOUND][k][j] = 2.*(double)BXW1END - mt->vS[V3S][k][j][IEND-1];
                b->b1Stag[V3S][X1ENDBOUND][k][j] = 2.*(double)BXW1END - td->vStaggered[V3S][k][j][IEND-1];
                #endif
            #else
            #if X1UBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V_1][X1BEGBOUND][k][j] = d->vCentroid[V1][k][j][IBEG];
                b->b1_2Stag[V_1][X1BEGBOUND][k][j] = mt->vC[V_1][k][j][IBEG];
                b->b1Stag[V_1][X1BEGBOUND][k][j] = td->vCentroid[V1][k][j][IBEG];
                #elif X1UBEGBOUNDARY == GRAD
                b->b0Stag[V_1][X1BEGBOUND][k][j] = d->vCentroid[V1][k][j][IBEG] - 
                    (double)BXU1BEG*g->dPoints[X1][IBEG];
                b->b1_2Stag[V_1][X1BEGBOUND][k][j] = mt->vC[V_1][k][j][IBEG] - 
                    (double)BXU1BEG*g->dPoints[X1][IBEG];
                b->b1Stag[V_1][X1BEGBOUND][k][j] = td->vCentroid[V1][k][j][IBEG] - 
                    (double)BXU1BEG*g->dPoints[X1][IBEG];
                #elif X1UBEGBOUNDARY == CONST
                b->b0Stag[V_1][X1BEGBOUND][k][j] = 2.*(double)BXU1BEG - d->vCentroid[V1][k][j][IBEG];
                b->b1_2Stag[V_1][X1BEGBOUND][k][j] = 2.*(double)BXU1BEG - mt->vC[V_1][k][j][IBEG];
                b->b1Stag[V_1][X1BEGBOUND][k][j] = 2.*(double)BXU1BEG - td->vCentroid[V1][k][j][IBEG];
                #endif

                #if X1VBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V_2][X1BEGBOUND][k][j] = d->vCentroid[V2][k][j][IBEG];
                b->b1_2Stag[V_2][X1BEGBOUND][k][j] = mt->vC[V_2][k][j][IBEG];
                b->b1Stag[V_2][X1BEGBOUND][k][j] = td->vCentroid[V2][k][j][IBEG];
                #elif X1VBEGBOUNDARY == GRAD
                b->b0Stag[V_2][X1BEGBOUND][k][j] = d->vCentroid[V2][k][j][IBEG] - 
                    (double)BXV1BEG*g->dPoints[X1][IBEG];
                b->b1_2Stag[V_2][X1BEGBOUND][k][j] = mt->vC[V_2][k][j][IBEG] - 
                    (double)BXV1BEG*g->dPoints[X1][IBEG];
                b->b1Stag[V_2][X1BEGBOUND][k][j] = td->vCentroid[V2][k][j][IBEG] - 
                    (double)BXV1BEG*g->dPoints[X1][IBEG];
                #elif X1VBEGBOUNDARY == CONST
                b->b0Stag[V_2][X1BEGBOUND][k][j] = 2.*(double)BXV1BEG - d->vCentroid[V2][k][j][IBEG];
                b->b1_2Stag[V_2][X1BEGBOUND][k][j] = 2.*(double)BXV1BEG - mt->vC[V_2][k][j][IBEG];
                b->b1Stag[V_2][X1BEGBOUND][k][j] = 2.*(double)BXV1BEG - td->vCentroid[V2][k][j][IBEG];
                #endif

                #if X1WBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V_3][X1BEGBOUND][k][j] = d->vCentroid[V3][k][j][IBEG];
                b->b1_2Stag[V_3][X1BEGBOUND][k][j] = mt->vC[V_3][k][j][IBEG];
                b->b1Stag[V_3][X1BEGBOUND][k][j] = td->vCentroid[V3][k][j][IBEG];
                #elif X1WBEGBOUNDARY == GRAD
                b->b0Stag[V_3][X1BEGBOUND][k][j] = d->vCentroid[V3][k][j][IBEG] - 
                    (double)BXW1BEG*g->dPoints[X1][IBEG];
                b->b1_2Stag[V_3][X1BEGBOUND][k][j] = mt->vC[V_3][k][j][IBEG] - 
                    (double)BXW1BEG*g->dPoints[X1][IBEG];
                b->b1Stag[V_3][X1BEGBOUND][k][j] = td->vCentroid[V3][k][j][IBEG] - 
                    (double)BXW1BEG*g->dPoints[X1][IBEG];
                #elif X1WBEGBOUNDARY == CONST
                b->b0Stag[V_3][X1BEGBOUND][k][j] = 2.*(double)BXW1BEG - d->vCentroid[V3][k][j][IEND];
                b->b1_2Stag[V_3][X1BEGBOUND][k][j] = 2.*(double)BXW1BEG - mt->vC[V_3][k][j][IEND];
                b->b1Stag[V_3][X1BEGBOUND][k][j] = 2.*(double)BXW1BEG - td->vCentroid[V3][k][j][IEND];
                #endif

                #if X1UENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V_1][X1ENDBOUND][k][j] = d->vCentroid[V1][k][j][IEND];
                b->b1_2Stag[V_1][X1ENDBOUND][k][j] = mt->vC[V_1][k][j][IEND];
                b->b1Stag[V_1][X1ENDBOUND][k][j] = td->vCentroid[V1][k][j][IEND];
                #elif X1UENDBOUNDARY == GRAD
                b->b0Stag[V_1][X1ENDBOUND][k][j] = d->vCentroid[V1][k][j][IEND] + 
                    (double)BXU1END*g->dPoints[X1][IEND];
                b->b1_2Stag[V_1][X1ENDBOUND][k][j] = mt->vC[V_1][k][j][IEND] + 
                    (double)BXU1END*g->dPoints[X1][IEND];
                b->b1Stag[V_1][X1ENDBOUND][k][j] = td->vCentroid[V1][k][j][IEND] + 
                    (double)BXU1END*g->dPoints[X1][IEND];
                #elif X1UENDBOUNDARY == CONST
                b->b0Stag[V_1][X1ENDBOUND][k][j] = 2.*(double)BXU1END - d->vCentroid[V1][k][j][IEND];
                b->b1_2Stag[V_1][X1ENDBOUND][k][j] = 2.*(double)BXU1END - mt->vC[V_1][k][j][IEND];
                b->b1Stag[V_1][X1ENDBOUND][k][j] = 2.*(double)BXU1END - td->vCentroid[V1][k][j][IEND];
                #endif

                #if X1VENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V_2][X1ENDBOUND][k][j] = d->vCentroid[V2][k][j][IEND];
                b->b1_2Stag[V_2][X1ENDBOUND][k][j] = mt->vC[V_2][k][j][IEND];
                b->b1Stag[V_2][X1ENDBOUND][k][j] = td->vCentroid[V2][k][j][IEND];
                #elif X1VENDBOUNDARY == GRAD
                b->b0Stag[V_2][X1ENDBOUND][k][j] = d->vCentroid[V2][k][j][IEND] + 
                    (double)BXV1END*g->dPoints[X1][IEND];
                b->b1_2Stag[V_2][X1ENDBOUND][k][j] = mt->vC[V_2][k][j][IEND] + 
                    (double)BXV1END*g->dPoints[X1][IEND];
                b->b1Stag[V_2][X1ENDBOUND][k][j] = td->vCentroid[V2][k][j][IEND] + 
                    (double)BXV1END*g->dPoints[X1][IEND];
                #elif X1VENDBOUNDARY == CONST
                b->b0Stag[V_2][X1ENDBOUND][k][j] = 2.*(double)BXV1END - d->vCentroid[V2][k][j][IEND];
                b->b1_2Stag[V_2][X1ENDBOUND][k][j] = 2.*(double)BXV1END - mt->vC[V_2][k][j][IEND];
                b->b1Stag[V_2][X1ENDBOUND][k][j] = 2.*(double)BXV1END - td->vCentroid[V2][k][j][IEND];
                #endif

                #if X1WENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V_3][X1ENDBOUND][k][j] = d->vCentroid[V3][k][j][IEND];
                b->b1_2Stag[V_3][X1ENDBOUND][k][j] = mt->vC[V_3][k][j][IEND];
                b->b1Stag[V_3][X1ENDBOUND][k][j] = td->vCentroid[V3][k][j][IEND];
                #elif X1WENDBOUNDARY == GRAD
                b->b0Stag[V_3][X1ENDBOUND][k][j] = d->vCentroid[V3][k][j][IEND] + 
                    (double)BXW1END*g->dPoints[X1][IEND];
                b->b1_2Stag[V_3][X1ENDBOUND][k][j] = mt->vC[V_3][k][j][IEND] + 
                    (double)BXW1END*g->dPoints[X1][IEND];
                b->b1Stag[V_3][X1ENDBOUND][k][j] = td->vCentroid[V3][k][j][IEND] + 
                    (double)BXW1END*g->dPoints[X1][IEND];
                #elif X1WENDBOUNDARY == CONST
                b->b0Stag[V_3][X1ENDBOUND][k][j] = 2.*(double)BXW1END - d->vCentroid[V3][k][j][IEND];
                b->b1_2Stag[V_3][X1ENDBOUND][k][j] = 2.*(double)BXW1END - mt->vC[V_3][k][j][IEND];
                b->b1Stag[V_3][X1ENDBOUND][k][j] = 2.*(double)BXW1END - td->vCentroid[V3][k][j][IEND];
                #endif
            #endif
        }
    }

    for (int i=0; i<=X1POINTS+1; i++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            #if X2PBEGBOUNDARY == ZERO_GRAD
            b->bCentroids[PRS][X2BEGBOUND][i][k] = td->vCentroid[PRS][k][JBEG][i];
            #elif X2PBEGBOUNDARY == GRAD
            b->bCentroids[PRS][X2BEGBOUND][i][k] = td->vCentroid[PRS][k][JBEG][i] - 
                (double)BXP2BEG*g->dPoints[X2][JBEG];
            #elif X2PBEGBOUNDARY == CONST
            b->bCentroids[PRS][X2BEGBOUND][i][k] = 2.*(double)BXP2BEG - td->vCentroid[PRS][k][JBEG][i];
            #endif

            #if X2PENDBOUNDARY == ZERO_GRAD
            b->bCentroids[PRS][X2ENDBOUND][i][k] = td->vCentroid[PRS][k][JEND][i];
            #elif X2PENDBOUNDARY == GRAD
            b->bCentroids[PRS][X2ENDBOUND][i][k] = td->vCentroid[PRS][k][JEND][i] + 
                (double)BXP2END*g->dPoints[X2][JEND];
            #elif X2PENDBOUNDARY == CONST
            b->bCentroids[PRS][X2ENDBOUND][i][k] = 2.*(double)BXP2END - td->vCentroid[PRS][k][JEND][i];
            #endif

            #if STAGGERED_GRID == YES
                #if X2UBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V1S][X2BEGBOUND][i][k] = d->vStaggered[V1][k][JBEG][i];
                b->b1_2Stag[V1S][X2BEGBOUND][i][k] = mt->vS[V1S][k][JBEG][i];
                b->b1Stag[V1S][X2BEGBOUND][i][k] = td->vStaggered[V1][k][JBEG][i];
                #elif X2UBEGBOUNDARY == GRAD
                b->b0Stag[V1S][X2BEGBOUND][i][k] = d->vStaggered[V1][k][JBEG][i] - 
                    (double)BXU2BEG*g->dPoints[X2][JBEG];
                b->b1_2Stag[V1S][X2BEGBOUND][i][k] = mt->vS[V1S][k][JBEG][i] - 
                    (double)BXU2BEG*g->dPoints[X2][JBEG];
                b->b1Stag[V1S][X2BEGBOUND][i][k] = td->vStaggered[V1][k][JBEG][i] - 
                    (double)BXU2BEG*g->dPoints[X2][JBEG];
                #elif X2UBEGBOUNDARY == CONST
                b->b0Stag[V1S][X2BEGBOUND][i][k] = 2*(double)BXU2BEG - d->vStaggered[V1][k][JBEG][i];
                b->b1_2Stag[V1S][X2BEGBOUND][i][k] = 2*(double)BXU2BEG - mt->vS[V1S][k][JBEG][i];
                b->b1Stag[V1S][X2BEGBOUND][i][k] = 2*(double)BXU2BEG - td->vStaggered[V1][k][JBEG][i];
                #endif

                #if X2VBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V2S][X2BEGBOUND][i][k] = d->vStaggered[V2][k][JBEG][i];
                b->b1_2Stag[V2S][X2BEGBOUND][i][k] = mt->vS[V2S][k][JBEG][i];
                b->b1Stag[V2S][X2BEGBOUND][i][k] = td->vStaggered[V2][k][JBEG][i];
                #elif X2VBEGBOUNDARY == GRAD
                b->b0Stag[V2S][X2BEGBOUND][i][k] = d->vStaggered[V2][k][JBEG][i] - 
                    (double)BXV2BEG*g->dPoints[X2][JBEG];
                b->b1_2Stag[V2S][X2BEGBOUND][i][k] = mt->vS[V2S][k][JBEG][i] - 
                    (double)BXV2BEG*g->dPoints[X2][JBEG];
                b->b1Stag[V2S][X2BEGBOUND][i][k] = td->vStaggered[V2][k][JBEG][i] - 
                    (double)BXV2BEG*g->dPoints[X2][JBEG];
                #elif X2VBEGBOUNDARY == CONST
                b->b0Stag[V2S][X2BEGBOUND][i][k] = (double)BXV2BEG;
                b->b1_2Stag[V2S][X2BEGBOUND][i][k] = (double)BXV2BEG;
                b->b1Stag[V2S][X2BEGBOUND][i][k] = (double)BXV2BEG;
                #endif

                #if X2WBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V3S][X2BEGBOUND][i][k] = d->vStaggered[V3][k][JBEG][i];
                b->b1_2Stag[V3S][X2BEGBOUND][i][k] = mt->vS[V3S][k][JBEG][i];
                b->b1Stag[V3S][X2BEGBOUND][i][k] = td->vStaggered[V3][k][JBEG][i];
                #elif X2WBEGBOUNDARY == GRAD
                b->b0Stag[V3S][X2BEGBOUND][i][k] = d->vStaggered[V3][k][JBEG][i] - 
                    (double)BXW2BEG*g->dPoints[X2][JBEG];
                b->b1_2Stag[V3S][X2BEGBOUND][i][k] = mt->vS[V3S][k][JBEG][i] - 
                    (double)BXW2BEG*g->dPoints[X2][JBEG];
                b->b1Stag[V3S][X2BEGBOUND][i][k] = td->vStaggered[V3][k][JBEG][i] - 
                    (double)BXW2BEG*g->dPoints[X2][JBEG];
                #elif X2WBEGBOUNDARY == CONST
                b->b0Stag[V3S][X2BEGBOUND][i][k] = (double)BXW2BEG;
                b->b1_2Stag[V3S][X2BEGBOUND][i][k] = (double)BXW2BEG;
                b->b1Stag[V3S][X2BEGBOUND][i][k] = (double)BXW2BEG;
                #endif

                #if X2UENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V1S][X2ENDBOUND][i][k] = d->vStaggered[V1][k][JEND][i];
                b->b1_2Stag[V1S][X2ENDBOUND][i][k] = mt->vS[V1S][k][JEND][i];
                b->b1Stag[V1S][X2ENDBOUND][i][k] = td->vStaggered[V1][V1S][k][JEND][i];
                #elif X2UENDBOUNDARY == GRAD
                b->b0Stag[V1S][X2ENDBOUND][i][k] = d->vStaggered[V1][k][JEND][i] + 
                    (double)BXU2END*g->dPoints[X2][JEND];
                b->b1_2Stag[V1S][X2ENDBOUND][i][k] = mt->vS[V1S][k][JEND][i] + 
                    (double)BXU2END*g->dPoints[X2][JEND];
                b->b1Stag[V1S][X2ENDBOUND][i][k] = td->vStaggered[V1][k][JEND][i] + 
                    (double)BXU2END*g->dPoints[X2][JEND];
                #elif X2UENDBOUNDARY == CONST
                b->b0Stag[V1S][X2ENDBOUND][i][k] = 2.*(double)BXU2END - d->vStaggered[V1][k][JEND][i];
                b->b1_2Stag[V1S][X2ENDBOUND][i][k] = 2.*(double)BXU2END - mt->vS[V1S][k][JEND][i];
                b->b1Stag[V1S][X2ENDBOUND][i][k] = 2.*(double)BXU2END - td->vStaggered[V1][k][JEND][i];
                #endif

                #if X2VENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V2S][X2ENDBOUND][i][k] = d->vStaggered[V2][k][JEND][i];
                b->b1_2Stag[V2S][X2ENDBOUND][i][k] = mt->vS[V2S][k][JEND][i];
                b->b1Stag[V2S][X2ENDBOUND][i][k] = td->vStaggered[V2][k][JEND][i];
                #elif X2VENDBOUNDARY == GRAD
                b->b0Stag[V2S][X2ENDBOUND][i][k] = d->vStaggered[V2][k][JEND][i] + 
                    (double)BXV2END*g->dPoints[X2][JEND];
                b->b1_2Stag[V2S][X2ENDBOUND][i][k] = mt->vS[V2S][k][JEND][i] + 
                    (double)BXV2END*g->dPoints[X2][JEND];
                b->b1Stag[V2S][X2ENDBOUND][i][k] = td->vStaggered[V2][k][JEND][i] + 
                    (double)BXV2END*g->dPoints[X2][JEND];
                #elif X2VENDBOUNDARY == CONST
                b->b0Stag[V2S][X2ENDBOUND][i][k] = 2.*(double)BXV2END - d->vStaggered[V2][k][JEND-1][i];
                b->b1_2Stag[V2S][X2ENDBOUND][i][k] = 2.*(double)BXV2END - mt->vS[V2S][k][JEND-1][i];
                b->b1Stag[V2S][X2ENDBOUND][i][k] = 2.*(double)BXV2END - td->vStaggered[V2][k][JEND-1][i];
                #endif

                #if X2WENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V3S][X2ENDBOUND][i][k] = d->vStaggered[V3][k][JEND][i];
                b->b1_2Stag[V3S][X2ENDBOUND][i][k] = mt->vS[V3S][k][JEND][i];
                b->b1Stag[V3S][X2ENDBOUND][i][k] = td->vStaggered[V3][k][JEND][i];
                #elif X2WENDBOUNDARY == GRAD
                b->b0Stag[V3S][X2ENDBOUND][i][k] = d->vStaggered[V3][k][JEND][i] + 
                    (double)BXW2END*g->dPoints[X2][JEND];
                b->b1_2Stag[V3S][X2ENDBOUND][i][k] = mt->vS[V3S][k][JEND][i] + 
                    (double)BXW2END*g->dPoints[X2][JEND];
                b->b1Stag[V3S][X2ENDBOUND][i][k] = td->vStaggered[V3][k][JEND][i] + 
                    (double)BXW2END*g->dPoints[X2][JEND];
                #elif X2WENDBOUNDARY == CONST
                b->b0Stag[V3S][X2ENDBOUND][i][k] = 2.*(double)BXW2END - d->vStaggered[V3][k][JEND-1][i];
                b->b1_2Stag[V3S][X2ENDBOUND][i][k] = 2.*(double)BXW2END - mt->vS[V3S][k][JEND-1][i];
                b->b1Stag[V3S][X2ENDBOUND][i][k] = 2.*(double)BXW2END - td->vStaggered[V3][k][JEND-1][i];
                #endif
            #else
            #if X2UBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V_1][X2BEGBOUND][i][k] = d->vCentroid[V1][k][JBEG][i];
                b->b1_2Stag[V_1][X2BEGBOUND][i][k] = mt->vC[V_1][k][JBEG][i];
                b->b1Stag[V_1][X2BEGBOUND][i][k] = td->vCentroid[V1][k][JBEG][i];
                #elif X2UBEGBOUNDARY == GRAD
                b->b0Stag[V_1][X2BEGBOUND][i][k] = d->vCentroid[V1][k][JBEG][i] - 
                    (double)BXU2BEG*g->dPoints[X2][JBEG];
                b->b1_2Stag[V_1][X2BEGBOUND][i][k] = mt->vC[V_1][k][JBEG][i] - 
                    (double)BXU2BEG*g->dPoints[X2][JBEG];
                b->b1Stag[V_1][X2BEGBOUND][i][k] = td->vCentroid[V1][k][JBEG][i] - 
                    (double)BXU2BEG*g->dPoints[X2][JBEG];
                #elif X2UBEGBOUNDARY == CONST
                b->b0Stag[V_1][X2BEGBOUND][i][k] = 2.*(double)BXU2BEG - d->vCentroid[V1][k][JBEG][i];
                b->b1_2Stag[V_1][X2BEGBOUND][i][k] = 2.*(double)BXU2BEG - mt->vC[V_1][k][JBEG][i];
                b->b1Stag[V_1][X2BEGBOUND][i][k] = 2.*(double)BXU2BEG - td->vCentroid[V1][k][JBEG][i];
                #endif

                #if X2VBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V_2][X2BEGBOUND][i][k] = d->vCentroid[V2][k][JBEG][i];
                b->b1_2Stag[V_2][X2BEGBOUND][i][k] = mt->vC[V_2][k][JBEG][i];
                b->b1Stag[V_2][X2BEGBOUND][i][k] = td->vCentroid[V2][k][JBEG][i];
                #elif X2VBEGBOUNDARY == GRAD
                b->b0Stag[V_2][X2BEGBOUND][i][k] = d->vCentroid[V2][k][JBEG][i] - 
                    (double)BXV2BEG*g->dPoints[X2][JBEG];
                b->b1_2Stag[V_2][X2BEGBOUND][i][k] = mt->vC[V_2][k][JBEG][i] - 
                    (double)BXV2BEG*g->dPoints[X2][JBEG];
                b->b1Stag[V_2][X2BEGBOUND][i][k] = td->vCentroid[V2][k][JBEG][i] - 
                    (double)BXV2BEG*g->dPoints[X2][JBEG];
                #elif X2VBEGBOUNDARY == CONST
                b->b0Stag[V_2][X2BEGBOUND][i][k] = 2.*(double)BXV2BEG - d->vCentroid[V2][k][JBEG][i];
                b->b1_2Stag[V_2][X2BEGBOUND][i][k] = 2.*(double)BXV2BEG - mt->vC[V_2][k][JBEG][i];
                b->b1Stag[V_2][X2BEGBOUND][i][k] = 2.*(double)BXV2BEG - td->vCentroid[V2][k][JBEG][i];
                #endif

                #if X2WBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V_3][X2BEGBOUND][i][k] = d->vCentroid[V3][k][JBEG][i];
                b->b1_2Stag[V_3][X2BEGBOUND][i][k] = mt->vC[V_3][k][JBEG][i];
                b->b1Stag[V_3][X2BEGBOUND][i][k] = td->vCentroid[V3][k][JBEG][i];
                #elif X2WBEGBOUNDARY == GRAD
                b->b0Stag[V_3][X2BEGBOUND][i][k] = d->vCentroid[V3][k][JBEG][i] - 
                    (double)BXW2BEG*g->dPoints[X2][JBEG];
                b->b1_2Stag[V_3][X2BEGBOUND][i][k] = mt->vC[V_3][k][JBEG][i] - 
                    (double)BXW2BEG*g->dPoints[X2][JBEG];
                b->b1Stag[V_3][X2BEGBOUND][i][k] = td->vCentroid[V3][k][JBEG][i] - 
                    (double)BXW2BEG*g->dPoints[X2][JBEG];
                #elif X2WBEGBOUNDARY == CONST
                b->b0Stag[V_3][X2BEGBOUND][i][k] = 2.*(double)BXW2BEG - d->vCentroid[V3][k][JBEG][i];
                b->b1_2Stag[V_3][X2BEGBOUND][i][k] = 2.*(double)BXW2BEG - mt->vC[V_3][k][JBEG][i];
                b->b1Stag[V_3][X2BEGBOUND][i][k] = 2.*(double)BXW2BEG - td->vCentroid[V3][k][JBEG][i];
                #endif

                #if X2UENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V_1][X2ENDBOUND][i][k] = d->vCentroid[V1][k][JEND][i];
                b->b1_2Stag[V_1][X2ENDBOUND][i][k] = mt->vC[V_1][k][JEND][i];
                b->b1Stag[V_1][X2ENDBOUND][i][k] = td->vCentroid[V1][V_1][k][JEND][i];
                #elif X2UENDBOUNDARY == GRAD
                b->b0Stag[V_1][X2ENDBOUND][i][k] = d->vCentroid[V1][k][JEND][i] + 
                    (double)BXU2END*g->dPoints[X2][JEND];
                b->b1_2Stag[V_1][X2ENDBOUND][i][k] = mt->vC[V_1][k][JEND][i] + 
                    (double)BXU2END*g->dPoints[X2][JEND];
                b->b1Stag[V_1][X2ENDBOUND][i][k] = td->vCentroid[V1][k][JEND][i] + 
                    (double)BXU2END*g->dPoints[X2][JEND];
                #elif X2UENDBOUNDARY == CONST
                b->b0Stag[V_1][X2ENDBOUND][i][k] = 2.*(double)BXU2END - d->vCentroid[V1][k][JEND][i];
                b->b1_2Stag[V_1][X2ENDBOUND][i][k] = 2.*(double)BXU2END - mt->vC[V_1][k][JEND][i];
                b->b1Stag[V_1][X2ENDBOUND][i][k] = 2.*(double)BXU2END - td->vCentroid[V1][k][JEND][i];
                #endif

                #if X2VENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V_2][X2ENDBOUND][i][k] = d->vCentroid[V2][k][JEND][i];
                b->b1_2Stag[V_2][X2ENDBOUND][i][k] = mt->vC[V_2][k][JEND][i];
                b->b1Stag[V_2][X2ENDBOUND][i][k] = td->vCentroid[V2][k][JEND][i];
                #elif X2VENDBOUNDARY == GRAD
                b->b0Stag[V_2][X2ENDBOUND][i][k] = d->vCentroid[V2][k][JEND][i] + 
                    (double)BXV2END*g->dPoints[X2][JEND];
                b->b1_2Stag[V_2][X2ENDBOUND][i][k] = mt->vC[V_2][k][JEND][i] + 
                    (double)BXV2END*g->dPoints[X2][JEND];
                b->b1Stag[V_2][X2ENDBOUND][i][k] = td->vCentroid[V2][k][JEND][i] + 
                    (double)BXV2END*g->dPoints[X2][JEND];
                #elif X2VENDBOUNDARY == CONST
                b->b0Stag[V_2][X2ENDBOUND][i][k] = 2.*(double)BXV2END - d->vCentroid[V2][k][JEND][i];
                b->b1_2Stag[V_2][X2ENDBOUND][i][k] = 2.*(double)BXV2END - mt->vC[V_2][k][JEND][i];
                b->b1Stag[V_2][X2ENDBOUND][i][k] = 2.*(double)BXV2END - td->vCentroid[V2][k][JEND][i];
                #endif

                #if X2WENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V_3][X2ENDBOUND][i][k] = d->vCentroid[V3][k][JEND][i];
                b->b1_2Stag[V_3][X2ENDBOUND][i][k] = mt->vC[V_3][k][JEND][i];
                b->b1Stag[V_3][X2ENDBOUND][i][k] = td->vCentroid[V3][k][JEND][i];
                #elif X2WENDBOUNDARY == GRAD
                b->b0Stag[V_3][X2ENDBOUND][i][k] = d->vCentroid[V3][k][JEND][i] + 
                    (double)BXW2END*g->dPoints[X2][JEND];
                b->b1_2Stag[V_3][X2ENDBOUND][i][k] = mt->vC[V_3][k][JEND][i] + 
                    (double)BXW2END*g->dPoints[X2][JEND];
                b->b1Stag[V_3][X2ENDBOUND][i][k] = td->vCentroid[V3][k][JEND][i] + 
                    (double)BXW2END*g->dPoints[X2][JEND];
                #elif X2WENDBOUNDARY == CONST
                b->b0Stag[V_3][X2ENDBOUND][i][k] = 2.*(double)BXW2END - d->vCentroid[V3][k][JEND][i];
                b->b1_2Stag[V_3][X2ENDBOUND][i][k] = 2.*(double)BXW2END - mt->vC[V_3][k][JEND][i];
                b->b1Stag[V_3][X2ENDBOUND][i][k] = 2.*(double)BXW2END - td->vCentroid[V3][k][JEND][i];
                #endif
            #endif
        }
    }
    
    for (int j=0; j<=X2POINTS+1; j++) {
        for (int i=0; i<=X1POINTS+1; i++) {
            #if X3PBEGBOUNDARY == ZERO_GRAD
            b->bCentroids[PRS][X3BEGBOUND][j][i] = td->vCentroid[PRS][KBEG][j][i];
            #elif X3PBEGBOUNDARY == GRAD
            b->bCentroids[PRS][X3BEGBOUND][j][i] = td->vCentroid[PRS][KBEG][j][i] + 
                (double)BXP3BEG*g->dPoints[X3][KBEG];
            #elif X3PBEGBOUNDARY == CONST
            b->bCentroids[PRS][X3BEGBOUND][j][i] = 2.*(double)BXP3BEG - td->vCentroid[PRS][KBEG][j][i];
            #endif

            #if X3PENDBOUNDARY == ZERO_GRAD
            b->bCentroids[PRS][X3ENDBOUND][j][i] = td->vCentroid[PRS][KEND][j][i];
            #elif X3PENDBOUNDARY == GRAD
            b->bCentroids[PRS][X3ENDBOUND][j][i] = td->vCentroid[PRS][KEND][j][i] - 
                (double)BXP3END*g->dPoints[X3][KEND];
            #elif X3PENDBOUNDARY == CONST
            b->bCentroids[PRS][X3ENDBOUND][j][i] = 2.*(double)BXP3END - td->vCentroid[PRS][KEND][j][i];
            #endif

            #if STAGGERED_GRID == YES
                #if X3UBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V1S][X3BEGBOUND][j][i] = d->vStaggered[V1S][KBEG][j][i];
                b->b1_2Stag[V1S][X3BEGBOUND][j][i] = mt->vS[V1S][KBEG][j][i];
                b->b1Stag[V1S][X3BEGBOUND][j][i] = td->vStaggered[V1S][KBEG][j][i];
                #elif X3UBEGBOUNDARY == GRAD
                b->b0Stag[V1S][X3BEGBOUND][j][i] = d->vStaggered[V1S][KBEG][j][i] - 
                    (double)BXU3BEG*g->dPoints[X3][KBEG];
                b->b1_2Stag[V1S][X3BEGBOUND][j][i] = mt->vS[V1S][KBEG][j][i] - 
                    (double)BXU3BEG*g->dPoints[X3][KBEG];
                b->b1Stag[V1S][X3BEGBOUND][j][i] = td->vStaggered[V1S][KBEG][j][i] - 
                    (double)BXU3BEG*g->dPoints[X3][KBEG];
                #elif X3UBEGBOUNDARY == CONST
                b->b0Stag[V1S][X3BEGBOUND][j][i] = (double)BXU3BEG;
                b->b1_2Stag[V1S][X3BEGBOUND][j][i] = (double)BXU3BEG;
                b->b1Stag[V1S][X3BEGBOUND][j][i] = (double)BXU3BEG;
                #endif

                #if X3VBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V2S][X3BEGBOUND][j][i] = d->vStaggered[V2S][KBEG][j][i];
                b->b1_2Stag[V2S][X3BEGBOUND][j][i] = mt->vS[V2S][KBEG][j][i];
                b->b1Stag[V2S][X3BEGBOUND][j][i] = td->vStaggered[V2S][KBEG][j][i];
                #elif X3VBEGBOUNDARY == GRAD
                b->b0Stag[V2S][X3BEGBOUND][j][i] = d->vStaggered[V2S][KBEG][j][i] - 
                    (double)BXV3BEG*g->dPoints[X3][KBEG];
                b->b1_2Stag[V2S][X3BEGBOUND][j][i] = mt->vS[V2S][KBEG][j][i] - 
                    (double)BXV3BEG*g->dPoints[X3][KBEG];
                b->b1Stag[V2S][X3BEGBOUND][j][i] = td->vStaggered[V2S][KBEG][j][i] - 
                    (double)BXV3BEG*g->dPoints[X3][KBEG];
                #elif X3VBEGBOUNDARY == CONST
                b->b0Stag[V2S][X3BEGBOUND][j][i] = (double)BXV3BEG;
                b->b1_2Stag[V2S][X3BEGBOUND][j][i] = (double)BXV3BEG;
                b->b1Stag[V2S][X3BEGBOUND][j][i] = (double)BXV3BEG;
                #endif

                #if X3WBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V3S][X3BEGBOUND][j][i] = d->vStaggered[V3S][KBEG][j][i];
                b->b1_2Stag[V3S][X3BEGBOUND][j][i] = mt->vS[V3S][KBEG][j][i];
                b->b1Stag[V3S][X3BEGBOUND][j][i] = td->vStaggered[V3S][KBEG][j][i];
                #elif X3WBEGBOUNDARY == GRAD
                b->b0Stag[V3S][X3BEGBOUND][j][i] = d->vStaggered[V3S][KBEG][j][i] - 
                    (double)BXW3BEG*g->dPoints[X3][KBEG];
                b->b1_2Stag[V3S][X3BEGBOUND][j][i] = mt->vS[V3S][KBEG][j][i] - 
                    (double)BXW3BEG*g->dPoints[X3][KBEG];
                b->b1Stag[V3S][X3BEGBOUND][j][i] = td->vStaggered[V3S][KBEG][j][i] - 
                    (double)BXW3BEG*g->dPoints[X3][KBEG];
                #elif X3WBEGBOUNDARY == CONST
                b->b0Stag[V3S][X3BEGBOUND][j][i] = (double)BXW3BEG;
                b->b1_2Stag[V3S][X3BEGBOUND][j][i] = (double)BXW3BEG;
                b->b1Stag[V3S][X3BEGBOUND][j][i] = (double)BXW3BEG;
                #endif

                #if X3UENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V1S][X3ENDBOUND][j][i] = d->vStaggered[V1S][KEND][j][i];
                b->b1_2Stag[V1S][X3ENDBOUND][j][i] = mt->vS[V1S][KEND][j][i];
                b->b1Stag[V1S][X3ENDBOUND][j][i] = td->vStaggered[V1S][KEND][j][i];
                #elif X3UENDBOUNDARY == GRAD
                b->b0Stag[V1S][X3ENDBOUND][j][i] = d->vStaggered[V1S][KEND][j][i] + 
                    (double)BXU3END*g->dPoints[X3][KEND];
                b->b1_2Stag[V1S][X3ENDBOUND][j][i] = mt->vS[V1S][KEND][j][i] + 
                    (double)BXU3END*g->dPoints[X3][KEND];
                b->b1Stag[V1S][X3ENDBOUND][j][i] = td->vStaggered[V1S][KEND][j][i] + 
                    (double)BXU3END*g->dPoints[X3][KEND];
                #elif X3UENDBOUNDARY == CONST
                b->b0Stag[V1S][X3ENDBOUND][j][i] = 2.*(double)BXU3END - d->vStaggered[V1S][KEND-1][j][i];
                b->b1_2Stag[V1S][X3ENDBOUND][j][i] = 2.*(double)BXU3END - mt->vS[V1S][KEND-1][j][i];
                b->b1Stag[V1S][X3ENDBOUND][j][i] = 2.*(double)BXU3END - td->vStaggered[V1S][KEND-1][j][i];
                #endif

                #if X3VENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V2S][X3ENDBOUND][j][i] = d->vStaggered[V2S][KEND][j][i];
                b->b1_2Stag[V2S][X3ENDBOUND][j][i] = mt->vS[V2S][KEND][j][i];
                b->b1Stag[V2S][X3ENDBOUND][j][i] = td->vStaggered[V2S][KEND][j][i];
                #elif X3VENDBOUNDARY == GRAD
                b->b0Stag[V2S][X3ENDBOUND][j][i] = d->vStaggered[V2S][KEND][j][i] + 
                    (double)BXV3END*g->dPoints[X3][KEND];
                b->b1_2Stag[V2S][X3ENDBOUND][j][i] = mt->vS[V2S][KEND][j][i] + 
                    (double)BXV3END*g->dPoints[X3][KEND];
                b->b1Stag[V2S][X3ENDBOUND][j][i] = td->vStaggered[V2S][KEND][j][i] + 
                    (double)BXV3END*g->dPoints[X3][KEND];
                #elif X3VENDBOUNDARY == CONST
                b->b0Stag[V2S][X3ENDBOUND][j][i] = 2.*(double)BXV3END - d->vStaggered[V2S][KEND-1][j][i];
                b->b1_2Stag[V2S][X3ENDBOUND][j][i] = 2.*(double)BXV3END - mt->vS[V2S][KEND-1][j][i];
                b->b1Stag[V2S][X3ENDBOUND][j][i] = 2.*(double)BXV3END - td->vStaggered[V2S][KEND-1][j][i];
                #endif

                #if X3WENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V3S][X3ENDBOUND][j][i] = d->vStaggered[V3S][KEND][j][i];
                b->b1_2Stag[V3S][X3ENDBOUND][j][i] = mt->vS[V3S][KEND][j][i];
                b->b1Stag[V3S][X3ENDBOUND][j][i] = td->vStaggered[V3S][KEND][j][i];
                #elif X3WENDBOUNDARY == GRAD
                b->b0Stag[V3S][X3ENDBOUND][j][i] = d->vStaggered[V3S][KEND][j][i] + 
                    (double)BXW3END*g->dPoints[X3][KEND];
                b->b1_2Stag[V3S][X3ENDBOUND][j][i] = mt->vS[V3S][KEND][j][i] + 
                    (double)BXW3END*g->dPoints[X3][KEND];
                b->b1Stag[V3S][X3ENDBOUND][j][i] = td->vStaggered[V3S][KEND][j][i] + 
                    (double)BXW3END*g->dPoints[X3][KEND];
                #elif X3WENDBOUNDARY == CONST
                b->b0Stag[V3S][X3ENDBOUND][j][i] = 2.*(double)BXW3END - d->vStaggered[V3S][KEND-1][j][i];
                b->b1_2Stag[V3S][X3ENDBOUND][j][i] = 2.*(double)BXW3END - mt->vS[V3S][KEND-1][j][i];
                b->b1Stag[V3S][X3ENDBOUND][j][i] = 2.*(double)BXW3END - td->vStaggered[V3S][KEND-1][j][i];
                #endif
            #else
            #if X3UBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V_1][X3BEGBOUND][j][i] = d->vCentroid[V1][KBEG][j][i];
                b->b1_2Stag[V_1][X3BEGBOUND][j][i] = mt->vC[V_1][KBEG][j][i];
                b->b1Stag[V_1][X3BEGBOUND][j][i] = td->vCentroid[V1][KBEG][j][i];
                #elif X3UBEGBOUNDARY == GRAD
                b->b0Stag[V_1][X3BEGBOUND][j][i] = d->vCentroid[V1][KBEG][j][i] - 
                    (double)BXU3BEG*g->dPoints[X3][KBEG];
                b->b1_2Stag[V_1][X3BEGBOUND][j][i] = mt->vC[V_1][KBEG][j][i] - 
                    (double)BXU3BEG*g->dPoints[X3][KBEG];
                b->b1Stag[V_1][X3BEGBOUND][j][i] = td->vCentroid[V1][KBEG][j][i] - 
                    (double)BXU3BEG*g->dPoints[X3][KBEG];
                #elif X3UBEGBOUNDARY == CONST
                b->b0Stag[V_1][X3BEGBOUND][j][i] = 2.*(double)BXU3BEG - d->vCentroid[V1][KBEG][j][i];
                b->b1_2Stag[V_1][X3BEGBOUND][j][i] = 2.*(double)BXU3BEG - mt->vC[V_1][KBEG][j][i];
                b->b1Stag[V_1][X3BEGBOUND][j][i] = 2.*(double)BXU3BEG - td->vCentroid[V1][KBEG][j][i];
                #endif

                #if X3VBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V_2][X3BEGBOUND][j][i] = d->vCentroid[V2][KBEG][j][i];
                b->b1_2Stag[V_2][X3BEGBOUND][j][i] = mt->vC[V_2][KBEG][j][i];
                b->b1Stag[V_2][X3BEGBOUND][j][i] = td->vCentroid[V2][KBEG][j][i];
                #elif X3VBEGBOUNDARY == GRAD
                b->b0Stag[V_2][X3BEGBOUND][j][i] = d->vCentroid[V2][KBEG][j][i] - 
                    (double)BXV3BEG*g->dPoints[X3][KBEG];
                b->b1_2Stag[V_2][X3BEGBOUND][j][i] = mt->vC[V_2][KBEG][j][i] - 
                    (double)BXV3BEG*g->dPoints[X3][KBEG];
                b->b1Stag[V_2][X3BEGBOUND][j][i] = td->vCentroid[V2][KBEG][j][i] - 
                    (double)BXV3BEG*g->dPoints[X3][KBEG];
                #elif X3VBEGBOUNDARY == CONST
                b->b0Stag[V_2][X3BEGBOUND][j][i] = 2.*(double)BXV3BEG - d->vCentroid[V2][KBEG][j][i];
                b->b1_2Stag[V_2][X3BEGBOUND][j][i] = 2.*(double)BXV3BEG - mt->vC[V_2][KBEG][j][i];
                b->b1Stag[V_2][X3BEGBOUND][j][i] = 2.*(double)BXV3BEG - td->vCentroid[V2][KBEG][j][i];
                #endif

                #if X3WBEGBOUNDARY == ZERO_GRAD
                b->b0Stag[V_3][X3BEGBOUND][j][i] = d->vCentroid[V3][KBEG][j][i];
                b->b1_2Stag[V_3][X3BEGBOUND][j][i] = mt->vC[V_3][KBEG][j][i];
                b->b1Stag[V_3][X3BEGBOUND][j][i] = td->vCentroid[V3][KBEG][j][i];
                #elif X3WBEGBOUNDARY == GRAD
                b->b0Stag[V_3][X3BEGBOUND][j][i] = d->vCentroid[V3][KBEG][j][i] - 
                    (double)BXW3BEG*g->dPoints[X3][KBEG];
                b->b1_2Stag[V_3][X3BEGBOUND][j][i] = mt->vC[V_3][KBEG][j][i] - 
                    (double)BXW3BEG*g->dPoints[X3][KBEG];
                b->b1Stag[V_3][X3BEGBOUND][j][i] = td->vCentroid[V3][KBEG][j][i] - 
                    (double)BXW3BEG*g->dPoints[X3][KBEG];
                #elif X3WBEGBOUNDARY == CONST
                b->b0Stag[V_3][X3BEGBOUND][j][i] = 2.*(double)BXW3BEG - d->vCentroid[V3][KBEG][j][i];
                b->b1_2Stag[V_3][X3BEGBOUND][j][i] = 2.*(double)BXW3BEG - mt->vC[V_3][KBEG][j][i];
                b->b1Stag[V_3][X3BEGBOUND][j][i] = 2.*(double)BXW3BEG - td->vCentroid[V3][KBEG][j][i];
                #endif

                #if X3UENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V_1][X3ENDBOUND][j][i] = d->vCentroid[V1][KEND][j][i];
                b->b1_2Stag[V_1][X3ENDBOUND][j][i] = mt->vC[V_1][KEND][j][i];
                b->b1Stag[V_1][X3ENDBOUND][j][i] = td->vCentroid[V1][KEND][j][i];
                #elif X3UENDBOUNDARY == GRAD
                b->b0Stag[V_1][X3ENDBOUND][j][i] = d->vCentroid[V1][KEND][j][i] + 
                    (double)BXU3END*g->dPoints[X3][KEND];
                b->b1_2Stag[V_1][X3ENDBOUND][j][i] = mt->vC[V_1][KEND][j][i] + 
                    (double)BXU3END*g->dPoints[X3][KEND];
                b->b1Stag[V_1][X3ENDBOUND][j][i] = td->vCentroid[V1][KEND][j][i] + 
                    (double)BXU3END*g->dPoints[X3][KEND];
                #elif X3UENDBOUNDARY == CONST
                b->b0Stag[V_1][X3ENDBOUND][j][i] = 2.*(double)BXU3END - d->vCentroid[V1][KEND][j][i];
                b->b1_2Stag[V_1][X3ENDBOUND][j][i] = 2.*(double)BXU3END - mt->vC[V_1][KEND][j][i];
                b->b1Stag[V_1][X3ENDBOUND][j][i] = 2.*(double)BXU3END - td->vCentroid[V1][KEND][j][i];
                #endif

                #if X3VENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V_2][X3ENDBOUND][j][i] = d->vCentroid[V2][KEND][j][i];
                b->b1_2Stag[V_2][X3ENDBOUND][j][i] = mt->vC[V_2][KEND][j][i];
                b->b1Stag[V_2][X3ENDBOUND][j][i] = td->vCentroid[V2][KEND][j][i];
                #elif X3VENDBOUNDARY == GRAD
                b->b0Stag[V_2][X3ENDBOUND][j][i] = d->vCentroid[V2][KEND][j][i] + 
                    (double)BXV3END*g->dPoints[X3][KEND];
                b->b1_2Stag[V_2][X3ENDBOUND][j][i] = mt->vC[V_2][KEND][j][i] + 
                    (double)BXV3END*g->dPoints[X3][KEND];
                b->b1Stag[V_2][X3ENDBOUND][j][i] = td->vCentroid[V2][KEND][j][i] + 
                    (double)BXV3END*g->dPoints[X3][KEND];
                #elif X3VENDBOUNDARY == CONST
                b->b0Stag[V_2][X3ENDBOUND][j][i] = 2.*(double)BXV3END - d->vCentroid[V2][KEND][j][i];
                b->b1_2Stag[V_2][X3ENDBOUND][j][i] = 2.*(double)BXV3END - mt->vC[V_2][KEND][j][i];
                b->b1Stag[V_2][X3ENDBOUND][j][i] = 2.*(double)BXV3END - td->vCentroid[V2][KEND][j][i];
                #endif

                #if X3WENDBOUNDARY == ZERO_GRAD
                b->b0Stag[V_3][X3ENDBOUND][j][i] = d->vCentroid[V3][KEND][j][i];
                b->b1_2Stag[V_3][X3ENDBOUND][j][i] = mt->vC[V_3][KEND][j][i];
                b->b1Stag[V_3][X3ENDBOUND][j][i] = td->vCentroid[V3][KEND][j][i];
                #elif X3WENDBOUNDARY == GRAD
                b->b0Stag[V_3][X3ENDBOUND][j][i] = d->vCentroid[V3][KEND][j][i] + 
                    (double)BXW3END*g->dPoints[X3][KEND];
                b->b1_2Stag[V_3][X3ENDBOUND][j][i] = mt->vC[V_3][KEND][j][i] + 
                    (double)BXW3END*g->dPoints[X3][KEND];
                b->b1Stag[V_3][X3ENDBOUND][j][i] = td->vCentroid[V3][KEND][j][i] + 
                    (double)BXW3END*g->dPoints[X3][KEND];
                #elif X3WENDBOUNDARY == CONST
                b->b0Stag[V_3][X3ENDBOUND][j][i] = 2.*(double)BXW3END - d->vCentroid[V3][KEND][j][i];
                b->b1_2Stag[V_3][X3ENDBOUND][j][i] = 2.*(double)BXW3END - mt->vC[V_3][KEND][j][i];
                b->b1Stag[V_3][X3ENDBOUND][j][i] = 2.*(double)BXW3END - td->vCentroid[V3][KEND][j][i];
                #endif
            #endif
        }
    }

    for (int j=0; j<=X2POINTS+1; j++) {
        for (int i=0; i<=X1POINTS+1; i++) {
            #if STAGGERED_GRID == YES
            td->vCentroid[PRS][KBEG-1][j][i] = b->bCentroids[PRS][X3BEGBOUND][j][i];
            td->vCentroid[PRS][KEND+1][j][i] = b->bCentroids[PRS][X3ENDBOUND][j][i];

            td->vStaggered[V1S][KBEG-1][j][i] = b->b1Stag[V1S][X3BEGBOUND][j][i];
            td->vStaggered[V1S][KEND+1][j][i] = b->b1Stag[V1S][X3ENDBOUND][j][i];
            td->vStaggered[V2S][KBEG-1][j][i] = b->b1Stag[V2S][X3BEGBOUND][j][i];
            td->vStaggered[V2S][KEND+1][j][i] = b->b1Stag[V2S][X3ENDBOUND][j][i];
            td->vStaggered[V3S][KBEG-1][j][i] = b->b1Stag[V3S][X3BEGBOUND][j][i];
            td->vStaggered[V3S][KEND+1][j][i] = b->b1Stag[V3S][X3ENDBOUND][j][i];

            td->vStaggered[V3S][KEND][j][i] = BXW3END;

            d->vStaggered[V1S][KBEG-1][j][i] = b->b0Stag[V1S][X3BEGBOUND][j][i];
            d->vStaggered[V1S][KEND+1][j][i] = b->b0Stag[V1S][X3ENDBOUND][j][i];
            d->vStaggered[V2S][KBEG-1][j][i] = b->b0Stag[V2S][X3BEGBOUND][j][i];
            d->vStaggered[V2S][KEND+1][j][i] = b->b0Stag[V2S][X3ENDBOUND][j][i];
            d->vStaggered[V3S][KBEG-1][j][i] = b->b0Stag[V3S][X3BEGBOUND][j][i];
            d->vStaggered[V3S][KEND+1][j][i] = b->b0Stag[V3S][X3ENDBOUND][j][i];

            d->vStaggered[V3S][KEND][j][i] = BXW3END;

            mt->vS[V1S][KBEG-1][j][i] = b->b1_2Stag[V1S][X3BEGBOUND][j][i];
            mt->vS[V1S][KEND+1][j][i] = b->b1_2Stag[V1S][X3ENDBOUND][j][i];
            mt->vS[V2S][KBEG-1][j][i] = b->b1_2Stag[V2S][X3BEGBOUND][j][i];
            mt->vS[V2S][KEND+1][j][i] = b->b1_2Stag[V2S][X3ENDBOUND][j][i];
            mt->vS[V3S][KBEG-1][j][i] = b->b1_2Stag[V3S][X3BEGBOUND][j][i];
            mt->vS[V3S][KEND+1][j][i] = b->b1_2Stag[V3S][X3ENDBOUND][j][i];

            mt->vS[V3S][KEND][j][i] = BXW3END;
            #else
            td->vCentroid[PRS][KBEG-1][j][i] = b->bCentroids[PRS][X3BEGBOUND][j][i];
            td->vCentroid[PRS][KEND+1][j][i] = b->bCentroids[PRS][X3ENDBOUND][j][i];

            td->vCentroid[V1][KBEG-1][j][i] = b->b1Stag[V_1][X3BEGBOUND][j][i];
            td->vCentroid[V1][KEND+1][j][i] = b->b1Stag[V_1][X3ENDBOUND][j][i];
            td->vCentroid[V2][KBEG-1][j][i] = b->b1Stag[V_2][X3BEGBOUND][j][i];
            td->vCentroid[V2][KEND+1][j][i] = b->b1Stag[V_2][X3ENDBOUND][j][i];
            td->vCentroid[V3][KBEG-1][j][i] = b->b1Stag[V_3][X3BEGBOUND][j][i];
            td->vCentroid[V3][KEND+1][j][i] = b->b1Stag[V_3][X3ENDBOUND][j][i];

            d->vCentroid[V1][KBEG-1][j][i] = b->b0Stag[V_1][X3BEGBOUND][j][i];
            d->vCentroid[V1][KEND+1][j][i] = b->b0Stag[V_1][X3ENDBOUND][j][i];
            d->vCentroid[V2][KBEG-1][j][i] = b->b0Stag[V_2][X3BEGBOUND][j][i];
            d->vCentroid[V2][KEND+1][j][i] = b->b0Stag[V_2][X3ENDBOUND][j][i];
            d->vCentroid[V3][KBEG-1][j][i] = b->b0Stag[V_3][X3BEGBOUND][j][i];
            d->vCentroid[V3][KEND+1][j][i] = b->b0Stag[V_3][X3ENDBOUND][j][i];

            mt->vC[V_1][KBEG-1][j][i] = b->b1_2Stag[V_1][X3BEGBOUND][j][i];
            mt->vC[V_1][KEND+1][j][i] = b->b1_2Stag[V_1][X3ENDBOUND][j][i];
            mt->vC[V_2][KBEG-1][j][i] = b->b1_2Stag[V_2][X3BEGBOUND][j][i];
            mt->vC[V_2][KEND+1][j][i] = b->b1_2Stag[V_2][X3ENDBOUND][j][i];
            mt->vC[V_3][KBEG-1][j][i] = b->b1_2Stag[V_3][X3BEGBOUND][j][i];
            mt->vC[V_3][KEND+1][j][i] = b->b1_2Stag[V_3][X3ENDBOUND][j][i];
            #endif
        }
    }

    for (int i=0; i<=X1POINTS+1; i++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            #if STAGGERED_GRID == YES
            td->vCentroid[PRS][k][JBEG-1][i] = b->bCentroids[PRS][X2BEGBOUND][i][k];
            td->vCentroid[PRS][k][JEND+1][i] = b->bCentroids[PRS][X2ENDBOUND][i][k];

            td->vStaggered[V1S][k][JBEG-1][i] = b->b1Stag[V1S][X2BEGBOUND][i][k];
            td->vStaggered[V1S][k][JEND+1][i] = b->b1Stag[V1S][X2ENDBOUND][i][k];
            td->vStaggered[V2S][k][JBEG-1][i] = b->b1Stag[V2S][X2BEGBOUND][i][k];
            td->vStaggered[V2S][k][JEND+1][i] = b->b1Stag[V2S][X2ENDBOUND][i][k];
            td->vStaggered[V3S][k][JBEG-1][i] = b->b1Stag[V3S][X2BEGBOUND][i][k];
            td->vStaggered[V3S][k][JEND+1][i] = b->b1Stag[V3S][X2ENDBOUND][i][k];

            td->vStaggered[V2S][k][JEND][i] = BXV2END;

            d->vStaggered[V1S][k][JBEG-1][i] = b->b0Stag[V1S][X2BEGBOUND][i][k];
            d->vStaggered[V1S][k][JEND+1][i] = b->b0Stag[V1S][X2ENDBOUND][i][k];
            d->vStaggered[V2S][k][JBEG-1][i] = b->b0Stag[V2S][X2BEGBOUND][i][k];
            d->vStaggered[V2S][k][JEND+1][i] = b->b0Stag[V2S][X2ENDBOUND][i][k];
            d->vStaggered[V3S][k][JBEG-1][i] = b->b0Stag[V3S][X2BEGBOUND][i][k];
            d->vStaggered[V3S][k][JEND+1][i] = b->b0Stag[V3S][X2ENDBOUND][i][k];

            d->vStaggered[V2S][k][JEND][i] = BXV2END;

            mt->vS[V1S][k][JBEG-1][i] =  b->b1_2Stag[V1S][X2BEGBOUND][i][k];
            mt->vS[V1S][k][JEND+1][i] =  b->b1_2Stag[V1S][X2ENDBOUND][i][k];
            mt->vS[V2S][k][JBEG-1][i] =  b->b1_2Stag[V2S][X2BEGBOUND][i][k];
            mt->vS[V2S][k][JEND+1][i] =  b->b1_2Stag[V2S][X2ENDBOUND][i][k];
            mt->vS[V3S][k][JBEG-1][i] =  b->b1_2Stag[V3S][X2BEGBOUND][i][k];
            mt->vS[V3S][k][JEND+1][i] =  b->b1_2Stag[V3S][X2ENDBOUND][i][k];

            mt->vS[V2S][k][JEND][i] =  BXV2END;
            #else
            td->vCentroid[PRS][k][JBEG-1][i] = b->bCentroids[PRS][X2BEGBOUND][i][k];
            td->vCentroid[PRS][k][JEND+1][i] = b->bCentroids[PRS][X2ENDBOUND][i][k];

            td->vCentroid[V1][k][JBEG-1][i] = b->b1Stag[V_1][X2BEGBOUND][i][k];
            td->vCentroid[V1][k][JEND+1][i] = b->b1Stag[V_1][X2ENDBOUND][i][k];
            td->vCentroid[V2][k][JBEG-1][i] = b->b1Stag[V_2][X2BEGBOUND][i][k];
            td->vCentroid[V2][k][JEND+1][i] = b->b1Stag[V_2][X2ENDBOUND][i][k];
            td->vCentroid[V3][k][JBEG-1][i] = b->b1Stag[V_3][X2BEGBOUND][i][k];
            td->vCentroid[V3][k][JEND+1][i] = b->b1Stag[V_3][X2ENDBOUND][i][k];

            d->vCentroid[V1][k][JBEG-1][i] = b->b0Stag[V_1][X2BEGBOUND][i][k];
            d->vCentroid[V1][k][JEND+1][i] = b->b0Stag[V_1][X2ENDBOUND][i][k];
            d->vCentroid[V2][k][JBEG-1][i] = b->b0Stag[V_2][X2BEGBOUND][i][k];
            d->vCentroid[V2][k][JEND+1][i] = b->b0Stag[V_2][X2ENDBOUND][i][k];
            d->vCentroid[V3][k][JBEG-1][i] = b->b0Stag[V_3][X2BEGBOUND][i][k];
            d->vCentroid[V3][k][JEND+1][i] = b->b0Stag[V_3][X2ENDBOUND][i][k];

            mt->vC[V_1][k][JBEG-1][i] =  b->b1_2Stag[V_1][X2BEGBOUND][i][k];
            mt->vC[V_1][k][JEND+1][i] =  b->b1_2Stag[V_1][X2ENDBOUND][i][k];
            mt->vC[V_2][k][JBEG-1][i] =  b->b1_2Stag[V_2][X2BEGBOUND][i][k];
            mt->vC[V_2][k][JEND+1][i] =  b->b1_2Stag[V_2][X2ENDBOUND][i][k];
            mt->vC[V_3][k][JBEG-1][i] =  b->b1_2Stag[V_3][X2BEGBOUND][i][k];
            mt->vC[V_3][k][JEND+1][i] =  b->b1_2Stag[V_3][X2ENDBOUND][i][k];
            #endif
        }
    }

    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            #if STAGGERED_GRID == YES
            td->vCentroid[PRS][k][j][IBEG-1] = b->bCentroids[PRS][X1BEGBOUND][k][j];
            td->vCentroid[PRS][k][j][IEND+1] = b->bCentroids[PRS][X1ENDBOUND][k][j];

            td->vStaggered[V1S][k][j][IBEG-1] = b->b1Stag[V1S][X1BEGBOUND][k][j];
            td->vStaggered[V1S][k][j][IEND+1] = b->b1Stag[V1S][X1ENDBOUND][k][j];
            td->vStaggered[V2S][k][j][IBEG-1] = b->b1Stag[V2S][X1BEGBOUND][k][j];
            td->vStaggered[V2S][k][j][IEND+1] = b->b1Stag[V2S][X1ENDBOUND][k][j];
            td->vStaggered[V3S][k][j][IBEG-1] = b->b1Stag[V3S][X1BEGBOUND][k][j];
            td->vStaggered[V3S][k][j][IEND+1] = b->b1Stag[V3S][X1ENDBOUND][k][j];

            td->vStaggered[V1S][k][j][IEND] = BXU1END;

            d->vStaggered[V1S][k][j][IBEG-1] = b->b0Stag[V1S][X1BEGBOUND][k][j];
            d->vStaggered[V1S][k][j][IEND+1] = b->b0Stag[V1S][X1ENDBOUND][k][j];
            d->vStaggered[V2S][k][j][IBEG-1] = b->b0Stag[V2S][X1BEGBOUND][k][j];
            d->vStaggered[V2S][k][j][IEND+1] = b->b0Stag[V2S][X1ENDBOUND][k][j];
            d->vStaggered[V3S][k][j][IBEG-1] = b->b0Stag[V3S][X1BEGBOUND][k][j];
            d->vStaggered[V3S][k][j][IEND+1] = b->b0Stag[V3S][X1ENDBOUND][k][j];

            d->vStaggered[V1S][k][j][IEND] = BXU1END;

            mt->vS[V1S][k][j][IBEG-1] =  b->b1_2Stag[V1S][X1BEGBOUND][k][j];
            mt->vS[V1S][k][j][IEND+1] =  b->b1_2Stag[V1S][X1ENDBOUND][k][j];
            mt->vS[V2S][k][j][IBEG-1] =  b->b1_2Stag[V2S][X1BEGBOUND][k][j];
            mt->vS[V2S][k][j][IEND+1] =  b->b1_2Stag[V2S][X1ENDBOUND][k][j];
            mt->vS[V3S][k][j][IBEG-1] =  b->b1_2Stag[V3S][X1BEGBOUND][k][j];
            mt->vS[V3S][k][j][IEND+1] =  b->b1_2Stag[V3S][X1ENDBOUND][k][j];

            mt->vS[V1S][k][j][IEND] =  BXU1END;
            #else
            td->vCentroid[PRS][k][j][IBEG-1] = b->bCentroids[PRS][X1BEGBOUND][k][j];
            td->vCentroid[PRS][k][j][IEND+1] = b->bCentroids[PRS][X1ENDBOUND][k][j];

            td->vCentroid[V1][k][j][IBEG-1] = b->b1Stag[V_1][X1BEGBOUND][k][j];
            td->vCentroid[V1][k][j][IEND+1] = b->b1Stag[V_1][X1ENDBOUND][k][j];
            td->vCentroid[V2][k][j][IBEG-1] = b->b1Stag[V_2][X1BEGBOUND][k][j];
            td->vCentroid[V2][k][j][IEND+1] = b->b1Stag[V_2][X1ENDBOUND][k][j];
            td->vCentroid[V3][k][j][IBEG-1] = b->b1Stag[V_3][X1BEGBOUND][k][j];
            td->vCentroid[V3][k][j][IEND+1] = b->b1Stag[V_3][X1ENDBOUND][k][j];

            d->vCentroid[V1][k][j][IBEG-1] = b->b0Stag[V_1][X1BEGBOUND][k][j];
            d->vCentroid[V1][k][j][IEND+1] = b->b0Stag[V_1][X1ENDBOUND][k][j];
            d->vCentroid[V2][k][j][IBEG-1] = b->b0Stag[V_2][X1BEGBOUND][k][j];
            d->vCentroid[V2][k][j][IEND+1] = b->b0Stag[V_2][X1ENDBOUND][k][j];
            d->vCentroid[V3][k][j][IBEG-1] = b->b0Stag[V_3][X1BEGBOUND][k][j];
            d->vCentroid[V3][k][j][IEND+1] = b->b0Stag[V_3][X1ENDBOUND][k][j];

            mt->vC[V_1][k][j][IBEG-1] =  b->b1_2Stag[V_1][X1BEGBOUND][k][j];
            mt->vC[V_1][k][j][IEND+1] =  b->b1_2Stag[V_1][X1ENDBOUND][k][j];
            mt->vC[V_2][k][j][IBEG-1] =  b->b1_2Stag[V_2][X1BEGBOUND][k][j];
            mt->vC[V_2][k][j][IEND+1] =  b->b1_2Stag[V_2][X1ENDBOUND][k][j];
            mt->vC[V_3][k][j][IBEG-1] =  b->b1_2Stag[V_3][X1BEGBOUND][k][j];
            mt->vC[V_3][k][j][IEND+1] =  b->b1_2Stag[V_3][X1ENDBOUND][k][j];
            #endif
        }
    }
    #endif
}

void update_pressure(Boundary* b, Data* d, Grid* g, tData* td, midT* mt) {

    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            #if X1PBEGBOUNDARY == ZERO_GRAD
            b->bCentroids[PRS][X1BEGBOUND][k][j] = td->vCentroid[PRS][k][j][IBEG];
            #elif X1PBEGBOUNDARY == GRAD
            b->bCentroids[PRS][X1BEGBOUND][k][j] = td->vCentroid[PRS][k][j][IBEG] - 
                (double)BXP1BEG*g->dPoints[X1][IBEG];
            #elif X1PBEGBOUNDARY == CONST
            b->bCentroids[PRS][X1BEGBOUND][k][j] = 2.*(double)BXP1BEG - td->vCentroid[PRS][k][j][IBEG];
            #endif

            #if X1PENDBOUNDARY == ZERO_GRAD
            b->bCentroids[PRS][X1ENDBOUND][k][j] = td->vCentroid[PRS][k][j][IEND];
            #elif X1PENDBOUNDARY == GRAD
            b->bCentroids[PRS][X1ENDBOUND][k][j] = td->vCentroid[PRS][k][j][IEND] + 
                (double)BXP1END*g->dPoints[X1][IEND];
            #elif X1PENDBOUNDARY == CONST
            b->bCentroids[PRS][X1ENDBOUND][k][j] = 2.*(double)BXP1END - td->vCentroid[PRS][k][j][IEND];
            #endif
        }
    }

    for (int i=0; i<=X1POINTS+1; i++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            #if X2PBEGBOUNDARY == ZERO_GRAD
            b->bCentroids[PRS][X2BEGBOUND][i][k] = td->vCentroid[PRS][k][JBEG][i];
            #elif X2PBEGBOUNDARY == GRAD
            b->bCentroids[PRS][X2BEGBOUND][i][k] = td->vCentroid[PRS][k][JBEG][i] - 
                (double)BXP2BEG*g->dPoints[X2][JBEG];
            #elif X2PBEGBOUNDARY == CONST
            b->bCentroids[PRS][X2BEGBOUND][i][k] = 2.*(double)BXP2BEG - td->vCentroid[PRS][k][JBEG][i];
            #endif

            #if X2PENDBOUNDARY == ZERO_GRAD
            b->bCentroids[PRS][X2ENDBOUND][i][k] = td->vCentroid[PRS][k][JEND][i];
            #elif X2PENDBOUNDARY == GRAD
            b->bCentroids[PRS][X2ENDBOUND][i][k] = td->vCentroid[PRS][k][JEND][i] + 
                (double)BXP2END*g->dPoints[X2][JEND];
            #elif X2PENDBOUNDARY == CONST
            b->bCentroids[PRS][X2ENDBOUND][i][k] = 2.*(double)BXP2END - td->vCentroid[PRS][k][JEND][i];
            #endif
        }
    }
    
    for (int j=0; j<=X2POINTS+1; j++) {
        for (int i=0; i<=X1POINTS+1; i++) {
            #if X3PBEGBOUNDARY == ZERO_GRAD
            b->bCentroids[PRS][X3BEGBOUND][j][i] = td->vCentroid[PRS][KBEG][j][i];
            #elif X3PBEGBOUNDARY == GRAD
            b->bCentroids[PRS][X3BEGBOUND][j][i] = td->vCentroid[PRS][KBEG][j][i] + 
                (double)BXP3BEG*g->dPoints[X3][KBEG];
            #elif X3PBEGBOUNDARY == CONST
            b->bCentroids[PRS][X3BEGBOUND][j][i] = 2.*(double)BXP3BEG - td->vCentroid[PRS][KBEG][j][i];
            #endif

            #if X3PENDBOUNDARY == ZERO_GRAD
            b->bCentroids[PRS][X3ENDBOUND][j][i] = td->vCentroid[PRS][KEND][j][i];
            #elif X3PENDBOUNDARY == GRAD
            b->bCentroids[PRS][X3ENDBOUND][j][i] = td->vCentroid[PRS][KEND][j][i] - 
                (double)BXP3END*g->dPoints[X3][KEND];
            #elif X3PENDBOUNDARY == CONST
            b->bCentroids[PRS][X3ENDBOUND][j][i] = 2.*(double)BXP3END - td->vCentroid[PRS][KEND][j][i];
            #endif
        }
    }

    for (int j=0; j<=X2POINTS+1; j++) {
        for (int i=0; i<=X1POINTS+1; i++) {
            td->vCentroid[PRS][KBEG-1][j][i] = b->bCentroids[PRS][X3BEGBOUND][j][i];
            td->vCentroid[PRS][KEND+1][j][i] = b->bCentroids[PRS][X3ENDBOUND][j][i];
        }
    }

    for (int i=0; i<=X1POINTS+1; i++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            td->vCentroid[PRS][k][JBEG-1][i] = b->bCentroids[PRS][X2BEGBOUND][i][k];
            td->vCentroid[PRS][k][JEND+1][i] = b->bCentroids[PRS][X2ENDBOUND][i][k];
        }
    }

    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            td->vCentroid[PRS][k][j][IBEG-1] = b->bCentroids[PRS][X1BEGBOUND][k][j];
            td->vCentroid[PRS][k][j][IEND+1] = b->bCentroids[PRS][X1ENDBOUND][k][j];
        }
    }
}

void update_current_velocity(Boundary* b, Data* d, Grid* g, tData* td, midT* mt) {

}

void update_fractional_velocity(Boundary* b, Data* d, Grid* g, tData* td, midT* mt) {
    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            #if STAGGERED_GRID == YES
                #if X1UBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V1S][X1BEGBOUND][k][j] = mt->vS[V1S][k][j][IBEG];
                #elif X1UBEGBOUNDARY == GRAD
                b->b1_2Stag[V1S][X1BEGBOUND][k][j] = mt->vS[V1S][k][j][IBEG] - 
                    (double)BXU1BEG*g->dPoints[X1][IBEG];
                #elif X1UBEGBOUNDARY == CONST
                b->b1_2Stag[V1S][X1BEGBOUND][k][j] = (double)BXU1BEG;
                #endif

                #if X1VBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V2S][X1BEGBOUND][k][j] = mt->vS[V2S][k][j][IBEG];
                #elif X1VBEGBOUNDARY == GRAD
                b->b1_2Stag[V2S][X1BEGBOUND][k][j] = mt->vS[V2S][k][j][IBEG] - 
                    (double)BXV1BEG*g->dPoints[X1][IBEG];
                #elif X1VBEGBOUNDARY == CONST
                b->b1_2Stag[V2S][X1BEGBOUND][k][j] = 2.*(double)BXV1BEG - mt->vS[V2S][k][j][IBEG];
                #endif

                #if X1WBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V3S][X1BEGBOUND][k][j] = mt->vS[V3S][k][j][IBEG];
                #elif X1WBEGBOUNDARY == GRAD
                b->b1_2Stag[V3S][X1BEGBOUND][k][j] = mt->vS[V3S][k][j][IBEG] - 
                    (double)BXW1BEG*g->dPoints[X1][IBEG];
                #elif X1WBEGBOUNDARY == CONST
                b->b1_2Stag[V3S][X1BEGBOUND][k][j] = (double)BXW1BEG;
                #endif

                #if X1UENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V1S][X1ENDBOUND][k][j] = mt->vS[V1S][k][j][IEND];
                #elif X1UENDBOUNDARY == GRAD
                b->b1_2Stag[V1S][X1ENDBOUND][k][j] = mt->vS[V1S][k][j][IEND] + 
                    (double)BXU1END*g->dPoints[X1][IEND];
                #elif X1UENDBOUNDARY == CONST
                b->b1_2Stag[V1S][X1ENDBOUND][k][j] = 2.*(double)BXU1END - mt->vS[V1S][k][j][IEND-1];
                #endif

                #if X1VENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V2S][X1ENDBOUND][k][j] = mt->vS[V2S][k][j][IEND];
                #elif X1VENDBOUNDARY == GRAD
                b->b1_2Stag[V2S][X1ENDBOUND][k][j] = mt->vS[V2S][k][j][IEND] + 
                    (double)BXV1END*g->dPoints[X1][IEND];
                #elif X1VENDBOUNDARY == CONST
                b->b1_2Stag[V2S][X1ENDBOUND][k][j] = 2.*(double)BXV1END - mt->vS[V2S][k][j][IEND];
                #endif

                #if X1WENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V3S][X1ENDBOUND][k][j] = mt->vS[V3S][k][j][IEND];
                #elif X1WENDBOUNDARY == GRAD
                b->b1_2Stag[V3S][X1ENDBOUND][k][j] = mt->vS[V3S][k][j][IEND] + 
                    (double)BXW1END*g->dPoints[X1][IEND];
                #elif X1WENDBOUNDARY == CONST
                b->b1_2Stag[V3S][X1ENDBOUND][k][j] = 2.*(double)BXW1END - mt->vS[V3S][k][j][IEND-1];
                #endif
            #else
                #if X1UBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_1][X1BEGBOUND][k][j] = mt->vC[V_1][k][j][IBEG];
                #elif X1UBEGBOUNDARY == GRAD
                b->b1_2Stag[V_1][X1BEGBOUND][k][j] = mt->vC[V_1][k][j][IBEG] - 
                    (double)BXU1BEG*g->dPoints[X1][IBEG];
                #elif X1UBEGBOUNDARY == CONST
                b->b1_2Stag[V_1][X1BEGBOUND][k][j] = 2.*(double)BXU1BEG - mt->vC[V_1][k][j][IBEG];
                #endif

                #if X1VBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_2][X1BEGBOUND][k][j] = mt->vC[V_2][k][j][IBEG];
                #elif X1VBEGBOUNDARY == GRAD
                b->b1_2Stag[V_2][X1BEGBOUND][k][j] = mt->vC[V_2][k][j][IBEG] - 
                    (double)BXV1BEG*g->dPoints[X1][IBEG];
                #elif X1VBEGBOUNDARY == CONST
                b->b1_2Stag[V_2][X1BEGBOUND][k][j] = 2.*(double)BXV1BEG - mt->vC[V_2][k][j][IBEG];
                #endif

                #if X1WBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_3][X1BEGBOUND][k][j] = mt->vC[V_3][k][j][IBEG];
                #elif X1WBEGBOUNDARY == GRAD
                b->b1_2Stag[V_3][X1BEGBOUND][k][j] = mt->vC[V_3][k][j][IBEG] - 
                    (double)BXW1BEG*g->dPoints[X1][IBEG];
                #elif X1WBEGBOUNDARY == CONST
                b->b1_2Stag[V_3][X1BEGBOUND][k][j] = 2.*(double)BXW1BEG - mt->vC[V_3][k][j][IBEG];
                #endif

                #if X1UENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_1][X1ENDBOUND][k][j] = mt->vC[V_1][k][j][IEND];
                #elif X1UENDBOUNDARY == GRAD
                b->b1_2Stag[V_1][X1ENDBOUND][k][j] = mt->vC[V_1][k][j][IEND] + 
                    (double)BXU1END*g->dPoints[X1][IEND];
                #elif X1UENDBOUNDARY == CONST
                b->b1_2Stag[V_1][X1ENDBOUND][k][j] = 2.*(double)BXU1END - mt->vC[V_1][k][j][IEND];
                #endif

                #if X1VENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_2][X1ENDBOUND][k][j] = mt->vC[V_2][k][j][IEND];
                #elif X1VENDBOUNDARY == GRAD
                b->b1_2Stag[V_2][X1ENDBOUND][k][j] = mt->vC[V_2][k][j][IEND] + 
                    (double)BXV1END*g->dPoints[X1][IEND];
                #elif X1VENDBOUNDARY == CONST
                b->b1_2Stag[V_2][X1ENDBOUND][k][j] = 2.*(double)BXV1END - mt->vC[V_2][k][j][IEND];
                #endif

                #if X1WENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_3][X1ENDBOUND][k][j] = mt->vC[V_3][k][j][IEND];
                #elif X1WENDBOUNDARY == GRAD
                b->b1_2Stag[V_3][X1ENDBOUND][k][j] = mt->vC[V_3][k][j][IEND] + 
                    (double)BXW1END*g->dPoints[X1][IEND];
                #elif X1WENDBOUNDARY == CONST
                b->b1_2Stag[V_3][X1ENDBOUND][k][j] = 2.*(double)BXW1END - mt->vC[V_3][k][j][IEND];
                #endif
            #endif
        }
    }

    for (int i=0; i<=X1POINTS+1; i++) {
        for (int k=0; k<=X3POINTS+1; k++) {

            #if STAGGERED_GRID == YES
                #if X2UBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V1S][X2BEGBOUND][i][k] = mt->vS[V1S][k][JBEG][i];
                #elif X2UBEGBOUNDARY == GRAD
                b->b1_2Stag[V1S][X2BEGBOUND][i][k] = mt->vS[V1S][k][JBEG][i] - 
                    (double)BXU2BEG*g->dPoints[X2][JBEG];
                #elif X2UBEGBOUNDARY == CONST
                b->b1_2Stag[V1S][X2BEGBOUND][i][k] = 2*(double)BXU2BEG - mt->vS[V1S][k][JBEG][i];
                #endif

                #if X2VBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V2S][X2BEGBOUND][i][k] = mt->vS[V2S][k][JBEG][i];
                #elif X2VBEGBOUNDARY == GRAD
                b->b1_2Stag[V2S][X2BEGBOUND][i][k] = mt->vS[V2S][k][JBEG][i] - 
                    (double)BXV2BEG*g->dPoints[X2][JBEG];
                #elif X2VBEGBOUNDARY == CONST
                b->b1_2Stag[V2S][X2BEGBOUND][i][k] = (double)BXV2BEG;
                #endif

                #if X2WBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V3S][X2BEGBOUND][i][k] = mt->vS[V3S][k][JBEG][i];
                #elif X2WBEGBOUNDARY == GRAD
                b->b1_2Stag[V3S][X2BEGBOUND][i][k] = mt->vS[V3S][k][JBEG][i] - 
                    (double)BXW2BEG*g->dPoints[X2][JBEG];
                #elif X2WBEGBOUNDARY == CONST
                b->b1_2Stag[V3S][X2BEGBOUND][i][k] = (double)BXW2BEG;
                #endif

                #if X2UENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V1S][X2ENDBOUND][i][k] = mt->vS[V1S][k][JEND][i];
                #elif X2UENDBOUNDARY == GRAD
                b->b1_2Stag[V1S][X2ENDBOUND][i][k] = mt->vS[V1S][k][JEND][i] + 
                    (double)BXU2END*g->dPoints[X2][JEND];
                #elif X2UENDBOUNDARY == CONST
                b->b1_2Stag[V1S][X2ENDBOUND][i][k] = 2.*(double)BXU2END - mt->vS[V1S][k][JEND][i];
                #endif

                #if X2VENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V2S][X2ENDBOUND][i][k] = mt->vS[V2S][k][JEND][i];
                #elif X2VENDBOUNDARY == GRAD
                b->b1_2Stag[V2S][X2ENDBOUND][i][k] = mt->vS[V2S][k][JEND][i] + 
                    (double)BXV2END*g->dPoints[X2][JEND];
                #elif X2VENDBOUNDARY == CONST
                b->b1_2Stag[V2S][X2ENDBOUND][i][k] = 2.*(double)BXV2END - mt->vS[V2S][k][JEND-1][i];
                #endif

                #if X2WENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V3S][X2ENDBOUND][i][k] = mt->vS[V3S][k][JEND][i];
                #elif X2WENDBOUNDARY == GRAD
                b->b1_2Stag[V3S][X2ENDBOUND][i][k] = mt->vS[V3S][k][JEND][i] + 
                    (double)BXW2END*g->dPoints[X2][JEND];
                #elif X2WENDBOUNDARY == CONST
                b->b1_2Stag[V3S][X2ENDBOUND][i][k] = 2.*(double)BXW2END - mt->vS[V3S][k][JEND-1][i];
                #endif
            #else
            #if X2UBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_1][X2BEGBOUND][i][k] = mt->vC[V_1][k][JBEG][i];
                #elif X2UBEGBOUNDARY == GRAD
                b->b1_2Stag[V_1][X2BEGBOUND][i][k] = mt->vC[V_1][k][JBEG][i] - 
                    (double)BXU2BEG*g->dPoints[X2][JBEG];
                #elif X2UBEGBOUNDARY == CONST
                b->b1_2Stag[V_1][X2BEGBOUND][i][k] = 2*(double)BXU2BEG - mt->vC[V_1][k][JBEG][i];
                #endif

                #if X2VBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_2][X2BEGBOUND][i][k] = mt->vC[V_2][k][JBEG][i];
                #elif X2VBEGBOUNDARY == GRAD
                b->b1_2Stag[V_2][X2BEGBOUND][i][k] = mt->vC[V_2][k][JBEG][i] - 
                    (double)BXV2BEG*g->dPoints[X2][JBEG];
                #elif X2VBEGBOUNDARY == CONST
                b->b1_2Stag[V_2][X2BEGBOUND][i][k] = 2.*(double)BXV2BEG - mt->vC[V_2][k][JBEG][i];
                #endif

                #if X2WBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_3][X2BEGBOUND][i][k] = mt->vC[V_3][k][JBEG][i];
                #elif X2WBEGBOUNDARY == GRAD
                b->b1_2Stag[V_3][X2BEGBOUND][i][k] = mt->vC[V_3][k][JBEG][i] - 
                    (double)BXW2BEG*g->dPoints[X2][JBEG];
                #elif X2WBEGBOUNDARY == CONST
                b->b1_2Stag[V_3][X2BEGBOUND][i][k] = 2.*(double)BXW2BEG - mt->vC[V_3][k][JBEG][i];
                #endif

                #if X2UENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_1][X2ENDBOUND][i][k] = mt->vC[V_1][k][JEND][i];
                #elif X2UENDBOUNDARY == GRAD
                b->b1_2Stag[V_1][X2ENDBOUND][i][k] = mt->vC[V_1][k][JEND][i] + 
                    (double)BXU2END*g->dPoints[X2][JEND];
                #elif X2UENDBOUNDARY == CONST
                b->b1_2Stag[V_1][X2ENDBOUND][i][k] = 2.*(double)BXU2END - mt->vC[V_1][k][JEND][i];
                #endif

                #if X2VENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_2][X2ENDBOUND][i][k] = mt->vC[V_2][k][JEND][i];
                #elif X2VENDBOUNDARY == GRAD
                b->b1_2Stag[V_2][X2ENDBOUND][i][k] = mt->vC[V_2][k][JEND][i] + 
                    (double)BXV2END*g->dPoints[X2][JEND];
                #elif X2VENDBOUNDARY == CONST
                b->b1_2Stag[V_2][X2ENDBOUND][i][k] = 2.*(double)BXV2END - mt->vC[V_2][k][JEND][i];
                #endif

                #if X2WENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_3][X2ENDBOUND][i][k] = mt->vC[V_3][k][JEND][i];
                #elif X2WENDBOUNDARY == GRAD
                b->b1_2Stag[V_3][X2ENDBOUND][i][k] = mt->vC[V_3][k][JEND][i] + 
                    (double)BXW2END*g->dPoints[X2][JEND];
                #elif X2WENDBOUNDARY == CONST
                b->b1_2Stag[V_3][X2ENDBOUND][i][k] = 2.*(double)BXW2END - mt->vC[V_3][k][JEND][i];
                #endif
            #endif
        }
    }
    
    for (int j=0; j<=X2POINTS+1; j++) {
        for (int i=0; i<=X1POINTS+1; i++) {
            #if STAGGERED_GRID == YES
                #if X3UBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V1S][X3BEGBOUND][j][i] = mt->vS[V1S][KBEG][j][i];
                #elif X3UBEGBOUNDARY == GRAD
                b->b1_2Stag[V1S][X3BEGBOUND][j][i] = mt->vS[V1S][KBEG][j][i] - 
                    (double)BXU3BEG*g->dPoints[X3][KBEG];
                #elif X3UBEGBOUNDARY == CONST
                b->b1_2Stag[V1S][X3BEGBOUND][j][i] = (double)BXU3BEG;
                #endif

                #if X3VBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V2S][X3BEGBOUND][j][i] = mt->vS[V2S][KBEG][j][i];
                #elif X3VBEGBOUNDARY == GRAD
                b->b1_2Stag[V2S][X3BEGBOUND][j][i] = mt->vS[V2S][KBEG][j][i] - 
                    (double)BXV3BEG*g->dPoints[X3][KBEG];
                #elif X3VBEGBOUNDARY == CONST
                b->b1_2Stag[V2S][X3BEGBOUND][j][i] = (double)BXV3BEG;
                #endif

                #if X3WBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V3S][X3BEGBOUND][j][i] = mt->vS[V3S][KBEG][j][i];
                #elif X3WBEGBOUNDARY == GRAD
                b->b1_2Stag[V3S][X3BEGBOUND][j][i] = mt->vS[V3S][KBEG][j][i] - 
                    (double)BXW3BEG*g->dPoints[X3][KBEG];
                #elif X3WBEGBOUNDARY == CONST
                b->b1_2Stag[V3S][X3BEGBOUND][j][i] = (double)BXW3BEG;
                #endif

                #if X3UENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V1S][X3ENDBOUND][j][i] = mt->vS[V1S][KEND][j][i];
                #elif X3UENDBOUNDARY == GRAD
                b->b1_2Stag[V1S][X3ENDBOUND][j][i] = mt->vS[V1S][KEND][j][i] + 
                    (double)BXU3END*g->dPoints[X3][KEND];
                #elif X3UENDBOUNDARY == CONST
                b->b1_2Stag[V1S][X3ENDBOUND][j][i] = 2.*(double)BXU3END - mt->vS[V1S][KEND-1][j][i];
                #endif

                #if X3VENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V2S][X3ENDBOUND][j][i] = mt->vS[V2S][KEND][j][i];
                #elif X3VENDBOUNDARY == GRAD
                b->b1_2Stag[V2S][X3ENDBOUND][j][i] = mt->vS[V2S][KEND][j][i] + 
                    (double)BXV3END*g->dPoints[X3][KEND];
                #elif X3VENDBOUNDARY == CONST
                b->b1_2Stag[V2S][X3ENDBOUND][j][i] = 2.*(double)BXV3END - mt->vS[V2S][KEND-1][j][i];
                #endif

                #if X3WENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V3S][X3ENDBOUND][j][i] = mt->vS[V3S][KEND][j][i];
                #elif X3WENDBOUNDARY == GRAD
                b->b1_2Stag[V3S][X3ENDBOUND][j][i] = mt->vS[V3S][KEND][j][i] + 
                    (double)BXW3END*g->dPoints[X3][KEND];
                #elif X3WENDBOUNDARY == CONST
                b->b1_2Stag[V3S][X3ENDBOUND][j][i] = 2.*(double)BXW3END - mt->vS[V3S][KEND-1][j][i];
                #endif
            #else
            #if X3UBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_1][X3BEGBOUND][j][i] = mt->vC[V_1][KBEG][j][i];
                #elif X3UBEGBOUNDARY == GRAD
                b->b1_2Stag[V_1][X3BEGBOUND][j][i] = mt->vC[V_1][KBEG][j][i] - 
                    (double)BXU3BEG*g->dPoints[X3][KBEG];
                #elif X3UBEGBOUNDARY == CONST
                b->b1_2Stag[V_1][X3BEGBOUND][j][i] = 2.*(double)BXU3BEG - mt->vC[V_1][KBEG][j][i];
                #endif

                #if X3VBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_2][X3BEGBOUND][j][i] = mt->vC[V_2][KBEG][j][i];
                #elif X3VBEGBOUNDARY == GRAD
                b->b1_2Stag[V_2][X3BEGBOUND][j][i] = mt->vC[V_2][KBEG][j][i] - 
                    (double)BXV3BEG*g->dPoints[X3][KBEG];
                #elif X3VBEGBOUNDARY == CONST
                b->b1_2Stag[V_2][X3BEGBOUND][j][i] = 2.*(double)BXV3BEG - mt->vC[V_2][KBEG][j][i];
                #endif

                #if X3WBEGBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_3][X3BEGBOUND][j][i] = mt->vC[V_3][KBEG][j][i];
                #elif X3WBEGBOUNDARY == GRAD
                b->b1_2Stag[V_3][X3BEGBOUND][j][i] = mt->vC[V_3][KBEG][j][i] - 
                    (double)BXW3BEG*g->dPoints[X3][KBEG];
                #elif X3WBEGBOUNDARY == CONST
                b->b1_2Stag[V_3][X3BEGBOUND][j][i] = 2.*(double)BXW3BEG - mt->vC[V_3][KBEG][j][i];
                #endif

                #if X3UENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_1][X3ENDBOUND][j][i] = mt->vC[V_1][KEND][j][i];
                #elif X3UENDBOUNDARY == GRAD
                b->b1_2Stag[V_1][X3ENDBOUND][j][i] = mt->vC[V_1][KEND][j][i] + 
                    (double)BXU3END*g->dPoints[X3][KEND];
                #elif X3UENDBOUNDARY == CONST
                b->b1_2Stag[V_1][X3ENDBOUND][j][i] = 2.*(double)BXU3END - mt->vC[V_1][KEND][j][i];
                #endif

                #if X3VENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_2][X3ENDBOUND][j][i] = mt->vC[V_2][KEND][j][i];
                #elif X3VENDBOUNDARY == GRAD
                b->b1_2Stag[V_2][X3ENDBOUND][j][i] = mt->vC[V_2][KEND][j][i] + 
                    (double)BXV3END*g->dPoints[X3][KEND];
                #elif X3VENDBOUNDARY == CONST
                b->b1_2Stag[V_2][X3ENDBOUND][j][i] = 2.*(double)BXV3END - mt->vC[V_2][KEND][j][i];
                #endif

                #if X3WENDBOUNDARY == ZERO_GRAD
                b->b1_2Stag[V_3][X3ENDBOUND][j][i] = mt->vC[V_3][KEND][j][i];
                #elif X3WENDBOUNDARY == GRAD
                b->b1_2Stag[V_3][X3ENDBOUND][j][i] = mt->vC[V_3][KEND][j][i] + 
                    (double)BXW3END*g->dPoints[X3][KEND];
                #elif X3WENDBOUNDARY == CONST
                b->b1_2Stag[V_3][X3ENDBOUND][j][i] = 2.*(double)BXW3END - mt->vC[V_3][KEND][j][i];
                #endif
            #endif
        }
    }

    for (int j=0; j<=X2POINTS+1; j++) {
        for (int i=0; i<=X1POINTS+1; i++) {
            #if STAGGERED_GRID == YES
            mt->vS[V1S][KBEG-1][j][i] = b->b1_2Stag[V1S][X3BEGBOUND][j][i];
            mt->vS[V1S][KEND+1][j][i] = b->b1_2Stag[V1S][X3ENDBOUND][j][i];
            mt->vS[V2S][KBEG-1][j][i] = b->b1_2Stag[V2S][X3BEGBOUND][j][i];
            mt->vS[V2S][KEND+1][j][i] = b->b1_2Stag[V2S][X3ENDBOUND][j][i];
            mt->vS[V3S][KBEG-1][j][i] = b->b1_2Stag[V3S][X3BEGBOUND][j][i];
            mt->vS[V3S][KEND+1][j][i] = b->b1_2Stag[V3S][X3ENDBOUND][j][i];

            mt->vS[V3S][KEND][j][i] = BXW3END;
            #else
            mt->vC[V_1][KBEG-1][j][i] = b->b1_2Stag[V_1][X3BEGBOUND][j][i];
            mt->vC[V_1][KEND+1][j][i] = b->b1_2Stag[V_1][X3ENDBOUND][j][i];
            mt->vC[V_2][KBEG-1][j][i] = b->b1_2Stag[V_2][X3BEGBOUND][j][i];
            mt->vC[V_2][KEND+1][j][i] = b->b1_2Stag[V_2][X3ENDBOUND][j][i];
            mt->vC[V_3][KBEG-1][j][i] = b->b1_2Stag[V_3][X3BEGBOUND][j][i];
            mt->vC[V_3][KEND+1][j][i] = b->b1_2Stag[V_3][X3ENDBOUND][j][i];
            #endif
        }
    }

    for (int i=0; i<=X1POINTS+1; i++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            #if STAGGERED_GRID == YES
            mt->vS[V1S][k][JBEG-1][i] =  b->b1_2Stag[V1S][X2BEGBOUND][i][k];
            mt->vS[V1S][k][JEND+1][i] =  b->b1_2Stag[V1S][X2ENDBOUND][i][k];
            mt->vS[V2S][k][JBEG-1][i] =  b->b1_2Stag[V2S][X2BEGBOUND][i][k];
            mt->vS[V2S][k][JEND+1][i] =  b->b1_2Stag[V2S][X2ENDBOUND][i][k];
            mt->vS[V3S][k][JBEG-1][i] =  b->b1_2Stag[V3S][X2BEGBOUND][i][k];
            mt->vS[V3S][k][JEND+1][i] =  b->b1_2Stag[V3S][X2ENDBOUND][i][k];

            mt->vS[V2S][k][JEND][i] =  BXV2END;
            #else
            mt->vC[V_1][k][JBEG-1][i] =  b->b1_2Stag[V_1][X2BEGBOUND][i][k];
            mt->vC[V_1][k][JEND+1][i] =  b->b1_2Stag[V_1][X2ENDBOUND][i][k];
            mt->vC[V_2][k][JBEG-1][i] =  b->b1_2Stag[V_2][X2BEGBOUND][i][k];
            mt->vC[V_2][k][JEND+1][i] =  b->b1_2Stag[V_2][X2ENDBOUND][i][k];
            mt->vC[V_3][k][JBEG-1][i] =  b->b1_2Stag[V_3][X2BEGBOUND][i][k];
            mt->vC[V_3][k][JEND+1][i] =  b->b1_2Stag[V_3][X2ENDBOUND][i][k];
            #endif
        }
    }

    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            #if STAGGERED_GRID == YES
            mt->vS[V1S][k][j][IBEG-1] =  b->b1_2Stag[V1S][X1BEGBOUND][k][j];
            mt->vS[V1S][k][j][IEND+1] =  b->b1_2Stag[V1S][X1ENDBOUND][k][j];
            mt->vS[V2S][k][j][IBEG-1] =  b->b1_2Stag[V2S][X1BEGBOUND][k][j];
            mt->vS[V2S][k][j][IEND+1] =  b->b1_2Stag[V2S][X1ENDBOUND][k][j];
            mt->vS[V3S][k][j][IBEG-1] =  b->b1_2Stag[V3S][X1BEGBOUND][k][j];
            mt->vS[V3S][k][j][IEND+1] =  b->b1_2Stag[V3S][X1ENDBOUND][k][j];

            mt->vS[V1S][k][j][IEND] =  BXU1END;
            #else
            mt->vC[V_1][k][j][IBEG-1] =  b->b1_2Stag[V_1][X1BEGBOUND][k][j];
            mt->vC[V_1][k][j][IEND+1] =  b->b1_2Stag[V_1][X1ENDBOUND][k][j];
            mt->vC[V_2][k][j][IBEG-1] =  b->b1_2Stag[V_2][X1BEGBOUND][k][j];
            mt->vC[V_2][k][j][IEND+1] =  b->b1_2Stag[V_2][X1ENDBOUND][k][j];
            mt->vC[V_3][k][j][IBEG-1] =  b->b1_2Stag[V_3][X1BEGBOUND][k][j];
            mt->vC[V_3][k][j][IEND+1] =  b->b1_2Stag[V_3][X1ENDBOUND][k][j];
            #endif
        }
    }
}

void update_next_velocity(Boundary* b, Data* d, Grid* g, tData* td, midT* mt) {
    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            #if STAGGERED_GRID == YES
                #if X1UBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V1S][X1BEGBOUND][k][j] = td->vStaggered[V1S][k][j][IBEG];
                #elif X1UBEGBOUNDARY == GRAD
                b->b1Stag[V1S][X1BEGBOUND][k][j] = td->vStaggered[V1S][k][j][IBEG] - 
                    (double)BXU1BEG*g->dPoints[X1][IBEG];
                #elif X1UBEGBOUNDARY == CONST
                b->b1Stag[V1S][X1BEGBOUND][k][j] = (double)BXU1BEG;
                #endif

                #if X1VBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V2S][X1BEGBOUND][k][j] = td->vStaggered[V2S][k][j][IBEG];
                #elif X1VBEGBOUNDARY == GRAD
                b->b1Stag[V2S][X1BEGBOUND][k][j] = td->vStaggered[V2S][k][j][IBEG] - 
                    (double)BXV1BEG*g->dPoints[X1][IBEG];
                #elif X1VBEGBOUNDARY == CONST
                b->b1Stag[V2S][X1BEGBOUND][k][j] = 2.*(double)BXV1BEG - td->vStaggered[V2S][k][j][IBEG];
                #endif

                #if X1WBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V3S][X1BEGBOUND][k][j] = td->vStaggered[V3S][k][j][IBEG];
                #elif X1WBEGBOUNDARY == GRAD
                b->b1Stag[V3S][X1BEGBOUND][k][j] = td->vStaggered[V3S][k][j][IBEG] - 
                    (double)BXW1BEG*g->dPoints[X1][IBEG];
                #elif X1WBEGBOUNDARY == CONST
                b->b1Stag[V3S][X1BEGBOUND][k][j] = (double)BXW1BEG;
                #endif

                #if X1UENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V1S][X1ENDBOUND][k][j] = td->vStaggered[V1S][k][j][IEND];
                #elif X1UENDBOUNDARY == GRAD
                b->b1Stag[V1S][X1ENDBOUND][k][j] = td->vStaggered[V1S][k][j][IEND] + 
                    (double)BXU1END*g->dPoints[X1][IEND];
                #elif X1UENDBOUNDARY == CONST
                b->b1Stag[V1S][X1ENDBOUND][k][j] = 2.*(double)BXU1END - td->vStaggered[V1S][k][j][IEND-1];
                #endif

                #if X1VENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V2S][X1ENDBOUND][k][j] = td->vStaggered[V2S][k][j][IEND];
                #elif X1VENDBOUNDARY == GRAD
                b->b1Stag[V2S][X1ENDBOUND][k][j] = td->vStaggered[V2S][k][j][IEND] + 
                    (double)BXV1END*g->dPoints[X1][IEND];
                #elif X1VENDBOUNDARY == CONST
                b->b1Stag[V2S][X1ENDBOUND][k][j] = 2.*(double)BXV1END - td->vStaggered[V2S][k][j][IEND];
                #endif

                #if X1WENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V3S][X1ENDBOUND][k][j] = td->vStaggered[V3S][k][j][IEND];
                #elif X1WENDBOUNDARY == GRAD
                b->b1Stag[V3S][X1ENDBOUND][k][j] = td->vStaggered[V3S][k][j][IEND] + 
                    (double)BXW1END*g->dPoints[X1][IEND];
                #elif X1WENDBOUNDARY == CONST
                b->b1Stag[V3S][X1ENDBOUND][k][j] = 2.*(double)BXW1END - td->vStaggered[V3S][k][j][IEND-1];
                #endif
            #else
            #if X1UBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V_1][X1BEGBOUND][k][j] = td->vCentroid[V1][k][j][IBEG];
                #elif X1UBEGBOUNDARY == GRAD
                b->b1Stag[V_1][X1BEGBOUND][k][j] = td->vCentroid[V1][k][j][IBEG] - 
                    (double)BXU1BEG*g->dPoints[X1][IBEG];
                #elif X1UBEGBOUNDARY == CONST
                b->b1Stag[V_1][X1BEGBOUND][k][j] = 2.*(double)BXU1BEG - td->vCentroid[V1][k][j][IBEG];
                #endif

                #if X1VBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V_2][X1BEGBOUND][k][j] = td->vCentroid[V2][k][j][IBEG];
                #elif X1VBEGBOUNDARY == GRAD
                b->b1Stag[V_2][X1BEGBOUND][k][j] = td->vCentroid[V2][k][j][IBEG] - 
                    (double)BXV1BEG*g->dPoints[X1][IBEG];
                #elif X1VBEGBOUNDARY == CONST
                b->b1Stag[V_2][X1BEGBOUND][k][j] = 2.*(double)BXV1BEG - td->vCentroid[V2][k][j][IBEG];
                #endif

                #if X1WBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V_3][X1BEGBOUND][k][j] = td->vCentroid[V3][k][j][IBEG];
                #elif X1WBEGBOUNDARY == GRAD
                b->b1Stag[V_3][X1BEGBOUND][k][j] = td->vCentroid[V3][k][j][IBEG] - 
                    (double)BXW1BEG*g->dPoints[X1][IBEG];
                #elif X1WBEGBOUNDARY == CONST
                b->b1Stag[V_3][X1BEGBOUND][k][j] = 2.*(double)BXW1BEG - td->vCentroid[V3][k][j][IBEG];
                #endif

                #if X1UENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V_1][X1ENDBOUND][k][j] = td->vCentroid[V1][k][j][IEND];
                #elif X1UENDBOUNDARY == GRAD
                b->b1Stag[V_1][X1ENDBOUND][k][j] = td->vCentroid[V1][k][j][IEND] + 
                    (double)BXU1END*g->dPoints[X1][IEND];
                #elif X1UENDBOUNDARY == CONST
                b->b1Stag[V_1][X1ENDBOUND][k][j] = 2.*(double)BXU1END - td->vCentroid[V1][k][j][IEND];
                #endif

                #if X1VENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V_2][X1ENDBOUND][k][j] = td->vCentroid[V2][k][j][IEND];
                #elif X1VENDBOUNDARY == GRAD
                b->b1Stag[V_2][X1ENDBOUND][k][j] = td->vCentroid[V2][k][j][IEND] + 
                    (double)BXV1END*g->dPoints[X1][IEND];
                #elif X1VENDBOUNDARY == CONST
                b->b1Stag[V_2][X1ENDBOUND][k][j] = 2.*(double)BXV1END - td->vCentroid[V2][k][j][IEND];
                #endif

                #if X1WENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V_3][X1ENDBOUND][k][j] = td->vCentroid[V3][k][j][IEND];
                #elif X1WENDBOUNDARY == GRAD
                b->b1Stag[V_3][X1ENDBOUND][k][j] = td->vCentroid[V3][k][j][IEND] + 
                    (double)BXW1END*g->dPoints[X1][IEND];
                #elif X1WENDBOUNDARY == CONST
                b->b1Stag[V_3][X1ENDBOUND][k][j] = 2.*(double)BXW1END - td->vCentroid[V3][k][j][IEND];
                #endif
            #endif
        }
    }

    for (int i=0; i<=X1POINTS+1; i++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            #if STAGGERED_GRID == YES
                #if X2UBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V1S][X2BEGBOUND][i][k] = td->vStaggered[V1S][k][JBEG][i];
                #elif X2UBEGBOUNDARY == GRAD
                b->b1Stag[V1S][X2BEGBOUND][i][k] = td->vStaggered[V1S][k][JBEG][i] - 
                    (double)BXU2BEG*g->dPoints[X2][JBEG];
                #elif X2UBEGBOUNDARY == CONST
                b->b1Stag[V1S][X2BEGBOUND][i][k] = 2*(double)BXU2BEG - td->vStaggered[V1S][k][JBEG][i];
                #endif

                #if X2VBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V2S][X2BEGBOUND][i][k] = td->vStaggered[V2S][k][JBEG][i];
                #elif X2VBEGBOUNDARY == GRAD
                b->b1Stag[V2S][X2BEGBOUND][i][k] = td->vStaggered[V2S][k][JBEG][i] - 
                    (double)BXV2BEG*g->dPoints[X2][JBEG];
                #elif X2VBEGBOUNDARY == CONST
                b->b1Stag[V2S][X2BEGBOUND][i][k] = (double)BXV2BEG;
                #endif

                #if X2WBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V3S][X2BEGBOUND][i][k] = td->vStaggered[V3S][k][JBEG][i];
                #elif X2WBEGBOUNDARY == GRAD
                b->b1Stag[V3S][X2BEGBOUND][i][k] = td->vStaggered[V3S][k][JBEG][i] - 
                    (double)BXW2BEG*g->dPoints[X2][JBEG];
                #elif X2WBEGBOUNDARY == CONST
                b->b1Stag[V3S][X2BEGBOUND][i][k] = (double)BXW2BEG;
                #endif

                #if X2UENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V1S][X2ENDBOUND][i][k] = td->vStaggered[V1S][V1S][k][JEND][i];
                #elif X2UENDBOUNDARY == GRAD
                b->b1Stag[V1S][X2ENDBOUND][i][k] = td->vStaggered[V1S][k][JEND][i] + 
                    (double)BXU2END*g->dPoints[X2][JEND];
                #elif X2UENDBOUNDARY == CONST
                b->b1Stag[V1S][X2ENDBOUND][i][k] = 2.*(double)BXU2END - td->vStaggered[V1S][k][JEND][i];
                #endif

                #if X2VENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V2S][X2ENDBOUND][i][k] = td->vStaggered[V2S][k][JEND][i];
                #elif X2VENDBOUNDARY == GRAD
                b->b1Stag[V2S][X2ENDBOUND][i][k] = td->vStaggered[V2S][k][JEND][i] + 
                    (double)BXV2END*g->dPoints[X2][JEND];
                #elif X2VENDBOUNDARY == CONST
                b->b1Stag[V2S][X2ENDBOUND][i][k] = 2.*(double)BXV2END - td->vStaggered[V2S][k][JEND-1][i];
                #endif

                #if X2WENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V3S][X2ENDBOUND][i][k] = td->vStaggered[V3S][k][JEND][i];
                #elif X2WENDBOUNDARY == GRAD
                b->b1Stag[V3S][X2ENDBOUND][i][k] = td->vStaggered[V3S][k][JEND][i] + 
                    (double)BXW2END*g->dPoints[X2][JEND];
                #elif X2WENDBOUNDARY == CONST
                b->b1Stag[V3S][X2ENDBOUND][i][k] = 2.*(double)BXW2END - td->vStaggered[V3S][k][JEND-1][i];
                #endif
            #else
            #if X2UBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V_1][X2BEGBOUND][i][k] = td->vCentroid[V1][k][JBEG][i];
                #elif X2UBEGBOUNDARY == GRAD
                b->b1Stag[V_1][X2BEGBOUND][i][k] = td->vCentroid[V1][k][JBEG][i] - 
                    (double)BXU2BEG*g->dPoints[X2][JBEG];
                #elif X2UBEGBOUNDARY == CONST
                b->b1Stag[V_1][X2BEGBOUND][i][k] = 2*(double)BXU2BEG - td->vCentroid[V1][k][JBEG][i];
                #endif

                #if X2VBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V_2][X2BEGBOUND][i][k] = td->vCentroid[V2][k][JBEG][i];
                #elif X2VBEGBOUNDARY == GRAD
                b->b1Stag[V_2][X2BEGBOUND][i][k] = td->vCentroid[V2][k][JBEG][i] - 
                    (double)BXV2BEG*g->dPoints[X2][JBEG];
                #elif X2VBEGBOUNDARY == CONST
                b->b1Stag[V_2][X2BEGBOUND][i][k] = 2.*(double)BXV2BEG - td->vCentroid[V2][k][JBEG][i];
                #endif

                #if X2WBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V_3][X2BEGBOUND][i][k] = td->vCentroid[V3][k][JBEG][i];
                #elif X2WBEGBOUNDARY == GRAD
                b->b1Stag[V_3][X2BEGBOUND][i][k] = td->vCentroid[V3][k][JBEG][i] - 
                    (double)BXW2BEG*g->dPoints[X2][JBEG];
                #elif X2WBEGBOUNDARY == CONST
                b->b1Stag[V_3][X2BEGBOUND][i][k] = 2.*(double)BXW2BEG - td->vCentroid[V3][k][JBEG][i];
                #endif

                #if X2UENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V_1][X2ENDBOUND][i][k] = td->vCentroid[V1][V_1][k][JEND][i];
                #elif X2UENDBOUNDARY == GRAD
                b->b1Stag[V_1][X2ENDBOUND][i][k] = td->vCentroid[V1][k][JEND][i] + 
                    (double)BXU2END*g->dPoints[X2][JEND];
                #elif X2UENDBOUNDARY == CONST
                b->b1Stag[V_1][X2ENDBOUND][i][k] = 2.*(double)BXU2END - td->vCentroid[V1][k][JEND][i];
                #endif

                #if X2VENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V_2][X2ENDBOUND][i][k] = td->vCentroid[V2][k][JEND][i];
                #elif X2VENDBOUNDARY == GRAD
                b->b1Stag[V_2][X2ENDBOUND][i][k] = td->vCentroid[V2][k][JEND][i] + 
                    (double)BXV2END*g->dPoints[X2][JEND];
                #elif X2VENDBOUNDARY == CONST
                b->b1Stag[V_2][X2ENDBOUND][i][k] = 2.*(double)BXV2END - td->vCentroid[V2][k][JEND][i];
                #endif

                #if X2WENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V_3][X2ENDBOUND][i][k] = td->vCentroid[V3][k][JEND][i];
                #elif X2WENDBOUNDARY == GRAD
                b->b1Stag[V_3][X2ENDBOUND][i][k] = td->vCentroid[V3][k][JEND][i] + 
                    (double)BXW2END*g->dPoints[X2][JEND];
                #elif X2WENDBOUNDARY == CONST
                b->b1Stag[V_3][X2ENDBOUND][i][k] = 2.*(double)BXW2END - td->vCentroid[V3][k][JEND][i];
                #endif
            #endif
        }
    }
    
    for (int j=0; j<=X2POINTS+1; j++) {
        for (int i=0; i<=X1POINTS+1; i++) {
            #if STAGGERED_GRID == YES
                #if X3UBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V1S][X3BEGBOUND][j][i] = td->vStaggered[V1S][KBEG][j][i];
                #elif X3UBEGBOUNDARY == GRAD
                b->b1Stag[V1S][X3BEGBOUND][j][i] = td->vStaggered[V1S][KBEG][j][i] - 
                    (double)BXU3BEG*g->dPoints[X3][KBEG];
                #elif X3UBEGBOUNDARY == CONST
                b->b1Stag[V1S][X3BEGBOUND][j][i] = (double)BXU3BEG;
                #endif

                #if X3VBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V2S][X3BEGBOUND][j][i] = td->vStaggered[V2S][KBEG][j][i];
                #elif X3VBEGBOUNDARY == GRAD
                b->b1Stag[V2S][X3BEGBOUND][j][i] = td->vStaggered[V2S][KBEG][j][i] - 
                    (double)BXV3BEG*g->dPoints[X3][KBEG];
                #elif X3VBEGBOUNDARY == CONST
                b->b1Stag[V2S][X3BEGBOUND][j][i] = (double)BXV3BEG;
                #endif

                #if X3WBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V3S][X3BEGBOUND][j][i] = td->vStaggered[V3S][KBEG][j][i];
                #elif X3WBEGBOUNDARY == GRAD
                b->b1Stag[V3S][X3BEGBOUND][j][i] = td->vStaggered[V3S][KBEG][j][i] - 
                    (double)BXW3BEG*g->dPoints[X3][KBEG];
                #elif X3WBEGBOUNDARY == CONST
                b->b1Stag[V3S][X3BEGBOUND][j][i] = (double)BXW3BEG;
                #endif

                #if X3UENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V1S][X3ENDBOUND][j][i] = td->vStaggered[V1S][KEND][j][i];
                #elif X3UENDBOUNDARY == GRAD
                b->b1Stag[V1S][X3ENDBOUND][j][i] = td->vStaggered[V1S][KEND][j][i] + 
                    (double)BXU3END*g->dPoints[X3][KEND];
                #elif X3UENDBOUNDARY == CONST
                b->b1Stag[V1S][X3ENDBOUND][j][i] = 2.*(double)BXU3END - td->vStaggered[V1S][KEND-1][j][i];
                #endif

                #if X3VENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V2S][X3ENDBOUND][j][i] = td->vStaggered[V2S][KEND][j][i];
                #elif X3VENDBOUNDARY == GRAD
                b->b1Stag[V2S][X3ENDBOUND][j][i] = td->vStaggered[V2S][KEND][j][i] + 
                    (double)BXV3END*g->dPoints[X3][KEND];
                #elif X3VENDBOUNDARY == CONST
                b->b1Stag[V2S][X3ENDBOUND][j][i] = 2.*(double)BXV3END - td->vStaggered[V2S][KEND-1][j][i];
                #endif

                #if X3WENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V3S][X3ENDBOUND][j][i] = td->vStaggered[V3S][KEND][j][i];
                #elif X3WENDBOUNDARY == GRAD
                b->b1Stag[V3S][X3ENDBOUND][j][i] = td->vStaggered[V3S][KEND][j][i] + 
                    (double)BXW3END*g->dPoints[X3][KEND];
                #elif X3WENDBOUNDARY == CONST
                b->b1Stag[V3S][X3ENDBOUND][j][i] = 2.*(double)BXW3END - td->vStaggered[V3S][KEND-1][j][i];
                #endif
            #else
            #if X3UBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V_1][X3BEGBOUND][j][i] = td->vCentroid[V1][KBEG][j][i];
                #elif X3UBEGBOUNDARY == GRAD
                b->b1Stag[V_1][X3BEGBOUND][j][i] = td->vCentroid[V1][KBEG][j][i] - 
                    (double)BXU3BEG*g->dPoints[X3][KBEG];
                #elif X3UBEGBOUNDARY == CONST
                b->b1Stag[V_1][X3BEGBOUND][j][i] = 2.*(double)BXU3BEG - td->vCentroid[V1][KBEG][j][i];
                #endif

                #if X3VBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V_2][X3BEGBOUND][j][i] = td->vCentroid[V2][KBEG][j][i];
                #elif X3VBEGBOUNDARY == GRAD
                b->b1Stag[V_2][X3BEGBOUND][j][i] = td->vCentroid[V2][KBEG][j][i] - 
                    (double)BXV3BEG*g->dPoints[X3][KBEG];
                #elif X3VBEGBOUNDARY == CONST
                b->b1Stag[V_2][X3BEGBOUND][j][i] = 2.*(double)BXV3BEG - td->vCentroid[V2][KBEG][j][i];
                #endif

                #if X3WBEGBOUNDARY == ZERO_GRAD
                b->b1Stag[V_3][X3BEGBOUND][j][i] = td->vCentroid[V3][KBEG][j][i];
                #elif X3WBEGBOUNDARY == GRAD
                b->b1Stag[V_3][X3BEGBOUND][j][i] = td->vCentroid[V3][KBEG][j][i] - 
                    (double)BXW3BEG*g->dPoints[X3][KBEG];
                #elif X3WBEGBOUNDARY == CONST
                b->b1Stag[V_3][X3BEGBOUND][j][i] = 2.*(double)BXW3BEG - td->vCentroid[V3][KBEG][j][i];
                #endif

                #if X3UENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V_1][X3ENDBOUND][j][i] = td->vCentroid[V1][KEND][j][i];
                #elif X3UENDBOUNDARY == GRAD
                b->b1Stag[V_1][X3ENDBOUND][j][i] = td->vCentroid[V1][KEND][j][i] + 
                    (double)BXU3END*g->dPoints[X3][KEND];
                #elif X3UENDBOUNDARY == CONST
                b->b1Stag[V_1][X3ENDBOUND][j][i] = 2.*(double)BXU3END - td->vCentroid[V1][KEND][j][i];
                #endif

                #if X3VENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V_2][X3ENDBOUND][j][i] = td->vCentroid[V2][KEND][j][i];
                #elif X3VENDBOUNDARY == GRAD
                b->b1Stag[V_2][X3ENDBOUND][j][i] = td->vCentroid[V2][KEND][j][i] + 
                    (double)BXV3END*g->dPoints[X3][KEND];
                #elif X3VENDBOUNDARY == CONST
                b->b1Stag[V_2][X3ENDBOUND][j][i] = 2.*(double)BXV3END - td->vCentroid[V2][KEND][j][i];
                #endif

                #if X3WENDBOUNDARY == ZERO_GRAD
                b->b1Stag[V_3][X3ENDBOUND][j][i] = td->vCentroid[V3][KEND][j][i];
                #elif X3WENDBOUNDARY == GRAD
                b->b1Stag[V_3][X3ENDBOUND][j][i] = td->vCentroid[V3][KEND][j][i] + 
                    (double)BXW3END*g->dPoints[X3][KEND];
                #elif X3WENDBOUNDARY == CONST
                b->b1Stag[V_3][X3ENDBOUND][j][i] = 2.*(double)BXW3END - td->vCentroid[V3][KEND][j][i];
                #endif
            #endif
        }
    }

    for (int j=0; j<=X2POINTS+1; j++) {
        for (int i=0; i<=X1POINTS+1; i++) {
            #if STAGGERED_GRID == YES
            td->vStaggered[V1S][KBEG-1][j][i] = b->b1Stag[V1S][X3BEGBOUND][j][i];
            td->vStaggered[V1S][KEND+1][j][i] = b->b1Stag[V1S][X3ENDBOUND][j][i];
            td->vStaggered[V2S][KBEG-1][j][i] = b->b1Stag[V2S][X3BEGBOUND][j][i];
            td->vStaggered[V2S][KEND+1][j][i] = b->b1Stag[V2S][X3ENDBOUND][j][i];
            td->vStaggered[V3S][KBEG-1][j][i] = b->b1Stag[V3S][X3BEGBOUND][j][i];
            td->vStaggered[V3S][KEND+1][j][i] = b->b1Stag[V3S][X3ENDBOUND][j][i];

            td->vStaggered[V3S][KEND][j][i] = BXW3END;
            #else
            td->vCentroid[V1][KBEG-1][j][i] = b->b1Stag[V_1][X3BEGBOUND][j][i];
            td->vCentroid[V1][KEND+1][j][i] = b->b1Stag[V_1][X3ENDBOUND][j][i];
            td->vCentroid[V2][KBEG-1][j][i] = b->b1Stag[V_2][X3BEGBOUND][j][i];
            td->vCentroid[V2][KEND+1][j][i] = b->b1Stag[V_2][X3ENDBOUND][j][i];
            td->vCentroid[V3][KBEG-1][j][i] = b->b1Stag[V_3][X3BEGBOUND][j][i];
            td->vCentroid[V3][KEND+1][j][i] = b->b1Stag[V_3][X3ENDBOUND][j][i];
            #endif
        }
    }

    for (int i=0; i<=X1POINTS+1; i++) {
        for (int k=0; k<=X3POINTS+1; k++) {
            #if STAGGERED_GRID == YES
            td->vStaggered[V1S][k][JBEG-1][i] = b->b1Stag[V1S][X2BEGBOUND][i][k];
            td->vStaggered[V1S][k][JEND+1][i] = b->b1Stag[V1S][X2ENDBOUND][i][k];
            td->vStaggered[V2S][k][JBEG-1][i] = b->b1Stag[V2S][X2BEGBOUND][i][k];
            td->vStaggered[V2S][k][JEND+1][i] = b->b1Stag[V2S][X2ENDBOUND][i][k];
            td->vStaggered[V3S][k][JBEG-1][i] = b->b1Stag[V3S][X2BEGBOUND][i][k];
            td->vStaggered[V3S][k][JEND+1][i] = b->b1Stag[V3S][X2ENDBOUND][i][k];

            td->vStaggered[V2S][k][JEND][i] = BXV2END;
            #else
            td->vCentroid[V1][k][JBEG-1][i] = b->b1Stag[V_1][X2BEGBOUND][i][k];
            td->vCentroid[V1][k][JEND+1][i] = b->b1Stag[V_1][X2ENDBOUND][i][k];
            td->vCentroid[V2][k][JBEG-1][i] = b->b1Stag[V_2][X2BEGBOUND][i][k];
            td->vCentroid[V2][k][JEND+1][i] = b->b1Stag[V_2][X2ENDBOUND][i][k];
            td->vCentroid[V3][k][JBEG-1][i] = b->b1Stag[V_3][X2BEGBOUND][i][k];
            td->vCentroid[V3][k][JEND+1][i] = b->b1Stag[V_3][X2ENDBOUND][i][k];
            #endif
        }
    }

    for (int k=0; k<=X3POINTS+1; k++) {
        for (int j=0; j<=X2POINTS+1; j++) {
            #if STAGGERED_GRID == YES
            td->vStaggered[V1S][k][j][IBEG-1] = b->b1Stag[V1S][X1BEGBOUND][k][j];
            td->vStaggered[V1S][k][j][IEND+1] = b->b1Stag[V1S][X1ENDBOUND][k][j];
            td->vStaggered[V2S][k][j][IBEG-1] = b->b1Stag[V2S][X1BEGBOUND][k][j];
            td->vStaggered[V2S][k][j][IEND+1] = b->b1Stag[V2S][X1ENDBOUND][k][j];
            td->vStaggered[V3S][k][j][IBEG-1] = b->b1Stag[V3S][X1BEGBOUND][k][j];
            td->vStaggered[V3S][k][j][IEND+1] = b->b1Stag[V3S][X1ENDBOUND][k][j];

            td->vStaggered[V1S][k][j][IEND] = BXU1END;
            #else
            td->vCentroid[V1][k][j][IBEG-1] = b->b1Stag[V_1][X1BEGBOUND][k][j];
            td->vCentroid[V1][k][j][IEND+1] = b->b1Stag[V_1][X1ENDBOUND][k][j];
            td->vCentroid[V2][k][j][IBEG-1] = b->b1Stag[V_2][X1BEGBOUND][k][j];
            td->vCentroid[V2][k][j][IEND+1] = b->b1Stag[V_2][X1ENDBOUND][k][j];
            td->vCentroid[V3][k][j][IBEG-1] = b->b1Stag[V_3][X1BEGBOUND][k][j];
            td->vCentroid[V3][k][j][IEND+1] = b->b1Stag[V_3][X1ENDBOUND][k][j];
            #endif
        }
    }
}