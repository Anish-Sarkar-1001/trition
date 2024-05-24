#include <stdio.h>
#include <math.h>
#include <string.h>

#include "struct.h"
#include "globals.h"
#include "code.h"
#include "main.h"

void write_output(Data *, Grid *, Boundary*, tData* td);
void data_output(Data *, tData *);
void grid_output(Grid *);
void boundary_output(Boundary *, Grid *);
void tostring(char [], int);

void write_output(Data* d, Grid* g, Boundary* b, tData* td) {
    data_output(d, td);
    grid_output(g);
    //boundary_output(b, g);
}

void data_output(Data* d, tData* td) {
    FILE* f1 = fopen("Data_v.txt", "w");
    fprintf(f1, "#\tData\n\n");
    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                fprintf(f1, "%g\t", td->vCentroid[V1][k][j][i]);
            } fprintf(f1, "\n");
        } fprintf(f1, "\n\n");
    }

    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                fprintf(f1, "%g\t", td->vCentroid[V2][k][j][i]);
            } fprintf(f1, "\n");
        } fprintf(f1, "\n\n");
    }

    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                fprintf(f1, "%g\t", td->vCentroid[V3][k][j][i]);
            } fprintf(f1, "\n");
        } fprintf(f1, "\n\n");
    }
    /*for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG-1; j<=JEND; j++) {
            for (int i=IBEG-1; i<=IEND; i++) {
                fprintf(f, "%g\t", td->vCentroid[PRS][k][j][i]);
            } fprintf(f, "\n");
        } fprintf(f, "\n\n");
    }*/
    fclose(f1);
    FILE* f = fopen("Data_p.txt", "w");
    fprintf(f, "#\tData\n\n");

    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                fprintf(f, "%g\t", td->vCentroid[PRS][k][j][i]);
            } fprintf(f, "\n");
        } fprintf(f, "\n\n");
    }
    fclose(f);
}

void grid_output(Grid* g) {
    FILE* f1 = fopen("Grid_v.txt", "w");
    #if PHYSICS == STEADY_HEAT_CONDUCTION
    /*fprintf(f, "#\tComputational Grid Information\n\n");
    fprintf(f, "\n#\tX1-Centroids\t\tX2-Centroids\t\tX3-Centroids\n");*/

    for (int i=0; i<MAX(MAX(X1POINTS, X2POINTS), X3POINTS); i++) {
            (i==0) ? fprintf(f, "%i", i) : fprintf(f, "\t%i", i);
        } fprintf(f, "\n");

    for (int i=0; i<=2; i++) {
        for (int j=0; j<g->npTotal[i]; j++) {
            (j==0) ? fprintf(f, "%g", g->xCentroids[i][j]) : fprintf(f, "\t%g", g->xCentroids[i][j]);
        } fprintf(f, "\n");
    }
    #elif PHYSICS == UNSTEADY_FLUID_FLOW

    for (int i=IBEG-1; i<=IEND-1; i++) {
        fprintf(f1, "%g\t", g->xCentroids[X1][i]);
    }
    fprintf(f1, "\n");

    for (int i=JBEG-1; i<=JEND-1; i++) {
        fprintf(f1, "%g\t", g->xCentroids[X2][i]);
    }
    #endif
    fclose(f1);
    FILE* f2 = fopen("Grid_p.txt", "w");
    #if PHYSICS == STEADY_HEAT_CONDUCTION
    /*fprintf(f, "#\tComputational Grid Information\n\n");
    fprintf(f, "\n#\tX1-Centroids\t\tX2-Centroids\t\tX3-Centroids\n");*/

    for (int i=0; i<MAX(MAX(X1POINTS, X2POINTS), X3POINTS); i++) {
            (i==0) ? fprintf(f, "%i", i) : fprintf(f, "\t%i", i);
        } fprintf(f, "\n");

    for (int i=0; i<=2; i++) {
        for (int j=0; j<g->npTotal[i]; j++) {
            (j==0) ? fprintf(f, "%g", g->xCentroids[i][j]) : fprintf(f, "\t%g", g->xCentroids[i][j]);
        } fprintf(f, "\n");
    }
    #elif PHYSICS == UNSTEADY_FLUID_FLOW

    for (int i=IBEG-1; i<=IEND-1; i++) {
        fprintf(f2, "%g\t", g->xCentroids[X1][i]);
    }
    fprintf(f2, "\n");

    for (int i=JBEG-1; i<=JEND-1; i++) {
        fprintf(f2, "%g\t", g->xCentroids[X2][i]);
    }
    #endif
    fclose(f2);
}

void boundary_output(Boundary* b, Grid* g) {
    #if PHYSICS == STEADY_HEAT_CONDUCTION
    for (int i=0; i<6; i++) {
        int check = (int)(floor(i/2))%3;

        FILE* f;
        char fileName[50] = "Boundary_";
        char index[1];
        char text[] = ".txt";
        tostring(index, check+1);
        strcat(fileName, index);
        strcat(fileName, text);
        printf("%s\t", fileName);
        if (i%2 == 0) {
            f = fopen(fileName, "w");
        }
        else {
            f = fopen(fileName, "a");
        }
        fprintf(f, "#\tExternal Boundary Information\n\n");
        fprintf(f, "\n#\tBEG\t\tEND\n\n");

        int outer_array = check-1<0 ? check+2 : check-1;
        int inner_array = check+1>=3 ? check-2 : check+1;
        for (int j=0; j<=g->npTotal[outer_array]+1; j++) {
            for (int k=0; k<=g->npTotal[inner_array]+1; k++) {
                fprintf(f, "%g\t", b->bCentroids[TMP][i][j][k]);
            } fprintf(f, "\n");
        } fclose(f);
    }
    #endif
}

void tostring(char str[], int num) {
    int i, rem, len = 0, n;
    n = num;
    while (n != 0) {
        len++;
        n /= 10;
    }
    for (i = 0; i < len; i++) {
        rem = num % 10;
        num = num / 10;
        str[len - (i + 1)] = rem + '0';
    }
    str[len] = '\0';
}
