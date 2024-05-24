#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdint.h>

#include "struct.h"

double* allocate_array_1D(double *, int);
double** allocate_array_2D(double **, int, int *);
double** allocate_array_uniform_2D(double **, int, int);
double*** allocate_array_3D(double ***, int, int *, int *);
double**** allocate_array_4D(double ****, int, int *, int *, int *);
double**** allocate_array_semi_uniform_4D(double ****, int, int, int *, int *);
double**** allocate_array_uniform_4D(double ****, int, int, int, int);
void generate_uniform_array(Grid*, int);
void generate_logminus_array(Grid*, int);

double* allocate_array_1D(double* arr, int a){
    
    arr = (double*)calloc(a, sizeof(double));
    return arr;
}

double** allocate_array_2D(double** arr, int a, int* b){
    
    arr = (double**)calloc(a, sizeof(double*));
    for (int i=0; i<a; i++){
        arr[i] = (double*)calloc(b[i], sizeof(double));
    }
    return arr;
}

double** allocate_array_uniform_2D(double** arr, int a, int b){
    
    arr = (double**)calloc(a, sizeof(double*));
    for (int i=0; i<a; i++){
        arr[i] = (double*)calloc(b, sizeof(double));
    }
    return arr;
}

double*** allocate_array_3D(double*** arr, int a, int* b, int* c){
    
    arr = (double***)calloc(a, sizeof(double**));
    for (int i=0; i<a; i++){
        arr[i] = (double**)calloc(b[i], sizeof(double*));
        for (int j=0; j<b[i]; j++){
            arr[i][j] = (double*)calloc(c[j], sizeof(double));
        }
    }
    return arr;
}

double**** allocate_array_4D(double**** arr, int a, int* b, int* c, int* d){
    
    arr = (double****)calloc(a, sizeof(double***));
    for (int i=0; i<a; i++){
        arr[i] = (double***)calloc(b[i], sizeof(double**));
        for (int j=0; j<b[i]; j++){
            arr[i][j] = (double**)calloc(c[j], sizeof(double*));
            for (int k=0; k<c[j]; k++){
                arr[i][j][k] = (double*)calloc(d[k], sizeof(double*));
            }
        }
    }
    return arr;
}

double**** allocate_array_semi_uniform_4D(double**** arr, int a, int b, int* c, int* d) {
    
    arr = (double****)calloc(a, sizeof(double***));
    for (int i=0; i<a; i++) {
        arr[i] = (double***)calloc(b, sizeof(double**));
        for (int j=0; j<b; j++) {
            arr[i][j] = (double**)calloc(c[j], sizeof(double*));
            for (int k=0; k<c[j]; k++) {
                arr[i][j][k] = (double*)calloc(d[j], sizeof(double*));
            }
        }
    }
    return arr;
}

double**** allocate_array_uniform_4D(double**** arr, int a, int b, int c, int d) {
    arr = (double****)calloc(a, sizeof(double***));
    for (int i=0; i<a; i++){
        arr[i] = (double***)calloc(b, sizeof(double**));
        for (int j=0; j<b; j++){
            arr[i][j] = (double**)calloc(c, sizeof(double*));
            for (int k=0; k<c; k++){
                arr[i][j][k] = (double*)calloc(d, sizeof(double*));
            }
        }
    }
    return arr;
}

void generate_uniform_array(Grid* g, int index){
    int N;
    double xL, xR;
    xL = g->xBeg[index];
    xR = g->xEnd[index];
    N = g->npTotal[index];
    for (int i=0; i<=g->npTotal[index]; i++){
            g->xPoints[index][i] = xL+(double)i*(xR-xL)/(double)N;
        }
        printf("Array generated along x%d\n\n", index+1);
}

void generate_logplus_array(Grid* g, int index){
    double dx, de, x, xL, xR;
        int N, counter;

        x = g->xBeg[index];
        xL = g->xBeg[index];
        xR = g->xEnd[index];
        dx = 0.0;
        N = (double)g->npTotal[index]+1;
        counter = 0;

        while(counter<=N){
            g->xPoints[index][counter] = x+dx;
            counter++;
            de = (1./N)*log10((xR+abs(xL)-xL)/abs(xL));
            dx = (x+abs(xL)-xL)*(pow(10.,de)-1.);
            x = x+dx;
        }
        printf("Array generated along x%d\n\n", index+1);
}

void generate_logminus_array(Grid* g, int index){
    double dx, de, x, xL, xR;
        int N, counter;

        x = g->xBeg[index];
        xL = g->xBeg[index];
        xR = g->xEnd[index];
        dx = 0.0;
        N = (double)g->npTotal[index]+1;
        counter = 0;

        while(counter<=N){
            g->xPoints[index][counter] = x+dx;
            counter++;
            de = -(1./N)*log10((xR+abs(xL)-xL)/abs(xL));
            dx = (x-abs(xL)-xR)*(pow(10.,de)-1.);
            x = x+dx;
        }
        printf("Array generated along x%d\n\n", index+1);
}
