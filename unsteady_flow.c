#include <math.h>
#include <stdio.h>

#include "main.h"
#include "code.h"
#include "struct.h"
#include "globals.h"

#if STAGGERED_GRID == YES
void calculate_flux(Data *, Grid *, tData *, midT *, Boundary *);
double advection_flux_prev_vsx1(int, int, int, Data*, Grid*);
double advection_flux_prev_vsx2(int, int, int, Data*, Grid*);
double advection_flux_prev_vsx3(int, int, int, Data*, Grid*);
double diffusion_flux_prev_vxs1(int, int, int, Data*, Grid*);
double diffusion_flux_prev_vxs2(int, int, int, Data*, Grid*);
double diffusion_flux_prev_vxs3(int, int, int, Data*, Grid*);
void pressure_poisson(Data*, Grid*, tData*, midT*);
void max_pressure(tData*);

#if TIMESTEPPING == EXPLICIT
double advection_flux_prev_vsx1(int i, int j, int k, Data* d, Grid* g) {
    double u_e, u_w, u_n, u_s, u_f, u_b;
    double v_s, v_n, w_f, w_b;
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double flux;

    S_e = g->v1StagArea[k][j][i][E];
    S_w = g->v1StagArea[k][j][i][W];
    S_n = g->v1StagArea[k][j][i][N];
    S_s = g->v1StagArea[k][j][i][S];
    S_f = g->v1StagArea[k][j][i][F];
    S_b = g->v1StagArea[k][j][i][S];

    u_e = (vel[V1S][k][j][i] + vel[V1S][k][j][i+1])/2.;
    u_w = (vel[V1S][k][j][i] + vel[V1S][k][j][i-1])/2.;
    u_n = (vel[V1S][k][j][i]/g->dPoints[X2][j] + vel[V1S][k][j+1][i]/g->dPoints[X2][j+1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j+1]);
    u_s = (vel[V1S][k][j][i]/g->dPoints[X2][j] + vel[V1S][k][j-1][i]/g->dPoints[X2][j-1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j-1]);
    u_f = (vel[V1S][k][j][i]/g->dPoints[X3][k] + vel[V1S][k+1][j][i]/g->dPoints[X3][k+1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k+1]);
    u_b = (vel[V1S][k][j][i]/g->dPoints[X3][k] + vel[V1S][k-1][j][i]/g->dPoints[X3][k-1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k-1]);

    v_n = (vel[V2S][k][j][i]/g->dPoints[X1][i] + vel[V2S][k][j][i+1]/g->dPoints[X1][i+1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i+1]);
    v_s = (vel[V2S][k][j-1][i]/g->dPoints[X1][i] + vel[V2S][k][j-1][i+1]/g->dPoints[X1][i+1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i+1]);
    w_f = (vel[V3S][k][j][i]/g->dPoints[X1][i] + vel[V3S][k][j][i+1]/g->dPoints[X1][i+1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i+1]);
    w_b = (vel[V3S][k-1][j][i]/g->dPoints[X1][i] + vel[V3S][k-1][j][i+1]/g->dPoints[X1][i+1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i+1]);

    flux = u_e*u_e*S_e - u_w*u_w*S_w + u_n*v_n*S_n - u_s*v_s*S_s + u_f*w_f*S_f - u_b*w_b*S_b;

    return flux;
}

double advection_flux_prev_vsx2(int i, int j, int k, Data* d, Grid* g) {
    double v_e, v_w, v_n, v_s, v_f, v_b;
    double u_e, u_w, w_f, w_b;
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double flux;

    S_e = g->v2StagArea[k][j][i][E];
    S_w = g->v2StagArea[k][j][i][W];
    S_n = g->v2StagArea[k][j][i][N];
    S_s = g->v2StagArea[k][j][i][S];
    S_f = g->v2StagArea[k][j][i][F];
    S_b = g->v2StagArea[k][j][i][B];

    v_e = (vel[V2S][k][j][i]/g->dPoints[X1][i] + vel[V2S][k][j][i+1]/g->dPoints[X1][i+1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i+1]);
    v_w = (vel[V2S][k][j][i]/g->dPoints[X1][i] + vel[V2S][k][j][i-1]/g->dPoints[X1][i-1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i-1]);
    v_n = (vel[V2S][k][j][i] + vel[V2S][k][j+1][i])/2.;
    v_s = (vel[V2S][k][j][i] + vel[V2S][k][j-1][i])/2.;
    v_f = (vel[V2S][k][j][i]/g->dPoints[X3][k] + vel[V2S][k+1][j][i]/g->dPoints[X3][k+1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k+1]);
    v_b = (vel[V2S][k][j][i]/g->dPoints[X3][k] + vel[V2S][k-1][j][i]/g->dPoints[X3][k-1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k-1]);

    u_e = (vel[V1S][k][j][i]/g->dPoints[X2][j] + vel[V1S][k][j+1][i]/g->dPoints[X2][j+1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j+1]);
    u_w = (vel[V1S][k][j][i-1]/g->dPoints[X2][j] + vel[V1S][k][j+1][i-1]/g->dPoints[X2][j+1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j+1]);
    w_f = (vel[V3S][k][j][i]/g->dPoints[X2][j] + vel[V3S][k][j+1][i]/g->dPoints[X2][j+1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j+1]);
    w_b = (vel[V3S][k-1][j][i]/g->dPoints[X2][j] + vel[V3S][k-1][j+1][i]/g->dPoints[X2][j+1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j+1]);

    flux = u_e*v_e*S_e - u_w*v_w*S_w + v_n*v_n*S_n - v_s*v_s*S_s + v_f*w_f*S_f - v_b*w_b*S_b;

    return flux;
}

double advection_flux_prev_vsx3(int i, int j, int k, Data* d, Grid* g) {
    double w_e, w_w, w_n, w_s, w_f, w_b;
    double u_e, u_w, v_n, v_s;
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double flux;

    S_e = g->v3StagArea[k][j][i][E];
    S_w = g->v3StagArea[k][j][i][W];
    S_n = g->v3StagArea[k][j][i][N];
    S_s = g->v3StagArea[k][j][i][S];
    S_f = g->v3StagArea[k][j][i][F];
    S_b = g->v3StagArea[k][j][i][B];

    w_e = (vel[V3S][k][j][i]/g->dPoints[X1][i] + vel[V3S][k][j][i+1]/g->dPoints[X1][i+1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i+1]);
    w_w = (vel[V3S][k][j][i]/g->dPoints[X1][i] + vel[V3S][k][j][i-1]/g->dPoints[X1][i-1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i-1]);
    w_n = (vel[V3S][k][j][i]/g->dPoints[X2][j] + vel[V3S][k][j+1][i]/g->dPoints[X2][j+1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j+1]);
    w_s = (vel[V3S][k][j][i]/g->dPoints[X2][j] + vel[V3S][k][j-1][i]/g->dPoints[X2][j-1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j-1]);
    w_f = (vel[V3S][k][j][i] + vel[V3S][k+1][j][i])/2.;
    w_b = (vel[V3S][k][j][i] + vel[V3S][k-1][j][i])/2.;

    u_e = (vel[V1S][k][j][i]/g->dPoints[X3][k] + vel[V1S][k+1][j][i]/g->dPoints[X3][k+1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k+1]);
    u_w = (vel[V1S][k][j][i-1]/g->dPoints[X3][k] + vel[V1S][k+1][j][i-1]/g->dPoints[X3][k+1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k+1]);
    v_n = (vel[V2S][k][j][i]/g->dPoints[X3][k] + vel[V2S][k+1][j][i]/g->dPoints[X3][k+1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k+1]);
    v_s = (vel[V2S][k][j-1][i]/g->dPoints[X3][k] + vel[V2S][k+1][j-1][i]/g->dPoints[X3][k+1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k+1]);

    flux = u_e*w_e*S_e - u_w*w_w*S_w + w_n*v_n*S_n - w_s*v_s*S_s + w_f*w_f*S_f - w_b*w_b*S_b;

    return flux;
}

double diffusion_flux_prev_vxs1(int i, int j, int k, Data* d, Grid* g) {
    double g_e, g_w, g_n, g_s, g_f, g_b;
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double flux;

    S_e = g->v1StagArea[k][j][i][E];
    S_w = g->v1StagArea[k][j][i][W];
    S_n = g->v1StagArea[k][j][i][N];
    S_s = g->v1StagArea[k][j][i][S];
    S_f = g->v1StagArea[k][j][i][F];
    S_b = g->v1StagArea[k][j][i][S];

    g_e = (vel[V1S][k][j][i+1] - vel[V1S][k][j][i])/g->dPoints[X1][i+1];
    g_w = (vel[V1S][k][j][i] - vel[V1S][k][j][i-1])/g->dPoints[X1][i];
    g_n = 2*(vel[V1S][k][j+1][i] - vel[V1S][k][j][i])/
        (g->dPoints[X2][j+1] + g->dPoints[X2][j]);
    g_s = 2*(vel[V1S][k][j][i] - vel[V1S][k][j-1][i])/
        (g->dPoints[X2][j] + g->dPoints[X2][j-1]);
    g_f = 2.*(vel[V1S][k+1][j][i] - vel[V1S][k][j][i])/
        (g->dPoints[X3][k+1] + g->dPoints[X3][k]);
    g_b = 2.*(vel[V1S][k][j][i] - vel[V1S][k-1][j][i])/
        (g->dPoints[X3][k] + g->dPoints[X3][k-1]);

    flux = (g_e*S_e - g_w*S_w + g_n*S_n - g_s*S_s + g_f*S_f - g_b*S_b)/g_Re;

    return flux;
}

double diffusion_flux_prev_vxs2(int i, int j, int k, Data* d, Grid* g) {
    double g_e, g_w, g_n, g_s, g_f, g_b;
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double flux;

    S_e = g->v2StagArea[k][j][i][E];
    S_w = g->v2StagArea[k][j][i][W];
    S_n = g->v2StagArea[k][j][i][N];
    S_s = g->v2StagArea[k][j][i][S];
    S_f = g->v2StagArea[k][j][i][F];
    S_b = g->v2StagArea[k][j][i][B];

    g_e = 2.*(vel[V2S][k][j][i+1] - vel[V2S][k][j][i])/
        (g->dPoints[X1][i+1] + g->dPoints[X1][i]);
    g_w = 2.*(vel[V2S][k][j][i] - vel[V2S][k][j][i-1])/
        (g->dPoints[X1][i] + g->dPoints[X1][i-1]);
    g_n = (vel[V2S][k][j+1][i] - vel[V2S][k][j][i])/g->dPoints[X2][j+1];
    g_s = (vel[V2S][k][j][i] - vel[V2S][k][j-1][i])/g->dPoints[X2][j];
    g_f = 2.*(vel[V2S][k+1][j][i] - vel[V2S][k][j][i])/
        (g->dPoints[X3][k+1] + g->dPoints[X3][k]);
    g_b = 2.*(vel[V2S][k][j][i] - vel[V2S][k-1][j][i])/
        (g->dPoints[X3][k] + g->dPoints[X3][k-1]);

    flux = (g_e*S_e - g_w*S_w + g_n*S_n - g_s*S_s + g_f*S_f - g_b*S_b)/g_Re;

    return flux;
}

double diffusion_flux_prev_vxs3(int i, int j, int k, Data* d, Grid* g) {
    double g_e, g_w, g_n, g_s, g_f, g_b;
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double flux;

    S_e = g->v3StagArea[k][j][i][E];
    S_w = g->v3StagArea[k][j][i][W];
    S_n = g->v3StagArea[k][j][i][N];
    S_s = g->v3StagArea[k][j][i][S];
    S_f = g->v3StagArea[k][j][i][F];
    S_b = g->v3StagArea[k][j][i][B];

    g_e = 2.*(vel[V3S][k][j][i+1] - vel[V3S][k][j][i])/
        (g->dPoints[X1][i+1] + g->dPoints[X1][i]);
    g_w = 2.*(vel[V3S][k][j][i] - vel[V3S][k][j][i-1])/
        (g->dPoints[X1][i] + g->dPoints[X1][i-1]);
    g_n = 2.*(vel[V3S][k][j+1][i] - vel[V3S][k][j][i])/
        (g->dPoints[X2][j+1] + g->dPoints[X2][j]);
    g_s = 2.*(vel[V3S][k][j][i] - vel[V3S][k][j-1][i])/
        (g->dPoints[X2][j] + g->dPoints[X2][j-1]);
    g_f = (vel[V3S][k+1][j][i] - vel[V3S][k][j][i])/g->dPoints[X3][k+1];
    g_b = (vel[V3S][k][j][i] - vel[V3S][k-1][j][i])/g->dPoints[X3][k];

    flux = (g_e*S_e - g_w*S_w + g_n*S_n - g_s*S_s + g_f*S_f - g_b*S_b)/g_Re;

    return flux;
}

void calculate_flux(Data* d, Grid* g, tData* td, midT* mt, Boundary* b) {
    
    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND-1; i++) {
                mt->vS[V1S][k][j][i] = 2.*(diffusion_flux_prev_vxs1(i,j,k,d,g) - advection_flux_prev_vsx1(i,j,k,d,g))*
                    g_dt/((g->dPoints[X1][i]+g->dPoints[X1][i+1])*g->dPoints[X2][j]*g->dPoints[X3][k])+
                    vel[V1S][k][j][i];
            }
        }
    }
    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND-1; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                mt->vS[V2S][k][j][i] = 2.*(diffusion_flux_prev_vxs2(i,j,k,d,g) - advection_flux_prev_vsx2(i,j,k,d,g))*
                    g_dt/((g->dPoints[X2][j]+g->dPoints[X2][j+1])*g->dPoints[X1][i]*g->dPoints[X3][k])+
                    vel[V2S][k][j][i];
            }
        }
    }
    for (int k=KBEG; k<=KEND-1; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                mt->vS[V3S][k][j][i] = 2.*(diffusion_flux_prev_vxs3(i,j,k,d,g) - advection_flux_prev_vsx3(i,j,k,d,g))*
                    g_dt/((g->dPoints[X3][k]+g->dPoints[X3][k+1])*g->dPoints[X2][j]*g->dPoints[X1][i])+
                    vel[V3S][k][j][i];
            }
        }
    }
    update_fractional_velocity(b, d, g, td, mt);
    pressure_poisson(d, g, td, mt);
    update_pressure(b, d, g, td, mt);
    max_pressure(td);
    double prs_prev = 0.0;
    double prs_new = td->vCentroid[PRS][td->kPos][td->jPos][td->iPos];

    for (int k=KBEG; k<=KEND; k++) {
            for (int j=JBEG; j<=JEND; j++) {
                for (int i=IBEG; i<=IEND; i++) {
                    td->vCentroid[PRS][k][j][i] = 0.0;
                }
            }
        }
    
    while (fabs(fabs(prs_new)-fabs(prs_prev)) > TOLLERANCE) {
        
        pressure_poisson(d, g, td, mt);
        update_pressure(b, d, g, td, mt);
        prs_prev = prs_new;
        prs_new = td->vCentroid[PRS][td->kPos][td->jPos][td->iPos];
    }

    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND-1; i++) {
                tvel[V1S][k][j][i] = mt->vS[V1S][k][j][i] - g_dt*2.*
                    (td->vCentroid[PRS][k][j][i+1] - td->vCentroid[PRS][k][j][i])/
                    (g->dPoints[X1][i] + g->dPoints[X1][i+1]);
            }
        }
    }
    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND-1; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                tvel[V2S][k][j][i] = mt->vS[V2S][k][j][i] - g_dt*2.*
                    (td->vCentroid[PRS][k][j+1][i] - td->vCentroid[PRS][k][j][i])/
                    (g->dPoints[X2][j] + g->dPoints[X2][j+1]);
            }
        }
    }
    for (int k=KBEG; k<=KEND-1; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                tvel[V3S][k][j][i] = mt->vS[V3S][k][j][i] - g_dt*2.*
                    (td->vCentroid[PRS][k+1][j][i] - td->vCentroid[PRS][k][j][i])/
                    (g->dPoints[X3][k] + g->dPoints[X3][k+1]);
            }
        }
    }
    update_next_velocity(b, d, g, td, mt);
}

void pressure_poisson(Data* d, Grid* g, tData* td, midT* mt) {
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double vel;
    double den, num;
    double dx1, dx2, dy1, dy2, dz1, dz2;

    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                S_e = g->dPoints[X2][j]*g->dPoints[X3][k];
                S_n = g->dPoints[X1][i]*g->dPoints[X3][k];
                S_f = g->dPoints[X1][i]*g->dPoints[X2][j];
                S_w = S_e;
                S_s = S_n;
                S_b = S_f;

                dx1 = (g->dPoints[X1][i] + g->dPoints[X1][i+1])/2.;
                dx2 = (g->dPoints[X1][i] + g->dPoints[X1][i-1])/2.;
                dy1 = (g->dPoints[X2][j] + g->dPoints[X2][j+1])/2.;
                dy2 = (g->dPoints[X2][j] + g->dPoints[X2][j-1])/2.;
                dz1 = (g->dPoints[X3][k] + g->dPoints[X3][k+1])/2.;
                dz2 = (g->dPoints[X3][k] + g->dPoints[X3][k-1])/2.;

                vel = mt->vS[V1S][k][j][i]*S_e - mt->vS[V1S][k][j][i-1]*S_w+
                    mt->vS[V2S][k][j][i]*S_n - mt->vS[V2S][k][j-1][i]*S_s+
                    mt->vS[V3S][k][j][i]*S_f - mt->vS[V3S][k-1][j][i]*S_b;

                den = S_e/dx1 + S_w/dx2 + S_n/dy1 + S_s/dy2 + S_f/dz1 + S_b/dz2;

                num = S_e*td->vCentroid[PRS][k][j][i+1]/dx1 + S_w*td->vCentroid[PRS][k][j][i-1]/dx2+
                    S_n*td->vCentroid[PRS][k][j+1][i]/dy1 + S_s*td->vCentroid[PRS][k][j-1][i]/dy2+
                    S_f*td->vCentroid[PRS][k+1][j][i]/dz1 + S_b*td->vCentroid[PRS][k-1][j][i]/dz2 - vel/g_dt;
                
                td->vCentroid[PRS][k][j][i] = num/den;
            }
        }
    }
}

void max_pressure(tData* td) {
    double max_P = 0.0;
    int index[3];
    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                max_P = MAX(td->vCentroid[PRS][k][j][i], max_P);
                if (max_P == td->vCentroid[PRS][k][j][i]) {
                    td->iPos = i;
                    td->jPos = j;
                    td->kPos = k;
                }
            }
        }
    }
}
#endif

#else
void calculate_flux(Data *, Grid *, tData *, midT *, Boundary *);
double advection_flux_prev_vsx1(int, int, int, Data*, Grid*, tData*, midT*);
double advection_flux_prev_vsx2(int, int, int, Data*, Grid*, tData*, midT*);
double advection_flux_prev_vsx3(int, int, int, Data*, Grid*, tData*, midT*);
double diffusion_flux_prev_vxs1(int, int, int, Data*, Grid*);
double diffusion_flux_prev_vxs2(int, int, int, Data*, Grid*);
double diffusion_flux_prev_vxs3(int, int, int, Data*, Grid*);
void pressure_poisson(Data*, Grid*, tData*, midT*);
void max_pressure(tData*);

#if TIMESTEPPING == EXPLICIT
double advection_flux_prev_vsx1(int i, int j, int k, Data* d, Grid* g, tData* td, midT* mt) {
    double u_ep, u_wp, u_np, u_sp, u_fp, u_bp;
    double u_e, u_w, v_s, v_n, w_f, w_b;
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double flux;

    S_e = g->Area[k][j][i][E];
    S_w = g->Area[k][j][i][W];
    S_n = g->Area[k][j][i][N];
    S_s = g->Area[k][j][i][S];
    S_f = g->Area[k][j][i][F];
    S_b = g->Area[k][j][i][S];

    u_ep = (td->vCentroid[V1][k][j][i]/g->dPoints[X1][i] + td->vCentroid[V1][k][j][i+1]/g->dPoints[X1][i+1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i+1]);
    u_wp = (td->vCentroid[V1][k][j][i]/g->dPoints[X1][i] + td->vCentroid[V1][k][j][i-1]/g->dPoints[X1][i-1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i-1]);
    u_np = (td->vCentroid[V1][k][j][i]/g->dPoints[X2][j] + td->vCentroid[V1][k][j+1][i]/g->dPoints[X2][j+1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j+1]);
    u_sp = (td->vCentroid[V1][k][j][i]/g->dPoints[X2][j] + td->vCentroid[V1][k][j-1][i]/g->dPoints[X2][j-1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j-1]);
    u_fp = (td->vCentroid[V1][k][j][i]/g->dPoints[X3][k] + td->vCentroid[V1][k+1][j][i]/g->dPoints[X3][k+1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k+1]);
    u_bp = (td->vCentroid[V1][k][j][i]/g->dPoints[X3][k] + td->vCentroid[V1][k-1][j][i]/g->dPoints[X3][k-1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k-1]);

    /*u_ep = (td->vCentroid[V1][k][j][i] + td->vCentroid[V1][k][j][i+1])/2.;
    u_wp = (td->vCentroid[V1][k][j][i] + td->vCentroid[V1][k][j][i-1])/2.;
    u_np = (td->vCentroid[V1][k][j][i] + td->vCentroid[V1][k][j+1][i])/2.;
    u_sp = (td->vCentroid[V1][k][j][i] + td->vCentroid[V1][k][j-1][i])/2.;
    u_fp = (td->vCentroid[V1][k][j][i] + td->vCentroid[V1][k+1][j][i])/2.;
    u_bp = (td->vCentroid[V1][k][j][i] + td->vCentroid[V1][k-1][j][i])/2.;*/

    /*u_e = (mt->vC[V_1][k][j][i] + mt->vC[V_1][k][j][i+1])/2. - 
        g_dt*(td->vCentroid[PRS][k][j][i+1] - td->vCentroid[PRS][k][j][i])/g->dPoints[X1][i];
    u_w = (mt->vC[V_1][k][j][i] + mt->vC[V_1][k][j][i-1])/2. -
        g_dt*(td->vCentroid[PRS][k][j][i] - td->vCentroid[PRS][k][j][i-1])/g->dPoints[X1][i];
    v_n = (mt->vC[V_2][k][j][i] + mt->vC[V_2][k][j+1][i])/2. -
        g_dt*(td->vCentroid[PRS][k][j+1][i] - td->vCentroid[PRS][k][j][i])/g->dPoints[X2][j];
    v_s = (mt->vC[V_2][k][j][i] + mt->vC[V_2][k][j-1][i])/2. - 
        g_dt*(td->vCentroid[PRS][k][j][i] - td->vCentroid[PRS][k][j-1][i])/g->dPoints[X2][j];
    w_f = 0.0;
    w_b = 0.0;*/

    u_e = (mt->vC[V_1][k][j][i]/g->dPoints[X1][i] + mt->vC[V_1][k][j][i+1]/g->dPoints[X1][i+1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i+1]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j][i+1] - td->vCentroid[PRS][k][j][i])/
        (g->dPoints[X1][i] + g->dPoints[X1][i+1]);
    u_w = (mt->vC[V_1][k][j][i]/g->dPoints[X1][i] + mt->vC[V_1][k][j][i-1]/g->dPoints[X1][i-1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i-1]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j][i] - td->vCentroid[PRS][k][j][i-1])/
        (g->dPoints[X1][i] + g->dPoints[X1][i-1]);
    v_n = (mt->vC[V_2][k][j][i]/g->dPoints[X2][j] + mt->vC[V_2][k][j+1][i]/g->dPoints[X2][j+1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j+1]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j+1][i] - td->vCentroid[PRS][k][j][i])/
        (g->dPoints[X2][j] + g->dPoints[X2][j+1]);
    v_s = (mt->vC[V_2][k][j-1][i]/g->dPoints[X2][j-1] + mt->vC[V_2][k][j][i]/g->dPoints[X2][j])/
        (1./g->dPoints[X2][j-1] + 1./g->dPoints[X2][j]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j][i] - td->vCentroid[PRS][k][j-1][i])/
        (g->dPoints[X2][j] + g->dPoints[X2][j-1]);
    w_f = (mt->vC[V_3][k][j][i]/g->dPoints[X3][k] + mt->vC[V_3][k+1][j][i]/g->dPoints[X3][k+1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k+1]) - 
        g_dt*2.*(td->vCentroid[PRS][k+1][j][i] - td->vCentroid[PRS][k][j][i])/
        (g->dPoints[X3][k+1] + g->dPoints[X3][k]);
    w_b = (mt->vC[V_3][k-1][j][i]/g->dPoints[X3][k-1] + mt->vC[V_3][k][j][i]/g->dPoints[X3][k])/
        (1./g->dPoints[X3][k-1] + 1./g->dPoints[X3][k]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j][i] - td->vCentroid[PRS][k-1][j][i])/
        (g->dPoints[X3][k] + g->dPoints[X3][k-1]);

    flux = u_ep*u_e*S_e - u_wp*u_w*S_w + u_np*v_n*S_n - u_sp*v_s*S_s + u_fp*w_f*S_f - u_bp*w_b*S_b;

    return flux;
}

double advection_flux_prev_vsx2(int i, int j, int k, Data* d, Grid* g, tData* td, midT* mt) {
    double v_ep, v_wp, v_np, v_sp, v_fp, v_bp;
    double u_e, u_w, v_s, v_n, w_f, w_b;
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double flux;

    S_e = g->Area[k][j][i][E];
    S_w = g->Area[k][j][i][W];
    S_n = g->Area[k][j][i][N];
    S_s = g->Area[k][j][i][S];
    S_f = g->Area[k][j][i][F];
    S_b = g->Area[k][j][i][B];
    
    
    v_ep = (td->vCentroid[V2][k][j][i]/g->dPoints[X1][i] + td->vCentroid[V2][k][j][i+1]/g->dPoints[X1][i+1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i+1]);
    v_wp = (td->vCentroid[V2][k][j][i]/g->dPoints[X1][i] + td->vCentroid[V2][k][j][i-1]/g->dPoints[X1][i-1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i-1]);
    v_np = (td->vCentroid[V2][k][j][i]/g->dPoints[X2][j] + td->vCentroid[V2][k][j+1][i]/g->dPoints[X2][j+1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j+1]);
    v_sp = (td->vCentroid[V2][k][j][i]/g->dPoints[X2][j] + td->vCentroid[V2][k][j-1][i]/g->dPoints[X2][j-1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j-1]);
    v_fp = (td->vCentroid[V2][k][j][i]/g->dPoints[X3][k] + td->vCentroid[V2][k+1][j][i]/g->dPoints[X3][k+1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k+1]);
    v_bp = (td->vCentroid[V2][k][j][i]/g->dPoints[X3][k] + td->vCentroid[V2][k-1][j][i]/g->dPoints[X3][k-1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k-1]);

    /*v_ep = (td->vCentroid[V2][k][j][i] + td->vCentroid[V2][k][j][i+1])/2.;
    v_wp = (td->vCentroid[V2][k][j][i] + td->vCentroid[V2][k][j][i-1])/2.;
    v_np = (td->vCentroid[V2][k][j][i] + td->vCentroid[V2][k][j+1][i])/2.;
    v_sp = (td->vCentroid[V2][k][j][i] + td->vCentroid[V2][k][j-1][i])/2.;
    v_fp = (td->vCentroid[V2][k][j][i] + td->vCentroid[V2][k+1][j][i])/2.;
    v_bp = (td->vCentroid[V2][k][j][i] + td->vCentroid[V2][k-1][j][i])/2.;*/

    /*u_e = (mt->vC[V_1][k][j][i] + mt->vC[V_1][k][j][i+1])/2. - 
        g_dt*(td->vCentroid[PRS][k][j][i+1] - td->vCentroid[PRS][k][j][i])/g->dPoints[X1][i];
    u_w = (mt->vC[V_1][k][j][i] + mt->vC[V_1][k][j][i-1])/2. - 
        g_dt*(td->vCentroid[PRS][k][j][i] - td->vCentroid[PRS][k][j][i-1])/g->dPoints[X1][i];
    v_n = (mt->vC[V_2][k][j][i] + mt->vC[V_2][k][j+1][i])/2. -
        g_dt*(td->vCentroid[PRS][k][j+1][i] - td->vCentroid[PRS][k][j][i])/g->dPoints[X2][j];
    v_s = (mt->vC[V_2][k][j][i] + mt->vC[V_2][k][j-1][i])/2. -
        g_dt*(td->vCentroid[PRS][k][j][i] - td->vCentroid[PRS][k][j-1][i])/g->dPoints[X2][j];
    w_f = 0.0;
    w_b = 0.0;*/
    u_e = (mt->vC[V_1][k][j][i]/g->dPoints[X1][i] + mt->vC[V_1][k][j][i+1]/g->dPoints[X1][i+1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i+1]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j][i+1] - td->vCentroid[PRS][k][j][i])/
        (g->dPoints[X1][i] + g->dPoints[X1][i+1]);
    u_w = (mt->vC[V_1][k][j][i]/g->dPoints[X1][i] + mt->vC[V_1][k][j][i-1]/g->dPoints[X1][i-1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i-1]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j][i] - td->vCentroid[PRS][k][j][i-1])/
        (g->dPoints[X1][i] + g->dPoints[X1][i-1]);
    v_n = (mt->vC[V_2][k][j][i]/g->dPoints[X2][j] + mt->vC[V_2][k][j+1][i]/g->dPoints[X2][j+1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j+1]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j+1][i] - td->vCentroid[PRS][k][j][i])/
        (g->dPoints[X2][j] + g->dPoints[X2][j+1]);
    v_s = (mt->vC[V_2][k][j-1][i]/g->dPoints[X2][j-1] + mt->vC[V_2][k][j][i]/g->dPoints[X2][j])/
        (1./g->dPoints[X2][j-1] + 1./g->dPoints[X2][j]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j][i] - td->vCentroid[PRS][k][j-1][i])/
        (g->dPoints[X2][j] + g->dPoints[X2][j-1]);
    w_f = (mt->vC[V_3][k][j][i]/g->dPoints[X3][k] + mt->vC[V_3][k+1][j][i]/g->dPoints[X3][k+1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k+1]) - 
        g_dt*2.*(td->vCentroid[PRS][k+1][j][i] - td->vCentroid[PRS][k][j][i])/
        (g->dPoints[X3][k+1] + g->dPoints[X3][k]);
    w_b = (mt->vC[V_3][k-1][j][i]/g->dPoints[X3][k-1] + mt->vC[V_3][k][j][i]/g->dPoints[X3][k])/
        (1./g->dPoints[X3][k-1] + 1./g->dPoints[X3][k]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j][i] - td->vCentroid[PRS][k-1][j][i])/
        (g->dPoints[X3][k] + g->dPoints[X3][k-1]);

    flux = u_e*v_ep*S_e - u_w*v_wp*S_w + v_n*v_np*S_n - v_s*v_sp*S_s + v_fp*w_f*S_f - v_bp*w_b*S_b;

    return flux;
}

double advection_flux_prev_vsx3(int i, int j, int k, Data* d, Grid* g, tData* td, midT* mt) {
    double w_ep, w_wp, w_np, w_sp, w_fp, w_bp;
    double u_e, u_w, v_s, v_n, w_f, w_b;
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double flux;

    S_e = g->Area[k][j][i][E];
    S_w = g->Area[k][j][i][W];
    S_n = g->Area[k][j][i][N];
    S_s = g->Area[k][j][i][S];
    S_f = g->Area[k][j][i][F];
    S_b = g->Area[k][j][i][B];

    w_ep = (td->vCentroid[V3][k][j][i]/g->dPoints[X1][i] + td->vCentroid[V3][k][j][i+1]/g->dPoints[X1][i+1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i+1]);
    w_wp = (td->vCentroid[V3][k][j][i]/g->dPoints[X1][i] + td->vCentroid[V3][k][j][i-1]/g->dPoints[X1][i-1])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i-1]);
    w_np = (td->vCentroid[V3][k][j][i]/g->dPoints[X2][j] + td->vCentroid[V3][k][j+1][i]/g->dPoints[X2][j+1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j+1]);
    w_sp = (td->vCentroid[V3][k][j][i]/g->dPoints[X2][j] + td->vCentroid[V3][k][j-1][i]/g->dPoints[X2][j-1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j-1]);
    w_fp = (td->vCentroid[V3][k][j][i]/g->dPoints[X3][k] + td->vCentroid[V3][k+1][j][i]/g->dPoints[X3][k+1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k+1]);
    w_bp = (td->vCentroid[V3][k-1][j][i]/g->dPoints[X3][k] + td->vCentroid[V3][k][j][i]/g->dPoints[X3][k])/
        (1./g->dPoints[X3][k-1] + 1./g->dPoints[X3][k]);

    u_e = (mt->vC[V_1][k][j][i]/g->dPoints[X1][i] + mt->vC[V_1][k][j][i+1]/g->dPoints[X1][i])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i+1]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j][i+1] - td->vCentroid[PRS][k][j][i])/
        (g->dPoints[X1][i] + g->dPoints[X1][i+1]);
    u_w = (mt->vC[V_1][k][j][i]/g->dPoints[X1][i] + mt->vC[V_1][k][j][i-1]/g->dPoints[X1][i])/
        (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i-1]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j][i] - td->vCentroid[PRS][k][j][i-1])/
        (g->dPoints[X1][i] + g->dPoints[X1][i-1]);
    v_n = (mt->vC[V_2][k][j][i]/g->dPoints[X2][j] + mt->vC[V_2][k][j+1][i]/g->dPoints[X2][j+1])/
        (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j+1]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j+1][i] - td->vCentroid[PRS][k][j][i])/
        (g->dPoints[X2][j] + g->dPoints[X2][j+1]);
    v_s = (mt->vC[V_2][k][j-1][i]/g->dPoints[X2][j-1] + mt->vC[V_2][k][j][i]/g->dPoints[X1][j])/
        (1./g->dPoints[X2][j-1] + 1./g->dPoints[X2][j]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j][i] - td->vCentroid[PRS][k][j-1][i])/
        (g->dPoints[X2][j] + g->dPoints[X2][j-1]);
    w_f = (mt->vC[V_3][k][j][i]/g->dPoints[X3][k] + mt->vC[V_3][k+1][j][i]/g->dPoints[X3][k+1])/
        (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k+1]) - 
        g_dt*2.*(td->vCentroid[PRS][k+1][j][i] - td->vCentroid[PRS][k][j][i])/
        (g->dPoints[X3][k+1] + g->dPoints[X3][k]);
    w_b = (mt->vC[V_3][k-1][j][i]/g->dPoints[X3][k] + mt->vC[V_3][k][j][i]/g->dPoints[X3][k])/
        (1./g->dPoints[X3][k-1] + 1./g->dPoints[X3][k]) - 
        g_dt*2.*(td->vCentroid[PRS][k][j][i] - td->vCentroid[PRS][k-1][j][i])/
        (g->dPoints[X3][k] + g->dPoints[X3][k-1]);

    flux = u_e*w_ep*S_e - u_w*w_wp*S_w + w_np*v_n*S_n - w_sp*v_s*S_s + w_fp*w_f*S_f - w_bp*w_b*S_b;

    return flux;
}

double diffusion_flux_prev_vxs1(int i, int j, int k, Data* d, Grid* g) {
    double g_e, g_w, g_n, g_s, g_f, g_b;
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double flux;

    S_e = g->Area[k][j][i][E];
    S_w = g->Area[k][j][i][W];
    S_n = g->Area[k][j][i][N];
    S_s = g->Area[k][j][i][S];
    S_f = g->Area[k][j][i][F];
    S_b = g->Area[k][j][i][S];

    /*g_e = (d->vCentroid[V1][k][j][i+1] - d->vCentroid[V1][k][j][i])/g->dPoints[X1][i];
    g_w = (d->vCentroid[V1][k][j][i] - d->vCentroid[V1][k][j][i-1])/g->dPoints[X1][i];
    g_n = (d->vCentroid[V1][k][j+1][i] - d->vCentroid[V1][k][j][i])/g->dPoints[X2][j];
    g_s = (d->vCentroid[V1][k][j][i] - d->vCentroid[V1][k][j-1][i])/g->dPoints[X2][j];
    g_f = 0.0;
    g_b = 0.0;*/

    g_e = 2.*(d->vCentroid[V1][k][j][i+1] - d->vCentroid[V1][k][j][i])/
        (g->dPoints[X1][i+1] + g->dPoints[X1][i]);
    g_w = 2.*(d->vCentroid[V1][k][j][i] - d->vCentroid[V1][k][j][i-1])/
        (g->dPoints[X1][i-1] + g->dPoints[X1][i]);;
    g_n = 2.*(d->vCentroid[V1][k][j+1][i] - d->vCentroid[V1][k][j][i])/
        (g->dPoints[X2][j+1] + g->dPoints[X2][j]);
    g_s = 2.*(d->vCentroid[V1][k][j][i] - d->vCentroid[V1][k][j-1][i])/
        (g->dPoints[X2][j] + g->dPoints[X2][j-1]);
    g_f = 2.*(d->vCentroid[V1][k+1][j][i] - d->vCentroid[V1][k][j][i])/
        (g->dPoints[X3][k+1] + g->dPoints[X3][k]);
    g_b = 2.*(d->vCentroid[V1][k][j][i] - d->vCentroid[V1][k-1][j][i])/
        (g->dPoints[X3][k] + g->dPoints[X3][k-1]);

    flux = (g_e*S_e - g_w*S_w + g_n*S_n - g_s*S_s + g_f*S_f - g_b*S_b)/g_Re;

    return flux;
}

double diffusion_flux_prev_vxs2(int i, int j, int k, Data* d, Grid* g) {
    double g_e, g_w, g_n, g_s, g_f, g_b;
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double flux;

    S_e = g->Area[k][j][i][E];
    S_w = g->Area[k][j][i][W];
    S_n = g->Area[k][j][i][N];
    S_s = g->Area[k][j][i][S];
    S_f = g->Area[k][j][i][F];
    S_b = g->Area[k][j][i][B];

    /*g_e = (d->vCentroid[V2][k][j][i+1] - d->vCentroid[V2][k][j][i])/g->dPoints[X1][i];
    g_w = (d->vCentroid[V2][k][j][i] - d->vCentroid[V2][k][j][i-1])/g->dPoints[X1][i];
    g_n = (d->vCentroid[V2][k][j+1][i] - d->vCentroid[V2][k][j][i])/g->dPoints[X2][j];
    g_s = (d->vCentroid[V2][k][j][i] - d->vCentroid[V2][k][j-1][i])/g->dPoints[X2][j];
    g_f = 0.0;
    g_b = 0.0;*/

    g_e = 2.*(d->vCentroid[V2][k][j][i+1] - d->vCentroid[V2][k][j][i])/
        (g->dPoints[X1][i+1] + g->dPoints[X1][i]);
    g_w = 2.*(d->vCentroid[V2][k][j][i] - d->vCentroid[V2][k][j][i-1])/
        (g->dPoints[X1][i] + g->dPoints[X1][i-1]);
    g_n = 2.*(d->vCentroid[V2][k][j+1][i] - d->vCentroid[V2][k][j][i])/
        (g->dPoints[X2][j+1] + g->dPoints[X2][j]);
    g_s = 2.*(d->vCentroid[V2][k][j][i] - d->vCentroid[V2][k][j-1][i])/
        (g->dPoints[X2][j-1] + g->dPoints[X2][j]);;
    g_f = 2.*(d->vCentroid[V2][k+1][j][i] - d->vCentroid[V2][k][j][i])/
        (g->dPoints[X3][k+1] + g->dPoints[X3][k]);
    g_b = 2.*(d->vCentroid[V2][k][j][i] - d->vCentroid[V2][k-1][j][i])/
        (g->dPoints[X3][k] + g->dPoints[X3][k-1]);

    flux = (g_e*S_e - g_w*S_w + g_n*S_n - g_s*S_s + g_f*S_f - g_b*S_b)/g_Re;

    return flux;
}

double diffusion_flux_prev_vxs3(int i, int j, int k, Data* d, Grid* g) {
    double g_e, g_w, g_n, g_s, g_f, g_b;
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double flux;

    S_e = g->Area[k][j][i][E];
    S_w = g->Area[k][j][i][W];
    S_n = g->Area[k][j][i][N];
    S_s = g->Area[k][j][i][S];
    S_f = g->Area[k][j][i][F];
    S_b = g->Area[k][j][i][B];

    g_e = 2.*(d->vCentroid[V3][k][j][i+1] - d->vCentroid[V3][k][j][i])/
        (g->dPoints[X1][i+1] + g->dPoints[X1][i]);
    g_w = 2.*(d->vCentroid[V3][k][j][i] - d->vCentroid[V3][k][j][i-1])/
        (g->dPoints[X1][i] + g->dPoints[X1][i-1]);
    g_n = 2.*(d->vCentroid[V3][k][j+1][i] - d->vCentroid[V3][k][j][i])/
        (g->dPoints[X2][j+1] + g->dPoints[X2][j]);
    g_s = 2.*(d->vCentroid[V3][k][j][i] - d->vCentroid[V3][k][j-1][i])/
        (g->dPoints[X2][j] + g->dPoints[X2][j-1]);
    g_f = 2.*(d->vCentroid[V3][k+1][j][i] - d->vCentroid[V3][k][j][i])/
        (g->dPoints[X3][k+1] + g->dPoints[X3][k]);
    g_b = 2.*(d->vCentroid[V3][k][j][i] - d->vCentroid[V3][k-1][j][i])/
        (g->dPoints[X3][k] + g->dPoints[X3][k-1]);

    flux = (g_e*S_e - g_w*S_w + g_n*S_n - g_s*S_s + g_f*S_f - g_b*S_b)/g_Re;

    return flux;
}

void calculate_flux(Data* d, Grid* g, tData* td, midT* mt, Boundary* b) {
    if (g_step == 1) {
        for (int k=KBEG; k<=KEND; k++) {
            for (int j=JBEG; j<=JEND; j++) {
                for (int i=IBEG; i<=IEND; i++) {
                    mt->vC[V_1][k][j][i] = (diffusion_flux_prev_vxs1(i,j,k,d,g) - advection_flux_prev_vsx1(i,j,k,d,g,td,mt))*
                        g_dt/(g->dPoints[X1][i]*g->dPoints[X2][j]*g->dPoints[X3][k]) +
                        d->vCentroid[V1][k][j][i];
                    mt->vC[V_2][k][j][i] = (diffusion_flux_prev_vxs2(i,j,k,d,g) - advection_flux_prev_vsx2(i,j,k,d,g,td,mt))*
                        g_dt/(g->dPoints[X2][j]*g->dPoints[X1][i]*g->dPoints[X3][k]) +
                        d->vCentroid[V2][k][j][i];
                    mt->vC[V_3][k][j][i] = (diffusion_flux_prev_vxs3(i,j,k,d,g) - advection_flux_prev_vsx3(i,j,k,d,g,td,mt))*
                        g_dt/(g->dPoints[X3][k]*g->dPoints[X2][j]*g->dPoints[X1][i]) +
                        d->vCentroid[V3][k][j][i];
                }
            }
        }
    }
    else {
        for (int k=KBEG; k<=KEND; k++) {
            for (int j=JBEG; j<=JEND; j++) {
                for (int i=IBEG; i<=IEND; i++) {
                    mt->vC[V_1][k][j][i] = (diffusion_flux_prev_vxs1(i,j,k,d,g) - d->advecFlux[V_1][k][j][i])*
                        g_dt/(g->dPoints[X1][i]*g->dPoints[X2][j]*g->dPoints[X3][k]) +
                        d->vCentroid[V1][k][j][i];
                    mt->vC[V_2][k][j][i] = (diffusion_flux_prev_vxs2(i,j,k,d,g) - d->advecFlux[V_2][k][j][i])*
                        g_dt/(g->dPoints[X2][j]*g->dPoints[X1][i]*g->dPoints[X3][k]) +
                        d->vCentroid[V2][k][j][i];
                    mt->vC[V_3][k][j][i] = (diffusion_flux_prev_vxs3(i,j,k,d,g) - d->advecFlux[V_3][k][j][i])*
                        g_dt/(g->dPoints[X3][k]*g->dPoints[X2][j]*g->dPoints[X1][i]) +
                        d->vCentroid[V3][k][j][i];
                }
            }
        }
    }
    update_fractional_velocity(b, d, g, td, mt);
    /*for (int k=KBEG; k<=KEND; k++) {
            for (int j=JBEG-1; j<=JEND+1; j++) {
                for (int i=IBEG-1; i<=IEND+1; i++) {
                    printf("%g,", td->vCentroid[PRS][k][j][i]);
                } printf("\n");
            } printf("\n\n");
        }*/
    pressure_poisson(d, g, td, mt);
    update_pressure(b, d, g, td, mt);
    max_pressure(td);
    double prs_prev = 0.0;
    double prs_new = td->vCentroid[PRS][td->kPos][td->jPos][td->iPos];

    /*for (int k=KBEG; k<=KEND; k++) {
            for (int j=JBEG; j<=JEND; j++) {
                for (int i=IBEG; i<=IEND; i++) {
                    td->vCentroid[PRS][k][j][i] = 0.0;
                }
            }
        }*/
    
    while (fabs(fabs(prs_new)-fabs(prs_prev)) > TOLLERANCE) {
        
        pressure_poisson(d, g, td, mt);
        update_pressure(b, d, g, td, mt);
        prs_prev = prs_new;
        prs_new = td->vCentroid[PRS][td->kPos][td->jPos][td->iPos];
    }

    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                td->vCentroid[V1][k][j][i] = mt->vC[V_1][k][j][i] - g_dt*2.*
                    (td->vCentroid[PRS][k][j][i+1] - td->vCentroid[PRS][k][j][i-1])/
                    (g->dPoints[X1][i] + g->dPoints[X1][i+1]);
                td->vCentroid[V2][k][j][i] = mt->vC[V_2][k][j][i] - g_dt*2.*
                    (td->vCentroid[PRS][k][j+1][i] - td->vCentroid[PRS][k][j-1][i])/
                    (g->dPoints[X2][j] + g->dPoints[X2][j+1]);
                td->vCentroid[V3][k][j][i] = mt->vC[V_3][k][j][i] - g_dt*2.*
                    (td->vCentroid[PRS][k+1][j][i] - td->vCentroid[PRS][k-1][j][i])/
                    (g->dPoints[X3][k] + g->dPoints[X3][k+1]);
            }
        }
    }
    update_next_velocity(b, d, g, td, mt);

    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                d->advecFlux[V_1][k][j][i] = advection_flux_prev_vsx1(i,j,k,d,g,td,mt);
                d->advecFlux[V_2][k][j][i] = advection_flux_prev_vsx2(i,j,k,d,g,td,mt);
                d->advecFlux[V_3][k][j][i] = advection_flux_prev_vsx3(i,j,k,d,g,td,mt);
            }
        }
    }
}

void pressure_poisson(Data* d, Grid* g, tData* td, midT* mt) {
    double S_e, S_w, S_n, S_s, S_f, S_b;
    double u_e, u_w, v_n, v_s, w_f, w_b;
    double vel;
    double den, num;
    double dx1, dx2, dy1, dy2, dz1, dz2;

    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                S_e = g->dPoints[X2][j]*g->dPoints[X3][k];
                S_n = g->dPoints[X1][i]*g->dPoints[X3][k];
                S_f = g->dPoints[X1][i]*g->dPoints[X2][j];
                S_w = S_e;
                S_s = S_n;
                S_b = S_f;

                dx1 = (g->dPoints[X1][i] + g->dPoints[X1][i+1])/2.;
                dx2 = (g->dPoints[X1][i] + g->dPoints[X1][i-1])/2.;
                dy1 = (g->dPoints[X2][j] + g->dPoints[X2][j+1])/2.;
                dy2 = (g->dPoints[X2][j] + g->dPoints[X2][j-1])/2.;
                dz1 = (g->dPoints[X3][k] + g->dPoints[X3][k+1])/2.;
                dz2 = (g->dPoints[X3][k] + g->dPoints[X3][k-1])/2.;

                u_e = (mt->vC[V_1][k][j][i]/g->dPoints[X1][i] + mt->vC[V_1][k][j][i+1]/g->dPoints[X1][i+1])/
                    (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i+1]);
                u_w = (mt->vC[V_1][k][j][i]/g->dPoints[X1][i] + mt->vC[V_1][k][j][i-1]/g->dPoints[X1][i-1])/
                    (1./g->dPoints[X1][i] + 1./g->dPoints[X1][i-1]);
                v_n = (mt->vC[V_2][k][j][i]/g->dPoints[X2][j] + mt->vC[V_2][k][j+1][i]/g->dPoints[X2][j+1])/
                    (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j+1]);
                v_s = (mt->vC[V_2][k][j][i]/g->dPoints[X2][i] + mt->vC[V_2][k][j-1][i]/g->dPoints[X2][j-1])/
                    (1./g->dPoints[X2][j] + 1./g->dPoints[X2][j-1]);
                w_f = (mt->vC[V_3][k][j][i]/g->dPoints[X3][k] + mt->vC[V_3][k+1][j][i]/g->dPoints[X3][k+1])/
                    (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k+1]);
                w_b = (mt->vC[V_3][k][j][i]/g->dPoints[X2][i] + mt->vC[V_3][k-1][j][i]/g->dPoints[X3][k-1])/
                    (1./g->dPoints[X3][k] + 1./g->dPoints[X3][k-1]);

                vel = u_e*S_e - u_w*S_w + v_n*S_n - v_s*S_s + w_f*S_f - w_b*S_b;

                den = S_e/dx1 + S_w/dx2 + S_n/dy1 + S_s/dy2 + S_f/dz1 + S_b/dz2;

                num = S_e*td->vCentroid[PRS][k][j][i+1]/dx1 + S_w*td->vCentroid[PRS][k][j][i-1]/dx2+
                    S_n*td->vCentroid[PRS][k][j+1][i]/dy1 + S_s*td->vCentroid[PRS][k][j-1][i]/dy2+
                    S_f*td->vCentroid[PRS][k+1][j][i]/dz1 + S_b*td->vCentroid[PRS][k-1][j][i]/dz2 - vel/g_dt;
                
                td->vCentroid[PRS][k][j][i] = num/den;
            }
        }
    }
}

void max_pressure(tData* td) {
    double max_P = 0.0;
    int index[3];
    for (int k=KBEG; k<=KEND; k++) {
        for (int j=JBEG; j<=JEND; j++) {
            for (int i=IBEG; i<=IEND; i++) {
                max_P = MAX(td->vCentroid[PRS][k][j][i], max_P);
                if (max_P == td->vCentroid[PRS][k][j][i]) {
                    td->iPos = i;
                    td->jPos = j;
                    td->kPos = k;
                }
            }
        }
    }
}
#endif

#endif
