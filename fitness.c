/* Tumour Parameters */
#include <values.h>
#include <math.h>
#include "RKF78.h"

#define lambda_par 0.336
#define lambdaLogTheta_par 9.284023094951992774680543277273261947643
#define NZero_par 20000000
#define CCUMj_par 127
#define NMax_par 9.5e+11
#define m_par 4
#define d_par 10
#define n_par 11

//static float k_j[d_par] = {0.3,0.1255,0.15925,0.33675,0.2255,0.1365,0.19175,0.28025,0.24275,0.10075};
static float k_j[d_par] = {0.12, 0.0502, 0.0637, 0.1347, 0.0902, 0.0546, 0.0767, 0.1121, 0.0971, 0.0403};
static float t_i[n_par] = {0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 33.0};
/*static float eta_kj[m_par][d_par] = {
        {0.0036, 0.0098, 0.0061, 0.0009, 0.0003, 0.0108, 0.0045, 0.0021, 0.0096, 0.0125},
        {0.0063, 0.0082, 0.0062, 0.0062, 0.0083, 0.013, 0.0039, 0.0019, 0.0015, 0.005},
        {0.0129, 0.0018, 0.0116, 0.0021, 0.009, 0.0129, 0.0054, 0.0049, 0.0093, 0.0066},
        {0.0053, 0.0086, 0.0067, 0.0029, 0.0089, 0.0054, 0.0042, 0.0095, 0.0112, 0.0092}
};//*/

static float eta_kj[m_par][d_par] = {
        {0.00612, 0.01666, 0.01037, 0.00153, 0.00051, 0.01836, 0.00765, 0.00357, 0.01632, 0.02125},
        {0.01071, 0.01394, 0.01054, 0.01054, 0.01411, 0.0221, 0.00663, 0.00323, 0.00255, 0.0085},
        {0.02193, 0.00306, 0.01972, 0.00357, 0.0153, 0.02193, 0.00918, 0.00833, 0.01581, 0.01122},
        {0.00901, 0.01462, 0.01139, 0.00493, 0.01513, 0.00918, 0.00714, 0.01615, 0.01904, 0.01564}
};//*/


typedef struct { double drift_i; } ODE_Parameters;

void Gompertz(double t, double N, double *der, void *Params){
    *der = ((N < 1.e-16) ? 0.0 :
            N*(lambdaLogTheta_par - lambda_par * log(N) - ((ODE_Parameters *) Params)->drift_i));
}

unsigned char TestIfConstraints2and3AreVerified(unsigned char *Cij){
    register unsigned char i, j, k;
    unsigned char npar = n_par - 1, *Cijofi;

    for(j=0; j < d_par; j++){
        unsigned int ccumj = 0U;
        for(i=0; i < npar; i++) ccumj += *(Cij + i*d_par + j);
        if(ccumj > CCUMj_par) return 0U;
    }

    for(i=0, Cijofi=Cij; i < npar; i++, Cijofi += d_par) {
        for(k=0; k < m_par; k++){
            double Cseffk = 0.0;
            for(j=0; j < d_par; j++) Cseffk += eta_kj[k][j] * Cijofi[j];
            if(Cseffk > 1.0) return 0U;
        }
    }

    return 1U;
}

double Curative_Fitness(unsigned char *Cij){
    register unsigned char i, j;
    ODE_Parameters GompertzParams;
    double N = NZero_par, t = t_i[0];
    double hmin = 1.e-8, hmax = 1.0, h = 1.e-3, tol = 1.e-8;
    unsigned char curative_counter = 0U, npar = n_par - 1;
    double integral = 0.0, lastt = t, lastN = N;

    //printf("Check Constraint 1 and 2\n");
    if(!TestIfConstraints2and3AreVerified(Cij)) return MAXDOUBLE;
    //printf("Passed!\n");

    for(i=0; i < npar; i++){ // Implementing treatment i
        double tfin = t_i[i+1];
        GompertzParams.drift_i = 0.0;
        for(j=0; j < d_par; j++) GompertzParams.drift_i += k_j[j] * *(Cij++);

        while(t+h < tfin) {
            RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
            if(N > NMax_par) return MAXDOUBLE;
            integral += (lastN + N)*(t - lastt);
            lastt = t; lastN = N;
        }

        do {
            h = tfin - t;
            RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
            if(N > NMax_par) return MAXDOUBLE;
            integral += (lastN + N)*(t - lastt);
            lastt = t; lastN = N;
        } while (t < tfin);

        /**
         if(N < 1000) {
            curative_counter++;
            if(curative_counter > 2) return N;
        }
        else curative_counter = 0U;//*/
    }
    // return MAXDOUBLE;
    return N;
}