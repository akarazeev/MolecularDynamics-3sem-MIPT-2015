//
//  main.c
//  MolecDynam
//
//  Created by Anton Wetret on 17/09/15.
//  Copyright (c) 2015 Anton Karazeev. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#define trace(x) printf("%f\n", x);
#define dtrace(x, y) printf("%f %f\n", x, y);

#define PRINT_TO_FILE 1
#define READ_INIT 0
#define USE_BERENDSEN 0

const char* READ_FROM = "MolecDynam/init_coord_N=100.xyz";

const double Temp0 = 1.3;
const int N = 64;
const double dt = 0.001;
const double iterations = 20000;
const double density = 0.1;

double rcut2 = 9;
double mAr = 1;
double K = 0;
double Temp = 0;
double utot = 0;

double v[N][3];
double f[N][3];
double r[N][3];
double rn[N][3];
double L2[3];
double L[3];
double m[N];

int flag = 0;

static inline double Potential(double x) {
    if (x < rcut2) {
        double res = 4.0 * ( (1.0/(float)powf(x,6.0)) - (1.0/(float)powf(x,3.0)) );
        return res;
    } else {
        double res = 0;
        return res;
    }
}

static inline double ForceDevByRange(double x) {
    double res = 48.0 * ( (1.0/(float)powf(x,7.0)) - (1.0/((float)powf(x,4.0) * 2.0)) );
    return res;
}

static inline void nearest_image() {
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            if (r[i][k] > 0) {
                rn[i][k] = fmod(r[i][k] + L2[k], L[k]) - L2[k];
            } else {
                rn[i][k] = fmod(r[i][k] - L2[k], L[k]) + L2[k];
            }
        }
    }
}

static inline void EqMotion() {
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            //            printf("F %f\n", f[i][k]);
            v[i][k] += ((float)f[i][k] * (float)dt) / (float)m[i];
            //            printf("V %f\n", v[i][k]);
            r[i][k] += (float)v[i][k] * (float)dt;
            //            printf("R %f\n", r[i][k]);
            assert(r[i][k] == r[i][k]);
        }
    }
}

static inline void ClearForces() {
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            f[i][k] = 0;
        }
    }
}

static inline void CalcForces() {
    double rij[3];
    utot = 0;
    double r2;
    for (int i = 1; i < N; ++i) {
        for (int j = 0; j < i; ++j) {
            r2 = 0;
            for (int k = 0; k < 3; ++k) {
                //                printf("Rn %f %f %d %d\n", rn[i][k], rn[j][k], i, j);
                rij[k] = rn[i][k] - rn[j][k];
                if (rij[k] > L2[k]) {
                    rij[k] -= L[k];
                } else if (rij[k] < -L2[k]) {
                    rij[k] += L[k];
                }
                r2 += rij[k] * rij[k];
                //                printf("R2 %f\n", r2);
            }
            double f_r = 0;
            utot += Potential(r2);
            f_r = ForceDevByRange(r2);
            for (int k = 0; k < 3; ++k) {
                //                printf("F/R %f\n", f_r);
                //                printf("Rij %f\n", rij[k]);
                assert(f_r == f_r);
                f[i][k] += (float)f_r * (float)rij[k];
                f[j][k] -= (float)f_r * (float)rij[k];
            }
        }
    }
}

static inline void CalcEnergy() {
    K = 0;
    double v2 = 0;
    for (int i = 0; i < N; ++i) {
        v2 = 0;
        for (int k = 0; k < 3; ++k) {
            //            printf("V %f\n", v[i][k]);
            v2 += (float)v[i][k] * (float)v[i][k];
            //            printf("V2 %f\n", v2);
        }
        K += (float)m[i] * ((float)v2 / 2.0);
        assert(K != INFINITY);
    }
}

static inline void CalcTemp() {
    Temp = 2.0 * ((K / N) / 3.0);
}

static inline void Thermostat(int cur_iter) {
    if (USE_BERENDSEN) {
        double tau = 0.001;
        double lambda = sqrt(1 + ( (dt/tau) * ((Temp0/Temp) - 1) ));
        if (fabs(Temp - Temp0) < 0.1 && cur_iter > 300) {
            flag = 1;
            printf("%d\n", cur_iter);
        }
                if (!flag) {
//        if (cur_iter < iterations/2) {
            for (int i = 0; i < N; ++i) {
                for (int k = 0; k < 3; ++k) {
                    v[i][k] *= lambda;
                }
            }
        }
    }
}

int main() {
    srand((unsigned int)time(NULL));
    double length = powf(N/density, 1.0/3.0);
    
    for (int k = 0; k < 3; ++k) {
        L[k] = length;
        L2[k] = L[k]/2.0;
        assert(L[k] > sqrt(rcut2));
    }
    for (int i = 0; i < N; ++i) {
        m[i] = mAr;
    }
    FILE* f_en;
    FILE* f_kin;
    FILE* f_poten;
    FILE* f_temp;
    FILE* f_xyz;
    FILE* f_velocity;
    FILE* f_init_coord;
    FILE* f_init_coord_r;
    if (PRINT_TO_FILE || READ_INIT) {
        f_init_coord = fopen("MolecDynam/init_coord.xyz", "w");
        f_init_coord_r = fopen(READ_FROM, "r");
        f_xyz = fopen("MolecDynam/t.xyz", "w");
        f_temp = fopen("molec_dynam_r/temp.csv", "w");
        f_velocity = fopen("molec_dynam_r/velocity.csv", "w");
        f_en = fopen("molec_dynam_r/energy.csv", "w");
        f_kin = fopen("molec_dynam_r/kinetic.csv", "w");
        f_poten = fopen("molec_dynam_r/poten.csv", "w");
    }
    if (!READ_INIT) {
        // Make Initial Coordinates
        double init[3];
        
        // Quantity of atoms pro line
        int quant = powf(N, 1.0/3.0);
        double step = L[0] / quant;
        double ic = (step/2.0) - L2[0];
        
        for (int i = 0; i < 3; ++i) {
            init[i] = ic;
        }
        
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < 3; ++k) {
                r[i][k] = init[k];
            }
            init[0] += step;
            if ((i+1) % (quant*quant) == 0) {
                init[0] = ic;
                init[1] = ic;
                init[2] += step;
            } else if ((i+1) % quant == 0) {
                init[0] = ic;
                init[1] += step;
            }
        }
        // Make Initial Velocities
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < 3; ++k) {
                v[i][k] = 1 - (((float)rand()/(float)RAND_MAX)*2);
            }
        }
        
    }
    if (PRINT_TO_FILE || READ_INIT) {
        if (READ_INIT) {
            // Read Initial Coordinates
            double tmp;
            fscanf(f_init_coord_r, "%lf", &tmp);
            assert(tmp == N);
            fscanf(f_init_coord_r, "%lf", &tmp);
            for (int i = 0; i < N; ++i) {
                char c;
                fscanf(f_init_coord_r, "%c", &c);
                for (int k = 0; k < 3; ++k) {
                    fscanf(f_init_coord_r, "%lf", &tmp);
                    r[i][k] = tmp;
                }
                fscanf(f_init_coord_r, "%c", &c);
                fscanf(f_init_coord_r, "%c", &c);
            }
        } else {
            // Print Coordinates to File
            fprintf(f_init_coord, "%d\n\n", N);
            for (int i = 0; i < N; ++i) {
                fprintf(f_init_coord, "%c ", (char) (97+(i%26)));
                for (int k = 0; k < 3; ++k) {
                    fprintf(f_init_coord, "%f ", r[i][k]);
                }
                fprintf(f_init_coord, "\n");
            }
            fclose(f_init_coord);
        }
    }
    
    for (int i = 0; i < iterations; ++i) {
        nearest_image();
        ClearForces();
        CalcForces();
        EqMotion();
        CalcEnergy();
        CalcTemp();
        Thermostat(i);
        
        if (PRINT_TO_FILE) {
            //            if (i == iterations/2) {
            if (i > iterations/2) {
                for (int j = 0; j < N; ++j) {
                    double v2 = 0;
                    for (int k = 0; k < 3; ++k) {
                        v2 += (float)v[j][k] * (float)v[j][k];
                    }
                    fprintf(f_velocity, "%f\n", sqrt(v2));
                }
                //                fclose(f_velocity);
            }
            //            if (i > iterations-500) {
            if (i < 500) {
                fprintf(f_xyz, "%d\n\n", N);
                for (int i = 0; i < N; ++i) {
                    fprintf(f_xyz, "%c ", (char) (97+(i%26)));
                    for (int k = 0; k < 3; ++k) {
                        fprintf(f_xyz, "%f ", rn[i][k]);
                    }
                    fprintf(f_xyz, "\n");
                }
            }
            //            if (i > 4) {
            fprintf(f_en, "%f,%i\n", K + utot, i);
            fprintf(f_kin, "%f,%i\n", K, i);
            fprintf(f_poten, "%f,%i\n", utot, i);
            fprintf(f_temp, "%f,%i\n", Temp, i);
            //            }
        }
        
        if (i*2 == iterations) {
            trace(Temp)
            dtrace(K, utot)
        }
        //        trace(Temp)
        //        dtrace(K, utot)
        assert(K != INFINITY);
    }
    if (PRINT_TO_FILE) {
        fclose(f_velocity);
        fclose(f_xyz);
        fclose(f_en);
        fclose(f_poten);
        fclose(f_kin);
        fclose(f_init_coord_r);
        fclose(f_temp);
        printf("Print to file: Done!\n");
    }
    printf("Result: Done!\n");
    return 0;
}
