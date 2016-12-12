//
//  main.c
//  MolecDynam
//
//  Created by Anton Karazeev on 17/09/15.
//  Copyright (c) 2015 Anton Karazeev. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#define PRINT_TO_FILE 1
#define USE_BERENDSEN 0

const int N = 343;
const int iterations = 15000;

const double Temp0 = 3;
const double tau = 1;
const double dt = 0.001;
const double iter_to_write = 4000;
const double density = 0.6;

const double rcut2 = 9;
const double mAr = 1;

/* current values of each iteration */
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
    if (rcut2 > 0) {
        if (x < rcut2) {
            double res = 4.0 * ( (1.0/(float)powf(x,6.0)) - (1.0/(float)powf(x,3.0)) );
            return res;
        } else {
            double res = 0;
            return res;
        }
    } else {
        double res = 4.0 * ( (1.0/(float)powf(x,6.0)) - (1.0/(float)powf(x,3.0)) );
        return res;
    }
}

static inline double ForceDivByRange(double x) {
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
            v[i][k] += ((float)f[i][k] * (float)dt) / (float)m[i];
            r[i][k] += (float)v[i][k] * (float)dt;
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
                rij[k] = rn[i][k] - rn[j][k];
                if (rij[k] > L2[k]) {
                    rij[k] -= L[k];
                } else if (rij[k] < -L2[k]) {
                    rij[k] += L[k];
                }
                r2 += rij[k] * rij[k];
            }
            double f_r = 0;
            utot += Potential(r2);
            f_r = ForceDivByRange(r2);
            for (int k = 0; k < 3; ++k) {
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
            v2 += (float)v[i][k] * (float)v[i][k];
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
        double lambda = sqrt(1 + ( (dt/tau) * ((Temp0/Temp) - 1) ));
        if (fabs(Temp - Temp0) < 0.1 && cur_iter > 300) {
            flag = 1;
        }
        if (!flag) {
            for (int i = 0; i < N; ++i) {
                for (int k = 0; k < 3; ++k) {
                    v[i][k] *= lambda;
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    srand((unsigned int)time(NULL));
    // Calculate the length of the box side
    double length = powf(N/density, 1.0/3.0);
    // Set boundaries
    for (int k = 0; k < 3; ++k) {
        L[k] = length;
        L2[k] = L[k]/2.0;
        if (rcut2 > 0) {
            assert(L[k] > sqrt(rcut2));
        }
    }
    // Set masses for every particle
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
    FILE* f_len;
    if (PRINT_TO_FILE) {
        f_init_coord = fopen("data/init_coord.xyz", "w");
        f_xyz = fopen("data/coordinates.xyz", "w");
        f_temp = fopen("data/temp.csv", "w");
        f_velocity = fopen("data/velocity.csv", "w");
        f_en = fopen("data/energy.csv", "w");
        f_kin = fopen("data/kinetic.csv", "w");
        f_poten = fopen("data/poten.csv", "w");
        f_len = fopen("data/len.csv", "w");
    }
    // Store length to file
    fprintf(f_len, "%f", length);
    
    /* Set Initial Coordinates */
    double init[3];
    // Quantity of atoms per line
    int quant = powf(N, 1.0/3.0);
    // Step of the mesh
    double step = L[0] / quant;
    double ic = (step/2.0) - L2[0];

    for (int i = 0; i < 3; ++i) {
        init[i] = ic;
    }
    // Set initial coordinates
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
    // Set initial velocities
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            v[i][k] = 1 - (((float)rand()/(float)RAND_MAX)*2);
        }
    }
    if (PRINT_TO_FILE) {
            // Print initial coordinates to file
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
    
    ///////////////
    // Main Part //
    ///////////////
    
    for (int i = 0; i < iterations; ++i) {
        nearest_image();
        ClearForces();
        CalcForces();
        EqMotion();
        CalcEnergy();
        CalcTemp();
        Thermostat(i);

        if (PRINT_TO_FILE) {
            // Store energy and temperature to file
            fprintf(f_en, "%f,%i\n", K + utot, i);
            fprintf(f_kin, "%f,%i\n", K, i);
            fprintf(f_poten, "%f,%i\n", utot, i);
            fprintf(f_temp, "%f,%i\n", Temp, i);
            // Store velocities to file
            if (i > iterations/2) {
                for (int j = 0; j < N; ++j) {
                    double v2 = 0;
                    for (int k = 0; k < 3; ++k) {
                        v2 += (float)v[j][k] * (float)v[j][k];
                    }
                    fprintf(f_velocity, "%f\n", sqrt(v2));
                }
            }
            // Store coordinates to file
            if (i < iter_to_write) {
                fprintf(f_xyz, "%d\n\n", N);
                for (int i = 0; i < N; ++i) {
                    fprintf(f_xyz, "%c ", (char) (97+(i%26)));
                    for (int k = 0; k < 3; ++k) {
                        fprintf(f_xyz, "%f ", rn[i][k]);
                    }
                    fprintf(f_xyz, "\n");
                }
            }
        }
        
        // Print to console
        if (i % (iterations / 10) == 0) {
            printf(">  iter: %d\n", i);
            printf("-> temp: %f\n", Temp);
            printf("-> utot: %f\n", utot);
            printf("-> K: %f\n", K);
            printf("\n");
        }
        assert(K != INFINITY);
    }
    
    if (PRINT_TO_FILE) {
        fclose(f_velocity);
        fclose(f_xyz);
        fclose(f_en);
        fclose(f_poten);
        fclose(f_kin);
        fclose(f_temp);
        fclose(f_len);
        printf("Print to file: Done!\n");
    }
    
    printf("Result: Done!\n");
    
    return 0;
}
