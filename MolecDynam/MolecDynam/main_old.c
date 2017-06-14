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

#define USE_BERENDSEN 0

const int npart = 125;
const int iterations = 3000;

// const double Temp0 = 0.1;
// const double dt = 0.001;
// const double density = 0.316;

// const double Temp0 = 1.05;
// const double dt = 0.005;
// const double density = 0.6;

// const double Temp0 = 0.728;
const double Temp0 = 0.728;
const double dt = 0.001;
const double density = 0.8442;

const double eps = 1.0;

const double tau = 1;
const double iter_to_write = iterations;

double rcut2 = 9;
const double max_vel = 1;
const double mAr = 1;

/* current values of each iteration */
double K = 0;
double Temp = 0;
double utot = 0;
double virial = 0;
double press = 0;
double volume = 0;
double length = 0;
double ecut = 0;

double v[npart][3];
double f[npart][3];
double r[npart][3];
double rm[npart][3];
double rn[npart][3];
double L2[3];
double L[3];
double m[npart];

/* to switch on/off the thermostat */
int flag = 0;

static inline double Potential(double r2) {
    if (r2 < rcut2 || rcut2 == 0) {
        double res = 4.0 * ( (1.0/(float)powf(r2, 6.0)) - (1.0/(float)powf(r2, 3.0)) );
        // double res = exp(-sqrt(r2)) - (4.0/(float)r2) + (5.0/powf(r2, 3.0/2.0)) - (10.0/powf(r2, 5.0));
        return res;
    } else {
        double res = 0;
        return res;
    }
}

static inline double ForceDivByRange(double r2) {
    double res = 48.0 * ( (1.0/(float)powf(r2, 7.0)) - (1.0/((float)powf(r2, 4.0) * 2.0)) );
    // double res = -exp(-sqrt(r2)) + (8.0/powf(r2, 3.0/2.0)) - (15.0/powf(r2, 2.0)) + (100.0/powf(r2, 11.0/2.0));
    return res;
}

static inline void NearestImage() {
    for (int i = 0; i < npart; ++i) {
        for (int k = 0; k < 3; ++k) {
            if (r[i][k] > 0) {
                rn[i][k] = fmod(r[i][k] + L2[k], L[k]) - L2[k];
            } else {
                rn[i][k] = fmod(r[i][k] - L2[k], L[k]) + L2[k];
            }
        }
    }
}

static inline void ClearForces() {
    for (int i = 0; i < npart; ++i) {
        for (int k = 0; k < 3; ++k) {
            f[i][k] = 0;
        }
    }
}

static inline void CalcForces() {
    double rij[3];
    utot = 0;
    double r2 = 0;
    virial = 0;
    for (int i = 0; i < npart-1; ++i) {
        for (int j = i+1; j < npart; ++j) {
            r2 = 0;
            for (int k = 0; k < 3; ++k) {
                rij[k] = r[i][k] - r[j][k];
                rij[k] -= (float)length * round(rij[k]/length);
                r2 += powf(rij[k], 2.0);
            }
            double f_r = 0;
            double r2i = 0;
            double r6i = 0;
            double ff = 0;
            if (r2 < rcut2) {
                r2i = 1.0/(float)r2;
                r6i = powf(r2i, 3.0);
                ff = 48.0 * r2i * r6i * (r6i - 0.5);
                utot += 4.0 * r6i * (r6i - 1.0);
                utot -= ecut;
            }
            for (int k = 0; k < 3; ++k) {
                assert(ff == ff);
                f[i][k] += (float)ff * (float)rij[k];
                f[j][k] -= (float)ff * (float)rij[k];
            }
        }
    }
    for (int i = 0; i < npart; ++i) {
        for (int k = 0; k < 3; ++k) {
            virial += (float)f[i][k] * (float)r[i][k];
        }
    }
    virial = (float)virial / (float)npart;
    // virial = (float)virial / 3.0;
    // double rij[3];
    // utot = 0;
    // double r2 = 0;
    // virial = 0;
    // for (int i = 1; i < npart; ++i) {
    //     for (int j = 0; j < i; ++j) {
    //         r2 = 0;
    //         for (int k = 0; k < 3; ++k) {
    //             rij[k] = rn[i][k] - rn[j][k];
    //             if (rij[k] > L2[k]) {
    //                 rij[k] -= L[k];
    //             } else if (rij[k] < -L2[k]) {
    //                 rij[k] += L[k];
    //             }
    //             r2 += powf(rij[k], 2.0);
    //         }
    //         double f_r = 0;
    //         utot += Potential(r2);
    //         // printf("Potential: %f\n", Potential(r2));
    //         f_r = ForceDivByRange(r2);
    //         for (int k = 0; k < 3; ++k) {
    //             assert(f_r == f_r);
    //             f[i][k] += (float)f_r * (float)rij[k];
    //             f[j][k] -= (float)f_r * (float)rij[k];
    //         }
    //     }
    // }
    // for (int i = 0; i < npart; ++i) {
    //     for (int k = 0; k < 3; ++k) {
    //         virial += (float)f[i][k] * (float)r[i][k];
    //     }
    // }
    // virial = (float)virial / (float)npart;
    // // virial = (float)virial / 3.0;
}

static inline void EqMotion() {
    double xx = 0;
    for (int i = 0; i < npart; ++i) {
        for (int k = 0; k < 3; ++k) {
            xx = (2.0 * (float)r[i][k]) - (float)rm[i][k] + ((float)powf(dt, 2.0) * ((float)f[i][k] / (float)m[i]));
            v[i][k] = ((float)xx - (float)rm[i][k]) / (2.0 * (float)dt);
            rm[i][k] = r[i][k];
            r[i][k] = xx;
        }
    }
}

static inline void CalcEnergy() {
    K = 0;
    double v2 = 0;
    for (int i = 0; i < npart; ++i) {
        v2 = 0;
        for (int k = 0; k < 3; ++k) {
            v2 += powf((float)v[i][k], 2.0);
        }
        K += (float)m[i] * ((float)v2 / 2.0);
        assert(K != INFINITY);
    }
}

static inline void CalcTempPress() {
    Temp = (2.0 * K) / (npart * 3.0);
    // Calculate Pressure
    // press = ((float)density * (float)Temp) + ((float)virial / (3.0 * (float)volume));
    press = ((float)density * (float)Temp) + ((float)virial / ((float)volume));
}

static inline void Thermostat(int cur_iter) {
    if (USE_BERENDSEN) {
        double lambda = sqrt(1 + ( (dt/tau) * ((Temp0/Temp) - 1) ));
        if (fabs(Temp - Temp0) < 0.1 && cur_iter > 1000) {
            // flag = 1;
        }
        if (!flag) {
            // puts("lol");
            for (int i = 0; i < npart; ++i) {
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
    length = powf((float)npart/(float)density, 1.0/3.0);
    volume = powf(length, 3.0);
    rcut2 = powf(length / 2.0, 2.0);
    double rcut2i = 1.0 / rcut2;
    ecut = 4.0 * (powf(rcut2i, 6) - powf(rcut2i, 3));
    int width = 11;
    puts("+--------------------");
    printf("| %-*s %d\n", width, "npart:", npart);
    printf("| %-*s %d\n", width, "iterations:", iterations);
    printf("| %-*s %f\n", width, "density:", density);
    printf("| %-*s %f\n", width, "length:", length);
    printf("| %-*s %f\n", width, "rcut1:", sqrt(rcut2));
    // Set boundaries
    for (int k = 0; k < 3; ++k) {
        L[k] = length;
        L2[k] = L[k]/2.0;
        if (rcut2 > 0) {
            assert(L[k] > sqrt(rcut2));
        }
    }
    // Set masses for every particle
    for (int i = 0; i < npart; ++i) {
        m[i] = mAr;
    }

    FILE* f_en;
    FILE* f_kin;
    FILE* f_poten;
    FILE* f_temp;
    FILE* f_xyz;
    FILE* f_vel;
    FILE* f_init_vel;
    FILE* f_init_coord;
    FILE* f_init_vel_nosq;
    FILE* f_len;
    FILE* f_press;
    FILE* f_vir;

    f_init_coord = fopen("data/init_coord.xyz", "w");
    f_xyz = fopen("data/coordinates.xyz", "w");
    f_temp = fopen("data/temp.csv", "w");
    f_vel = fopen("data/velocity.csv", "w");
    f_init_vel = fopen("data/init_velocity.csv", "w");
    f_init_vel_nosq = fopen("data/init_velocity_nosq.csv", "w");
    f_en = fopen("data/energy.csv", "w");
    f_kin = fopen("data/kinetic.csv", "w");
    f_poten = fopen("data/poten.csv", "w");
    f_len = fopen("data/len.csv", "w");
    f_press = fopen("data/press.csv", "w");
    f_vir = fopen("data/vir.csv", "w");

    // Print length to file
    fprintf(f_len, "%f", length);

    // Print Initial Coordinates ot file
    double init[3];
    // Quantity of atoms per line
    int quant = powf(npart, 1.0/3.0);
    // Step of the mesh
    double step = L[0] / quant;
    printf("| %-*s %f\n", width, "step:", step);
    puts("+--------------------");
    double ic = (step/2.0) - L2[0];

    for (int i = 0; i < 3; ++i) {
        init[i] = ic;
    }

    // Set initial coordinates
    for (int i = 0; i < npart; ++i) {
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
    // Print initial coordinates to file
    fprintf(f_init_coord, "%d\n\n", npart);
    for (int i = 0; i < npart; ++i) {
        fprintf(f_init_coord, "%c ", (char) (97+(i%26)));
        for (int k = 0; k < 3; ++k) {
            fprintf(f_init_coord, "%f ", r[i][k]);
        }
        fprintf(f_init_coord, "\n");
    }
    fclose(f_init_coord);

    // Set initial velocities and Print initial velocities to file
    double sumv = 0;
    double sumv2 = 0;
    for (int i = 0; i < npart; ++i) {
        for (int k = 0; k < 3; ++k) {
            v[i][k] = (1 - (((float)rand()/(float)RAND_MAX) * 2)) * max_vel;
            sumv2 += powf(v[i][k], 2.0);
            sumv += v[i][k];
        }
    }
    sumv = sumv / (3.0 * (float)npart);
    sumv2 = sumv2 / ((float)npart);
    double scale_factor = sqrt(3.0 * Temp0 / (float)sumv2);
    sumv2 = 0;
    for (int i = 0; i < npart; ++i) {
        double v2 = 0;
        for (int k = 0; k < 3; ++k) {
            v[i][k] = (v[i][k] - sumv) * scale_factor;
            rm[i][k] = r[i][k] - (v[i][k] * dt);
            v2 += powf((float)v[i][k], 2);
            fprintf(f_init_vel_nosq, "%f\n", v[i][k]);
        }
        sumv2 += v2;
        fprintf(f_init_vel, "%f\n", v2);
    }
    printf("v2: %f\n", sumv2);
    printf("temp: %f\n", sumv2/(npart * 3.0));
    fclose(f_init_vel_nosq);
    fclose(f_init_vel);

    ///////////////
    // Main Part //
    ///////////////

    for (int i = 0; i < iterations; ++i) {
        // NearestImage();
        ClearForces();
        CalcForces();
        EqMotion();
        CalcEnergy();
        CalcTempPress();
        // Thermostat(i);
        // Print energy and temperature to file
        fprintf(f_en, "%f,%i\n", (float)(K + utot)/(float)npart, i);
        fprintf(f_kin, "%f,%i\n", (float)K/(float)npart, i);
        fprintf(f_poten, "%f,%i\n", (float)utot/(float)npart, i);
        fprintf(f_temp, "%f,%i\n", Temp, i);
        fprintf(f_press, "%f,%i\n", press, i);
        fprintf(f_vir, "%f,%i\n", virial, i);
        // Print velocities to file
        if (i > iterations - 11) {
            for (int j = 0; j < npart; ++j) {
                double v2 = 0;
                for (int k = 0; k < 3; ++k) {
                    v2 += (float)v[j][k] * (float)v[j][k];
                }
                fprintf(f_vel, "%f\n", sqrt(v2));
            }
        }
        // Print coordinates to file
        if (i < iter_to_write) {
            fprintf(f_xyz, "%d\n\n", npart);
            for (int i = 0; i < npart; ++i) {
                if (i < npart/3.0) {
                    fprintf(f_xyz, "%c ", (char) 122);
                } else if (i > 2.0*npart/3.0) {
                    fprintf(f_xyz, "%c ", (char) 101);
                } else {
                    fprintf(f_xyz, "%c ", (char) 114);
                }
                for (int k = 0; k < 3; ++k) {
                    fprintf(f_xyz, "%f ", rn[i][k]);
                }
                fprintf(f_xyz, "\n");
            }
        }

        // Print to console
        if (i % (iterations / 5) == 0) {
            printf(">  iter: %d\n", i);
            printf("-> temp: %f\n", Temp);
            printf("-> utot: %f\n", utot);
            printf("-> K: %f\n", K);
            printf("-> press: %f\n", press);
            printf("\n");
        }
        assert(K != INFINITY);
    }

    fclose(f_vel);
    fclose(f_xyz);
    fclose(f_en);
    fclose(f_poten);
    fclose(f_kin);
    fclose(f_temp);
    fclose(f_len);
    fclose(f_press);
    fclose(f_vir);

    printf("Print to file: Done!\n");

    printf("Result: Done!\n");

    return 0;
}
