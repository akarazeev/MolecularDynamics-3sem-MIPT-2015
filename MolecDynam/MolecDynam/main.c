//
//  main.c
//  MolecDynam
//
//  Created by Anton Karazeev on 23/01/17.
//  Copyright (c) 2017 Anton Karazeev. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

const int npart = 512;
const int iterations = 10000;
const double iter_to_write = iterations;

const int N = 64;
const int iterations = 10000;

const double Temp0 = 3;
const double tau = 1;
const double dt = 0.001;
const double iter_to_write = iterations;
const double density = 0.6;

const double max_vel = 1;
const double mAr = 1;
const int nhis = 50; // number of bins
double g[nhis];
int ngr;

double rcut2 = -1; // init in main()

/* current values of each iteration */
double K = 0;
double Temp = 0;
double etot = 0;
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

static inline void CalcForces() {
    for (int i = 0; i < npart; ++i) {
        for (int k = 0; k < 3; ++k) {
            f[i][k] = 0;
        }
    }
    double rij[3];
    etot = 0;
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
            double r2i = 0;
            double r6i = 0;
            double ff = 0;
            if (r2 < rcut2) {
                r2i = 1.0/(float)r2;
                r6i = powf(r2i, 3.0);
                ff = 48.0 * r2i * r6i * (r6i - 0.5);
                virial += r2 * ff;
                etot += 4.0 * r6i * (r6i - 1.0);
                etot -= ecut;
            }
            for (int k = 0; k < 3; ++k) {
                assert(ff == ff);
                f[i][k] += (float)ff * (float)rij[k];
                f[j][k] -= (float)ff * (float)rij[k];
            }
        }
    }
}

static inline void EqMotion() {
    // double xx = 0;
    // double sumv2 = 0;
    // for (int i = 0; i < npart; ++i) {
    //     for (int k = 0; k < 3; ++k) {
    //         xx = (2.0 * (float)r[i][k]) - (float)rm[i][k] + ((float)powf(dt, 2.0) * ((float)f[i][k] / (float)m[i]));
    //         v[i][k] = ((float)xx - (float)rm[i][k]) / (2.0 * (float)dt);
    //         sumv2 += powf(v[i][k], 2.0);
    //     }
    // }
    // Temp = (float)sumv2 / ((float)npart * 3.0);
    // double alpha = sqrt((float)Temp / (float)Temp0);
    // sumv2 = 0;
    // for (int i = 0; i < npart; ++i) {
    //     for (int k = 0; k < 3; ++k) {
    //         v[i][k] = (float)v[i][k] / (float)alpha;
    //         xx = (2.0 * dt * v[i][k]) + rm[i][k];
    //         sumv2 += powf(v[i][k], 2.0);
    //         rm[i][k] = r[i][k];
    //         r[i][k] = xx;
    //     }
    // }
    // Temp = (float)sumv2 / ((float)npart * 3.0);
    // etot = (float)(etot + ((float)sumv2 / 2.0)) / (float)npart;
    // K = (float)sumv2 / ((float)npart * 2.0);
    // ^ constant Temperature
    double xx = 0;
    double sumv2 = 0;
    for (int i = 0; i < npart; ++i) {
        for (int k = 0; k < 3; ++k) {
            xx = (2.0 * (float)r[i][k]) - (float)rm[i][k] + ((float)powf(dt, 2.0) * ((float)f[i][k] / (float)m[i]));
            v[i][k] = ((float)xx - (float)rm[i][k]) / (2.0 * (float)dt);
            sumv2 += powf(v[i][k], 2.0);
            rm[i][k] = r[i][k];
            r[i][k] = xx;
        }
    }
    Temp = (float)sumv2 / ((float)npart * 3.0);
    etot = (float)(etot + ((float)sumv2 / 2.0)) / (float)npart;
    K = (float)sumv2 / ((float)npart * 2.0);
}

static inline void CalcPress() {
    virial = (float)virial / (float)npart;

    // press = ((float)density * (float)Temp) + ((float)virial / (3.0 * (float)volume));
    press = (density * (Temp + virial)) / 3.0;
    // printf("1: %f\n", (float)density * (float)Temp);
    // printf("2: %f\n", ((float)virial / (3.0 * (float)volume)));
    // printf("press: %f\n", press);
}

static inline void RadDistr() {

    double delg = (float)length / (2.0 * (float)nhis);

    ngr += 1;
    double rij[3];
    double r1 = 0;
    double r2 = 0;
    int ig = 0;
    for (int i = 0; i < npart-1; ++i) {
        for (int j = i+1; j < npart; ++j) {
            r2 = 0;
            for (int k = 0; k < 3; ++k) {
                rij[k] = r[i][k] - r[j][k];
                rij[k] -= (float)length * round(rij[k]/length);
                r2 += powf(rij[k], 2.0);
            }
            r1 = sqrt(r2);
            if (r1 < (float)length/2.0) {
                ig = round((float)r1/(float)delg);
                g[ig] += 2;
            }
        }
    }
    double vb;
    double nid;
    for (int i = 0; i < nhis; ++i) {
        r1 = delg * (i + 0.5);
        vb = ( powf((i+1), 3) - powf(i, 3)) * powf(delg, 3);
        nid = (4.0/3.0) * 3.14159 * vb * density;
        g[i] = (float)g[i] / (float)(ngr * npart * nid);
    }
}

int main(int argc, char** argv) {
    srand((unsigned int)time(NULL));
    // Calculate the length of the box side
    length = powf((float)npart/(float)density, 1.0/3.0);
    volume = powf((float)length, 3.0);
    if (rcut2 <= 0) {
        rcut2 = 0.9 * powf((float)length / 2.0, 2.0);
    }
    double rcut2i = 1.0 / rcut2;
    ecut = 4.0 * (powf(rcut2i, 6) - powf(rcut2i, 3));

    ngr = 0;
    for (int k = 0; k < nhis; ++k) {
        g[k] = 0;
    }

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
    FILE* f_raddistr;
    FILE* f_dens;

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
    f_dens = fopen("data/rho.csv", "w");
    f_raddistr = fopen("data/raddistr.csv", "w");
  
    // Print length to file
    fprintf(f_len, "%f", length);
    fclose(f_len);
    fprintf(f_dens, "%f", density);
    fclose(f_dens);


    /* Set Initial Coordinates */
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
    // double sumv = 0;
    double sumv[3];
    for (int k = 0; k < 3; ++k) {
        sumv[k] = 0;
    }
    double sumv2 = 0;
    for (int i = 0; i < npart; ++i) {
        for (int k = 0; k < 3; ++k) {
            v[i][k] = ((float)rand() / (float)RAND_MAX) - 0.5;
            sumv2 += powf(v[i][k], 2.0);
            sumv[k] += v[i][k];
        }
    }
    for (int k = 0; k < 3; ++k) {
        sumv[k] = (float)sumv[k] / (float)npart;
    }
    sumv2 = sumv2 / (float)npart;

    double scale_factor = sqrt(3.0 * Temp0 / (float)sumv2);
    sumv2 = 0;
    for (int i = 0; i < npart; ++i) {
        double v2 = 0;
        for (int k = 0; k < 3; ++k) {
            v[i][k] = (v[i][k] - sumv[k]) * scale_factor;
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

    /////////////////////////////////////////////////////////////////
    ////////////////////// Main Part ////////////////////////////////
    /////////////////////////////////////////////////////////////////

    for (int i = 0; i < iterations; ++i) {
        NearestImage();
        CalcForces();
        EqMotion();
        CalcPress();
        // Print energy and temperature to file
        fprintf(f_en, "%f,%i\n", etot, i);
        fprintf(f_kin, "%f,%i\n", K, i);
        fprintf(f_poten, "%f,%i\n", etot - K, i);
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
            printf("-> etot: %f\n", etot);
            printf("-> K: %f\n", K);
            printf("-> press: %f\n", press);
            printf("\n");
        }
        assert(K != INFINITY);
    }

    RadDistr();
    for (int i = 0; i < nhis; ++i) {
        fprintf(f_raddistr, "%f,%i\n", g[i], i);
    }

    fclose(f_raddistr);
    fclose(f_vel);
    fclose(f_xyz);
    fclose(f_en);
    fclose(f_poten);
    fclose(f_kin);
    fclose(f_temp);
    fclose(f_press);
    fclose(f_vir);

    printf("Result: Done!\n");

    return 0;
}
