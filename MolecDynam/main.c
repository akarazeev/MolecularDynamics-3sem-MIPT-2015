//
//  main.c
//  MolecDynam
//
//  Created by Anton Wetret on 17/09/15.
//  Copyright (c) 2015 Anton Karazeev. All rights reserved.
//

/*
 Element - Ar
 sigma = 3.40 angstrem = 3.40 * E-10 m
 epsilon = 1.040 * E-2 eV = 1.664 * E-21
 1u = 1.66 * E-27 kg
 mAr = 39.9u
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define trace(x) printf("%f\n", x);
#define dtrace(x, y) printf("%f %f\n", x, y);

const int N = 5;
const double dt = 0.001;

// Mass of 1 Ar atom
double mAr = 1;
double K = 0;
double utot = 0;
//double sigma = 1;
//double epsilon = 1;

double v[N][3];
double f[N][3];
double r[N][3];
double m[N];

double Potential(double x) {
//    if (x > 2.5) {
//        double res = (float)4 * ( (1/(float)pow(x,6)) - (1/(float)pow(x,3)) );
//        return res;
//    } else {
//        double res = -0.0163;
//        return res;
//    }
    double res = (float)4 * ( (1/(float)pow(x,6)) - (1/(float)pow(x,3)) );
    return res;
}

double ForceDevByRange(double x) {
    double res = (float)48 * ( (1/(float)pow(x,7)) - (1/((float)pow(x,4) * 2)) );
    return res;
}

void EqMotion() {
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            v[i][k] += ((float)f[i][k] * (float)dt) / (float)m[i];
            r[i][k] += (float)v[i][k] * (float)dt;
        }
    }
}

void ClearForces() {
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            f[i][k] = 0;
        }
    }
}

void CalcForces() {
    double rij[3];
    rij[0] = 0;
    rij[1] = 0;
    rij[2] = 0;
    utot = 0;
    double r2;
    for (int i = 1; i < N; ++i) {
        for (int j = 0; j < i; ++j) {
            // Squared range
            r2 = 0;
            for (int k = 0; k < 3; ++k) {
                rij[k] = r[i][k] - r[j][k];
                r2 += rij[k] * rij[k];
            }
            // Total Potential
            utot += Potential(r2);
            // Force/Range
            double f_r = ForceDevByRange(r2);
            for (int k = 0; k < 3; ++k) {
                f[i][k] += (float)f_r * (float)rij[k];
                f[j][k] -= (float)f_r * (float)rij[k];
//                trace(f[i][k])
//                trace(f[j][k])
            }
        }
    }
}

void CalcEnergy() {
    K = 0;
    double v2 = 0;
    for (int i = 0; i < N; ++i) {
        // Squared Velocity
        v2 = 0;
        for (int k = 0; k < 3; ++k) {
            v2 += (float)v[i][k] * (float)v[i][k];
        }
//        trace(v2)
        K += (float)m[i] * ((float)v2 / (float)2);
    }
}

int main() {
//    srand((unsigned int)time(NULL));
    FILE* f_en = fopen("/Users/AntonWetret/Documents/Rlang/molec_dynam/energy.csv", "w");
    FILE* f_coord0 = fopen("/Users/AntonWetret/Documents/Rlang/molec_dynam/data0.csv", "w");
    FILE* f_coord1 = fopen("/Users/AntonWetret/Documents/Rlang/molec_dynam/data1.csv", "w");
    FILE* f_coord2 = fopen("/Users/AntonWetret/Documents/Rlang/molec_dynam/data2.csv", "w");
    for (int i = 0; i < N; ++i) {
        m[i] = mAr;
    }
    for (int i = 0; i < N; ++i) {
        printf("%s\n", "Beginning coordinates");
        for (int k = 0; k < 3; ++k) {
            r[i][k] = ((float)rand()/(float)RAND_MAX) * 2;
            trace(r[i][k])
        }
    }
    for (int i = 0; i < 10000; ++i) {
        ClearForces();
        CalcForces();
        EqMotion();
        CalcEnergy();
//        for (int k = 0; k < 3; ++k) {
//            printf("%f", r[2][k]);
//            if (k != 2) {
//                printf(",");
//            }
//        }
//        printf("\n");

//        for (int k = 0; k < 3; ++k) {
//            fprintf(f_coord, "%f", r[4][k]);
//            if (k != 2) {
//                fprintf(f_coord, ",");
//            }
//        }
//        fprintf(f_coord, "\n");
        
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                if (j == 0) {
                    fprintf(f_coord0, "%f", r[j][k]);
                    if (k != 2) {
                        fprintf(f_coord0, ",");
                    }
                } else if (j == 1) {
                    fprintf(f_coord1, "%f", r[j][k]);
                    if (k != 2) {
                        fprintf(f_coord1, ",");
                    }
                } else if (j == 2) {
                    fprintf(f_coord2, "%f", r[j][k]);
                    if (k != 2) {
                        fprintf(f_coord2, ",");
                    }
                }
            }
            if (j == 0) {
                fprintf(f_coord0, "\n");
            } else if (j == 1) {
                fprintf(f_coord1, "\n");
            } else if (j == 2) {
                fprintf(f_coord2, "\n");
            }
        }
        
//        for (int j = 0; j < N; ++j) {
//            for (int k = 0; k < 3; ++k) {
//                printf("%f", r[j][k]);
//                if (k != 2) {
//                    printf(",");
//                }
//            }
//            printf("\n");
//        }
        
        dtrace(K, utot)
//        trace(K + utot)
//        printf("%f,%i\n", K + utot, i);
        fprintf(f_en, "%f,%i\n", K + utot, i);
    }
    fclose(f_en);
    fclose(f_coord0);
    fclose(f_coord1);
    fclose(f_coord2);
    return 0;
}
