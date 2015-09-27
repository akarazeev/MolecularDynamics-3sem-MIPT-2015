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
#include <assert.h>

#define trace(x) printf("%f\n", x);
#define dtrace(x, y) printf("%f %f\n", x, y);
#define PRINT_TO_FILE 0
#define USE_CUT_OFF 0

const int N = 5;
const double dt = 0.001;
const double iterations = 10000;

double rcut = 9;
double mAr = 1;
double K = 0;
double utot = 0;
//double sigma = 1;
//double epsilon = 1;

double v[N][3];
double f[N][3];
double r[N][3];
double rn[N][3];
double L2[3];
double L[3];
double m[N];

double Potential(double x) {
    trace(x)
    if (USE_CUT_OFF) {
        if (x < rcut) {
            double res = (float)4 * ( (1/(float)pow(x,6)) - (1/(float)pow(x,3)) );
            return res;
        } else {
            double res = 0;
            return res;
        }
    } else {
        double res = (float)4 * ( (1/(float)pow(x,6)) - (1/(float)pow(x,3)) );
        return res;
    }
}

double ForceDevByRange(double x) {
    double res = (float)48 * ( (1/(float)pow(x,7)) - (1/((float)pow(x,4) * 2)) );
    return res;
}

void nearest_image() {
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            if (r[i][k] > 0) {
                //???: What is L2 ???
                rn[i][k] = fmod(r[i][k] + L2[k], L[k]) - L2[k];
            } else {
                rn[i][k] = fmod(r[i][k] - L2[k], L[k]) + L2[k];
            }
        }
    }
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
    utot = 0;
    double r2;
    for (int i = 1; i < N; ++i) {
        for (int j = 0; j < i; ++j) {
            // Squared range
            r2 = 0;
            for (int k = 0; k < 3; ++k) {
//                rij[k] = r[i][k] - r[j][k];
                rij[k] = rn[i][k] - rn[j][k];
                //???: What's the difference between L and L2 ???
                if (rij[k] > L2[k]) {
                    rij[k] -= L[k];
                } else if (rij[k] < -L2[k]) {
                    rij[k] += L[k];
                }
                r2 += rij[k] * rij[k];
            }
            // Total Potential
            utot += Potential(r2);
            // Force/Range
            double f_r = ForceDevByRange(r2);
            for (int k = 0; k < 3; ++k) {
                f[i][k] += (float)f_r * (float)rij[k];
                f[j][k] -= (float)f_r * (float)rij[k];
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
        K += (float)m[i] * ((float)v2 / (float)2);
    }
}

int main() {
//    srand((unsigned int)time(NULL));
    for (int k = 0; k < 3; ++k) {
        L[k] = 15;
        assert(L[k] > rcut);
    }
    FILE* f_en;
    FILE* f_coord0;
    FILE* f_coord1;
    FILE* f_coord2;
    if (PRINT_TO_FILE) {
        f_en = fopen("molec_dynam_r/energy.csv", "w");
        f_coord0 = fopen("molec_dynam_r/data0.csv", "w");
        f_coord1 = fopen("molec_dynam_r/data1.csv", "w");
        f_coord2 = fopen("molec_dynam_r/data2.csv", "w");
    }
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
    for (int i = 0; i < iterations; ++i) {
        nearest_image();
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
        
        if (PRINT_TO_FILE) {
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
        if (PRINT_TO_FILE) {
            fprintf(f_en, "%f,%i\n", K + utot, i);
        }
    }
    if (PRINT_TO_FILE) {
        fclose(f_en);
        fclose(f_coord0);
        fclose(f_coord1);
        fclose(f_coord2);
    }
    return 0;
}
