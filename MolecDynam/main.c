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

#define PRINT_TO_FILE 1
#define USE_CUT_OFF 1
#define READ_INIT 1

const int N = 15;
const double dt = 0.001;
const double iterations = 10000;
double length = 10.0;

double rcut2 = 9;
double mAr = 1;
double K = 0;
double utot = 0;

double v[N][3];
double f[N][3];
double r[N][3];
double rn[N][3];
double L2[3];
double L[3];
double m[N];

double Potential(double x) {
    if (USE_CUT_OFF) {
        if (x < rcut2) {
            double res = 4.0 * ( (1/(float)pow(x,6)) - (1/(float)pow(x,3)) );
            return res;
        } else {
            double res = 0;
            return res;
        }
    } else {
        double res = 4.0 * ( (1/(float)pow(x,6)) - (1/(float)pow(x,3)) );
        return res;
    }
}

double ForceDevByRange(double x) {
    double res = 48.0 * ( (1/(float)pow(x,7)) - (1/((float)pow(x,4) * 2)) );
    return res;
}

void nearest_image() {
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

void EqMotion() {
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
            if (r2 != 0) {
                utot += Potential(r2);
                f_r = ForceDevByRange(r2);
            } else {
                utot += Potential(1);
                f_r = ForceDevByRange(1);
            }
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

void CalcEnergy() {
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

void CalcTemp() {
    
}

int main() {
    srand((unsigned int)time(NULL));
    for (int k = 0; k < 3; ++k) {
        L[k] = length;
        L2[k] = L[k]/2.0;
        assert(L[k] > sqrt(rcut2));
    }
    FILE* f_en;
    FILE* f_coord0;
    FILE* f_coord1;
    FILE* f_coord2;
    FILE* f_xyz;
    FILE* f_init_coord;
    FILE* f_init_coord_r;
    if (PRINT_TO_FILE || READ_INIT) {
        f_init_coord = fopen("MolecDynam/init_coord.xyz", "w");
        f_init_coord_r = fopen("MolecDynam/init_coord_r.xyz", "r");
        f_xyz = fopen("MolecDynam/t.xyz", "w");
        f_en = fopen("molec_dynam_r/energy.csv", "w");
        f_coord0 = fopen("molec_dynam_r/data0.csv", "w");
        f_coord1 = fopen("molec_dynam_r/data1.csv", "w");
        f_coord2 = fopen("molec_dynam_r/data2.csv", "w");
    }
    for (int i = 0; i < N; ++i) {
        m[i] = mAr;
    }
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            r[i][k] = ((float)rand()/(float)RAND_MAX) * 4;
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
            // Initial Coordinates
            fprintf(f_init_coord, "%d\n\n", N);
            for (int i = 0; i < N; ++i) {
                fprintf(f_init_coord, "%c ", (char) (97+i));
                for (int k = 0; k < 3; ++k) {
                    fprintf(f_init_coord, "%f ", r[i][k]);
                }
                fprintf(f_init_coord, "\n");
            }
        }
    }
    
    for (int i = 0; i < iterations; ++i) {
        nearest_image();
        ClearForces();
        CalcForces();
        EqMotion();
        CalcEnergy();
        
        if (PRINT_TO_FILE) {
            if (i > iterations-1000) {
                fprintf(f_xyz, "%d\n\n", N);
                for (int i = 0; i < N; ++i) {
                    fprintf(f_xyz, "%c ", (char) (97+i));
                    for (int k = 0; k < 3; ++k) {
                        fprintf(f_xyz, "%f ", rn[i][k]);
                    }
                    fprintf(f_xyz, "\n");
                }
            }
            for (int j = 2; j < 5; ++j) {
                for (int k = 0; k < 3; ++k) {
                    if (j == 2) {
                        fprintf(f_coord0, "%f", rn[j][k]);
                        if (k != 2) {
                            fprintf(f_coord0, ",");
                        }
                    } else if (j == 3) {
                        fprintf(f_coord1, "%f", rn[j][k]);
                        if (k != 2) {
                            fprintf(f_coord1, ",");
                        }
                    } else if (j == 4) {
                        fprintf(f_coord2, "%f", rn[j][k]);
                        if (k != 2) {
                            fprintf(f_coord2, ",");
                        }
                    }
                }
                if (j == 2) {
                    fprintf(f_coord0, "\n");
                } else if (j == 3) {
                    fprintf(f_coord1, "\n");
                } else if (j == 4) {
                    fprintf(f_coord2, "\n");
                }
            }
            fprintf(f_en, "%f,%i\n", K + utot, i);
        }
        
        dtrace(K, utot)
        assert(K != INFINITY);
    }
    if (PRINT_TO_FILE) {
        fclose(f_xyz);
        fclose(f_en);
        fclose(f_coord0);
        fclose(f_coord1);
        fclose(f_coord2);
        printf("Print to file: Done!\n");
    }
    printf("Result: Done!\n");
    return 0;
}
