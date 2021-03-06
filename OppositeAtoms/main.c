//
//  main.c
//  OppositeAtoms
//
//  Created by Anton Wetret on 13/10/15.
//  Copyright © 2015 Anton Karazeev. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#define trace(x) printf("%f\n", x);
#define dtrace(x, y) printf("%f %f\n", x, y);

#define PRINT_TO_FILE 1
#define USE_BERENDSEN 0

const double length = 30;
const double Temp0 = 103;
const int N = 2;
const double dt = 0.0006;
const double iterations = 10000;
const double density = 0.1;

double rcut2 = 9;
double mAr = 1;
double eps_dev_by_k_boltz = 121;
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

double Potential(double x) {
//    if (x < rcut2) {
        double res = 4.0 * ( (1.0/(float)powf(x,6.0)) - (1.0/(float)powf(x,3.0)) );
        return res;
//    } else {
//        double res = 0;
//        return res;
//    }
}

double ForceDevByRange(double x) {
    double res = -48.0 * ( (1.0/(float)powf(x,7.0)) - (1.0/((float)powf(x,4.0) * 2.0)) );
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
//            assert(r[i][k] == r[i][k]);
        }
    }
    r[0][0] = 0;
    r[0][1] = 0;
    r[0][2] = 0;
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
            //            trace(r2)
//            if (r2 != 0) {
                utot += Potential(r2);
                f_r = ForceDevByRange(r2);
//            } else {
//                utot += Potential(1.0);
//                f_r = ForceDevByRange(1.0);
//            }
            for (int k = 0; k < 3; ++k) {
                //                printf("F/R %f\n", f_r);
                //                printf("Rij %f\n", rij[k]);
//                assert(f_r == f_r);
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
    Temp = ((K / N) / 3.0) * (2.0 * eps_dev_by_k_boltz);
}

void Thermostat() {
    if (USE_BERENDSEN) {
        double tau = 1;
        double lambda = sqrt(1 + ( (dt/tau) * ((Temp0/Temp) - 1) ));
        if (fabs(Temp - Temp0) < 0.1) {
            flag = 1;
        }
        if (!flag) {
            //            trace(lambda)
            for (int i = 0; i < N; ++i) {
                for (int k = 0; k < 3; ++k) {
                    v[i][k] *= lambda;
                }
            }
        }
    }
}

int main() {
    for (int k = 0; k < 3; ++k) {
        L[k] = length;
        L2[k] = L[k]/2.0;
//        assert(L[k] > sqrt(rcut2));
    }

    m[0] = 1000000;
    m[1] = 1;

//    for (int i = 0; i < N; ++i) {
//        m[i] = mAr;
//    }
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
        fprintf(f_len, "%f", length);
        fclose(f_len);
    }
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            r[i][k] = 0;
        }
    }

    r[1][0] = 1.2;

    // Make Initial Velocities
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            v[i][k] = 0;
        }
    }

    v[1][0] = -100;
    v[1][1] = -66.876615;

    if (PRINT_TO_FILE ) {
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

    for (int i = 0; i < iterations; ++i) {
        nearest_image();
        ClearForces();
        CalcForces();
        EqMotion();
//        CalcEnergy();
        CalcTemp();
        Thermostat();

        if (PRINT_TO_FILE) {
            if (i == iterations/2) {
                for (int j = 0; j < N; ++j) {
                    double v2 = 0;
                    for (int k = 0; k < 3; ++k) {
                        v2 += (float)v[j][k] * (float)v[j][k];
                    }
                    fprintf(f_velocity, "%f\n", sqrt(v2));
                }
                fclose(f_velocity);
            }
            //            if (i > iterations-500) {
            if (i < 300) {
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
        fclose(f_xyz);
        fclose(f_en);
        fclose(f_poten);
        fclose(f_kin);
        fclose(f_temp);
        printf("Print to file: Done!\n");
    }
    printf("Result: Done!\n");
    return 0;
}
