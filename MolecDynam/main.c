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

const int N = 3;
const double dt = 0.1;

// Mass of 1 Ar atom
double mAr = 39.9;
double K = 0;
double momentum = 0;
double utot = 0;
double sigma = 3.4;
double epsilon = 1.04;

double v[N][3];
double f[N][3];
double r[N][3];
double m[N];

double Potential(double x) {
    //    if (x > 4 * sigma) {
    //        double res = 4 * epsilon * (pow(sigma/x, 12) - pow(sigma/x, 6));
    //        return res;
    //    } else {
    //        double res = 0;
    //        return res;
    //    }
    double res = 4 * epsilon * (pow(sigma/x, 12) - pow(sigma/x, 6));
    return res;
}

void EqMotion() {
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            v[i][k] += (f[i][k] * dt) / m[i];
            r[i][k] += v[i][k] * dt;
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
    for (int i = 1; i < N; ++i) {
        for (int j = 0; j < i; ++j) {
            double r2 = 0; // Squared range
            for (int k = 0; k < 3; ++k) {
                rij[k] = r[i][k] - r[j][k];
                r2 += rij[k] * rij[k];
            }
            double r1 = sqrt(r2);
            //            trace(r1)
            utot += Potential(r1); // Total Potential
//            dtrace(r1 ,Potential(r1))
            //            trace(utot)
            //            f_r = fr(r); // Force/Range
            //TODO: Force!!!
            double f_r = Potential(r1)/r1;
            for (int k = 0; k < 3; ++k) {
                f[i][k] += f_r * rij[k];
                f[j][k] -= f_r * rij[k];
                //                trace(f[i][k])
                //                trace(f[j][k])
            }
        }
    }
}

//FIXME: Works wrong because of Potentials
void CalcEnergy() {
    K = 0;
    double v2 = 0;
    for (int i = 0; i < N; ++i) {
        v2 = 0; // Squared Velocity
        for (int k = 0; k < 3; ++k) {
            v2 += v[i][k] * v[i][k];
        }
        //        trace(v2)
        K += m[i] * v2 / 2;
    }
}

void CalcMomentum() {
    momentum = 0;
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            momentum += m[i] * v[i][k];
        }
    }
}

int main() {
    srand(time(NULL));
    m[0] = mAr;
    m[1] = mAr;
    m[2] = mAr;
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            r[i][k] = 1 + rand() % 3;
            //            trace(r[i][k])
        }
    }
    for (int i = 0; i < 10; ++i) {
        ClearForces();
        CalcForces();
        EqMotion();
        CalcEnergy();
        CalcMomentum();
        //        for (int k = 0; k < 3; ++k) {
        //            printf("%f %f ", v[1][k], r[1][k]);
        //        }
        //        printf("\n");
        //        printf("%f %f\n", K, utot);
        dtrace(K, utot);
    }
    return 0;
}
