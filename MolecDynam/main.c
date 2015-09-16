//
//  main.c
//  MolecDynam
//
//  Created by Anton Wetret on 17/09/15.
//  Copyright (c) 2015 Anton Karazeev. All rights reserved.
//

#include <stdio.h>

const int N = 3;
const double dt = 0.1;

double v[N][3];
double f[N][3];
double r[N][3];
double m[N];

void EqMotion() {
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) {
            v[i][k] += f[i][k] * dt / m[i];
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
    double utot = 0;
    for (int i = 1; i < N; ++i) { // Why we start from 1 ?
        for (int j = 0; j < i; ++j) {
            double r2 = 0; // Squared range
            for (int k = 0; k < 3; ++k) {
                rij[k] = r[i][k] - r[j][k];
                r2 += rij[k] * rij[k];
            }
            double r1 = sqrt(r2);
//            utot += u(r1); // Potential
//            f_r = fr(r); // Force/Range
            for (int k = 0; k < 3; ++k) {
//                f[i][k] += f_r * rij[k];
//                f[j][k] -= f_r * rij[k];
            }
        }
    }
}

int main() {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 3; ++j) {
            printf("%f ", v[i][j]);
        }
        printf("\n");
    }
    for (int i = 0; i < 100; ++i) {
        ClearForces();
        CalcForces();
        EqMotion();
    }
    return 0;
}
