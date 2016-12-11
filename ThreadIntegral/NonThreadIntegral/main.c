//
//  main.c
//  ThreadIntegral
//
//  Created by Anton Wetret on 29/11/15.
//  Copyright © 2015 Anton Karazeev. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int     n, thread_numb, int_pro_thr;
double  a, b, step;

double  total = 0.0;

double MyFunc(double x) {
    double tmp = 1.0;
    double res = 1.0;
    for (int i = 1; tmp > 0.000001; ++i) {
        tmp = (tmp * x)/(double)i;
        res += tmp;
    }
    return res;
}

double IntegralTmp(double tmp_a, double tmp_b, int int_pro_thr, double step) {
    double res, x;
    
    res = (MyFunc(tmp_a) + MyFunc(tmp_b))/2.0;
    x = tmp_a;
    for (int i = 1; i < int_pro_thr; ++i) {
        x = tmp_a + (i * step);
        res += MyFunc(x);
    }
    res = res * step;
    
    return res;
}

void Calc(int rank) {
    double tmp_a, tmp_b;
    
    tmp_a = a + ((int)rank * int_pro_thr * step);
    tmp_b = tmp_a + (int_pro_thr * step);
    
    total += IntegralTmp(tmp_a, tmp_b, int_pro_thr, step);
}

int main() {
//    scanf("%d %d %lf %lf", &thread_numb, &n, &a, &b);
    thread_numb = 10000;
    n = 10000;
    a = 52;
    b = 56;
    assert(a >= 0 && b >= 0);
    
    /* Legth of intervals */
    step = (b-a)/n;
    
    /* Intervals pro thread */
    int_pro_thr = n/thread_numb;
    
    for (int i = 0; i < thread_numb; ++i) {
        Calc(i);
    }
    
    printf("Integral: %f\n", total);
    
    return 0;
}
