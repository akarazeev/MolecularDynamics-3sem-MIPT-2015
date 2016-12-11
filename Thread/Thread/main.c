#include <stdio.h>
#include <math.h>
#include <libiomp/omp.h>

long double myfunc(long double x, int steps) {
    long double tmp = 1.0;
    long double Res = 0.0;
    
    x = -x*x;
    
    int iter = 1000;
    for (int i = 1; i < iter; ++i) {
        x *= i;
    }
    
    for (int i = iter-1; i >= 1; --i) {
        x /= i;
    }
    
    for (int i = 1; i < 1000; ++i) {
        tmp = (tmp * x) / (long double)i;
        long double tmp2 = tmp/(long double)steps;
        Res += tmp2;
    }
    
    return Res;
}

int main() {
    long double res = 1.0;
    
    int steps = 10000;
    
//    #pragma omp parallel for num_threads(6)
//    #pragma omp parallel for reduction(+: res) num_threads(15)
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < steps; ++j) {
            res += myfunc((long double)j/(long double)steps, steps);
        }
        printf("%d\n", omp_get_num_threads());
//        printf("%d\n", omp_get_thread_num());
    }

    printf("%Lf\n", res);
    return 0;
}








////
////  main.c
////  Thread
////
////  Created by Anton Wetret on 27/10/15.
////  Copyright © 2015 Anton Karazeev. All rights reserved.
////
//
//#include <libiomp/omp.h>
//#include <stdio.h>
//
//double f(double x, double m) {
//    if ((int)m % 2 == 0) {
//        double tmp = f(x, m / 2);
//        return tmp * tmp;
//    }
//    else if (m == 1) {
//        return x;
//    }
//    else if (m == 0) {
//        return 1;
//    }
//    else {
//        return f(x, m/2) * f(x, 1 + (m / 2));
//    }
//}
//
//int main() {
////    #pragma omp parallel
////    printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
////
//    
//#pragma omp parallel for
//    for (int i = 1; i < 100000; ++i) {
//#pragma omp parallel for
//        for (int j = 1; j < i; ++j) {
//            if (i % j == 0) {
//                printf("%d\n", i);
//            }
//        }
//    }
//    
//    
////    int y = 0;
////    #pragma omp parallel for reduction (+: y)
////        for (int i = 0; i < 10; ++i) {
////            ++y;
////        }
////    
////    printf("%i\n", y);
////    
////    double x;
////    #pragma omp parallel shared (x)
////    {
////    #pragma omp critical
////        {
////        x = x + 1.0;
////            printf("%f\n", x);
////        }
////    }
//
////    printf("%f", f(3, 4));
//    
//    return 0;
//}

//
//  main.c
//  Integral
//
//  Created by Anton Wetret on 31/10/15.
//  Copyright © 2015 Anton Karazeev. All rights reserved.
//