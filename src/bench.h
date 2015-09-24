/**********************************************************************
 * Copyright (c) 2014 Pieter Wuille                                   *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _SECP256K1_BENCH_H_
#define _SECP256K1_BENCH_H_

#include <stdio.h>
#include <math.h>
#include "sys/time.h"

#define B_SCANNING (1UL << 0)
#define B_COUNT    (1UL << 1)
#define B_ITERS    (1UL << 2)
#define B_WINDOWG  (1UL << 3)
#define B_DO_TYPE  (1UL << 4)
#define B_FINISHED (1UL << 5)

#define B_SCALAR   (1UL << 0)
#define B_FIELD    (1UL << 1)
#define B_GROUP    (1UL << 2)
#define B_ECMULT   (1UL << 3)
#define B_HASH     (1UL << 4)
#define B_CONTEXT  (1UL << 5)
#define B_ECDH     (1UL << 6)
#define B_ECDSA    (1UL << 7)
#define B_SCHNORR  (1UL << 8)
#define B_ALL      (  ~0UL  )

typedef struct {
    char  *name;
    void (*func)(void *, int);
    void (*setup)(void *);
    void (*teardown)(void *);
    int    flags;
    int    default_runs;
    int    default_iters;
    short  run;
} bench_member;

int find_bench_member(bench_member hay[], char *needle, long unsigned int type) {
    int i=0;

    for (i=0; (*(hay+i)).name != NULL; i++){
        if (needle==NULL) {
            continue;
        }
        if(!strcmp(hay[i].name, needle)) {
            if (type == 0UL) {
                return i;
            } else if (hay[i].flags & type) {  
                return i;
            }       
        }
    }
    if (needle == NULL) {
        return i;
    }
    return 0;
}

static double gettimedouble(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_usec * 0.000001 + tv.tv_sec;
}

void print_number(double x) {
    double y = x;
    int c = 0;
    if (y < 0.0) y = -y;
    while (y < 100.0) {
        y *= 10.0;
        c++;
    }
    printf("%.*f", c, x);
}

void run_benchmark(char *name, void (*benchmark)(void*, int), void (*setup)(void*), void (*teardown)(void*), void* data, int count, int iters) {
    int i;
    double min = HUGE_VAL;
    double sum = 0.0;
    double max = 0.0;
    for (i = 0; i < count; i++) {
        double begin, total;
        if (setup) setup(data);
        begin = gettimedouble();
        benchmark(data, iters);
        total = gettimedouble() - begin;
        if (teardown) teardown(data);
        if (total < min) min = total;
        if (total > max) max = total;
        sum += total;
    }
    printf("%s: min ", name);
    print_number(min * 1000000.0 / iters);
    printf("us / avg ");
    print_number((sum / count) * 1000000.0 / iters);
    printf("us / max ");
    print_number(max * 1000000.0 / iters);
    printf("us\n");
}

#endif
