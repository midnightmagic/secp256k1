/**********************************************************************
 * Copyright (c) 2014-2015 Pieter Wuille                              *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/
#include <stdio.h>

#include "include/secp256k1.h"

#include "util.h"
#include "hash_impl.h"
#include "num_impl.h"
#include "field_impl.h"
#include "group_impl.h"
#include "scalar_impl.h"
#include "ecmult_const_impl.h"
#include "ecmult_impl.h"
#include "bench.h"
#include "secp256k1.c"

typedef struct {
    secp256k1_scalar_t scalar_x, scalar_y;
    secp256k1_fe_t fe_x, fe_y;
    secp256k1_ge_t ge_x, ge_y;
    secp256k1_gej_t gej_x, gej_y;
    unsigned char data[64];
    int wnaf[256];
} bench_inv_t;

void bench_setup(void* arg) {
    bench_inv_t *data = (bench_inv_t*)arg;

    static const unsigned char init_x[32] = {
        0x02, 0x03, 0x05, 0x07, 0x0b, 0x0d, 0x11, 0x13,
        0x17, 0x1d, 0x1f, 0x25, 0x29, 0x2b, 0x2f, 0x35,
        0x3b, 0x3d, 0x43, 0x47, 0x49, 0x4f, 0x53, 0x59,
        0x61, 0x65, 0x67, 0x6b, 0x6d, 0x71, 0x7f, 0x83
    };

    static const unsigned char init_y[32] = {
        0x82, 0x83, 0x85, 0x87, 0x8b, 0x8d, 0x81, 0x83,
        0x97, 0xad, 0xaf, 0xb5, 0xb9, 0xbb, 0xbf, 0xc5,
        0xdb, 0xdd, 0xe3, 0xe7, 0xe9, 0xef, 0xf3, 0xf9,
        0x11, 0x15, 0x17, 0x1b, 0x1d, 0xb1, 0xbf, 0xd3
    };

    secp256k1_scalar_set_b32(&data->scalar_x, init_x, NULL);
    secp256k1_scalar_set_b32(&data->scalar_y, init_y, NULL);
    secp256k1_fe_set_b32(&data->fe_x, init_x);
    secp256k1_fe_set_b32(&data->fe_y, init_y);
    CHECK(secp256k1_ge_set_xo_var(&data->ge_x, &data->fe_x, 0));
    CHECK(secp256k1_ge_set_xo_var(&data->ge_y, &data->fe_y, 1));
    secp256k1_gej_set_ge(&data->gej_x, &data->ge_x);
    secp256k1_gej_set_ge(&data->gej_y, &data->ge_y);
    memcpy(data->data, init_x, 32);
    memcpy(data->data + 32, init_y, 32);
}

void bench_scalar_add(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 2000000; i++) {
        secp256k1_scalar_add(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}

void bench_scalar_negate(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 2000000; i++) {
        secp256k1_scalar_negate(&data->scalar_x, &data->scalar_x);
    }
}

void bench_scalar_sqr(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 200000; i++) {
        secp256k1_scalar_sqr(&data->scalar_x, &data->scalar_x);
    }
}

void bench_scalar_mul(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 200000; i++) {
        secp256k1_scalar_mul(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}

#ifdef USE_ENDOMORPHISM
void bench_scalar_split(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 20000; i++) {
        secp256k1_scalar_t l, r;
        secp256k1_scalar_split_lambda(&l, &r, &data->scalar_x);
        secp256k1_scalar_add(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}
#endif

void bench_scalar_inverse(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 2000; i++) {
        secp256k1_scalar_inverse(&data->scalar_x, &data->scalar_x);
        secp256k1_scalar_add(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}

void bench_scalar_inverse_var(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 2000; i++) {
        secp256k1_scalar_inverse_var(&data->scalar_x, &data->scalar_x);
        secp256k1_scalar_add(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}

void bench_field_normalize(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 2000000; i++) {
        secp256k1_fe_normalize(&data->fe_x);
    }
}

void bench_field_normalize_weak(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 2000000; i++) {
        secp256k1_fe_normalize_weak(&data->fe_x);
    }
}

void bench_field_mul(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 200000; i++) {
        secp256k1_fe_mul(&data->fe_x, &data->fe_x, &data->fe_y);
    }
}

void bench_field_sqr(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 200000; i++) {
        secp256k1_fe_sqr(&data->fe_x, &data->fe_x);
    }
}

void bench_field_inverse(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 20000; i++) {
        secp256k1_fe_inv(&data->fe_x, &data->fe_x);
        secp256k1_fe_add(&data->fe_x, &data->fe_y);
    }
}

void bench_field_inverse_var(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 20000; i++) {
        secp256k1_fe_inv_var(&data->fe_x, &data->fe_x);
        secp256k1_fe_add(&data->fe_x, &data->fe_y);
    }
}

void bench_field_sqrt_var(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 20000; i++) {
        secp256k1_fe_sqrt_var(&data->fe_x, &data->fe_x);
        secp256k1_fe_add(&data->fe_x, &data->fe_y);
    }
}

void bench_group_double_var(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 200000; i++) {
        secp256k1_gej_double_var(&data->gej_x, &data->gej_x, NULL);
    }
}

void bench_group_add_var(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 200000; i++) {
        secp256k1_gej_add_var(&data->gej_x, &data->gej_x, &data->gej_y, NULL);
    }
}

void bench_group_add_affine(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 200000; i++) {
        secp256k1_gej_add_ge(&data->gej_x, &data->gej_x, &data->ge_y);
    }
}

void bench_group_add_affine_var(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 200000; i++) {
        secp256k1_gej_add_ge_var(&data->gej_x, &data->gej_x, &data->ge_y, NULL);
    }
}

void bench_ecmult_wnaf(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 20000; i++) {
        secp256k1_ecmult_wnaf(data->wnaf, 256, &data->scalar_x, WINDOW_A);
        secp256k1_scalar_add(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}

void bench_wnaf_const(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;

    for (i = 0; i < 20000; i++) {
        secp256k1_wnaf_const(data->wnaf, data->scalar_x, WINDOW_A);
        secp256k1_scalar_add(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}


void bench_sha256(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;
    secp256k1_sha256_t sha;

    for (i = 0; i < 20000; i++) {
        secp256k1_sha256_initialize(&sha);
        secp256k1_sha256_write(&sha, data->data, 32);
        secp256k1_sha256_finalize(&sha, data->data);
    }
}

void bench_hmac_sha256(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;
    secp256k1_hmac_sha256_t hmac;

    for (i = 0; i < 20000; i++) {
        secp256k1_hmac_sha256_initialize(&hmac, data->data, 32);
        secp256k1_hmac_sha256_write(&hmac, data->data, 32);
        secp256k1_hmac_sha256_finalize(&hmac, data->data);
    }
}

void bench_rfc6979_hmac_sha256(void* arg) {
    int i;
    bench_inv_t *data = (bench_inv_t*)arg;
    secp256k1_rfc6979_hmac_sha256_t rng;

    for (i = 0; i < 20000; i++) {
        secp256k1_rfc6979_hmac_sha256_initialize(&rng, data->data, 64);
        secp256k1_rfc6979_hmac_sha256_generate(&rng, data->data, 32);
    }
}

void bench_context_verify(void* arg) {
    int i;
    (void)arg;
    for (i = 0; i < 20; i++) {
        secp256k1_context_destroy(secp256k1_context_create(SECP256K1_CONTEXT_VERIFY));
    }
}

void bench_context_sign(void* arg) {
    int i;
    (void)arg;
    for (i = 0; i < 200; i++) {
        secp256k1_context_destroy(secp256k1_context_create(SECP256K1_CONTEXT_SIGN));
    }
}

int find_bench_member(bench_member hay[], char *needle, long unsigned int type) {
    int i=0;
    if (needle==NULL) {
        needle="";
    }
    for (i=0; hay+i != (bench_member *)NULL; i++){
        if(hay[i].name == NULL) {
            return i;
        }
        if(!strcmp(hay[i].name, needle)) {
            return i;
        }
    }
    return i;
}

int main(int argc, char **argv) {
    bench_inv_t data;
    int i=0, runs=0, iters=0, windowg=0;
    int whichbench=0, total_benchmarks=0;
    bench_flags bflags;
    bench_mode emode;

    if (argc < 2) {
        printf("No arguments.\n");
        return -1;
    }

    for (i=1; i<argc; i++) {
        switch(emode) {
            case B_SCANNING:
                if (!strcmp(argv[i], "scalar")) {
                    bflags=B_SCALAR;
                } elseif (!strcmp(argv[i], "field")) {
                    bflags=B_FIELD;
                } elseif (!strcmp(argv[i], "group")) {
                    bflags=B_GROUP;
                } elseif (!strcmp(argv[i], "ecmult")) {
                    bflags=B_ECMULT;
                } elseif (!strcmp(argv[i], "hash")) {
                    bflags=B_HASH;
                } elseif (!strcmp(argv[i], "context")) {
                    bflags=B_CONTEXT;
                } elseif (!strcmp(argv[i], "-c")) {
                    emode=B_COUNT;
                } elseif (!strcmp(argv[i], "-i")) {
                    emode=B_ITERS;
                } elseif (!strcmp(argv[i], "-w")) {
                    emode=B_WINDOWG;
                } elseif (!strcmp(argv[i], "go")) {
                    emode=B_DO_TYPE;
                } elseif (!strcmp(argv[i], "all")) {
                    bflags=B_ALL;
                    emode=B_DO_TYPE;
                } elseif (whichbench=find_bench_member(bench_table, argv[i])) {
                    run_benchmark(bench_table[whichbench].name, bench_table[whichbench].func, bench_table[whichbench].setup, \
                        NULL, &data, runs?runs:bench_table[whichbench].default_runs, iters?iters:bench_table[whichbench].default_iters);
                } else {
                    printf("Unrecognised keyword: %s\n", argv[i]);
                    return -1;
                }
                break;
            case B_COUNT:
                runs=atoi(argv[i]);
                emode=B_SCANNING;
                break;
            case B_ITERS:
                iters=atoi(argv[i]);
                emode=B_SCANNING;
                break;
            case B_WINDOWG:
                windowg=atoi(argv[i]);
                emode=B_SCANNING;
                break;
            case B_DO_TYPE:
                totalbenchmarks=find_bench_member(bench_table, NULL);
                for (whichbench=0; whichbench<totalbenchmarks; whichbench++) {
                    if (bench_table[whichbench].flags & bflags) {
                        run_benchmark(bench_table[whichbench].name, bench_table[whichbench].func, bench_table[whichbench].setup, \
                            NULL, &data, runs?runs:bench_table[whichbench].default_runs, iters?iters:bench_table[whichbench].default_iters);
                    }
                }
                break;
            default:
                break;
        } // End of main scanning loop thing
    return 0;
}
