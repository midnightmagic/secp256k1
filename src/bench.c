/**********************************************************************
 * Copyright (c) 2014-2015 Pieter Wuille                              *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/
#include <stdio.h>

#include "include/secp256k1.h"
#include "include/secp256k1_ecdh.h"
#include "include/secp256k1_recovery.h"
#include "include/secp256k1_schnorr.h"
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

#include "modules/ecdh/main_impl.h"
#include "modules/recovery/main_impl.h"
#include "modules/schnorr/schnorr.h"
#include "modules/schnorr/main_impl.h"

typedef struct {
    unsigned char key[32];
    unsigned char sig[64];
    unsigned char pubkey[33];
    int pubkeylen;
} bench_schnorr_sig_t;

typedef struct {
    secp256k1_scalar_t scalar_x, scalar_y;
    secp256k1_fe_t fe_x, fe_y;
    secp256k1_ge_t ge_x, ge_y;
    secp256k1_gej_t gej_x, gej_y;
    unsigned char data[64];
    int wnaf[256];

    secp256k1_context_t *ctx;
    secp256k1_pubkey_t point;
    unsigned char scalar[32];

    unsigned char msg[32];
    unsigned char sig[64];

    bench_schnorr_sig_t sigs[64];
    int numsigs;

    unsigned char key[32];

    int siglen;
    unsigned char ecdsa_sig[72];
    unsigned char pubkey[33];
    int pubkeylen;

    int windowG_override;
} bench_t;

void bench_setup(void* arg) {
    bench_t *data = (bench_t*)arg;

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

void bench_scalar_add(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_scalar_add(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}

void bench_scalar_negate(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_scalar_negate(&data->scalar_x, &data->scalar_x);
    }
}

void bench_scalar_sqr(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_scalar_sqr(&data->scalar_x, &data->scalar_x);
    }
}

void bench_scalar_mul(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_scalar_mul(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}

#ifdef USE_ENDOMORPHISM
void bench_scalar_split(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_scalar_t l, r;
        secp256k1_scalar_split_lambda(&l, &r, &data->scalar_x);
        secp256k1_scalar_add(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}
#endif

void bench_scalar_inverse(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_scalar_inverse(&data->scalar_x, &data->scalar_x);
        secp256k1_scalar_add(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}

void bench_scalar_inverse_var(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_scalar_inverse_var(&data->scalar_x, &data->scalar_x);
        secp256k1_scalar_add(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}

void bench_field_normalize(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_fe_normalize(&data->fe_x);
    }
}

void bench_field_normalize_weak(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_fe_normalize_weak(&data->fe_x);
    }
}

void bench_field_mul(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_fe_mul(&data->fe_x, &data->fe_x, &data->fe_y);
    }
}

void bench_field_sqr(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_fe_sqr(&data->fe_x, &data->fe_x);
    }
}

void bench_field_inverse(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_fe_inv(&data->fe_x, &data->fe_x);
        secp256k1_fe_add(&data->fe_x, &data->fe_y);
    }
}

void bench_field_inverse_var(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_fe_inv_var(&data->fe_x, &data->fe_x);
        secp256k1_fe_add(&data->fe_x, &data->fe_y);
    }
}

void bench_field_sqrt_var(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_fe_sqrt_var(&data->fe_x, &data->fe_x);
        secp256k1_fe_add(&data->fe_x, &data->fe_y);
    }
}

void bench_group_double_var(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_gej_double_var(&data->gej_x, &data->gej_x, NULL);
    }
}

void bench_group_add_var(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_gej_add_var(&data->gej_x, &data->gej_x, &data->gej_y, NULL);
    }
}

void bench_group_add_affine(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_gej_add_ge(&data->gej_x, &data->gej_x, &data->ge_y);
    }
}

void bench_group_add_affine_var(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_gej_add_ge_var(&data->gej_x, &data->gej_x, &data->ge_y, NULL);
    }
}

void bench_ecmult_wnaf(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_ecmult_wnaf(data->wnaf, 256, &data->scalar_x, WINDOW_A);
        secp256k1_scalar_add(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}

void bench_wnaf_const(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_wnaf_const(data->wnaf, data->scalar_x, WINDOW_A);
        secp256k1_scalar_add(&data->scalar_x, &data->scalar_x, &data->scalar_y);
    }
}


void bench_sha256(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;
    secp256k1_sha256_t sha;

    for (i = 0; i < iters; i++) {
        secp256k1_sha256_initialize(&sha);
        secp256k1_sha256_write(&sha, data->data, 32);
        secp256k1_sha256_finalize(&sha, data->data);
    }
}

void bench_hmac_sha256(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;
    secp256k1_hmac_sha256_t hmac;

    for (i = 0; i < iters; i++) {
        secp256k1_hmac_sha256_initialize(&hmac, data->data, 32);
        secp256k1_hmac_sha256_write(&hmac, data->data, 32);
        secp256k1_hmac_sha256_finalize(&hmac, data->data);
    }
}

void bench_rfc6979_hmac_sha256(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;
    secp256k1_rfc6979_hmac_sha256_t rng;

    for (i = 0; i < iters; i++) {
        secp256k1_rfc6979_hmac_sha256_initialize(&rng, data->data, 64);
        secp256k1_rfc6979_hmac_sha256_generate(&rng, data->data, 32);
    }
}

void bench_context_verify(void* arg, int iters) {
    int i;
    bench_t* data = (bench_t *)arg;
    for (i = 0; i < iters; i++) {
        secp256k1_context_destroy(secp256k1_context_create(SECP256K1_CONTEXT_VERIFY | SECP256K1_CONTEXT_WINDOWG_PACK(data->windowG_override) ));
    }
}

void bench_context_sign(void* arg, int iters) {
    int i;
    (void)arg;
    for (i = 0; i < iters; i++) {
        secp256k1_context_destroy(secp256k1_context_create(SECP256K1_CONTEXT_SIGN));
    }
}

void bench_ecdh_setup(void* arg) {
    int i;
    bench_t *data = (bench_t*)arg;

    const unsigned char point[] = {
        0x03,
        0x54, 0x94, 0xc1, 0x5d, 0x32, 0x09, 0x97, 0x06,
        0xc2, 0x39, 0x5f, 0x94, 0x34, 0x87, 0x45, 0xfd,
        0x75, 0x7c, 0xe3, 0x0e, 0x4e, 0x8c, 0x90, 0xfb,
        0xa2, 0xba, 0xd1, 0x84, 0xf8, 0x83, 0xc6, 0x9f
    };

    data->ctx = secp256k1_context_create(0);
    for (i = 0; i < 32; i++) data->scalar[i] = i + 1;
    CHECK(secp256k1_ec_pubkey_parse(data->ctx, &data->point, point, sizeof(point)) == 1);
}

void bench_ecdh(void* arg, int iters) {
    int i;
    unsigned char res[32];
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        CHECK(secp256k1_ecdh(data->ctx, res, &data->point, data->scalar) == 1);
    }
}

void bench_recover(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;
    secp256k1_pubkey_t pubkey;
    unsigned char pubkeyc[33];

    for (i = 0; i < iters; i++) {
        int j;
        int pubkeylen = 33;
        secp256k1_ecdsa_recoverable_signature_t sig;
        CHECK(secp256k1_ecdsa_recoverable_signature_parse_compact(data->ctx, &sig, data->sig, i % 2));
        CHECK(secp256k1_ecdsa_recover(data->ctx, &pubkey, &sig, data->msg));
        CHECK(secp256k1_ec_pubkey_serialize(data->ctx, pubkeyc, &pubkeylen, &pubkey, 1));
        for (j = 0; j < 32; j++) {
            data->sig[j + 32] = data->msg[j];    /* Move former message to S. */
            data->msg[j] = data->sig[j];         /* Move former R to message. */
            data->sig[j] = pubkeyc[j + 1];       /* Move recovered pubkey X coordinate to R (which must be a valid X coordinate). */
        }
    }
}

void bench_recover_setup(void* arg) {
    int i;
    bench_t *data = (bench_t*)arg;

    data->ctx = secp256k1_context_create(SECP256K1_CONTEXT_VERIFY | SECP256K1_CONTEXT_WINDOWG_PACK(data->windowG_override) );

    for (i = 0; i < 32; i++) data->msg[i] = 1 + i;
    for (i = 0; i < 64; i++) data->sig[i] = 65 + i;
}

void bench_recover_teardown(void* arg) {

    bench_t *data = (bench_t*)arg;

    secp256k1_context_destroy(data->ctx);
}

static void bench_schnorr_init(void* arg) {
    int i, k;
    bench_t* data = (bench_t*)arg;

    data->ctx = secp256k1_context_create(SECP256K1_CONTEXT_SIGN | SECP256K1_CONTEXT_VERIFY | SECP256K1_CONTEXT_WINDOWG_PACK(data->windowG_override) );

    data->numsigs = 1; /* XXX: Excisable */

    for (i = 0; i < 32; i++) data->msg[i] = 1 + i;
    for (k = 0; k < data->numsigs; k++) {
        secp256k1_pubkey_t pubkey;
        for (i = 0; i < 32; i++) data->sigs[k].key[i] = 33 + i + k;
        secp256k1_schnorr_sign(data->ctx, data->sigs[k].sig, data->msg, data->sigs[k].key, NULL, NULL);
        data->sigs[k].pubkeylen = 33;
        CHECK(secp256k1_ec_pubkey_create(data->ctx, &pubkey, data->sigs[k].key));
        CHECK(secp256k1_ec_pubkey_serialize(data->ctx, data->sigs[k].pubkey, &data->sigs[k].pubkeylen, &pubkey, 1));
    }
}

void bench_schnorr_teardown(void* arg) {
    bench_t* data = (bench_t*)arg;

    secp256k1_context_destroy(data->ctx);
}

static void bench_schnorr_verify(void* arg, int iters) {
    int i;
    bench_t* data = (bench_t*)arg;

    for (i = 0; i < iters / data->numsigs; i++) {
        secp256k1_pubkey_t pubkey;
        data->sigs[0].sig[(i >> 8) % 64] ^= (i & 0xFF);
        CHECK(secp256k1_ec_pubkey_parse(data->ctx, &pubkey, data->sigs[0].pubkey, data->sigs[0].pubkeylen));
        CHECK(secp256k1_schnorr_verify(data->ctx, data->sigs[0].sig, data->msg, &pubkey) == ((i & 0xFF) == 0));
        data->sigs[0].sig[(i >> 8) % 64] ^= (i & 0xFF);
    }
}

static void bench_sign_setup(void* arg) {
    int i;
    bench_t *data = (bench_t*)arg;

    for (i = 0; i < 32; i++) data->msg[i] = i + 1;
    for (i = 0; i < 32; i++) data->key[i] = i + 65;

    data->ctx = secp256k1_context_create(SECP256K1_CONTEXT_SIGN);
}

void bench_sign_teardown(void* arg) {
    bench_t *data = (bench_t*)arg;

    secp256k1_context_destroy(data->ctx);
}

static void bench_sign(void* arg, int iters) {
    int i;
    bench_t *data = (bench_t*)arg;

    unsigned char sig[74];
    for (i = 0; i < iters; i++) {
        int siglen = 74;
        int j;
        secp256k1_ecdsa_signature_t signature;
        CHECK(secp256k1_ecdsa_sign(data->ctx, &signature, data->msg, data->key, NULL, NULL));
        CHECK(secp256k1_ecdsa_signature_serialize_der(data->ctx, sig, &siglen, &signature));
        for (j = 0; j < 32; j++) {
            data->msg[j] = sig[j];
            data->key[j] = sig[j + 32];
        }
    }
}

void bench_verify_setup(void *arg) {
    bench_t* data = (bench_t*)arg;
    int i;
    secp256k1_pubkey_t pubkey;
    secp256k1_ecdsa_signature_t sig;

    data->ctx = secp256k1_context_create(SECP256K1_CONTEXT_SIGN | SECP256K1_CONTEXT_VERIFY | SECP256K1_CONTEXT_WINDOWG_PACK(data->windowG_override) );

    for (i = 0; i < 32; i++) data->msg[i] = 1 + i;
    for (i = 0; i < 32; i++) data->key[i] = 33 + i;
    data->siglen = 72;
    CHECK(secp256k1_ecdsa_sign(data->ctx, &sig, data->msg, data->key, NULL, NULL));
    CHECK(secp256k1_ecdsa_signature_serialize_der(data->ctx, data->sig, &data->siglen, &sig));
    CHECK(secp256k1_ec_pubkey_create(data->ctx, &pubkey, data->key));
    CHECK(secp256k1_ec_pubkey_serialize(data->ctx, data->pubkey, &data->pubkeylen, &pubkey, 1) == 1);
}

void bench_verify_teardown(void* arg) {
    bench_t* data = (bench_t*)arg;

    secp256k1_context_destroy(data->ctx);
}

static void bench_verify(void* arg, int iters) {
    int i;
    bench_t* data = (bench_t*)arg;

    for (i = 0; i < iters; i++) {
        secp256k1_pubkey_t pubkey;
        secp256k1_ecdsa_signature_t sig;
        data->sig[data->siglen - 1] ^= (i & 0xFF);
        data->sig[data->siglen - 2] ^= ((i >> 8) & 0xFF);
        data->sig[data->siglen - 3] ^= ((i >> 16) & 0xFF);
        CHECK(secp256k1_ec_pubkey_parse(data->ctx, &pubkey, data->pubkey, data->pubkeylen) == 1);
        CHECK(secp256k1_ecdsa_signature_parse_der(data->ctx, &sig, data->sig, data->siglen) == 1);
        CHECK(secp256k1_ecdsa_verify(data->ctx, &sig, data->msg, &pubkey) == (i == 0));
        data->sig[data->siglen - 1] ^= (i & 0xFF);
        data->sig[data->siglen - 2] ^= ((i >> 8) & 0xFF);
        data->sig[data->siglen - 3] ^= ((i >> 16) & 0xFF);
    }
}

bench_member bench_table[] = {
    { "scalar_add",               &bench_scalar_add,           &bench_setup, NULL, B_SCALAR,  10, 2000000, 0 },
    { "scalar_negate",            &bench_scalar_negate,        &bench_setup, NULL, B_SCALAR,  10, 2000000, 0 },
    { "scalar_sqr",               &bench_scalar_sqr,           &bench_setup, NULL, B_SCALAR,  10,  200000, 0 },
    { "scalar_mul",               &bench_scalar_mul,           &bench_setup, NULL, B_SCALAR,  10,  200000, 0 },
#ifdef USE_ENDOMORPHISM
    { "scalar_split",             &bench_scalar_split,         &bench_setup, NULL, B_SCALAR,  10,   20000, 0 },
#endif
    { "scalar_inverse",           &bench_scalar_inverse,       &bench_setup, NULL, B_SCALAR,  10,    2000, 0 },
    { "scalar_inverse_var",       &bench_scalar_inverse_var,   &bench_setup, NULL, B_SCALAR,  10,    2000, 0 },
    { "field_normalize",          &bench_field_normalize,      &bench_setup, NULL, B_FIELD,   10, 2000000, 0 },
    { "field_normalize_weak",     &bench_field_normalize_weak, &bench_setup, NULL, B_FIELD,   10, 2000000, 0 },
    { "field_sqr",                &bench_field_sqr,            &bench_setup, NULL, B_FIELD,   10,  200000, 0 },
    { "field_mul",                &bench_field_mul,            &bench_setup, NULL, B_FIELD,   10,  200000, 0 },
    { "field_inverse",            &bench_field_inverse,        &bench_setup, NULL, B_FIELD,   10,   20000, 0 },
    { "field_inverse_var",        &bench_field_inverse_var,    &bench_setup, NULL, B_FIELD,   10,   20000, 0 },
    { "field_sqrt_var",           &bench_field_sqrt_var,       &bench_setup, NULL, B_FIELD,   10,   20000, 0 },
    { "group_double_var",         &bench_group_double_var,     &bench_setup, NULL, B_GROUP,   10,  200000, 0 },
    { "group_add_var",            &bench_group_add_var,        &bench_setup, NULL, B_GROUP,   10,  200000, 0 },
    { "group_add_affine",         &bench_group_add_affine,     &bench_setup, NULL, B_GROUP,   10,  200000, 0 },
    { "group_add_affine_var",     &bench_group_add_affine_var, &bench_setup, NULL, B_GROUP,   10,  200000, 0 },
    { "wnaf_const",               &bench_wnaf_const,           &bench_setup, NULL, B_ECMULT,  10,   20000, 0 },
    { "ecmult_wnaf",              &bench_ecmult_wnaf,          &bench_setup, NULL, B_ECMULT,  10,   20000, 0 },
    { "hash_sha256",              &bench_sha256,               &bench_setup, NULL, B_HASH,    10,   20000, 0 },
    { "hash_hmac_sha256",         &bench_hmac_sha256,          &bench_setup, NULL, B_HASH,    10,   20000, 0 },
    { "hash_rfc6979_hmac_sha256", &bench_rfc6979_hmac_sha256,  &bench_setup, NULL, B_HASH,    10,   20000, 0 },
    { "context_verify",           &bench_context_verify,       &bench_setup, NULL, B_CONTEXT, 10,      20, 0 },
    { "context_sign",             &bench_context_sign,         &bench_setup, NULL, B_CONTEXT, 10,     200, 0 },
    { "ecdh",                     &bench_ecdh,                 &bench_ecdh_setup, NULL, B_ECDH,    10, 20000, 0 },
    { "ecdsa_recover",            &bench_recover,              &bench_recover_setup, &bench_recover_teardown, B_ECDSA,   10, 20000, 0 },
    { "schnorr_verify",           &bench_schnorr_verify,       &bench_schnorr_init,  &bench_schnorr_teardown, B_SCHNORR, 10, 20000, 0 },
    { "ecdsa_sign",               &bench_sign,                 &bench_sign_setup,    &bench_sign_teardown,    B_ECDSA,   10, 20000, 0 },
    { "ecdsa_verify",             &bench_verify,               &bench_verify_setup,  &bench_verify_teardown,  B_ECDSA,   10, 20000, 0 },
    { NULL }
};

void usage(char *av) {
    int i, t;

    t=find_bench_member(bench_table, NULL, 0UL);

    printf("%s [[scalar|field|group|ecmult|hash|context|{name}|-c n|-i n|-w n|-h|all|clear|go] ...]\n", av);
    printf("scalar  - queue all scalar-related benchmarks\n");
    printf("field   - queue all field-related benchmarks\n");
    printf("group   - queue all group-related benchmarks\n");
    printf("ecmult  - queue all ecmult-related benchmarks\n");
    printf("hash    - queue all hash-related benchmarks\n");
    printf("context - queue all context-related benchmarks\n");
    printf("{name}  - queue specific benchmark by name (see list below)\n");
    printf("-c n    - set the number of times the benchmark member should be run\n");
    printf("-i n    - set the number of iterations internal to the benchmark member\n");
    printf("-w n    - configure runtime WINDOW_G table optimization size\n");
    printf("all     - queue all benchmarks available\n");
    printf("clear   - clear the benchmark queue\n");
    printf("go      - perform all runnable benchmarks and print results\n\n");
    printf("Benchmarks available:\n");

    for (i=0; i<t; i++) {
        if (i>0) {
            printf(" ");
        }
        printf("%s", bench_table[i].name);
    }
    printf("\n");
}

void bench_set_type_runnable(unsigned int colour) {
    int i, t;
    t=find_bench_member(bench_table, NULL, 0UL);

    for(i=0; i<t; i++) {
        if (bench_table[i].flags & colour) {
            bench_table[i].run = 1;
        }
    }
}

int main(int argc, char **argv) {
    bench_t data;
    int i=0, runs=0, iters=0, windowg=WINDOW_G;
    int whichbench=0, total_benchmarks=0;
    long unsigned int bflags=0;
    int emode=B_SCANNING;

    data.windowG_override=WINDOW_G;

    if (argc < 2) {
        printf("No arguments.\n");
        return -1;
    }

    for (i=1; i<argc; i++) {
        switch(emode) {
            case B_SCANNING:
                if (!strcmp(argv[i], "scalar")) {
                    bench_set_type_runnable(B_SCALAR);
                } else if (!strcmp(argv[i], "field")) {
                    bench_set_type_runnable(B_FIELD);
                } else if (!strcmp(argv[i], "group")) {
                    bench_set_type_runnable(B_GROUP);
                } else if (!strcmp(argv[i], "ecmult")) {
                    bench_set_type_runnable(B_ECMULT);
                } else if (!strcmp(argv[i], "hash")) {
                    bench_set_type_runnable(B_HASH);
                } else if (!strcmp(argv[i], "context")) {
                    bench_set_type_runnable(B_CONTEXT);
                } else if (!strcmp(argv[i], "-c")) {
                    emode=B_COUNT;
                } else if (!strcmp(argv[i], "-i")) {
                    emode=B_ITERS;
                } else if (!strcmp(argv[i], "-w")) {
                    emode=B_WINDOWG;
                } else if (!strcmp(argv[i], "-h")) {
                    usage(argv[0]);
                    return 0;
                } else if (!strcmp(argv[i], "go")) {
                    total_benchmarks=find_bench_member(bench_table, NULL, 0UL);

                    for (whichbench=0; whichbench<total_benchmarks; whichbench++) {
                        if (bench_table[whichbench].flags & bflags || bench_table[whichbench].run == 1) {
                            run_benchmark(bench_table[whichbench].name, bench_table[whichbench].func, bench_table[whichbench].setup, \
                                bench_table[whichbench].teardown, &data, runs?runs:bench_table[whichbench].default_runs, \
                                iters?iters:bench_table[whichbench].default_iters);
                        }
                    }
                } else if (!strcmp(argv[i], "all")) {
                    bflags=B_ALL;
                } else if (!strcmp(argv[i], "clear")) {
                    bflags=0L;
                    for (whichbench=find_bench_member(bench_table, NULL, 0UL); whichbench>0; whichbench--) {
                        bench_table[whichbench].run=0;
                    }
                } else if ((whichbench=find_bench_member(bench_table, argv[i], 0UL)) > 0) {
                    bench_table[whichbench].run=1;
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
                if (windowg < 2 || windowg > 19) {
                    printf("windowg value out of bounds (2-19) resetting to 16\n");
                    windowg=16;
                }
                data.windowG_override=windowg;
                emode=B_SCANNING;
                break;
            default:
                break;
        }
    }
    return 0;
}
