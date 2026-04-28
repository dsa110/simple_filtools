/*
 * test_unpack.c - sanity tests for the n-bit unpacker.
 *
 * Verifies:
 *   1. unpack_2bit() matches the per-byte char2fourints reference for every
 *      possible byte value (exhaustive: all 256 bytes).
 *   2. unpack_2bit() matches the per-byte reference on a random byte stream
 *      of length 1<<20.
 *   3. unpack_to_uint8() dispatches correctly for nbits 1, 2, 4, 8.
 *   4. unpack_1bit / unpack_4bit produce the expected sigproc bit-order.
 */
#include "unpack.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

static int fails = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { fprintf(stderr, "FAIL: %s\n", msg); fails++; } \
} while (0)

/* Reference: per-byte 2-bit unpack adapted from dsa110-sigproc/src/pack_unpack.c.
 * (Same masks, same semantics, but writes to a buffer instead of int*.) */
static void ref_unpack_2bit_byte(uint8_t b, uint8_t *out4) {
    out4[0] = b & 0x3u;
    out4[1] = (b >> 2) & 0x3u;
    out4[2] = (b >> 4) & 0x3u;
    out4[3] = (b >> 6) & 0x3u;
}

static void test_2bit_exhaustive(void) {
    uint8_t expected[4], got[4];
    for (int b = 0; b < 256; b++) {
        uint8_t ub = (uint8_t)b;
        ref_unpack_2bit_byte(ub, expected);
        unpack_2bit(&ub, 1, got);
        CHECK(memcmp(expected, got, 4) == 0, "2-bit byte mismatch (exhaustive)");
    }
}

static void test_2bit_random(void) {
    const size_t N = 1 << 20;   /* 1 MiB */
    uint8_t *in  = (uint8_t *)malloc(N);
    uint8_t *out = (uint8_t *)malloc(4 * N);
    uint8_t *ref = (uint8_t *)malloc(4 * N);
    if (!in || !out || !ref) { fprintf(stderr, "OOM\n"); exit(2); }

    srand(42);
    for (size_t i = 0; i < N; i++) in[i] = (uint8_t)rand();

    unpack_2bit(in, N, out);
    for (size_t i = 0; i < N; i++) ref_unpack_2bit_byte(in[i], ref + 4*i);

    CHECK(memcmp(out, ref, 4 * N) == 0, "2-bit bulk unpack matches reference");
    free(in); free(out); free(ref);
}

static void test_dispatch_2bit(void) {
    uint8_t in[2] = {0xE4, 0x1B};   /* 0xE4 = 11100100b -> {0,1,2,3} ; 0x1B -> {3,2,1,0} */
    uint8_t out[8];
    long n = unpack_to_uint8(in, 2, 2, out);
    CHECK(n == 8, "unpack_to_uint8(nbits=2) returns 4*nin");
    uint8_t expect[8] = {0,1,2,3, 3,2,1,0};
    CHECK(memcmp(out, expect, 8) == 0, "dispatch nbits=2 ordering");
}

static void test_dispatch_1bit(void) {
    /* sigproc 1-bit: bit k is sample k (LSB first) */
    uint8_t in[1] = {0xA5};   /* 10100101b */
    uint8_t out[8];
    long n = unpack_to_uint8(in, 1, 1, out);
    CHECK(n == 8, "unpack_to_uint8(nbits=1) returns 8*nin");
    uint8_t expect[8] = {1,0,1,0, 0,1,0,1};
    CHECK(memcmp(out, expect, 8) == 0, "dispatch nbits=1 ordering");
}

static void test_dispatch_4bit(void) {
    uint8_t in[2] = {0x12, 0xFE};   /* 0x12 -> {2, 1}, 0xFE -> {0xE, 0xF} */
    uint8_t out[4];
    long n = unpack_to_uint8(in, 2, 4, out);
    CHECK(n == 4, "unpack_to_uint8(nbits=4) returns 2*nin");
    uint8_t expect[4] = {0x2, 0x1, 0xE, 0xF};
    CHECK(memcmp(out, expect, 4) == 0, "dispatch nbits=4 ordering");
}

static void test_dispatch_8bit(void) {
    uint8_t in[4]  = {1,2,3,4};
    uint8_t out[4] = {0};
    long n = unpack_to_uint8(in, 4, 8, out);
    CHECK(n == 4, "unpack_to_uint8(nbits=8) returns nin");
    CHECK(memcmp(in, out, 4) == 0, "dispatch nbits=8 passthrough");
}

static void test_dispatch_unsupported(void) {
    uint8_t in[1] = {0}, out[8] = {0};
    CHECK(unpack_to_uint8(in, 1, 3,  out) < 0, "unsupported nbits=3 returns -1");
    CHECK(unpack_to_uint8(in, 1, 16, out) < 0, "unsupported nbits=16 returns -1");
}

int main(void) {
    test_2bit_exhaustive();
    test_2bit_random();
    test_dispatch_1bit();
    test_dispatch_2bit();
    test_dispatch_4bit();
    test_dispatch_8bit();
    test_dispatch_unsupported();

    if (fails) {
        fprintf(stderr, "%d test(s) failed.\n", fails);
        return 1;
    }
    printf("test_unpack: OK (all 2-bit packings exhaustive + 1Mi random + n-bit dispatch).\n");
    return 0;
}
