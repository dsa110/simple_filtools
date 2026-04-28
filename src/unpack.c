/*
 * unpack.c - n-bit -> uint8 unpacking for SIGPROC filterbank data.
 */
#include "unpack.h"

#include <string.h>

void unpack_2bit(const uint8_t *in, size_t n_in, uint8_t *out) {
    /* Hot path. Compilers happily unroll this; keep it explicit so the
     * intent matches sigproc's char2fourints exactly. */
    for (size_t i = 0; i < n_in; i++) {
        uint8_t b = in[i];
        out[4*i + 0] = (uint8_t)( b       & 0x3u);
        out[4*i + 1] = (uint8_t)((b >> 2) & 0x3u);
        out[4*i + 2] = (uint8_t)((b >> 4) & 0x3u);
        out[4*i + 3] = (uint8_t)((b >> 6) & 0x3u);
    }
}

static void unpack_1bit(const uint8_t *in, size_t n_in, uint8_t *out) {
    for (size_t i = 0; i < n_in; i++) {
        uint8_t b = in[i];
        for (int k = 0; k < 8; k++) {
            out[8*i + k] = (uint8_t)((b >> k) & 0x1u);
        }
    }
}

static void unpack_4bit(const uint8_t *in, size_t n_in, uint8_t *out) {
    for (size_t i = 0; i < n_in; i++) {
        uint8_t b = in[i];
        out[2*i + 0] = (uint8_t)( b       & 0xFu);
        out[2*i + 1] = (uint8_t)((b >> 4) & 0xFu);
    }
}

long unpack_to_uint8(const uint8_t *in, size_t n_in, int nbits, uint8_t *out) {
    switch (nbits) {
    case 1:
        unpack_1bit(in, n_in, out);
        return (long)(8 * n_in);
    case 2:
        unpack_2bit(in, n_in, out);
        return (long)(4 * n_in);
    case 4:
        unpack_4bit(in, n_in, out);
        return (long)(2 * n_in);
    case 8:
        memcpy(out, in, n_in);
        return (long)n_in;
    default:
        return -1;
    }
}
