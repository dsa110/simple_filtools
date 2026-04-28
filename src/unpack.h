/*
 * unpack.h - n-bit sample unpacking for SIGPROC filterbank data.
 *
 * SIGPROC packs sub-byte samples LSB-first within each byte:
 *   nbits=1:  8 samples/byte, sample k -> bit k
 *   nbits=2:  4 samples/byte, sample k -> bits (2k .. 2k+1)
 *   nbits=4:  2 samples/byte, sample 0 = low nibble, sample 1 = high nibble
 *   nbits=8:  1 sample/byte (passthrough)
 * Reference: dsa110-sigproc/src/pack_unpack.c, read_block.c.
 *
 * For 2-bit data (this tool's primary case) we provide a hot bulk loop.
 */
#ifndef SIMPLE_FILTOOLS_UNPACK_H
#define SIMPLE_FILTOOLS_UNPACK_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Unpack `n_in` packed bytes into `out`. The number of output samples is
 * 8*n_in / nbits, which the caller is responsible for sizing.
 * Returns the number of samples written, or -1 on unsupported nbits. */
long unpack_to_uint8(const uint8_t *in, size_t n_in, int nbits, uint8_t *out);

/* Specialized 2-bit unpacker. Writes 4*n_in uint8 samples. */
void unpack_2bit(const uint8_t *in, size_t n_in, uint8_t *out);

/* Per-byte reference unpack for nbits=2 (matches sigproc char2fourints). */
static inline void unpack_2bit_byte(uint8_t b,
                                    uint8_t *s0, uint8_t *s1,
                                    uint8_t *s2, uint8_t *s3) {
    *s0 = (uint8_t)( b       & 0x3u);
    *s1 = (uint8_t)((b >> 2) & 0x3u);
    *s2 = (uint8_t)((b >> 4) & 0x3u);
    *s3 = (uint8_t)((b >> 6) & 0x3u);
}

#ifdef __cplusplus
}
#endif

#endif /* SIMPLE_FILTOOLS_UNPACK_H */
