/*
 * diag.h - on-disk format for rfidiag output (.diag) files.
 *
 * Layout:
 *   [DiagHeader]                                  // fixed-size struct
 *   { tag[4], len_u64, payload[len] } * N         // tagged sections
 *
 * Tags currently emitted (always little-endian on the host platforms we run on):
 *   "SUMC"  uint64[nchans]               per-channel sample sum
 *   "SUMQ"  uint64[nchans]               per-channel sample sumsq
 *   "HIST"  uint64[nchans * 2^nbits]     per-channel value histograms
 *   "CMEA"  float32[nchunks * nchans]    per-chunk per-channel mean
 *   "CRMS"  float32[nchunks * nchans]    per-chunk per-channel RMS
 *   "Z0DM"  float32[nsamples]            full-rate zero-DM time series
 *
 * The reader parses tag-by-tag and skips unknown tags, so the format is
 * forward-compatible.
 */
#ifndef SIMPLE_FILTOOLS_DIAG_H
#define SIMPLE_FILTOOLS_DIAG_H

#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define DIAG_MAGIC   "FILDIAG\0"
#define DIAG_VERSION 1u

#pragma pack(push, 1)
typedef struct {
    char     magic[8];        /* "FILDIAG\0" */
    uint32_t version;         /* DIAG_VERSION */
    uint32_t nchans;
    uint32_t nbits;
    uint32_t nchunks;
    uint64_t nsamples;        /* number of full-rate time samples processed */
    uint64_t nsamples_per_chunk;
    double   tsamp;           /* seconds */
    double   fch1;            /* MHz */
    double   foff;            /* MHz (signed) */
    double   tstart;          /* MJD */
    double   chunk_sec;       /* nominal chunk duration in seconds */
    double   t_offset;        /* seconds: data start offset relative to file tstart */
    uint32_t nhist;           /* histogram bins per channel = 2^nbits */
    uint32_t reserved;
} DiagHeader;
#pragma pack(pop)

/* Open output file and write the fixed DiagHeader. Returns FILE* or NULL. */
FILE *diag_open_write(const char *path, const DiagHeader *hdr);

/* Write a tagged section. tag must be exactly 4 chars (null padding allowed). */
int diag_write_section(FILE *fp, const char tag[4],
                       const void *data, uint64_t nbytes);

#ifdef __cplusplus
}
#endif

#endif /* SIMPLE_FILTOOLS_DIAG_H */
