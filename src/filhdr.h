/*
 * filhdr.h - minimal SIGPROC filterbank header reader.
 *
 * Ported from dsa110-sigproc (read_header.c, header.h, strings_equal.c) and
 * trimmed to a self-contained, globals-free API. Supports the standard
 * SIGPROC HEADER_START / HEADER_END text-tag binary header used by
 * filterbank, dedisperse, prepfold, etc.
 */
#ifndef SIMPLE_FILTOOLS_FILHDR_H
#define SIMPLE_FILTOOLS_FILHDR_H

#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define FIL_NAME_LEN 256

typedef struct {
    int      machine_id;
    int      telescope_id;
    int      data_type;
    int      nchans;
    int      nbits;
    int      nifs;
    int      nbeams;
    int      ibeam;
    int      barycentric;
    int      pulsarcentric;
    double   tstart;
    double   tsamp;
    double   fch1;
    double   foff;
    double   refdm;
    double   src_raj;
    double   src_dej;
    double   az_start;
    double   za_start;
    char     source_name[FIL_NAME_LEN];
    char     rawdatafile[FIL_NAME_LEN];

    /* file/io bookkeeping */
    long     header_bytes;     /* bytes occupied by the header */
    long long file_bytes;      /* total file size in bytes */
    long long data_bytes;      /* file_bytes - header_bytes */
    long long nsamples;        /* number of complete time samples in file */
} fil_header;

/*
 * Open a SIGPROC filterbank file and parse its header.
 * On success, returns a valid FILE* positioned at the start of the data
 * payload, populates *hdr, and returns 0 in *err.
 * On failure, returns NULL and sets *err to a non-zero value; if errmsg
 * is non-NULL, a human-readable description is written there (errmsg_len
 * bytes max).
 */
FILE *fil_open(const char *path, fil_header *hdr,
               int *err, char *errmsg, size_t errmsg_len);

/* Pretty-print key header fields to fp (for diagnostics). */
void fil_header_print(FILE *fp, const fil_header *hdr);

/* Number of bytes per complete time sample (over all channels and IFs). */
static inline long fil_bytes_per_sample(const fil_header *hdr) {
    /* nbits * nchans * nifs / 8, rounded up if not multiple of 8 */
    long bits = (long)hdr->nbits * hdr->nchans * hdr->nifs;
    return (bits + 7) / 8;
}

#ifdef __cplusplus
}
#endif

#endif /* SIMPLE_FILTOOLS_FILHDR_H */
