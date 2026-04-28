/*
 * filhdr.c - minimal SIGPROC filterbank header reader.
 *
 * Ported from dsa110-sigproc/src/{read_header.c, strings_equal.c, header.h}.
 * The original sigproc reader uses globals; this version returns everything
 * in a fil_header struct and never aborts the process on bad input.
 */

/* POSIX feature flags: keeps stat() and friends visible under -std=c99
 * on Linux/glibc. Must come before any system header. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif

#include "filhdr.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static int streq(const char *a, const char *b) {
    return strcmp(a, b) == 0;
}

/* read a length-prefixed sigproc string into out (must be >= 81 bytes).
 * returns total bytes consumed (4 + nchar) on success, -1 on EOF/error. */
static int get_string(FILE *fp, char *out, size_t outlen) {
    int nchar = 0;
    if (fread(&nchar, sizeof(int), 1, fp) != 1) return -1;
    if (nchar < 1 || nchar > 80) return -1;
    if ((size_t)nchar + 1 > outlen) return -1;
    if (fread(out, 1, (size_t)nchar, fp) != (size_t)nchar) return -1;
    out[nchar] = '\0';
    return (int)sizeof(int) + nchar;
}

/* read an int into *out from fp; advance bytes counter. */
static int read_int(FILE *fp, int *out, long *bytes) {
    if (fread(out, sizeof(int), 1, fp) != 1) return -1;
    *bytes += (long)sizeof(int);
    return 0;
}
static int read_double(FILE *fp, double *out, long *bytes) {
    if (fread(out, sizeof(double), 1, fp) != 1) return -1;
    *bytes += (long)sizeof(double);
    return 0;
}

static long long file_size_bytes(const char *path) {
    struct stat st;
    if (stat(path, &st) != 0) return -1;
    return (long long)st.st_size;
}

static void seterr(char *errmsg, size_t n, const char *fmt, ...) {
    if (!errmsg || n == 0) return;
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(errmsg, n, fmt, ap);
    va_end(ap);
}

FILE *fil_open(const char *path, fil_header *hdr,
               int *err, char *errmsg, size_t errmsg_len)
{
    if (err) *err = 0;
    if (!hdr || !path) {
        if (err) *err = 1;
        seterr(errmsg, errmsg_len, "fil_open: NULL argument");
        return NULL;
    }

    memset(hdr, 0, sizeof(*hdr));
    hdr->nifs = 1;          /* sigproc default */
    hdr->nbeams = 1;
    hdr->ibeam = 1;
    hdr->machine_id = 0;
    hdr->telescope_id = 0;

    long long fbytes = file_size_bytes(path);
    if (fbytes < 0) {
        if (err) *err = 1;
        seterr(errmsg, errmsg_len, "cannot stat '%s'", path);
        return NULL;
    }
    hdr->file_bytes = fbytes;

    FILE *fp = fopen(path, "rb");
    if (!fp) {
        if (err) *err = 1;
        seterr(errmsg, errmsg_len, "cannot open '%s'", path);
        return NULL;
    }

    long total = 0;
    char tag[81];
    int n = get_string(fp, tag, sizeof(tag));
    if (n < 0 || !streq(tag, "HEADER_START")) {
        if (err) *err = 2;
        seterr(errmsg, errmsg_len,
               "no HEADER_START in '%s' (not a sigproc filterbank file?)", path);
        fclose(fp);
        return NULL;
    }
    total += n;

    int expecting_rawdatafile = 0, expecting_source_name = 0;
    /* dummies for sigproc tags we read but ignore */
    int    tmp_i;
    double tmp_d;

    for (;;) {
        n = get_string(fp, tag, sizeof(tag));
        if (n < 0) {
            if (err) *err = 3;
            seterr(errmsg, errmsg_len, "truncated header in '%s'", path);
            fclose(fp);
            return NULL;
        }
        total += n;
        if (streq(tag, "HEADER_END")) break;

        if (streq(tag, "rawdatafile")) {
            expecting_rawdatafile = 1;
        } else if (streq(tag, "source_name")) {
            expecting_source_name = 1;
        } else if (streq(tag, "machine_id"))     { if (read_int   (fp, &hdr->machine_id,    &total) < 0) goto trunc; }
        else if (streq(tag, "telescope_id"))     { if (read_int   (fp, &hdr->telescope_id,  &total) < 0) goto trunc; }
        else if (streq(tag, "data_type"))        { if (read_int   (fp, &hdr->data_type,     &total) < 0) goto trunc; }
        else if (streq(tag, "nchans"))           { if (read_int   (fp, &hdr->nchans,        &total) < 0) goto trunc; }
        else if (streq(tag, "nbits"))            { if (read_int   (fp, &hdr->nbits,         &total) < 0) goto trunc; }
        else if (streq(tag, "nifs"))             { if (read_int   (fp, &hdr->nifs,          &total) < 0) goto trunc; }
        else if (streq(tag, "nbeams"))           { if (read_int   (fp, &hdr->nbeams,        &total) < 0) goto trunc; }
        else if (streq(tag, "ibeam"))            { if (read_int   (fp, &hdr->ibeam,         &total) < 0) goto trunc; }
        else if (streq(tag, "barycentric"))      { if (read_int   (fp, &hdr->barycentric,   &total) < 0) goto trunc; }
        else if (streq(tag, "pulsarcentric"))    { if (read_int   (fp, &hdr->pulsarcentric, &total) < 0) goto trunc; }
        else if (streq(tag, "tstart"))           { if (read_double(fp, &hdr->tstart,        &total) < 0) goto trunc; }
        else if (streq(tag, "tsamp"))            { if (read_double(fp, &hdr->tsamp,         &total) < 0) goto trunc; }
        else if (streq(tag, "fch1"))             { if (read_double(fp, &hdr->fch1,          &total) < 0) goto trunc; }
        else if (streq(tag, "foff"))             { if (read_double(fp, &hdr->foff,          &total) < 0) goto trunc; }
        else if (streq(tag, "refdm"))            { if (read_double(fp, &hdr->refdm,         &total) < 0) goto trunc; }
        else if (streq(tag, "src_raj"))          { if (read_double(fp, &hdr->src_raj,       &total) < 0) goto trunc; }
        else if (streq(tag, "src_dej"))          { if (read_double(fp, &hdr->src_dej,       &total) < 0) goto trunc; }
        else if (streq(tag, "az_start"))         { if (read_double(fp, &hdr->az_start,      &total) < 0) goto trunc; }
        else if (streq(tag, "za_start"))         { if (read_double(fp, &hdr->za_start,      &total) < 0) goto trunc; }
        else if (streq(tag, "period"))           { if (read_double(fp, &tmp_d,              &total) < 0) goto trunc; }
        else if (streq(tag, "nbins"))            { if (read_int   (fp, &tmp_i,              &total) < 0) goto trunc; }
        else if (streq(tag, "npuls"))            { /* int64 */
            int64_t v;
            if (fread(&v, sizeof(v), 1, fp) != 1) goto trunc;
            total += (long)sizeof(v);
        }
        else if (streq(tag, "nsamples"))         { if (read_int   (fp, &tmp_i,              &total) < 0) goto trunc; }
        else if (streq(tag, "signed"))           { /* one byte */
            char c;
            if (fread(&c, 1, 1, fp) != 1) goto trunc;
            total += 1;
        }
        else if (streq(tag, "FREQUENCY_START") || streq(tag, "FREQUENCY_END")) {
            /* not supported; sigproc uses these for non-contiguous channel tables.
             * we just skip over them by re-reading FREQUENCY_END below. */
        }
        else if (streq(tag, "fchannel")) {
            if (read_double(fp, &tmp_d, &total) < 0) goto trunc;
        }
        else if (expecting_rawdatafile) {
            strncpy(hdr->rawdatafile, tag, FIL_NAME_LEN - 1);
            expecting_rawdatafile = 0;
        }
        else if (expecting_source_name) {
            strncpy(hdr->source_name, tag, FIL_NAME_LEN - 1);
            expecting_source_name = 0;
        }
        else {
            if (err) *err = 4;
            seterr(errmsg, errmsg_len,
                   "unknown sigproc header tag '%s' in '%s'", tag, path);
            fclose(fp);
            return NULL;
        }
    }
    (void)tmp_i; (void)tmp_d;

    /* sanity check that we are positioned where the running byte counter says */
    long pos = ftell(fp);
    if (pos != total) {
        if (err) *err = 5;
        seterr(errmsg, errmsg_len,
               "header byte mismatch in '%s' (counter=%ld, ftell=%ld)",
               path, total, pos);
        fclose(fp);
        return NULL;
    }

    if (hdr->nchans <= 0 || hdr->nbits <= 0 || hdr->tsamp <= 0.0) {
        if (err) *err = 6;
        seterr(errmsg, errmsg_len,
               "invalid header in '%s' (nchans=%d nbits=%d tsamp=%g)",
               path, hdr->nchans, hdr->nbits, hdr->tsamp);
        fclose(fp);
        return NULL;
    }
    if (hdr->nifs <= 0) hdr->nifs = 1;

    hdr->header_bytes = total;
    hdr->data_bytes = hdr->file_bytes - total;
    long bps = fil_bytes_per_sample(hdr);
    hdr->nsamples = (bps > 0) ? hdr->data_bytes / bps : 0;

    return fp;

trunc:
    if (err) *err = 3;
    seterr(errmsg, errmsg_len, "truncated header value in '%s'", path);
    fclose(fp);
    return NULL;
}

void fil_header_print(FILE *fp, const fil_header *hdr) {
    fprintf(fp, "  source_name : %s\n", hdr->source_name[0] ? hdr->source_name : "(unset)");
    fprintf(fp, "  telescope_id: %d\n", hdr->telescope_id);
    fprintf(fp, "  machine_id  : %d\n", hdr->machine_id);
    fprintf(fp, "  nchans      : %d\n", hdr->nchans);
    fprintf(fp, "  nbits       : %d\n", hdr->nbits);
    fprintf(fp, "  nifs        : %d\n", hdr->nifs);
    fprintf(fp, "  fch1 (MHz)  : %.6f\n", hdr->fch1);
    fprintf(fp, "  foff (MHz)  : %.6f\n", hdr->foff);
    fprintf(fp, "  tstart (MJD): %.10f\n", hdr->tstart);
    fprintf(fp, "  tsamp (s)   : %.9g\n", hdr->tsamp);
    fprintf(fp, "  header_bytes: %ld\n", hdr->header_bytes);
    fprintf(fp, "  data_bytes  : %lld\n", hdr->data_bytes);
    fprintf(fp, "  nsamples    : %lld  (~%.3f s)\n",
            hdr->nsamples, hdr->nsamples * hdr->tsamp);
}
