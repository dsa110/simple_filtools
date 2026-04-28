/*
 * chop_fil - cut a time slice out of a SIGPROC filterbank file.
 *
 *   chop_fil <input.fil> -s <start_sec> -d <dur_sec> [-o <output.fil>]
 *
 * Copies the original header verbatim (preserving any exotic tags), then
 * patches the `tstart` value in the new header to reflect the new start
 * time:
 *
 *     tstart_new = tstart_old + start_samp * tsamp / 86400
 *
 * where start_samp = round(start_sec / tsamp). After the header, copies
 * `dur_samp = round(dur_sec / tsamp)` time samples of raw data.
 *
 * For nbits<8 files, sample-aligned cuts are byte-aligned iff
 * nchans*nbits is a multiple of 8 -- true for all standard SIGPROC
 * layouts (e.g. DSA-2000's 2-bit / 2048-channel files).
 */

/* POSIX feature flags: needed on Linux/glibc to expose fseeko / off_t
 * when compiling with -std=c99. Must come before any system header. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif

#include "filhdr.h"

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

static void die(const char *msg) {
    fprintf(stderr, "chop_fil: %s\n", msg);
    exit(1);
}

static void usage(FILE *fp, const char *argv0) {
    fprintf(fp,
        "Usage: %s <input.fil> -s <start_sec> -d <dur_sec> [-o <output.fil>]\n"
        "Options:\n"
        "  -s <start_sec>       Start of the slice, in seconds from file start\n"
        "  -d <dur_sec>         Duration of the slice, in seconds\n"
        "  -o <output.fil>      Output file (default: <input>.cut.fil)\n"
        "  -h, --help           Show this help\n",
        argv0);
}

/* Locate the 8-byte `tstart` value inside a copy of the raw header.
 * SIGPROC stores each tag as int32 length + ASCII bytes; the value
 * follows the tag with no separator. Returns the byte offset of the
 * value within `buf`, or -1 if not found. */
static long find_tstart_value_offset(const uint8_t *buf, long n) {
    const char *needle = "tstart";
    long nlen = (long)strlen(needle);
    for (long i = 0; i + nlen + 8 <= n; i++) {
        /* The tag is preceded by its int32 length (= 6 for "tstart"). */
        if (i >= 4 && memcmp(buf + i, needle, (size_t)nlen) == 0) {
            int32_t plen;
            memcpy(&plen, buf + i - 4, sizeof(plen));
            if (plen == nlen) return i + nlen;
        }
    }
    return -1;
}

int main(int argc, char **argv) {
    const char *in_path = NULL;
    char out_path[1024] = {0};
    double start_sec = -1.0;
    double dur_sec   = -1.0;

    if (argc < 2) { usage(stderr, argv[0]); return 1; }
    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        usage(stdout, argv[0]); return 0;
    }
    in_path = argv[1];

    for (int i = 2; i < argc; i++) {
        const char *a = argv[i];
        if      (!strcmp(a, "-s") && i + 1 < argc) { start_sec = atof(argv[++i]); }
        else if (!strcmp(a, "-d") && i + 1 < argc) { dur_sec   = atof(argv[++i]); }
        else if (!strcmp(a, "-o") && i + 1 < argc) { strncpy(out_path, argv[++i], sizeof(out_path) - 1); }
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { usage(stdout, argv[0]); return 0; }
        else {
            fprintf(stderr, "chop_fil: unknown option '%s'\n", a);
            usage(stderr, argv[0]);
            return 1;
        }
    }
    if (start_sec < 0.0) die("missing or negative -s start_sec");
    if (dur_sec  <= 0.0) die("missing or non-positive -d dur_sec");
    if (out_path[0] == '\0') {
        snprintf(out_path, sizeof(out_path), "%s.cut.fil", in_path);
    }

    fil_header h;
    int err = 0;
    char errmsg[256];
    FILE *fin = fil_open(in_path, &h, &err, errmsg, sizeof(errmsg));
    if (!fin) { fprintf(stderr, "chop_fil: %s\n", errmsg); return 2; }
    if (h.nbits != 1 && h.nbits != 2 && h.nbits != 4 && h.nbits != 8 && h.nbits != 16 && h.nbits != 32) {
        fprintf(stderr, "chop_fil: unsupported nbits=%d\n", h.nbits);
        fclose(fin); return 2;
    }
    if ((long long)h.nbits * h.nchans % 8 != 0) {
        fprintf(stderr,
                "chop_fil: nbits*nchans=%d not a multiple of 8; cannot cleanly cut on sample boundaries\n",
                h.nbits * h.nchans);
        fclose(fin); return 2;
    }

    long bps = fil_bytes_per_sample(&h);
    long long start_samp = (long long)floor(start_sec / h.tsamp + 0.5);
    long long dur_samp   = (long long)floor(dur_sec   / h.tsamp + 0.5);
    if (start_samp >= h.nsamples) {
        fprintf(stderr,
                "chop_fil: start_sec=%g (sample %lld) is past end of file (%lld samples, %.3f s)\n",
                start_sec, start_samp, h.nsamples, h.nsamples * h.tsamp);
        fclose(fin); return 2;
    }
    long long avail = h.nsamples - start_samp;
    if (dur_samp > avail) {
        fprintf(stderr,
                "chop_fil: warning, requested %lld samples (%g s) but only %lld available; truncating\n",
                dur_samp, dur_sec, avail);
        dur_samp = avail;
    }

    /* Read the original header bytes verbatim. */
    if (fseek(fin, 0L, SEEK_SET) != 0) die("seek to header start");
    uint8_t *hdr_buf = (uint8_t *)malloc((size_t)h.header_bytes);
    if (!hdr_buf) die("out of memory (header buffer)");
    if (fread(hdr_buf, 1, (size_t)h.header_bytes, fin) != (size_t)h.header_bytes)
        die("short read of input header");

    /* Patch tstart in place. */
    double tstart_new = h.tstart + (double)start_samp * h.tsamp / 86400.0;
    long off = find_tstart_value_offset(hdr_buf, h.header_bytes);
    if (off < 0) {
        fprintf(stderr,
                "chop_fil: warning, could not find 'tstart' tag in header; "
                "writing header verbatim with original tstart=%.10f\n",
                h.tstart);
    } else {
        memcpy(hdr_buf + off, &tstart_new, sizeof(double));
    }

    /* Open output and write header. */
    FILE *fout = fopen(out_path, "wb");
    if (!fout) {
        fprintf(stderr, "chop_fil: cannot open '%s' for writing: %s\n",
                out_path, strerror(errno));
        free(hdr_buf); fclose(fin); return 3;
    }
    if (fwrite(hdr_buf, 1, (size_t)h.header_bytes, fout) != (size_t)h.header_bytes)
        die("short write of output header");
    free(hdr_buf);

    /* Seek input to the first byte of the requested data slice. */
    long long byte_offset = (long long)h.header_bytes + start_samp * bps;
    if (fseeko(fin, (off_t)byte_offset, SEEK_SET) != 0) {
        fprintf(stderr, "chop_fil: seek to byte %lld failed: %s\n",
                byte_offset, strerror(errno));
        fclose(fin); fclose(fout); return 3;
    }

    /* Stream `dur_samp * bps` bytes from input to output in 1 MiB chunks. */
    long long bytes_to_copy = dur_samp * bps;
    enum { BUFSZ = 1 << 20 };
    uint8_t *buf = (uint8_t *)malloc(BUFSZ);
    if (!buf) die("out of memory (copy buffer)");
    long long copied = 0;
    while (copied < bytes_to_copy) {
        long long want = bytes_to_copy - copied;
        if (want > BUFSZ) want = BUFSZ;
        size_t got = fread(buf, 1, (size_t)want, fin);
        if (got == 0) {
            fprintf(stderr,
                    "chop_fil: short read at offset %lld (copied %lld / %lld bytes)\n",
                    byte_offset + copied, copied, bytes_to_copy);
            break;
        }
        if (fwrite(buf, 1, got, fout) != got) die("write to output failed");
        copied += (long long)got;
    }
    free(buf);
    fclose(fin);
    fclose(fout);

    long long actual_samp = copied / bps;
    fprintf(stderr,
            "chop_fil: %s -> %s\n"
            "  tstart : %.10f -> %.10f  (+%.6f s)\n"
            "  samples: [%lld .. %lld)  (%lld samp, %.6f s)\n"
            "  bytes  : header=%ld + data=%lld  (= %lld total)\n",
            in_path, out_path,
            h.tstart, tstart_new, (double)start_samp * h.tsamp,
            start_samp, start_samp + actual_samp, actual_samp, actual_samp * h.tsamp,
            h.header_bytes, copied, h.header_bytes + copied);
    return 0;
}
