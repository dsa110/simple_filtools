/*
 * rfidiag - streaming RFI diagnostics for SIGPROC filterbank data.
 *
 * Treats samples as raw integer intensities (e.g. {0,1,2,3} for nbits=2);
 * data are positive-definite / chi-squared (NOT Gaussian). Per-channel and
 * zero-DM statistics use the mean and direct RMS, not median/MAD.
 *
 * Output: a single binary .diag file (see diag.h). A companion python
 * script (python/plot_diagnostics.py) turns the .diag into PNG plots.
 *
 * Memory:
 *   - per-channel accumulators ........ O(nchans * nhist)         (small)
 *   - per-chunk arrays  (cmean, crms) . O(nchunks * nchans * 4)   (~38 MB for 0.5s/40min)
 *   - full-rate zero-DM ............... O(nsamples * 4)           (~192 MB for 40 min @ 50us)
 *   - one chunk's unpacked samples .... O(nchans * chunk_nsamp)   (~20 MB)
 * Total at default settings on 40-min/2-bit/2048ch input: ~250 MB.
 */

#include "filhdr.h"
#include "unpack.h"
#include "diag.h"

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>

static void die(const char *msg) {
    fprintf(stderr, "rfidiag: %s\n", msg);
    exit(1);
}

static void usage(FILE *fp, const char *argv0) {
    fprintf(fp,
        "Usage: %s <input.fil> [options]\n"
        "Options:\n"
        "  -o <out.diag>        Output diagnostics file (default: <input>.diag)\n"
        "  -s <start_sec>       Skip the first <start_sec> seconds (default: 0)\n"
        "  -d <dur_sec>         Process only <dur_sec> seconds (default: all)\n"
        "  -c <chunk_sec>       Chunk size for chunk-averaged stats and waterfall\n"
        "                       (default: 0.5)\n"
        "  -Z                   Do not store full-rate zero-DM in the output\n"
        "                       (saves ~4 bytes/sample, e.g. ~192 MB for 40 min @ 50us)\n"
        "  -q                   Quiet (no progress)\n"
        "  -h                   Show this help\n",
        argv0);
}

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    const char *in_path = NULL;
    char out_path[1024] = {0};
    double start_sec = 0.0;
    double dur_sec = -1.0;       /* <0 => to end of file */
    double chunk_sec = 0.5;
    int store_zerodm = 1;
    int quiet = 0;

    int i = 1;
    if (argc < 2) { usage(stderr, argv[0]); return 1; }
    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        usage(stdout, argv[0]); return 0;
    }
    in_path = argv[i++];
    while (i < argc) {
        const char *a = argv[i];
        if      (!strcmp(a, "-o") && i + 1 < argc) { strncpy(out_path, argv[++i], sizeof(out_path)-1); }
        else if (!strcmp(a, "-s") && i + 1 < argc) { start_sec = atof(argv[++i]); }
        else if (!strcmp(a, "-d") && i + 1 < argc) { dur_sec   = atof(argv[++i]); }
        else if (!strcmp(a, "-c") && i + 1 < argc) { chunk_sec = atof(argv[++i]); }
        else if (!strcmp(a, "-Z"))                 { store_zerodm = 0; }
        else if (!strcmp(a, "-q"))                 { quiet = 1; }
        else if (!strcmp(a, "-h"))                 { usage(stdout, argv[0]); return 0; }
        else {
            fprintf(stderr, "rfidiag: unknown option '%s'\n", a);
            usage(stderr, argv[0]);
            return 1;
        }
        i++;
    }
    if (out_path[0] == '\0') {
        snprintf(out_path, sizeof(out_path), "%s.diag", in_path);
    }
    if (chunk_sec <= 0) die("-c chunk_sec must be > 0");
    if (start_sec < 0)  die("-s start_sec must be >= 0");

    /* --- Open + parse header --- */
    fil_header h;
    int err = 0;
    char errmsg[256];
    FILE *fin = fil_open(in_path, &h, &err, errmsg, sizeof(errmsg));
    if (!fin) {
        fprintf(stderr, "rfidiag: %s\n", errmsg);
        return 2;
    }
    if (h.nifs != 1) {
        fprintf(stderr, "rfidiag: only nifs=1 is supported (got nifs=%d)\n", h.nifs);
        fclose(fin);
        return 2;
    }
    if (h.nbits != 1 && h.nbits != 2 && h.nbits != 4 && h.nbits != 8) {
        fprintf(stderr, "rfidiag: unsupported nbits=%d (supported: 1,2,4,8)\n", h.nbits);
        fclose(fin);
        return 2;
    }

    if (!quiet) {
        fprintf(stderr, "Input:   %s\n", in_path);
        fil_header_print(stderr, &h);
        fprintf(stderr, "Output:  %s\n", out_path);
    }

    /* --- Determine sample range --- */
    long long start_samp = (long long)floor(start_sec / h.tsamp + 0.5);
    if (start_samp > h.nsamples) start_samp = h.nsamples;
    long long end_samp = h.nsamples;
    if (dur_sec >= 0) {
        long long n = (long long)floor(dur_sec / h.tsamp + 0.5);
        if (start_samp + n < end_samp) end_samp = start_samp + n;
    }
    long long nsamples = end_samp - start_samp;
    if (nsamples <= 0) {
        fprintf(stderr, "rfidiag: no samples to process (start_samp=%lld nsamples_total=%lld)\n",
                start_samp, h.nsamples);
        fclose(fin);
        return 2;
    }

    /* fseek to start of selected data */
    long bps = fil_bytes_per_sample(&h);
    long long byte_offset = (long long)h.header_bytes + start_samp * bps;
    if (fseeko(fin, (off_t)byte_offset, SEEK_SET) != 0) {
        fprintf(stderr, "rfidiag: seek to byte %lld failed: %s\n", byte_offset, strerror(errno));
        fclose(fin);
        return 2;
    }

    /* --- Chunking --- */
    long long chunk_nsamp = (long long)floor(chunk_sec / h.tsamp + 0.5);
    if (chunk_nsamp < 1) chunk_nsamp = 1;
    /* For nbits<8 we need chunk_nsamp * nchans * nbits to be a whole number of
     * bytes. Round chunk_nsamp up to a multiple of (8/nbits) so that one chunk
     * always lands on byte boundaries. */
    if (h.nbits < 8) {
        long long spb = 8 / h.nbits;     /* samples per byte (at fixed channel) */
        /* actually we want (chunk_nsamp * nchans) to be divisible by spb */
        long long need = spb;
        long long g = 1;
        long long a = (long long)h.nchans, b = need;
        while (b) { long long t = a % b; a = b; b = t; } g = a;     /* gcd */
        long long m = need / g;          /* required factor for chunk_nsamp */
        if (chunk_nsamp % m != 0) chunk_nsamp = (chunk_nsamp / m + 1) * m;
    }
    long long nchunks = (nsamples + chunk_nsamp - 1) / chunk_nsamp;

    long bytes_per_chunk = (long)((chunk_nsamp * (long long)h.nchans * h.nbits + 7) / 8);
    long samples_per_chunk = (long)chunk_nsamp;
    int  nhist = 1 << h.nbits;

    if (!quiet) {
        fprintf(stderr,
                "Range:   samples [%lld..%lld) = %lld samples (%.3f s)\n"
                "Chunks:  %lld chunks of %lld samples (%.3f s each, %ld bytes packed)\n"
                "Hist:    %d bins per channel (nbits=%d)\n"
                "Zero-DM stored: %s\n",
                start_samp, end_samp, nsamples, nsamples * h.tsamp,
                nchunks, chunk_nsamp, chunk_nsamp * h.tsamp, bytes_per_chunk,
                nhist, h.nbits,
                store_zerodm ? "yes" : "no");
    }

    /* --- Allocations --- */
    uint8_t *packed   = (uint8_t *)malloc((size_t)bytes_per_chunk);
    uint8_t *unpacked = (uint8_t *)malloc((size_t)samples_per_chunk * h.nchans);
    if (!packed || !unpacked) die("out of memory (chunk buffers)");

    /* per-channel running totals */
    uint64_t *tot_sum   = (uint64_t *)calloc((size_t)h.nchans, sizeof(uint64_t));
    uint64_t *tot_sumsq = (uint64_t *)calloc((size_t)h.nchans, sizeof(uint64_t));
    uint64_t *tot_hist  = (uint64_t *)calloc((size_t)h.nchans * nhist, sizeof(uint64_t));
    if (!tot_sum || !tot_sumsq || !tot_hist) die("out of memory (totals)");

    /* per-chunk per-channel float outputs */
    float *cmean = (float *)malloc((size_t)nchunks * h.nchans * sizeof(float));
    float *crms  = (float *)malloc((size_t)nchunks * h.nchans * sizeof(float));
    if (!cmean || !crms) die("out of memory (chunk arrays)");

    /* full-rate zero-DM */
    float *zdm = NULL;
    if (store_zerodm) {
        zdm = (float *)malloc((size_t)nsamples * sizeof(float));
        if (!zdm) die("out of memory (zero-DM); consider -Z");
    }

    /* per-chunk scratch */
    uint64_t *cs_sum   = (uint64_t *)malloc((size_t)h.nchans * sizeof(uint64_t));
    uint64_t *cs_sumsq = (uint64_t *)malloc((size_t)h.nchans * sizeof(uint64_t));
    if (!cs_sum || !cs_sumsq) die("out of memory (scratch)");

    /* --- Main streaming loop --- */
    long long samples_done = 0;
    long long t_global = 0;
    double t0 = now_sec();
    double t_progress = t0;

    for (long long k = 0; k < nchunks; k++) {
        long long want = chunk_nsamp;
        if (samples_done + want > nsamples) want = nsamples - samples_done;
        long long want_bytes = (want * (long long)h.nchans * h.nbits + 7) / 8;

        size_t got = fread(packed, 1, (size_t)want_bytes, fin);
        if (got != (size_t)want_bytes) {
            fprintf(stderr, "rfidiag: short read at chunk %lld (got %zu wanted %lld)\n",
                    k, got, want_bytes);
            break;
        }

        long unpacked_n = unpack_to_uint8(packed, (size_t)got, h.nbits, unpacked);
        if (unpacked_n < 0) die("unpack failed (unsupported nbits)");
        long n_samp = (long)(unpacked_n / h.nchans);
        if (n_samp < want) {
            fprintf(stderr, "rfidiag: warning, partial chunk at k=%lld (%ld < %lld)\n",
                    k, n_samp, want);
            want = n_samp;
        }

        /* Reset per-chunk scratch */
        memset(cs_sum,   0, (size_t)h.nchans * sizeof(uint64_t));
        memset(cs_sumsq, 0, (size_t)h.nchans * sizeof(uint64_t));

        /* Hot loop. Index order: t outer, c inner -> linear access. */
        const int nchans = h.nchans;
        for (long t = 0; t < n_samp; t++) {
            const uint8_t *row = unpacked + (size_t)t * nchans;
            uint64_t z = 0;
            for (int c = 0; c < nchans; c++) {
                uint8_t v = row[c];
                cs_sum[c]   += v;
                cs_sumsq[c] += (uint32_t)v * (uint32_t)v;
                tot_hist[(size_t)c * nhist + v]++;
                z += v;
            }
            if (zdm) zdm[t_global + t] = (float)z;
        }

        /* Per-chunk reductions and accumulate into totals */
        const double inv = 1.0 / (double)n_samp;
        for (int c = 0; c < nchans; c++) {
            double mean   = (double)cs_sum[c]   * inv;
            double meansq = (double)cs_sumsq[c] * inv;
            double var    = meansq - mean * mean;
            if (var < 0.0) var = 0.0;
            cmean[(size_t)k * nchans + c] = (float)mean;
            crms [(size_t)k * nchans + c] = (float)sqrt(var);
            tot_sum[c]   += cs_sum[c];
            tot_sumsq[c] += cs_sumsq[c];
        }

        samples_done += n_samp;
        t_global     += n_samp;

        if (!quiet) {
            double now = now_sec();
            if (now - t_progress > 1.0 || k == nchunks - 1) {
                double frac = (double)samples_done / (double)nsamples;
                double elapsed = now - t0;
                double eta = (frac > 0) ? elapsed * (1.0 / frac - 1.0) : 0.0;
                fprintf(stderr, "\r  %5.1f%%  %lld/%lld chunks   %.1f s elapsed   ~%.0f s remaining   ",
                        100.0 * frac, k + 1, nchunks, elapsed, eta);
                fflush(stderr);
                t_progress = now;
            }
        }
    }
    if (!quiet) fprintf(stderr, "\n");

    /* --- Write output --- */
    DiagHeader dh;
    memset(&dh, 0, sizeof(dh));
    memcpy(dh.magic, DIAG_MAGIC, 8);
    dh.version            = DIAG_VERSION;
    dh.nchans             = (uint32_t)h.nchans;
    dh.nbits              = (uint32_t)h.nbits;
    dh.nchunks            = (uint32_t)nchunks;
    dh.nsamples           = (uint64_t)samples_done;
    dh.nsamples_per_chunk = (uint64_t)chunk_nsamp;
    dh.tsamp              = h.tsamp;
    dh.fch1               = h.fch1;
    dh.foff               = h.foff;
    dh.tstart             = h.tstart;
    dh.chunk_sec          = chunk_nsamp * h.tsamp;
    dh.t_offset           = start_samp * h.tsamp;
    dh.nhist              = (uint32_t)nhist;

    FILE *fout = diag_open_write(out_path, &dh);
    if (!fout) {
        fprintf(stderr, "rfidiag: cannot open output '%s': %s\n", out_path, strerror(errno));
        return 3;
    }

    if (diag_write_section(fout, "SUMC", tot_sum,   (uint64_t)h.nchans * sizeof(uint64_t))   < 0) die("write SUMC");
    if (diag_write_section(fout, "SUMQ", tot_sumsq, (uint64_t)h.nchans * sizeof(uint64_t))   < 0) die("write SUMQ");
    if (diag_write_section(fout, "HIST", tot_hist,  (uint64_t)h.nchans * nhist * sizeof(uint64_t)) < 0) die("write HIST");
    if (diag_write_section(fout, "CMEA", cmean,     (uint64_t)nchunks * h.nchans * sizeof(float))  < 0) die("write CMEA");
    if (diag_write_section(fout, "CRMS", crms,      (uint64_t)nchunks * h.nchans * sizeof(float))  < 0) die("write CRMS");
    if (zdm) {
        if (diag_write_section(fout, "Z0DM", zdm, (uint64_t)samples_done * sizeof(float)) < 0) die("write Z0DM");
    }

    fclose(fout);
    fclose(fin);

    /* --- Cleanup --- */
    free(packed); free(unpacked);
    free(tot_sum); free(tot_sumsq); free(tot_hist);
    free(cmean); free(crms);
    free(zdm);
    free(cs_sum); free(cs_sumsq);

    if (!quiet) {
        double dt = now_sec() - t0;
        double mb = (double)samples_done * h.nchans * h.nbits / 8.0 / (1024.0*1024.0);
        fprintf(stderr, "Done. Processed %.1f MB in %.2f s (%.1f MB/s).\n",
                mb, dt, mb / (dt > 0 ? dt : 1.0));
    }
    return 0;
}
