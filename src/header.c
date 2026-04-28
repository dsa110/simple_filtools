/*
 * header - print parameters of a SIGPROC filterbank file.
 *
 * Standalone reimplementation of dsa110-sigproc's `header` tool, on top
 * of simple_filtools' fil_header reader. With no flags, prints a
 * full human-readable summary; with one or more flags, prints the
 * requested values (one per line) for use in shell scripts.
 *
 *   Usage: header <file.fil> [-tstart ...]
 *
 * Notable differences from sigproc:
 *   - we do not implement the obsolete WAPP / PSPM / BPP raw formats,
 *   - we do not track per-channel `frequency_table`; `-frequencies`
 *     prints fch1 + i*foff for i in 0..nchans-1,
 *   - we do not track `signed`, `gal_l/b`, `scan_number`, `header_tobs`,
 *     `raw_fch1/foff`, `period`, `nbins`, so `-obsdb` and `-scan_number`
 *     are not supported.
 */

#include "filhdr.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------- helpers ported from dsa110-sigproc/{aliases,angle_split,mjd}.c ---------- */

static const char *telescope_name(int telescope_id) {
    switch (telescope_id) {
        case  0: return "Fake";
        case  1: return "Arecibo";
        case  2: return "Ooty";
        case  3: return "Nancay";
        case  4: return "Parkes";
        case  5: return "Jodrell";
        case  6: return "GBT";
        case  7: return "GMRT";
        case  8: return "Effelsberg";
        case  9: return "140ft";
        case 10: return "SRT";
        case 64: return "MeerKAT";
        case 65: return "KAT-7";
        case 66: return "DSA-110";   /* dsa110-sigproc local extension */
        case 82: return "DSA-2000";  /* placeholder; harmless if unused */
        default: return "???????";
    }
}

static const char *backend_name(int machine_id) {
    switch (machine_id) {
        case  0: return "FAKE";
        case  1: return "PSPM";
        case  2: return "WAPP";
        case  3: return "AOFTM";
        case  4: return "BPP";
        case  5: return "OOTY";
        case  6: return "SCAMP";
        case  7: return "GMRTFB";
        case  8: return "PULSAR2000";
        case  9: return "PARSPEC";
        case 10: return "BPSR";
        case 14: return "GMRTNEW";
        case 64: return "KAT";
        case 65: return "KAT-DC2";
        default: return "?????";
    }
}

static const char *data_category(int data_type) {
    switch (data_type) {
        case 0: return "raw data";
        case 1: return "filterbank";
        case 2: return "time series";
        case 3: return "pulse profiles";
        case 4: return "amplitude spectrum";
        case 5: return "complex spectrum";
        case 6: return "dedispersed subbands";
        default: return "unknown!";
    }
}

/* SIGPROC stores RA/Dec as DDMMSS.s packed in a single double:
 *   123456.78  ->  dd=12, mm=34, ss=56.78  */
static void angle_split(double angle, int *dd, int *mm, double *ss) {
    int negative = 0;
    if (angle < 0.0) { angle = -angle; negative = 1; }
    *dd = (int)(angle / 10000.0);
    angle -= (double)(*dd) * 10000.0;
    *mm = (int)(angle / 100.0);
    *ss = angle - 100.0 * (*mm);
    if (negative) *dd = -(*dd);
}

/* MJD -> Gregorian (year, month, day). Fliegel & van Flandern (CACM 1968). */
static void mjd_to_date(double mjd, int *year, int *month, int *day) {
    long jd = (long)floor(mjd + 2400001.0);   /* JD at start of UT day */
    long L = jd + 68569;
    long N = (4 * L) / 146097;
    L = L - (146097 * N + 3) / 4;
    long I = (4000 * (L + 1)) / 1461001;
    L = L - (1461 * I) / 4 + 31;
    long J = (80 * L) / 2447;
    *day   = (int)(L - (2447 * J) / 80);
    L = J / 11;
    *month = (int)(J + 2 - 12 * L);
    *year  = (int)(100 * (N - 49) + I + L);
}

/* ---------- cli ---------- */

static void usage(FILE *fp, const char *argv0) {
    fprintf(fp,
        "Usage: %s <file.fil> [flag ...]\n"
        "\n"
        "With no flags, prints a full human-readable summary. Otherwise prints\n"
        "the requested value(s), one per line. Supported flags:\n"
        "\n"
        "  -source_name      -telescope        -machine          -datatype\n"
        "  -data_type        -frame            -barycentric      -pulsarcentric\n"
        "  -headersize       -datasize         -nsamples         -tobs\n"
        "  -fch1             -foff             -bandwidth        -fmid\n"
        "  -nchans           -nbits            -nifs             -nbeam\n"
        "  -ibeam            -tstart           -tsamp            -mjd\n"
        "  -utstart          -date             -src_raj          -src_dej\n"
        "  -ra_deg           -dec_deg          -az_start         -za_start\n"
        "  -refdm   (alias -dm)\n"
        "  -frequencies      -h, --help\n",
        argv0);
}

static void print_full(FILE *fp, const char *path, const fil_header *h) {
    int year, month, day;
    int rah, ram, ded, dem;
    double ras, des;
    char sra[8], sde[8], decsign;

    angle_split(h->src_raj, &rah, &ram, &ras);
    angle_split(h->src_dej, &ded, &dem, &des);
    decsign = (h->src_dej >= 0.0) ? '+' : '-';
    mjd_to_date(h->tstart, &year, &month, &day);

    if (ras < 10.0) snprintf(sra, sizeof(sra), "0%.1f", ras);
    else            snprintf(sra, sizeof(sra), "%.1f",  ras);
    if (des < 10.0) snprintf(sde, sizeof(sde), "0%.1f", des);
    else            snprintf(sde, sizeof(sde), "%.1f",  des);

    fprintf(fp, "Data file                        : %s\n", path);
    fprintf(fp, "Header size (bytes)              : %ld\n", h->header_bytes);
    if (h->data_bytes > 0)
        fprintf(fp, "Data size (bytes)                : %lld\n", h->data_bytes);

    const char *frame =
        h->pulsarcentric ? "(pulsarcentric)" :
        h->barycentric   ? "(barycentric)"   : "(topocentric)";
    fprintf(fp, "Data type                        : %s %s\n",
            data_category(h->data_type), frame);

    fprintf(fp, "Telescope                        : %s\n", telescope_name(h->telescope_id));
    fprintf(fp, "Datataking Machine               : %s\n", backend_name(h->machine_id));

    if (h->source_name[0])
        fprintf(fp, "Source Name                      : %s\n", h->source_name);
    if (h->src_raj != 0.0)
        fprintf(fp, "Source RA (J2000)                : %02d:%02d:%s\n", rah, ram, sra);
    if (h->src_dej != 0.0)
        fprintf(fp, "Source DEC (J2000)               : %c%02d:%02d:%s\n",
                decsign, abs(ded), dem, sde);
    if (h->az_start != 0.0 && h->az_start != -1.0)
        fprintf(fp, "Start AZ (deg)                   : %f\n", h->az_start);
    if (h->za_start != 0.0 && h->za_start != -1.0)
        fprintf(fp, "Start ZA (deg)                   : %f\n", h->za_start);

    switch (h->data_type) {
        case 0:
        case 1:
            fprintf(fp, "Frequency of channel 1 (MHz)     : %f\n", h->fch1);
            fprintf(fp, "Channel bandwidth      (MHz)     : %f\n", h->foff);
            fprintf(fp, "Number of channels               : %d\n", h->nchans);
            fprintf(fp, "Number of beams                  : %d\n", h->nbeams);
            fprintf(fp, "Beam number                      : %d\n", h->ibeam);
            break;
        case 2:
            fprintf(fp, "Reference DM (pc/cc)             : %f\n", h->refdm);
            fprintf(fp, "Reference frequency    (MHz)     : %f\n", h->fch1);
            break;
        case 6:
            fprintf(fp, "Reference DM (pc/cc)             : %f\n", h->refdm);
            fprintf(fp, "Frequency of channel 1 (MHz)     : %f\n", h->fch1);
            fprintf(fp, "Channel bandwidth      (MHz)     : %f\n", h->foff);
            fprintf(fp, "Number of channels               : %d\n", h->nchans);
            break;
    }

    fprintf(fp, "Time stamp of first sample (MJD) : %.12f\n", h->tstart);
    fprintf(fp, "Gregorian date (YYYY/MM/DD)      : %4d/%02d/%02d\n", year, month, day);
    fprintf(fp, "Sample time (us)                 : %.5f\n", h->tsamp * 1.0e6);

    if (h->nsamples > 0) {
        fprintf(fp, "Number of samples                : %lld\n", h->nsamples);
        double tobs = (double)h->nsamples * h->tsamp;
        const char *unit = "(seconds)   ";
        if (tobs > 60.0) { tobs /= 60.0; unit = "(minutes)   ";
            if (tobs > 60.0) { tobs /= 60.0; unit = "(hours)     ";
                if (tobs > 24.0) { tobs /= 24.0; unit = "(days)      "; } } }
        fprintf(fp, "Observation length %s  : %.1f\n", unit, tobs);
    }
    fprintf(fp, "Number of bits per sample        : %d\n", h->nbits);
    fprintf(fp, "Number of IFs                    : %d\n", h->nifs);
}

static int handle_flag(const char *arg, const fil_header *h) {
    int rah, ram, ded, dem;
    double ras, des;
    char sra[8], sde[8], decsign;

    if (!strcmp(arg, "-source_name"))   { puts(h->source_name); return 0; }
    if (!strcmp(arg, "-telescope"))     { puts(telescope_name(h->telescope_id)); return 0; }
    if (!strcmp(arg, "-machine"))       { puts(backend_name(h->machine_id)); return 0; }
    if (!strcmp(arg, "-datatype"))      { puts(data_category(h->data_type)); return 0; }
    if (!strcmp(arg, "-data_type"))     { printf("%d\n", h->data_type); return 0; }
    if (!strcmp(arg, "-frame")) {
        puts(h->pulsarcentric ? "pulsarcentric" :
             h->barycentric   ? "barycentric"   : "topocentric");
        return 0;
    }
    if (!strcmp(arg, "-barycentric"))   { printf("%d\n", h->barycentric); return 0; }
    if (!strcmp(arg, "-pulsarcentric")) { printf("%d\n", h->pulsarcentric); return 0; }
    if (!strcmp(arg, "-headersize"))    { printf("%ld\n", h->header_bytes); return 0; }
    if (!strcmp(arg, "-datasize"))      { printf("%lld\n", h->data_bytes); return 0; }
    if (!strcmp(arg, "-nsamples"))      { printf("%lld\n", h->nsamples); return 0; }
    if (!strcmp(arg, "-tobs"))          { printf("%f\n", (double)h->nsamples * h->tsamp); return 0; }
    if (!strcmp(arg, "-fch1"))          { printf("%.3f\n", h->fch1); return 0; }
    if (!strcmp(arg, "-foff"))          { printf("%f\n", h->foff); return 0; }
    if (!strcmp(arg, "-bandwidth"))     { printf("%.3f\n", fabs(h->foff) * (double)h->nchans); return 0; }
    if (!strcmp(arg, "-fmid"))          { printf("%.3f\n", h->fch1 + h->foff * h->nchans / 2.0); return 0; }
    if (!strcmp(arg, "-nchans"))        { printf("%d\n", h->nchans); return 0; }
    if (!strcmp(arg, "-nbits"))         { printf("%d\n", h->nbits); return 0; }
    if (!strcmp(arg, "-nifs"))          { printf("%d\n", h->nifs); return 0; }
    if (!strcmp(arg, "-nbeam"))         { printf("%d\n", h->nbeams); return 0; }
    if (!strcmp(arg, "-ibeam"))         { printf("%d\n", h->ibeam); return 0; }
    if (!strcmp(arg, "-tstart"))        { printf("%.12f\n", h->tstart); return 0; }
    if (!strcmp(arg, "-tsamp"))         { printf("%.5f\n", h->tsamp * 1.0e6); return 0; }
    if (!strcmp(arg, "-mjd"))           { printf("%d\n", (int)floor(h->tstart)); return 0; }
    if (!strcmp(arg, "-az_start"))      { printf("%f\n", h->az_start); return 0; }
    if (!strcmp(arg, "-za_start"))      { printf("%f\n", h->za_start); return 0; }
    if (!strcmp(arg, "-refdm") || !strcmp(arg, "-dm")) { printf("%f\n", h->refdm); return 0; }
    if (!strcmp(arg, "-utstart")) {
        double frac = h->tstart - floor(h->tstart);
        int uth = (int)floor(24.0 * frac);  frac -= (double)uth / 24.0;
        int utm = (int)floor(1440.0 * frac); frac -= (double)utm / 1440.0;
        int uts = (int)floor(86400.0 * frac);
        printf("%02d:%02d:%02d\n", uth, utm, uts);
        return 0;
    }
    if (!strcmp(arg, "-date")) {
        int year, month, day;
        mjd_to_date(h->tstart, &year, &month, &day);
        printf("%4d/%02d/%02d\n", year, month, day);
        return 0;
    }
    if (!strcmp(arg, "-src_raj")) {
        angle_split(h->src_raj, &rah, &ram, &ras);
        if (ras < 10.0) snprintf(sra, sizeof(sra), "0%.1f", ras);
        else            snprintf(sra, sizeof(sra), "%.1f",  ras);
        printf("%02d:%02d:%s\n", rah, ram, sra);
        return 0;
    }
    if (!strcmp(arg, "-src_dej")) {
        angle_split(h->src_dej, &ded, &dem, &des);
        decsign = (h->src_dej >= 0.0) ? '+' : '-';
        if (des < 10.0) snprintf(sde, sizeof(sde), "0%.1f", des);
        else            snprintf(sde, sizeof(sde), "%.1f",  des);
        printf("%c%02d:%02d:%s\n", decsign, abs(ded), dem, sde);
        return 0;
    }
    if (!strcmp(arg, "-ra_deg")) {
        angle_split(h->src_raj, &rah, &ram, &ras);
        printf("%f\n", rah * 15.0 + ram / 4.0 + ras / 240.0);
        return 0;
    }
    if (!strcmp(arg, "-dec_deg")) {
        angle_split(h->src_dej, &ded, &dem, &des);
        decsign = (h->src_dej >= 0.0) ? '+' : '-';
        printf("%c%f\n", decsign, abs(ded) + dem / 60.0 + des / 3600.0);
        return 0;
    }
    if (!strcmp(arg, "-frequencies")) {
        for (int j = 0; j < h->nchans; j++)
            printf("%f\n", h->fch1 + j * h->foff);
        return 0;
    }
    return -1;
}

int main(int argc, char **argv) {
    if (argc < 2) { usage(stderr, argv[0]); return 1; }
    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        usage(stdout, argv[0]); return 0;
    }
    const char *path = argv[1];

    fil_header h;
    int err = 0;
    char errmsg[256];
    FILE *fp = fil_open(path, &h, &err, errmsg, sizeof(errmsg));
    if (!fp) { fprintf(stderr, "header: %s\n", errmsg); return 2; }
    fclose(fp);

    if (argc == 2) {
        print_full(stdout, path, &h);
        return 0;
    }

    int rc = 0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            usage(stdout, argv[0]);
            continue;
        }
        if (handle_flag(argv[i], &h) < 0) {
            fprintf(stderr, "header: unknown flag '%s'\n", argv[i]);
            usage(stderr, argv[0]);
            rc = 1;
        }
    }
    return rc;
}
