/*
 * diag.c - writer for .diag files emitted by rfidiag.
 */
#include "diag.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

FILE *diag_open_write(const char *path, const DiagHeader *hdr) {
    FILE *fp = fopen(path, "wb");
    if (!fp) return NULL;
    if (fwrite(hdr, sizeof(*hdr), 1, fp) != 1) {
        fclose(fp);
        return NULL;
    }
    return fp;
}

int diag_write_section(FILE *fp, const char tag[4],
                       const void *data, uint64_t nbytes) {
    char tbuf[4] = {0,0,0,0};
    memcpy(tbuf, tag, 4);
    if (fwrite(tbuf, 1, 4, fp) != 4) return -1;
    if (fwrite(&nbytes, sizeof(nbytes), 1, fp) != 1) return -1;
    if (nbytes > 0) {
        if (fwrite(data, 1, nbytes, fp) != nbytes) return -1;
    }
    return 0;
}
