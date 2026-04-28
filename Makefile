# simple_filtools - standalone RFI diagnostics for SIGPROC filterbank data
#
# Build:
#   make              # builds rfidiag
#   make test         # builds + runs unit and end-to-end tests (needs python+numpy)
#   make clean

CC      ?= cc
CFLAGS  ?= -O3 -Wall -Wextra -std=c99
LDLIBS  ?= -lm

HDRS    := src/filhdr.h src/unpack.h src/diag.h

RFIDIAG_SRC   := src/filhdr.c src/unpack.c src/diag.c src/rfidiag.c
CHOP_SRC      := src/filhdr.c src/chop_fil.c
HEADER_SRC    := src/filhdr.c src/header.c

OBJ     := $(sort $(RFIDIAG_SRC:.c=.o) $(CHOP_SRC:.c=.o) $(HEADER_SRC:.c=.o))

BINS    := rfidiag chop_fil header

PYTHON  ?= python3

.PHONY: all clean test test-unpack test-numpy

all: $(BINS)

rfidiag: $(RFIDIAG_SRC:.c=.o)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

chop_fil: $(CHOP_SRC:.c=.o)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

header: $(HEADER_SRC:.c=.o)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

src/%.o: src/%.c $(HDRS)
	$(CC) $(CFLAGS) -Isrc -c $< -o $@

tests/test_unpack: tests/test_unpack.c src/unpack.c src/unpack.h
	$(CC) $(CFLAGS) -Isrc -o $@ tests/test_unpack.c src/unpack.c $(LDLIBS)

test-unpack: tests/test_unpack
	./tests/test_unpack

test-numpy: rfidiag
	$(PYTHON) tests/check_against_numpy.py

test: test-unpack test-numpy
	@echo "All tests passed."

clean:
	rm -f $(OBJ) $(BINS) tests/test_unpack
	rm -rf tests/_tmp
