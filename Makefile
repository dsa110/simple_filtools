# simple_filtools - standalone RFI diagnostics for SIGPROC filterbank data
#
# Build:
#   make              # builds rfidiag
#   make test         # builds + runs unit and end-to-end tests (needs python+numpy)
#   make clean

CC      ?= cc
CFLAGS  ?= -O3 -Wall -Wextra -std=c99
LDLIBS  ?= -lm

SRC     := src/filhdr.c src/unpack.c src/diag.c src/rfidiag.c
HDRS    := src/filhdr.h src/unpack.h src/diag.h
OBJ     := $(SRC:.c=.o)
BIN     := rfidiag

PYTHON  ?= python3

.PHONY: all clean test test-unpack test-numpy

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDLIBS)

src/%.o: src/%.c $(HDRS)
	$(CC) $(CFLAGS) -Isrc -c $< -o $@

tests/test_unpack: tests/test_unpack.c src/unpack.c src/unpack.h
	$(CC) $(CFLAGS) -Isrc -o $@ tests/test_unpack.c src/unpack.c $(LDLIBS)

test-unpack: tests/test_unpack
	./tests/test_unpack

test-numpy: $(BIN)
	$(PYTHON) tests/check_against_numpy.py

test: test-unpack test-numpy
	@echo "All tests passed."

clean:
	rm -f $(OBJ) $(BIN) tests/test_unpack
	rm -rf tests/_tmp
