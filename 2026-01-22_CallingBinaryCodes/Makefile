# Detect OS
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
	# macOS
	CC ?= gcc
	CFLAGS ?= -Wall -Wextra -O2 -fPIC -std=c99 -Wno-nullability-completeness
	LIBFLAGS = -dynamiclib
else
	# Linux (and others)
	CC ?= gcc
	CFLAGS ?= -Wall -Wextra -O2 -fPIC -std=c99
	LIBFLAGS = -shared
endif

LDLIBS ?= -ldl

TARGET = binary_playground
SRC = binary_playground.c
HDR = binary_playground.h
LIBNAME = libbinary_playground.dylib

.PHONY: all lib clean run runlib

all: $(TARGET)

$(TARGET): $(SRC) $(HDR)
	$(CC) $(CFLAGS) -o $@ $(SRC) $(LDLIBS)

lib: $(SRC) $(HDR)
    # Build a dynamic library (macOS .dylib) and a Linux-style .so
    # Use $(CFLAGS) so warning-suppression flags are applied to library builds.
	$(CC) $(CFLAGS) -dynamiclib -o $(LIBNAME) $(SRC) $(LDLIBS) || true
    # Also build a shared object for Linux compatibility. If the platform
    # does not support `-shared`, ignore the failure of this step.
	$(CC) $(CFLAGS) -shared -fPIC -o libbinary_playground.so $(SRC) $(LDLIBS) || true

run: $(TARGET)
	./$(TARGET)

runlib: lib $(TARGET)
	# Runs the program using the dynamic library via --use-lib
	./$(TARGET) --use-lib --mode scalar --a 1 --b 2 --c 3

clean:
	rm -f $(TARGET) $(LIBNAME) libbinary_playground.so