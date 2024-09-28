# MIT License
# Copyright (c) 2024 Lauri Lorenzo Fiestas
# https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

CC = gcc
CFLAGS = -Wall -Wextra -Werror -Wpedantic -std=c99 -Iinclude -Itests/include
SRCS = src/iir.c src/fft.c

.PHONY: all release perf debug tests clean

all: release

release: build/pnstest.x64
build/pnstest.x64: $(SRCS) tests/pnstest.c include/dsp/iir.h
	mkdir -p build
	$(CC) -o $@ $(SRCS) tests/pnstest.c $(CFLAGS) -O3 -flto -s -shared -Wl,--subsystem,windows

perf: CFLAGS += -DPERFORMANCE_TEST
perf: build/pnstest.x64

debug: build/pnstestd.x64
build/pnstestd.x64: $(SRCS) tests/pnstest.c include/dsp/iir.h
	mkdir -p build
	$(CC) -o $@ $(SRCS) tests/pnstest.c $(CFLAGS) -O0 -ggdb3 -gdwarf -shared -Wl,--subsystem,windows

tests: build/tests.exe
build/tests.exe: $(SRCS) tests/tests.c include/dsp/fft.h
	mkdir -p build
	$(CC) -o $@ tests/tests.c $(CFLAGS) -lgpc -lm -lpthread -ggdb3 -gdwarf

tests:
	./build/tests.exe

clean:
	rm -rf build
