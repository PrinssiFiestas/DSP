# MIT License
# Copyright (c) 2024 Lauri Lorenzo Fiestas
# https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

CC = gcc
CFLAGS = -Wall -Wextra -Werror -Wpedantic -std=c99 -Iinclude -Itests/include
SRCS = tests/pnstest.c src/iir.c

.PHONY: all release perf debug clean

all: release

release: build/pnstest.x64
build/pnstest.x64: $(SRCS) include/dsp/iir.h
	mkdir -p build
	$(CC) -o $@ $(SRCS) $(CFLAGS) -O3 -s -shared -Wl,--subsystem,windows

perf: CFLAGS += -DPERFORMANCE_TEST
perf: build/pnstest.x64

debug: build/pnstestd.x64
build/pnstestd.x64: $(SRCS) include/dsp/iir.h
	mkdir -p build
	$(CC) -o $@ $(SRCS) $(CFLAGS) -O0 -ggdb3 -gdwarf -shared -Wl,--subsystem,windows

clean:
	rm -rf build
