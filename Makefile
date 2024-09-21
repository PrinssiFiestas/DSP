# MIT License
# Copyright (c) 2024 Lauri Lorenzo Fiestas
# https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

CC = gcc
CFLAGS = -Wall -Wextra -Werror -Wpedantic -std=c99 -Iinclude -Itests/include -O3
SRCS = tests/pnstest.c src/iir.c

build/pnstest.x64: $(SRCS) include/dsp/iir.h
	mkdir -p build
	$(CC) -o $@ $(SRCS) $(CFLAGS) -s -shared -Wl,--subsystem,windows

clean:
	rm -rf build
