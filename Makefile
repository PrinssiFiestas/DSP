# MIT License
# Copyright (c) 2024 Lauri Lorenzo Fiestas
# https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

CC = gcc
CFLAGS = -Wall -Wextra -Werror -Iinclude -Itests/include -O3

build/pnstest.x64: tests/pnstest.c
	mkdir -p build
	$(CC) -o $@ $< $(CFLAGS) -s -shared -Wl,--subsystem,windows
