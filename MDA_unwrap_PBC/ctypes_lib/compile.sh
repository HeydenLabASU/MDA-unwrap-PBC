#!/bin/bash

gcc -fpic -c unwrap.c
gcc -shared unwrap.o -o unwrap.so
