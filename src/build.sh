#!/bin/bash
IMPL_FILES="gaussian.cpp integral_en.cpp"
icc -lgsl sweep.cpp -o sweep -O3 $IMPL_FILES
icc -lgsl   min.cpp -o min   -O3 $IMPL_FILES