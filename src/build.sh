#!/bin/bash
IMPL_FILES="gaussian.cpp integral_en.cpp"
icc -lgsl mgsc.cpp -o mgsc -O3 $IMPL_FILES