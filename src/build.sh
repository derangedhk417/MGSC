#!/bin/bash
IMPL_FILES="gaussian.cpp"
icc -lgsl mgsc.cpp -o mgsc -O3 $IMPL_FILES