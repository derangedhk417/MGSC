#!/bin/bash
icc -lgsl mgsc.cpp -o mgsc -O3
icc -lgsl sweep.cpp -o sweep -O3
icc -lgsl grid.cpp -o grid -O3