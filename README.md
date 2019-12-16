# MGSC
Molecular Ground State Energy Calculator

This repository contains three programs of interest. All of them will be built when running `build.sh` from the `src/` directory. You'll need the intel compiler for this script to work. You'll also need the GNU Scientific Library (GSL) installed. 
  
You'll want to run all of the following while in the `src/` directory.
   
**sweep**  
This program computes the energy of a single term wavefunction with radial symmetry for a wide array of values. You'll want to run it with something like

```bash
./sweep > ../results/sweep
```

Once it has finished, you can graph the results with

```bash
python3 graph.py ../results/sweep
```
This program filters out some of the extra points that make it hard to see the interesting parts of the curve.

**min**  
This program performs the conjugate gradient minimization that was mentioned in my paper and in the speech. Run it with

```bash
./min > ../results/min
```

I reduced the grid of initial values from 64 by 64 to 16 by 16. This should make it run in a reasonable amount of time. Graph the results with

```bash
python3 grid_graph.py ../results/min
```

**mgsc**  
This is a very basic monte carlo algorithm for attempting to find the real ground state energy. It works by guessing random terms to add to the wavefunction and accepting or rejecting them based on whether or not they lower the energy. The best result I've obtained so far using this method is `-13.2525 eV`. 

```bash
./mgsc
```

## Documentation

The `doc/` directory contains several latex files that have documentation on the math implemented here, as well as some ideas for the future. Not everything is implemented in the code exactly the same way it is written in the documents. For example, I kept the electron-nucleus integrals in cartesian coordinates, but they are in polar coordinates in their final form in the documentation.
