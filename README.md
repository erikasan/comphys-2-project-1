# Project 1 in FYS4411 by Erik Alexander Sandvik & Simon Elias Schrader

## How to run our code
 In the folder code/cpp/, all executables can be compiled the following way
 ```bash
# Create build-directory
mkdir build

# Move into the build-directory
cd build

# Run CMake to create a Makefile
cmake ../

# Make the Makefile using two threads
make -j2
```
this will produce 5 executables, *test*, *vmc*, *vmc_parallel*, *gradientdescent* and *gradientdescent_parallel*.
These are run the following way (except for test which does not take arguments):
 ```bash
 #Parallel files:
 ./filename numberOfDimensions numberOfParticles numberOfMCCycles alpha stepLength numberOfEquilibrationSteps seed wF_type sampler_type Filename_Blocking numberOfThreads (path)
 
# Nonparallel files:
  ./filename numberOfDimensions numberOfParticles numberOfMCCycles alpha stepLength numberOfEquilibrationSteps seed wF_type sampler_type Filename_Blocking (path)
```
### Path has a default, which is "../../../output/". If that folder does not exist, one has either to create it or use a different path!
where
- wF_type is either "HO" (Harmonic Oscillator), "NHO" (Numerical Harmonic Oscillator) or "EO" (Elliptic Oscillator), where only Elliptical Oscillator contains the repulsion term
- sampler_type is "VMC" (Brute force) or "IMP" (Importance Sampling)
- Filename_Blocking is the filename of the blocking data written to file (one file per thread, gets an extra number after the file name). Use "no" if no output is desired. There will still be output to "sympleharmonic.csv". 
- alpha is either the initial guess (for gradient descent) or the actual value (for VMC). 
- For the parallel simulations, each thread gets the seed (input_seed+thread_id) to keep the process sufficiently random.

If any other parameter should be tuned, such as writing out the position distribution (which can only be done in "vmc") or changing other parameters, this needs to be done directly in the source code for the corresponding main files.
There is a bunch of python files that automatize runs and have several runs run at once, found in the code/python folder.
