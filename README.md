# hpc-assignment1
The first assignment for HPC course

Compile using the following command:

`gcc -fopenmp -o galsim galsim.c file_operations.c graphics.c -lm -lX11`

Run using the following command:

`./galsim [number of planets] [initial condition filename] [nsteps] [timestep] [isWithGraphics] [threadNum]`
