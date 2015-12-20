# SoftwareEngg_Fall15
Final Project code for Software Engg course I took in Fall 2015 (CS294-73) by Phil Colella. 
The software can perform molecular dynamics simulations of Na/Cl atoms. The PDF files in the directory, 
particularly in src/finalproj_writeup.pdf has information about the implementation. Also documentation 
created using doxygen can be found in src/documentation.pdf.

Running the Code

We have provided a Makefile (GNUMakefile) that can be used to compile the code on
any machine, provided the correct compiler and FFTW libraries are specified in the
Makefile. Default input files for running the simulation have been provided (runtime.cfg
and SodiumChloride.prm). The executable that is created is MD.exe. Once the code is
compiled, the user needs to type in ./MD.exe –p runtime.cfg. This should run the
simulation.
