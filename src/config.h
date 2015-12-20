#ifndef _CONFIG_H_
#define _CONFIG_H_
#define DIM     3

/* if more atoms are added atomtypes needs to be changed*/
#define atypes  2
#define Cl	1
#define Na	0

#define kB                  1.38064852 // Boltzmann constant (the exponent has been accounted for in the code)
#define Navo                6.023 // Avogadro's number (the exponent has been accounted for in the code)
#define k_e                 2.3*pow(10,-10) // Electrostatic constant for Coulomb's law
#define cubeFactor          2.0 // Factor scaling between cutoff and cubesize
#define gridFactor          2.0 // Factor scaling between cutoff and grid size
#define reassignCubesFreq   1000 // Frequency to reassign the cubes for neighbor lists
#define numKVectors         10 // Maximum number of k vectors to use in reciprocal space
#define alphaForGamma       1.0*pow(10,4) // Parameters to shape the error function which dampens realspace electrostatics, and shapes the Gaussian for reciprocal space electrostatics (nm-1)
#endif
