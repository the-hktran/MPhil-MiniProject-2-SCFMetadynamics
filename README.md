Self-Consistent Field Metadynamics
==============

Mini Project 2<br />
MPhil in Scientific Computing<br />
Churchill College<br />
University of Cambridge

*This work was produced as part of the miniproject requirements of the MPhil in Scientific Computing course I undertook as a student of Churchill College, University of Cambridge. This work was done under the supervision of Dr. Alex Thom and with funding from the Sir Winston Churchill Foundation of the USA.*

## Introduction

This is a program designed to calculate multiple SCF solutions using a technique inspired by metadynamics, following the original publication by Thom and Head-Gordon (Phys. Rev. Lett. 101, 193001 (2008)).

## Compiling

The SCF Metadynamics program uses the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library for linear algebra. The program also utilizes the C++11 standard. The following should allow for compilation of the program.

```
% g++ -std=c++11 -I /path/to/eigen Metadynamics.cpp SCF.cpp Fock.cpp ReadInput.cpp -O3 -o SCFMD
```

## Running the Program

### Single Point

The program takes three command line inputs. These are, in order, the input filename, the overlap matrix filename, the output filename. The overlap matrix filename is not yet implemented, as the program only takes output from Q-Chem where the overlap matrix is known to be the identity, so the second input does not matter. The first input is the file that contains all the settings and values of the integrals. The third input is the filename of the output. Alternatively the program can be simply ran without any command line inputs. A prompt will ask for these to be input individually.

### Scan

A scan can be ran by running the program with no command line inputs. Enter the input filename, overlap filename (not yet implemented), and output filename. When prompted about running a scan, enter "1" and a prompt will appear asking for parameters regarding the scan. These are, in order, the integer that labels the first step, the integer that labels the last step, the parameter associated with the first step, and the step size of the parameter. More details about the naming of the input files will be discussed in the following section.

## Format of the Input File

The input file is formatted as such. First, a string of input parameters should be listed. In order, they are:
- (Integer, Positive) The number of spacial orbitals.
- (Integer, Positive) The number of electrons.
- (Integer, Positive) The number of solutions desired.
- (Integer, 1 / 0) Option to use / not use DIIS error.
- (Integer, 1 / 0) Option to use / not use MOM.
- (Integer, 2 / 1 / 0) Option to change converged density matrix by Randomization / Randomization with unity trace / Rotation of orbitals.
- (Integer, Positive) Maximum number of SCF cycles.
- (Double, Positive) Starting Norm of the biasing potential.
- (Double, Positive) Starting Lambda of the biasing potential.

Next, the values for the two electron integrals (nm|kl) are listed in the following format
```
(nm|kl)     n     m     k     l
```
It should be noted that the nuclear repulsion has n, m, k, and l set to zero and the one electron integrals are labelled by n and m while k and l are set to zero. This is the format of Q-Chem.

### Scan Input File

More care has to be taken for the filename of the scan input files. The content of each input file is the same, but the filename should have the following format.
```
(BASEFILENAME)_(LABEL_INTEGER)
```
BASEFILENAME is the filename that is input into the program. LABEL\_INTEGER is an integer between the label of the first and the label of the last step of the scan. For example, if a series of input files have been generated with the filenames input\_1, ..., input\_5 corresponding to a bond length of 1.00, 1.50, 2.00, 2.50, 3.00, then the input filename is "input" with starting integer "1" and ending integer "5" where the starting value is "1.00" and the step size is "0.50."
