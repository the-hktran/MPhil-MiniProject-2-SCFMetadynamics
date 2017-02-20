#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include "ReadInput.h"

void InputObj::GetInputName()
{
    std::cout << "SCF MetaD: Enter integral filename:\nFCI: ";
    std::cin >> IntegralsInput;
    std::cout << "SCF MetaD: Enter overlap filename:\nFCI: ";
    std::cin >> OverlapInput;
    std::cout << "SCF MetaD: Enter output filename:\nFCI: ";
    std::cin >> OutputName;
}

void InputObj::SetNames(char* Int, char* Overlap, char* Out)
{
    IntegralsInput = Int;
    OverlapInput = Overlap;
    OutputName = Out;
}

/*****                                    FORMAT OF THE INPUT FILE                                 *****/
/* The input file is taken from the integrals computed in QChem (version 4.4). To obtain these integrals,
   the line "CORRELATION idump" must be added to the %rem section of the QChem input. This outputs a list
   of integrals in the format 
                (ij|kl)     i   j   k   l
   in the file titled "INTDUMP". Further modification must be made to this file in order to use it in my
   program. The first two lines must be deleted and the following entered (separated by spaces)
            Number of alpha electrons
            Number of orbitals to be occupied by the alpha electrons
            Number of beta electrons
            Number of orbitals to be occupied by the beta electrons
            Number of solutions desired from Davidson Diagonalization
   As an example, here is the first few lines of an input file for H2, which has one alpha and one beta 
   electron and a space of four orbitals and we are looking for 2 solutions.
            1 4 1 4 2
                0.64985185942031   1   1   1   1
                0.16712550470738   1   3   1   1
                0.080102886434995  1   2   1   2
                0.07936780580498   1   4   1   2
                (And the rest of the integrals)                                                        */
void InputObj::Set()
{
    std::ifstream IntegralsFile(IntegralsInput.c_str());
    double tmpDouble;
    int tmpInt1, tmpInt2, tmpInt3, tmpInt4;
    while(!IntegralsFile.eof())
    {
        /* We have to include all 8-fold permuation symmetries. This holds each integral in chemistry notation. We represent
        (ij|kl) as "i j k l". h_ij is "i j 0 0", as given in QChem. */
        IntegralsFile >> tmpDouble >> tmpInt1 >> tmpInt2 >> tmpInt3 >> tmpInt4;
        Integrals[std::to_string(tmpInt1) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt4)] = tmpDouble;
        Integrals[std::to_string(tmpInt3) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt2)] = tmpDouble;
        Integrals[std::to_string(tmpInt2) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt3)] = tmpDouble;
        Integrals[std::to_string(tmpInt4) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt1)] = tmpDouble;
        Integrals[std::to_string(tmpInt2) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt4)] = tmpDouble;
        Integrals[std::to_string(tmpInt4) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt2)] = tmpDouble;
        Integrals[std::to_string(tmpInt1) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt3)] = tmpDouble;
        Integrals[std::to_string(tmpInt3) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt1)] = tmpDouble;
    }

    std::ifstream OverlapFile(OverlapInput.c_str());
}