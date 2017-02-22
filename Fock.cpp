#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <string>

/* This function calculated <mn||kl> and takes as arguments, the orbital numbers, the spin of the orbitals, and the map 
   Don't know if this will be needed. I'm just leaving it here for now. */
double TwoElectronIntegral(unsigned short int m, unsigned short int n, unsigned short int k, unsigned short int l, bool m_isAlpha, bool n_isAlpha, bool k_isAlpha, bool l_isAlpha, std::map<std::string, double> &Integrals)
{
    double mknl = 0; // First term. (mk|nl)
    double mlnk = 0; // Second term. (ml|nk)
    
    /* Deal with first term first */
    if((m_isAlpha != k_isAlpha) || (n_isAlpha != l_isAlpha)) // Means spin component is different.
    {
        mknl = 0;
    }
    else
    {
        mknl = Integrals[std::to_string(m) + " " + std::to_string(k) + " " + std::to_string(n) + " " + std::to_string(l)];
    }
    /* Now, the second term */
    if((m_isAlpha != l_isAlpha) || (n_isAlpha != k_isAlpha))
    {
        mlnk = 0;
    }
    else
    {
        mlnk = Integrals[std::to_string(m) + " " + std::to_string(l) + " " + std::to_string(n) + " " + std::to_string(k)];
    }
    return mknl - mlnk;
}

double ExchangeTerm(int m, int n, Eigen::MatrixXd &DensityMatrix, std::map<std::string, double> &Integrals)
{
    double XTerm = 0; // Sum_ij D_ij [(mn|ij) - 0.5 (mj|in)]
    for(int i = 0; i < DensityMatrix.rows(); i++)
    {
        for(int j = 0; j < DensityMatrix.cols(); j++)
        {
            XTerm += DensityMatrix(i, j) * (2 * Integrals[std::to_string(m + 1) + " " + std::to_string(n + 1) + " " + std::to_string(i + 1) + " " + std::to_string(j + 1)]
                                              - Integrals[std::to_string(m + 1) + " " + std::to_string(i + 1) + " " + std::to_string(j + 1) + " " + std::to_string(n + 1)]);
        }
    }
    return XTerm;
}

/* Takes a matrix of doubles as a reference and stores FockMatrix into it. 

   Variables 
   FockMatrix - Stores FockMatrix. Passed by reference. 
   DensityMatrix - Density matrix.
   Integrals - Maps to list of two electron integrals from QChem. */
void BuildFockMatrix(Eigen::MatrixXd &FockMatrix, Eigen::MatrixXd &DensityMatrix, std::map<std::string, double> &Integrals)
{
    for(int m = 0; m < FockMatrix.rows(); m++)
    {
        for(int n = m; n < FockMatrix.cols(); n++)
        {
            FockMatrix(m, n) = Integrals[std::to_string(m + 1) + " " + std::to_string(n + 1) + " 0 0"] 
                             + ExchangeTerm(m, n, DensityMatrix, Integrals);
            FockMatrix(n, m) = FockMatrix(m, n);
        }
    }
}
