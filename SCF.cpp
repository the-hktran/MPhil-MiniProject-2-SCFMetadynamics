#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <string>
#include <map>
#include "ReadInput.h"
#include <Eigen/Sparse>

void BuildFockMatrix(Eigen::MatrixXd &FockMatrix, Eigen::MatrixXd &DensityMatrix, std::map<std::string, double> &Integrals);

double CalcDensityRMS(Eigen::MatrixXd &DensityMatrix, Eigen::MatrixXd &DensityMatrixPrev)
{
    double DensityRMS = 0;
    for(int i = 0; i < DensityMatrix.rows(); i++)
    {
        for(int j = 0; j < DensityMatrix.cols(); j++)
        {
            DensityRMS += (DensityMatrix(i, j) - DensityMatrixPrev(i, j)) * (DensityMatrix(i, j) - DensityMatrixPrev(i, j));
        }
    }
    DensityRMS = sqrt(DensityRMS);
    return DensityRMS;
}

double SCFIteration(Eigen::MatrixXd &FockMatrix, Eigen::MatrixXd &DensityMatrix, std::map<std::string, double> &Integrals, Eigen::MatrixXd &HCore, Eigen::MatrixXd &SOrtho)
{
    double Energy = 0;
    BuildFockMatrix(FockMatrix, DensityMatrix, Integrals);
    Eigen::MatrixXd FockOrtho = SOrtho.transpose() * FockMatrix * SOrtho;
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > EigensystemFockOrtho(FockOrtho);
    Eigen::MatrixXd CoeffMatrix = SOrtho * EigensystemFockOrtho.eigenvectors();
    DensityMatrix = CoeffMatrix * CoeffMatrix.transpose();

    for(int i = 0; i < FockMatrix.rows(); i++)
    {
        for(int j = 0; j < FockMatrix.cols(); j++)
        {
            Energy += 0.5 * (DensityMatrix(i, j) * (HCore(i, j) + FockMatrix(i, j)));
        }
    }
    return Energy;
}

int main(int argc, char* argv[])
{
    InputObj Input;
    if(argc == 4)
    {
        Input.SetNames(argv[1], argv[2], argv[3]);
    }
    else
    {
        Input.GetInputName();
    }
    Input.Set();

    unsigned int NumAO = Input.OverlapMatrix.rows();

    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > EigensystemS(Input.OverlapMatrix);
    Eigen::SparseMatrix< double > LambdaSOrtho; // Holds the inverse sqrt matrix of eigenvalues of S ( Lambda^-1/2 )
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(int i = 0; i < NumAO; i++)
    {
        tripletList.push_back(T(i, i, 1 / sqrt(EigensystemS.eigenvalues()[i])));
    }
    LambdaSOrtho.setFromTriplets(tripletList.begin(), tripletList.end());
    
    Eigen::MatrixXd SOrtho = EigensystemS.eigenvectors() * LambdaSOrtho * EigensystemS.eigenvectors().transpose();

    Eigen::MatrixXd FockMatrix(NumAO, NumAO);
    Eigen::MatrixXd DensityMatrix = Eigen::MatrixXd::Zero(NumAO, NumAO); // First guess of D set to zero.
    Eigen::MatrixXd HCore(NumAO, NumAO);
    BuildFockMatrix(HCore, DensityMatrix, Input.Integrals); // Form HCore (D is zero)

    /* We go through the SCF step once. */
    double Energy = SCFIteration(FockMatrix, DensityMatrix, Input.Integrals, HCore, SOrtho);
    Eigen::MatrixXd DensityMatrixPrev;
    double EnergyPrev;
    double DensityRMS = 1;

    while(DensityRMS > 10E-6)
    {
        EnergyPrev = Energy;
        DensityMatrixPrev = DensityMatrix;
        Energy = SCFIteration(FockMatrix, DensityMatrix, Input.Integrals, HCore, SOrtho);
        DensityRMS = CalcDensityRMS(DensityMatrix, DensityMatrixPrev);
    }
    return 0;
}
