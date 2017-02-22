#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <string>
#include <map>
#include "ReadInput.h"
#include <Eigen/Sparse>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <ctime>

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

double SCFIteration(Eigen::MatrixXd &DensityMatrix, std::map<std::string, double> &Integrals, Eigen::MatrixXd &HCore, Eigen::MatrixXd &SOrtho, unsigned short int NumOcc)
{
    Eigen::MatrixXd FockMatrix(DensityMatrix.rows(), DensityMatrix.cols());
    double Energy = 0;
    BuildFockMatrix(FockMatrix, DensityMatrix, Integrals);
    Eigen::MatrixXd FockOrtho = SOrtho.transpose() * FockMatrix * SOrtho;
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > EigensystemFockOrtho(FockOrtho);
    Eigen::MatrixXd CoeffMatrix = SOrtho * EigensystemFockOrtho.eigenvectors();

	/* Density matrix: C(occ) * C(occ)^T */
	for (int i = 0; i < DensityMatrix.rows(); i++)
	{
		for (int j = 0; j < DensityMatrix.cols(); j++)
		{
			double DensityElement = 0;
			for (int k = 0; k < NumOcc; k++)
			{
				DensityElement += CoeffMatrix(i, k) * CoeffMatrix(j, k);
			}
			DensityMatrix(i, j) = DensityElement;
		}
	}

	/* Now calculate the energy. */
    for(int i = 0; i < FockMatrix.rows(); i++)
    {
        for(int j = 0; j < FockMatrix.cols(); j++)
        {
            Energy += (DensityMatrix(i, j) * (HCore(i, j) + FockMatrix(i, j)));
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

	std::ofstream Output(Input.OutputName);

	Output << "Self-Consistent Field Metadynamics Calculation" << std::endl;

    unsigned int NumAO = Input.OverlapMatrix.rows();

    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > EigensystemS(Input.OverlapMatrix);
    Eigen::SparseMatrix< double > LambdaSOrtho(NumAO, NumAO); // Holds the inverse sqrt matrix of eigenvalues of S ( Lambda^-1/2 )
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(int i = 0; i < NumAO; i++)
    {
        tripletList.push_back(T(i, i, 1 / sqrt(EigensystemS.eigenvalues()[i])));
    }
    LambdaSOrtho.setFromTriplets(tripletList.begin(), tripletList.end());
    
    Eigen::MatrixXd SOrtho = EigensystemS.eigenvectors() * LambdaSOrtho * EigensystemS.eigenvectors().transpose();
    Eigen::MatrixXd DensityMatrix = Eigen::MatrixXd::Zero(NumAO, NumAO); // First guess of D set to zero.
    Eigen::MatrixXd HCore(NumAO, NumAO);
    BuildFockMatrix(HCore, DensityMatrix, Input.Integrals); // Form HCore (D is zero)

	Output << "Beginning search for Solution " << 1 << std::endl;
	Output << "Iteration\tEnergy" << std::endl;
	std::cout << "SCF MetaD: Beginning search for Solution " << 1 << std::endl;
	clock_t ClockStart = clock();

    /* We go through the SCF step once. */
	std::cout << "SCF MetaD: Iteration 1...";
    double Energy = SCFIteration(DensityMatrix, Input.Integrals, HCore, SOrtho, Input.NumOcc);
    Eigen::MatrixXd DensityMatrixPrev;
    double EnergyPrev = 1;
    double DensityRMS = 1;
	std::cout << " complete with an energy of " << Energy + Input.Integrals["0 0 0 0"] << std::endl;
	Output << 1 << "\t" << Energy + Input.Integrals["0 0 0 0"] << std::endl;

    unsigned short int Count = 2;

    while(fabs(DensityRMS) > 10E-12 || fabs(Energy - EnergyPrev) > 10E-12)
    {
        std::cout << "SCF MetaD: Iteration " << Count << "...";
        EnergyPrev = Energy;
        DensityMatrixPrev = DensityMatrix;
        Energy = SCFIteration(DensityMatrix, Input.Integrals, HCore, SOrtho, Input.NumOcc);
        DensityRMS = CalcDensityRMS(DensityMatrix, DensityMatrixPrev);
        std::cout << " complete with an energy of " << Energy + Input.Integrals["0 0 0 0"] << std::endl;
		Output << Count << "\t" << Energy + Input.Integrals["0 0 0 0"] << std::endl;
        Count++;
    }
	std::cout << "SCF MetaD: Solution " << 1 << " has converged with energy " << Energy + Input.Integrals["0 0 0 0"] << std::endl;
	std::cout << "SCF MetaD: This solution took " << (clock() - ClockStart) / CLOCKS_PER_SEC << " seconds." << std::endl;
	Output << "Solution " << 1 << " has converged with energy " << Energy + Input.Integrals["0 0 0 0"] << std::endl;
	Output << "This solution took " << (clock() - ClockStart) / CLOCKS_PER_SEC << " seconds." << std::endl;

    return 0;
}
