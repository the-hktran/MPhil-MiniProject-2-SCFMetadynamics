#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <tuple>
#include "ReadInput.h"
#include <fstream>
#include <map>
#include <stdlib.h> 

double SCF(std::vector< std::tuple< Eigen::MatrixXd, double, double > > &Bias, int SolnNum, Eigen::MatrixXd &DensityMatrix, InputObj &Input, std::ofstream &Output, Eigen::MatrixXd &SOrtho, Eigen::MatrixXd &HCore, std::vector< double > &AllEnergies, Eigen::MatrixXd &CoeffMatrix, std::vector<int> &OccupiedOrbitals, std::vector<int> &VirtualOrbitals);
void BuildFockMatrix(Eigen::MatrixXd &FockMatrix, Eigen::MatrixXd &DensityMatrix, std::map<std::string, double> &Integrals, std::vector< std::tuple< Eigen::MatrixXd, double, double > > &Bias, int NumElectrons);

/* This function makes a new density matrix to be used for the next metadynamics iteration. It switches an occupied orbitals with
   a random virtual orbital */
void NewDensityMatrix(Eigen::MatrixXd &DensityMatrix, Eigen::MatrixXd &CoeffMatrix, std::vector<int> OccupiedOrbitals, std::vector<int> VirtualOrbitals)
{
    int ExcludedOcc = rand() % OccupiedOrbitals.size();
    int IncludedVirt = rand() % VirtualOrbitals.size();

    double Cos = rand() / RAND_MAX; // A random value between 0 and 1;
    double Sin = sqrt(1 - Cos * Cos); // Corresponding sin value.

    Eigen::MatrixXd RotatedCoeff = CoeffMatrix;
    for (int i = 0; i < RotatedCoeff.rows(); i++)
    {
        RotatedCoeff(i, OccupiedOrbitals[ExcludedOcc]) = Cos * CoeffMatrix(i, OccupiedOrbitals[ExcludedOcc]) - Sin * CoeffMatrix(i, VirtualOrbitals[IncludedVirt]);
    }
    for (int i = 0; i < DensityMatrix.rows(); i++)
	{
		for (int j = 0; j < DensityMatrix.cols(); j++)
		{
			double DensityElement = 0;
			for (int k = 0; k < OccupiedOrbitals.size(); k++)
			{
				DensityElement += RotatedCoeff(i, OccupiedOrbitals[k]) * RotatedCoeff(j, OccupiedOrbitals[k]);
			}
			DensityMatrix(i, j) = DensityElement;
		}
	}
    // for (int i = 0; i < DensityMatrix.rows(); i++)
	// {
	// 	for (int j = 0; j < DensityMatrix.cols(); j++)
	// 	{
	// 		double DensityElement = 0;
	// 		for (int k = 0; k < OccupiedOrbitals.size(); k++)
	// 		{
    //             if(k == ExcludedOcc)
    //             {
    //                 DensityElement += CoeffMatrix(i, VirtualOrbitals[IncludedVirt]) * CoeffMatrix(j, VirtualOrbitals[IncludedVirt]);
    //             }
    //             else
    //             {
	// 			    DensityElement += CoeffMatrix(i, OccupiedOrbitals[k]) * CoeffMatrix(j, OccupiedOrbitals[k]);
    //             }
	// 		}
	// 		DensityMatrix(i, j) = DensityElement;
	// 	}
	// }

}

/* Increases N_x and decreases lambda_x when SCF converges to solution x */
void ModifyBias(std::vector< std::tuple< Eigen::MatrixXd, double, double > > &Bias, short int WhichSoln)
{
    if(WhichSoln == -1) // Means the solution was positive, not reconverged. Don't do anything.
    {
        return;
    }
    double BiasScale = 1.1; // Scale to increase and decrease parameters. Hard coded for now.
    double NewNorm = std::get<1>(Bias[WhichSoln]) * BiasScale; // Increase height of Gaussian.
    double NewLambda = std::get<2>(Bias[WhichSoln]) / BiasScale; // Increase width of Gaussian (lambda is the inverse variance).
    if(NewNorm > 100)
    {
        NewNorm = 100 * (rand() / RAND_MAX);
    }
    if(NewLambda < 1E-10)
    {
        NewLambda = 2 * (rand() / RAND_MAX) + 1E-10;
    }
    std::tuple< Eigen::MatrixXd, double, double > NewTuple = std::make_tuple(std::get<0>(Bias[WhichSoln]), NewNorm, NewLambda);
    Bias[WhichSoln] = NewTuple;
}

double Metric(int NumElectrons, Eigen::MatrixXd &FirstDensityMatrix, Eigen::MatrixXd &SecondDensityMatrix)
{
    double d = 0;
    for(int i = 0; i < FirstDensityMatrix.rows(); i++)
    {
        for(int j = 0; j < FirstDensityMatrix.cols(); j++)
        {
            d += FirstDensityMatrix(i, j) * SecondDensityMatrix(j, i);
        }
    }
    d = (double)NumElectrons * (1 - d);
    return d;
}

double BiasMatrixElement(int Row, int Col, std::vector< std::tuple< Eigen::MatrixXd, double, double > > &Bias, Eigen::MatrixXd &CurrentDensity, int NumElectrons)
{
    double BiasElement = 0;
    for(int i = 0; i < Bias.size(); i++)
    {
        BiasElement += std::get<0>(Bias[i])(Row, Col) * std::get<1>(Bias[i]) * std::get<2>(Bias[i]) * exp(-1 * std::get<2>(Bias[i]) * Metric(NumElectrons, CurrentDensity, std::get<0>(Bias[i])));
    }
    return BiasElement;
}

int main(int argc, char* argv[])
{
    std::vector< std::tuple< Eigen::MatrixXd, double, double > > Bias; // Tuple containing DensityMatrix, N_x, lambda_x

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
	Output << "\n" << Input.NumSoln << " solutions desired." << std::endl;

    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > EigensystemS(Input.OverlapMatrix);
    Eigen::SparseMatrix< double > LambdaSOrtho(Input.NumAO, Input.NumAO); // Holds the inverse sqrt matrix of eigenvalues of S ( Lambda^-1/2 )
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for(int i = 0; i < Input.NumAO; i++)
    {
        tripletList.push_back(T(i, i, 1 / sqrt(EigensystemS.eigenvalues()[i])));
    }
    LambdaSOrtho.setFromTriplets(tripletList.begin(), tripletList.end());
    
    Eigen::MatrixXd SOrtho = EigensystemS.eigenvectors() * LambdaSOrtho * EigensystemS.eigenvectors().transpose();

    /* Initialize the density matrix. We're going to be smart about it and use the correct ground state density
       corresponding to Q-Chem outputs. Q-Chem uses an MO basis for its output, so the density matrix has ones
       along the diagonal for occupied orbitals. */
    Eigen::MatrixXd DensityMatrix = Eigen::MatrixXd::Zero(Input.NumAO, Input.NumAO); //
    for(int i = 0; i < Input.NumOcc; i++)
    {
        DensityMatrix(i, i) = 1;
    }
    
    Eigen::MatrixXd HCore(Input.NumAO, Input.NumAO);
    Eigen::MatrixXd ZeroMatrix = Eigen::MatrixXd::Zero(Input.NumAO, Input.NumAO);
    BuildFockMatrix(HCore, ZeroMatrix, Input.Integrals, Bias, Input.NumElectrons); // Form HCore (D is zero)

    double Energy;

    std::vector< double > AllEnergies;
    Eigen::MatrixXd CoeffMatrix = Eigen::MatrixXd::Zero(Input.NumAO, Input.NumAO);
    std::vector<int> OccupiedOrbitals(Input.NumOcc);
    std::vector<int> VirtualOrbitals(Input.NumAO - Input.NumOcc);
    for(int i = 0; i < Input.NumAO; i++)
    {
        if(i < Input.NumOcc)
        {
            OccupiedOrbitals[i] = i;
        }
        else
        {
            VirtualOrbitals[i - Input.NumOcc] = i;
        }
    }

    for(int i = 0; i < Input.NumSoln; i++)
    {
        std::tuple< Eigen::MatrixXd, double, double > tmpTuple;
        NewDensityMatrix(DensityMatrix, CoeffMatrix, OccupiedOrbitals, VirtualOrbitals); // CoeffMatrix is zero so this doesn't do anything the  first time.
        Energy = SCF(Bias, i + 1, DensityMatrix, Input, Output, SOrtho, HCore, AllEnergies, CoeffMatrix, OccupiedOrbitals, VirtualOrbitals);
        tmpTuple = std::make_tuple(DensityMatrix, 0.1, 1);
        Bias.push_back(tmpTuple);
    }
    // for(int i = 0; i < Bias.size(); i++)
    // {
    //     std::cout << std::get<0>(Bias[i]) << "\n" << std::get<1>(Bias[i]) << "\t" << std::get<2>(Bias[i]) << std::endl;
    // }

    return 0;
}