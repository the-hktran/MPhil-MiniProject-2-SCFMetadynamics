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

double SCF(std::vector< std::tuple< Eigen::MatrixXd, double, double > > &Bias, int SolnNum, Eigen::MatrixXd &DensityMatrix, InputObj &Input, std::ofstream &Output, Eigen::MatrixXd SOrtho, Eigen::MatrixXd HCore);
void BuildFockMatrix(Eigen::MatrixXd &FockMatrix, Eigen::MatrixXd &DensityMatrix, std::map<std::string, double> &Integrals, std::vector< std::tuple< Eigen::MatrixXd, double, double > > &Bias, int NumElectrons);

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
    std::vector< std::tuple< Eigen::MatrixXd, double, double > > Bias; // Tuple containing DensityMatrix, N_x, lambda_x, d^2_0x

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
    Eigen::MatrixXd DensityMatrix = Eigen::MatrixXd::Zero(Input.NumAO, Input.NumAO); // First guess of D set to zero.
    Eigen::MatrixXd HCore(Input.NumAO, Input.NumAO);
    BuildFockMatrix(HCore, DensityMatrix, Input.Integrals, Bias, Input.NumElectrons); // Form HCore (D is zero)

    double CurrentEnergy;
    double PreviousEnergy;

    for(int i = 0; i < Input.NumSoln; i++)
    {
        PreviousEnergy = CurrentEnergy;
        std::tuple< Eigen::MatrixXd, double, double > tmpTuple;
        CurrentEnergy = SCF(Bias, i + 1, DensityMatrix, Input, Output, SOrtho, HCore);
        tmpTuple = std::make_tuple(DensityMatrix, 1, 1);
        Bias.push_back(tmpTuple);
    }

    return 0;
}