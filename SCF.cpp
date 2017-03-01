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
#include <queue>
#include <algorithm>

void BuildFockMatrix(Eigen::MatrixXd &FockMatrix, Eigen::MatrixXd &DensityMatrix, std::map<std::string, double> &Integrals, std::vector< std::tuple< Eigen::MatrixXd, double, double > > &Bias, int NumElectrons);
double Metric(int NumElectrons, Eigen::MatrixXd &FirstDensityMatrix, Eigen::MatrixXd &SecondDensityMatrix);
void ModifyBias(std::vector< std::tuple< Eigen::MatrixXd, double, double > > &Bias);
void NewDensityMatrix(Eigen::MatrixXd &DensityMatrix, Eigen::MatrixXd &CoeffMatrix, std::vector<int> OccupiedOrbitals, std::vector<int> VirtualOrbitals);

double MatrixDot(Eigen::MatrixXd FirstMatrix, Eigen::MatrixXd SecondMatrix)
{
    double Dot = 0;
    for(int i = 0; i < FirstMatrix.rows(); i++)
    {
        for(int j = 0; j < FirstMatrix.cols(); j++)
        {
            Dot += FirstMatrix(i, j) * SecondMatrix(i, j);
        }
    }
    return Dot;
}

void DIIS(Eigen::MatrixXd &FockMatrix, std::vector< Eigen::MatrixXd > &AllFockMatrices, std::vector< Eigen::MatrixXd > &AllErrorMatrices)
{
    Eigen::MatrixXd B(AllErrorMatrices.size() + 1, AllErrorMatrices.size() + 1); // This is the linear system we solve for the coefficients.
    for(int i = 0; i < AllErrorMatrices.size(); i++)
    {
        for(int j = i; j < AllErrorMatrices.size(); j++)
        {
            B(i, j) = MatrixDot(AllErrorMatrices[i], AllErrorMatrices[j]);
            B(j, i) = B(i, j);
        }
    }
    for(int i = 0; i < AllErrorMatrices.size(); i++)
    {
        B(i, AllErrorMatrices.size()) = -1;
        B(AllErrorMatrices.size(), i) = -1;
    }
    B(AllErrorMatrices.size(), AllErrorMatrices.size()) = 0;

    Eigen::VectorXd b = Eigen::VectorXd::Zero(AllErrorMatrices.size() + 1);
    b(AllErrorMatrices.size()) = -1;
    // We want to solve for c in Bc = b
    Eigen::VectorXd c = B.colPivHouseholderQr().solve(b);

    /* And now put it together to get the new FockMatrix. If we are on the first iteration, then c1 = 1 and we end up where we
       started. */
    FockMatrix = Eigen::MatrixXd::Zero(FockMatrix.rows(), FockMatrix.cols());
    for(int i = 0; i < AllFockMatrices.size(); i++)
    {
        FockMatrix += c[i] * AllFockMatrices[i];
    }
}

void MaximumOverlapMethod(Eigen::MatrixXd &DensityMatrix, Eigen::MatrixXd &CoeffMatrix, Eigen::MatrixXd &CoeffMatrixPrev, Eigen::MatrixXd &OverlapMatrix, int NumOcc, int NumAO, std::vector<int> &OccupiedOrbitals, std::vector<int> &VirtualOrbitals) // MoM 
{
    std::priority_queue< std::pair<double, int> > PQueue;
    // std::cout << "\n" << CoeffMatrix << std::endl;
    // if(CoeffMatrixPrev(0, 0) == 1) // Means that this is the first iteration. This step is important because the signs of the eigenvector affect the overlap. 
    // {
    //     CoeffMatrixPrev = CoeffMatrix;
    // }
    for(int j = 0; j < CoeffMatrix.cols(); j++)
    {
        double Projection_j = 0;
        for(int i = 0; i < CoeffMatrixPrev.cols(); i++)
        {
            Projection_j += CoeffMatrix.col(j).dot(CoeffMatrixPrev.col(i));
        }
        for(int nu = 0; nu < CoeffMatrix.rows(); nu++)
        {
            for(int mu = 0; mu < CoeffMatrix.rows(); mu++)
            {
                for(int i = 0; i < CoeffMatrix.cols(); i++)
                {
                    Projection_j += (CoeffMatrixPrev.col(i)[mu] * OverlapMatrix(mu, nu) * CoeffMatrix.col(j)[nu]); // C_mu,i * S_mu,nu * C_nu,j
                }
            }
        }
        PQueue.push(std::pair<double, int>(Projection_j, j));
    }

    for(int i = 0; i < NumOcc; i++)
    {
        OccupiedOrbitals[NumOcc - 1 - i] = PQueue.top().second; // Insert largest indices into occupied orbital list, smallest first.
        PQueue.pop(); // Remove largest element.
    }
    int VirtIndex = 0;
    for(int i = 0; i < NumAO; i++)
    {
        if(std::find(OccupiedOrbitals.begin(), OccupiedOrbitals.end(), i) == OccupiedOrbitals.end()) // Means not in occupied orbital list.
        {
            VirtualOrbitals[VirtIndex] = i;
            VirtIndex++;
        }
    }

    /* Density matrix: C(occ) * C(occ)^T */
	for (int i = 0; i < DensityMatrix.rows(); i++)
	{
		for (int j = 0; j < DensityMatrix.cols(); j++)
		{
			double DensityElement = 0;
			for (int k = 0; k < NumOcc; k++)
			{
				DensityElement += CoeffMatrix(i, OccupiedOrbitals[k]) * CoeffMatrix(j, OccupiedOrbitals[k]);
			}
			DensityMatrix(i, j) = DensityElement;
		}
    }

    std::cout << "\nCoeff\n" << CoeffMatrix << "\nDensity\n" << DensityMatrix << std::endl;
    for(int i = 0; i < OccupiedOrbitals.size(); i++)
    {
        std::cout << OccupiedOrbitals[i] << std::endl;
    }
}

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

/// <summary>
/// One SCF iteration. Takes a density matrix, generates the corresponding Fock matrix, and then computes a new density
/// matrix and the corresponding energy. The energy is returned.
/// </summary>
/// <param name="DensityMatrix">
/// Density matrix of the current iteration.
/// </param>
/// <param name="Integrals">
/// Map to value of two electron integrals.
/// </param>
/// <param name="HCore">
/// Fock matrix generate from a zero density matrix. Stored because it can be reused for the generation of a new fock
/// matrix and of the energy.
/// </param>
/// <param name="SOrtho">
/// S^(-1/2), the symmetric orthogonalization matrix
/// </param>
/// <param name="NumOcc">
/// Number of occupied orbitals. Used to calculate the density matrix.
/// </param>
double SCFIteration(Eigen::MatrixXd &DensityMatrix, InputObj &Input, Eigen::MatrixXd &HCore, Eigen::MatrixXd &SOrtho, std::vector< std::tuple< Eigen::MatrixXd, double, double > > &Bias, Eigen::MatrixXd &CoeffMatrix, std::vector< Eigen::MatrixXd > &AllFockMatrices, std::vector< Eigen::MatrixXd > &AllErrorMatrices, Eigen::MatrixXd &CoeffMatrixPrev, std::vector<int> &OccupiedOrbitals, std::vector<int> &VirtualOrbitals)
{
    Eigen::MatrixXd FockMatrix(DensityMatrix.rows(), DensityMatrix.cols());
    double Energy = 0;
    BuildFockMatrix(FockMatrix, DensityMatrix, Input.Integrals, Bias, Input.NumElectrons); // Calculates and stores fock matrix. Includes bias.
    AllFockMatrices.push_back(FockMatrix); // Store this iteration's Fock matrix.
    
    Eigen::MatrixXd ErrorMatrix = FockMatrix * DensityMatrix * Input.OverlapMatrix - Input.OverlapMatrix * DensityMatrix * FockMatrix; // Current error matrix
    AllErrorMatrices.push_back(ErrorMatrix);
    DIIS(FockMatrix, AllFockMatrices, AllErrorMatrices); // Generates F' using DIIS and stores it in FockMatrix.

    Eigen::MatrixXd FockOrtho = SOrtho.transpose() * FockMatrix * SOrtho;
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > EigensystemFockOrtho(FockOrtho); // Eigenvectors ordered from lowest to highest eigenvalues
    CoeffMatrix = SOrtho * EigensystemFockOrtho.eigenvectors();

	/* Density matrix: C(occ) * C(occ)^T */
    if(Bias.empty())
    {
        for (int i = 0; i < DensityMatrix.rows(); i++)
        {
            for (int j = 0; j < DensityMatrix.cols(); j++)
            {
                double DensityElement = 0;
                for (int k = 0; k < Input.NumOcc; k++)
                {
                    DensityElement += CoeffMatrix(i, OccupiedOrbitals[k]) * CoeffMatrix(j, OccupiedOrbitals[k]);
                }
                DensityMatrix(i, j) = DensityElement;
            }
        }
        // for(int i = 0; i < Input.NumAO; i++)
        // {
        //     if(i < Input.NumOcc)
        //     {
        //         OccupiedOrbitals[i] = i;
        //     }
        //     else
        //     {
        //         VirtualOrbitals[Input.NumOcc - i] = i;
        //     }
        // }
        std::cout << "\nDensityMatrix\n" << DensityMatrix << std::endl;
    }
    else /* Use MoM to generate new density matrix, as well as determine the occupied and virtual orbitals */
    {
        MaximumOverlapMethod(DensityMatrix, CoeffMatrix, CoeffMatrixPrev, Input.OverlapMatrix, Input.NumOcc, Input.NumAO, OccupiedOrbitals, VirtualOrbitals); // MoM 
        CoeffMatrixPrev = CoeffMatrix; // Now that we finish the MoM iteration, set CoeffMatrixPrev.
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

double SCF(std::vector< std::tuple< Eigen::MatrixXd, double, double > > &Bias, int SolnNum, Eigen::MatrixXd &DensityMatrix, InputObj &Input, std::ofstream &Output, Eigen::MatrixXd &SOrtho, Eigen::MatrixXd &HCore, std::vector< double > &AllEnergies, Eigen::MatrixXd &CoeffMatrix, std::vector<int> &OccupiedOrbitals, std::vector<int> &VirtualOrbitals)
{
	double SCFTol = 10E-12;

    // DensityMatrix = Eigen::MatrixXd::Zero(Input.NumAO, Input.NumAO);

	Output << "Beginning search for Solution " << SolnNum << std::endl;
	Output << "Iteration\tEnergy" << std::endl;
	std::cout << "SCF MetaD: Beginning search for Solution " << SolnNum << std::endl;
	clock_t ClockStart = clock();

    double Energy = 1;
    Eigen::MatrixXd DensityMatrixPrev;
    double EnergyPrev = 1;
    double DensityRMS = 1;
    unsigned short int Count = 1; // Counts number of iterations.
    bool isUniqueSoln = false; // Will tell us if the solution is unique by checking against all previous energies.
    // CoeffMatrix = Eigen::MatrixXd::Zero(Input.NumAO, Input.NumAO); // Stores coefficient, we need this to generate the next density matrix.

    while(!isUniqueSoln)
    {
        std::vector< Eigen::MatrixXd > AllFockMatrices; // Holds previous fock matrices for DIIS procedure.
        std::vector< Eigen::MatrixXd > AllErrorMatrices; // Error matrices for DIIS
        Eigen::MatrixXd CoeffMatrixPrev = Eigen::MatrixXd::Identity(Input.NumAO, Input.NumAO); // For MoM
        DensityRMS = 1; // Reset loop.
        while(fabs(DensityRMS) > SCFTol || fabs(Energy - EnergyPrev) > SCFTol)
        {
            std::cout << "SCF MetaD: Iteration " << Count << "...";
            EnergyPrev = Energy;
            DensityMatrixPrev = DensityMatrix;
            Energy = SCFIteration(DensityMatrix, Input, HCore, SOrtho, Bias, CoeffMatrix, AllFockMatrices, AllErrorMatrices, CoeffMatrixPrev, OccupiedOrbitals, VirtualOrbitals);
            DensityRMS = CalcDensityRMS(DensityMatrix, DensityMatrixPrev);
            std::cout << " complete with a biased energy of " << Energy + Input.Integrals["0 0 0 0"] << std::endl;
            // Output << Count << "\t" << Energy + Input.Integrals["0 0 0 0"] << std::endl;
            Count++;
            if(Count > 500)
            {
            //     // NewDensityMatrix(DensityMatrix, CoeffMatrix, Input.NumOcc, Input.NumAO);
            //     // Count = 1;
            //     DensityMatrix = Eigen::MatrixXd::Random(Input.NumAO, Input.NumAO);
            //     Count = 1;
                AllFockMatrices.clear();
                AllErrorMatrices.clear();
                // NewDensityMatrix(DensityMatrix, CoeffMatrix, OccupiedOrbitals, VirtualOrbitals);
                Count = 1;
            }
            std::cout << "\nDensityRMS\n" << DensityRMS << std::endl;
            std::string pause;
            std::getline(std::cin, pause);
        } // Means we have converged with the bias. Now we remove the bias and converge to the minimum.
        Count = 1;
        std::vector< std::tuple< Eigen::MatrixXd, double, double > > EmptyBias;
        // CoeffMatrixPrev = CoeffMatrix;
        AllFockMatrices.clear();
        AllErrorMatrices.clear();
        DensityRMS = 1;
        while(fabs(DensityRMS) > SCFTol || fabs(Energy - EnergyPrev) > SCFTol)
        {
            std::cout << "SCF MetaD: Iteration " << Count << "...";
            EnergyPrev = Energy;
            DensityMatrixPrev = DensityMatrix;
            Energy = SCFIteration(DensityMatrix, Input, HCore, SOrtho, EmptyBias, CoeffMatrix, AllFockMatrices, AllErrorMatrices, CoeffMatrixPrev, OccupiedOrbitals, VirtualOrbitals);
            DensityRMS = CalcDensityRMS(DensityMatrix, DensityMatrixPrev);
            std::cout << " complete with an energy of " << Energy + Input.Integrals["0 0 0 0"] << std::endl;
            // Output << Count << "\t" << Energy + Input.Integrals["0 0 0 0"] << std::endl;
            Count++;
            if(Count > 500)
            {
            //     // NewDensityMatrix(DensityMatrix, CoeffMatrix, Input.NumOcc, Input.NumAO);
            //     // Count = 1;
            //     DensityMatrix = Eigen::MatrixXd::Random(Input.NumAO, Input.NumAO);
            //     Count = 1;
                AllFockMatrices.clear();
                AllErrorMatrices.clear();
                // NewDensityMatrix(DensityMatrix, CoeffMatrix, Input.NumOcc, Input.NumAO);
                Count = 1;
            }
            std::cout << "\nDensityRMS\n" << DensityRMS << std::endl;
            std::string pause;
            std::getline(std::cin, pause);
        }
        Count = 1;
        isUniqueSoln = true;
        for(int i = 0; i < AllEnergies.size(); i++)
        {
            if(fabs(Energy + Input.Integrals["0 0 0 0"] - AllEnergies[i]) < 10E-3)
            {
                isUniqueSoln = false;
            }
        }
        if(!isUniqueSoln)
        {
            std::cout << "SCF MetaD: Solution is not unique. Retrying solution " << SolnNum << "." << std::endl;
            ModifyBias(Bias);
            // NewDensityMatrix(DensityMatrix, CoeffMatrix, OccupiedOrbitals, VirtualOrbitals);
            // DensityMatrix = Eigen::MatrixXd::Random(DensityMatrix.rows(), DensityMatrix.cols());
        }
    }

    AllEnergies.push_back(Energy + Input.Integrals["0 0 0 0"]);

    // double EnergyBias = 0;
    // for(int i = 0; i < Bias.size(); i++)
    // {
    //     EnergyBias += std::get<1>(Bias[i]) * exp(-1 * std::get<2>(Bias[i]) * Metric(Input.NumElectrons, DensityMatrix, std::get<0>(Bias[i])));
    // }

	std::cout << "SCF MetaD: Solution " << SolnNum << " has converged with energy " << Energy + Input.Integrals["0 0 0 0"] << std::endl;
	std::cout << "SCF MetaD: This solution took " << (clock() - ClockStart) / CLOCKS_PER_SEC << " seconds." << std::endl;
	Output << "Solution " << SolnNum << " has converged with energy " << Energy + Input.Integrals["0 0 0 0"] << std::endl;
	Output << "This solution took " << (clock() - ClockStart) / CLOCKS_PER_SEC << " seconds." << std::endl;

    return Energy;
}
