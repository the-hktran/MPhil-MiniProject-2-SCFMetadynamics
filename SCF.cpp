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

void GenerateRandomDensity(Eigen::MatrixXd &DensityMatrix)
{
    for(int i = 0; i < DensityMatrix.rows(); i++)
    {
        for(int j = 0; j < DensityMatrix.cols(); j++)
        {
            DensityMatrix(i, j) = rand() / RAND_MAX;
        }
    }
    DensityMatrix / DensityMatrix.trace();
}
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
    /* This holds the set of total summed overlaps described below. It is a queue because this is an efficient way to choose largest elements in a set, at
       the cost of more complex insertion. The first in the pair is the value of the summed overlap and the second index is the orbital number, which ranges
       over all orbitals in the new coefficient matrix, denoted j below. */
    std::priority_queue< std::pair<double, int> > PQueue;

    /* Now we calculate the projections for the j'th new orbital in the new coefficient matrix. This is the j'th column in CoeffMatrix. 
       What we do is calculate the total summed overlap of EVERY new molecular orbital with every OCCUPIED molecular orbital in the
       old coefficient matrix, and then we take the NumOcc highest new orbitals to be the new overlap orbitals. */
    for(int j = 0; j < CoeffMatrix.cols(); j++) // Loop over EVERY orbital in the new coefficient matrix.
    {
        double NewOrbitalTotalOverlap = 0; // Holds the sum of < jthNewOrbital | OldOccupiedOrbital > over all OldOccupiedOrbital
        for(int nu = 0; nu < CoeffMatrix.rows(); nu++) // nu and mu index the basis functions, so they run down the coefficient matrix and label rows. This is needed since basis functions may not be orthogonal.
        {
            for(int mu = 0; mu < CoeffMatrix.rows(); mu++) // See previous comment.
            {
                for(int i = 0; i < OccupiedOrbitals.size(); i++) // i runs over all occupied orbitals in the old coefficient matrix.
                {
                    // Note that I am calling the column index first, so it may look like the matrix is transposed even though it's standard.
                    // This is C_mu,i(old) * S_mu,nu * C_nu,j(new) where i is an occupied orbital. For each j, we sum over mu and nu. If S is the identity
                    // matrix, then for each j we are calculating sum_mu C_mu,i(old) * C_nu,j(new), or C_i(old) dot C_j(new)
                    NewOrbitalTotalOverlap += (CoeffMatrixPrev.col(OccupiedOrbitals[i])[mu] * OverlapMatrix(mu, nu) * CoeffMatrix.col(j)[nu]); 
                }
            }
        }
        PQueue.push(std::pair<double, int>(NewOrbitalTotalOverlap, j)); 
        // Adds the value of the summed overlap to the queue.
        // The first value is the value of the overlap, which is used to order the list. The second value is the orbital number,
        // which is what we actually care about, and it will be the value we pull out of this queue.
    } // end loop over new orbitals j.
   
    /* Now we take the NumOcc highest summed overlap values and get their corresponding orbitals, adding them to the list of occupied
       orbitals. This is then used to calculate the density matrix and used in the MOM procedure for the next step. In order to this,
       we search for the highest element in PQueue using top(). Then we remove this element using pop(). Then top() gives the next 
       highest element the next time we call it. */
    for(int i = 0; i < NumOcc; i++)
    {
        OccupiedOrbitals[NumOcc - 1 - i] = PQueue.top().second; // Insert largest indices into occupied orbital list, smallest first.
        PQueue.pop(); // Remove largest element. The next time we call top(), we will get the next highest element.
    }
    int VirtIndex = 0;
    for(int i = 0; i < NumAO; i++) // Determine which orbitals are not in the occupied orbitals list and add them to the virtual orbitals list.
    {
        if(std::find(OccupiedOrbitals.begin(), OccupiedOrbitals.end(), i) == OccupiedOrbitals.end()) // Means not in occupied orbital list.
        {
            VirtualOrbitals[VirtIndex] = i;
            VirtIndex++;
        }
    }

    /* This calculates the density matrix (P) using P = C(occ) * C(occ)^T */
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
    // std::cout << "\nCoeff\n" << CoeffMatrix << "\nDensity\n" << DensityMatrix << std::endl;
    // for(int i = 0; i < OccupiedOrbitals.size(); i++)
    // {
    //     std::cout << OccupiedOrbitals[i] << std::endl;
    // }
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
    // DensityRMS = sqrt(DensityRMS);
    return DensityRMS;
}

/* Calculates the sum of the squared elements of Matrix. Note that this is not actually RMS, I just use it as an error. */
double CalcMatrixRMS(Eigen::MatrixXd &Matrix)
{
    double MatrixRMS = 0;
    for(int i = 0; i < Matrix.rows(); i++)
    {
        for(int j = 0; j < Matrix.cols(); j++)
        {
            MatrixRMS += Matrix(i, j) * Matrix(i, j);
        }
    }
    return MatrixRMS;
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
    Eigen::MatrixXd FockMatrix(DensityMatrix.rows(), DensityMatrix.cols()); // This will hold the FockMatrix.
    BuildFockMatrix(FockMatrix, DensityMatrix, Input.Integrals, Bias, Input.NumElectrons); // Calculates and stores fock matrix. Includes bias.
    AllFockMatrices.push_back(FockMatrix); // Store this iteration's Fock matrix for the DIIS procedure.
    
    Eigen::MatrixXd ErrorMatrix = FockMatrix * DensityMatrix * Input.OverlapMatrix - Input.OverlapMatrix * DensityMatrix * FockMatrix; // DIIS error matrix of the current iteration: FPS - SPF
    AllErrorMatrices.push_back(ErrorMatrix); // Save error matrix for DIIS.
    DIIS(FockMatrix, AllFockMatrices, AllErrorMatrices); // Generates F' using DIIS and stores it in FockMatrix.

    Eigen::MatrixXd FockOrtho = SOrtho.transpose() * FockMatrix * SOrtho; // Fock matrix in orthonormal basis.
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > EigensystemFockOrtho(FockOrtho); // Eigenvectors and eigenvalues ordered from lowest to highest eigenvalues
    CoeffMatrix = SOrtho * EigensystemFockOrtho.eigenvectors(); // Multiply the matrix of coefficients by S^-1/2 to get coefficients for nonorthonormal basis.

	/* Density matrix: C(occ) * C(occ)^T */
    if(Bias.empty())
    {
        // if(!Bias.empty())
        // {
        //     for(int i = 0; i < Input.NumAO; i++)
        //     {
        //         if(i < Input.NumOcc)
        //         {
        //             OccupiedOrbitals[i] = i;
        //         }
        //         else
        //         {
        //             VirtualOrbitals[Input.NumOcc - i] = i;
        //         }
        //     }
        // }
        
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
        std::cout << "\nDensityMatrix\n" << DensityMatrix << std::endl;
    }
    else /* Use MoM to generate new density matrix, as well as determine the occupied and virtual orbitals */
    {
        MaximumOverlapMethod(DensityMatrix, CoeffMatrix, CoeffMatrixPrev, Input.OverlapMatrix, Input.NumOcc, Input.NumAO, OccupiedOrbitals, VirtualOrbitals); // MoM 
        CoeffMatrixPrev = CoeffMatrix; // Now that we finish the MoM iteration, set CoeffMatrixPrev.
    }

	/* Now calculate the HF energy. E = sum_ij P_ij * (HCore_ij + F_ij) */
    double Energy = 0; // HF energy.
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
	double SCFTol = 10E-6; // SCF will terminate when the DIIS error is below this amount. 

    // DensityMatrix = Eigen::MatrixXd::Zero(Input.NumAO, Input.NumAO);

	Output << "Beginning search for Solution " << SolnNum << std::endl;
	Output << "Iteration\tEnergy" << std::endl;
	std::cout << "SCF MetaD: Beginning search for Solution " << SolnNum << std::endl;
	clock_t ClockStart = clock();

    double Energy = 1; // HF energy of SCF iteration.
    double DIISError = 1; // Square sum of DIIS error matrix of the current iteration. Used to test convergence.
    // Eigen::MatrixXd DensityMatrixPrev; // Stores density matrix of the previous iteration to test density matrix convergence.
    double EnergyPrev = 1; // Stores energy of previous iteration to test energy convergence.
    // double DensityRMS = 1; // Stores squared different between two sequential density matrices.
    unsigned short int Count = 1; // Counts number of iterations.
    bool isUniqueSoln = false; // Will tell us if the solution is unique by checking against all previous energies.
    // CoeffMatrix = Eigen::MatrixXd::Zero(Input.NumAO, Input.NumAO); // Stores coefficient, we need this to generate the next density matrix.

    while(!isUniqueSoln)
    {
        std::vector< Eigen::MatrixXd > AllFockMatrices; // Holds previous fock matrices for DIIS procedure.
        std::vector< Eigen::MatrixXd > AllErrorMatrices; // Error matrices for DIIS
        Eigen::MatrixXd CoeffMatrixPrev = Eigen::MatrixXd::Identity(Input.NumAO, Input.NumAO); // Two sequential coefficient matrices are stored for MOM.
        // DensityRMS = 1; // Reset loop.
        Count = 1;
        DIISError = 1; /// First guess of D set to zero. Reset loop.
        while(fabs(DIISError) > SCFTol * SCFTol && !Bias.empty()) // (fabs(DensityRMS) > SCFTol * SCFTol || fabs(Energy - EnergyPrev) > SCFTol * (fabs(Energy) + 1) || Count < 100)
        {
            std::cout << "SCF MetaD: Iteration " << Count << "...";
            EnergyPrev = Energy;
            // DensityMatrixPrev = DensityMatrix;
            Energy = SCFIteration(DensityMatrix, Input, HCore, SOrtho, Bias, CoeffMatrix, AllFockMatrices, AllErrorMatrices, CoeffMatrixPrev, OccupiedOrbitals, VirtualOrbitals);
            // DensityRMS = CalcDensityRMS(DensityMatrix, DensityMatrixPrev);
            DIISError = CalcMatrixRMS(AllErrorMatrices[AllErrorMatrices.size() - 1]);
            std::cout << " complete with a biased energy of " << Energy + Input.Integrals["0 0 0 0"] << std::endl;
            // Output << Count << "\t" << Energy + Input.Integrals["0 0 0 0"] << std::endl;
            Count++;

            /* This is a work-around that I put in. The first guess of the density is a zero matrix and this is not good. Unfortunately, DIIS
               rarely corrects this so I find that it helps to clear the Fock and Error matrices after a few iterations and we have a more reasonable
               guess of the coefficient, and thus density, matrices. Then DIIS converges to a reasonable solution. */
            if(Count == 5)
            {
                AllFockMatrices.clear();
                AllErrorMatrices.clear();
            }

            if(Count % 500 == 0) // Shouldn't take this long.
            {
            //     // NewDensityMatrix(DensityMatrix, CoeffMatrix, Input.NumOcc, Input.NumAO);
            //     // Count = 1;
            //     DensityMatrix = Eigen::MatrixXd::Random(Input.NumAO, Input.NumAO);
            //     Count = 1;
                AllFockMatrices.clear();
                AllErrorMatrices.clear();
                NewDensityMatrix(DensityMatrix, CoeffMatrix, OccupiedOrbitals, VirtualOrbitals);
                // GenerateRandomDensity(DensityMatrix);
                // Count = 1;
            }
            std::cout << "\nDIIS Error\n" << DIISError << std::endl;
            // std::string pause;
            // std::getline(std::cin, pause);
        } // Means we have converged with the bias. Now we remove the bias and converge to the minimum
        Count = 1;
        std::vector< std::tuple< Eigen::MatrixXd, double, double > > EmptyBias; // Same type as Bias, but it's empty so it's the same as having no bias.
        // CoeffMatrixPrev = CoeffMatrix;
        AllFockMatrices.clear();
        AllErrorMatrices.clear();
        DIISError = 1; // Reset for the next loop to start.
        while(fabs(DIISError) > SCFTol * SCFTol || Count < 10) // (fabs(DensityRMS) > SCFTol * SCFTol || fabs(Energy - EnergyPrev) > SCFTol * (fabs(Energy) + 1) || Count < 100)
        {
            std::cout << "SCF MetaD: Iteration " << Count << "...";
            EnergyPrev = Energy;
            // DensityMatrixPrev = DensityMatrix;
            Energy = SCFIteration(DensityMatrix, Input, HCore, SOrtho, EmptyBias, CoeffMatrix, AllFockMatrices, AllErrorMatrices, CoeffMatrixPrev, OccupiedOrbitals, VirtualOrbitals);
            // DensityRMS = CalcDensityRMS(DensityMatrix, DensityMatrixPrev);
            DIISError = CalcMatrixRMS(AllErrorMatrices[AllErrorMatrices.size() - 1]);
            std::cout << " complete with an energy of " << Energy + Input.Integrals["0 0 0 0"] << std::endl;
            // Output << Count << "\t" << Energy + Input.Integrals["0 0 0 0"] << std::endl;
            Count++;

            // if(Count == 5)
            // {
            //     AllFockMatrices.clear();
            //     AllErrorMatrices.clear();
            // }

            if(Count % 500 == 0)
            {
            //     // NewDensityMatrix(DensityMatrix, CoeffMatrix, Input.NumOcc, Input.NumAO);
            //     // Count = 1;
            //     DensityMatrix = Eigen::MatrixXd::Random(Input.NumAO, Input.NumAO);
            //     Count = 1;
                AllFockMatrices.clear();
                AllErrorMatrices.clear();
                // NewDensityMatrix(DensityMatrix, CoeffMatrix, OccupiedOrbitals, VirtualOrbitals);
                // Count = 1;
            }
            std::cout << "\nDIIS Error\n" << DIISError << std::endl;
            // std::string pause;
            // std::getline(std::cin, pause);
        }
        isUniqueSoln = true;
        
        if(Energy + Input.Integrals["0 0 0 0"] > 0)
        {
            isUniqueSoln = false;
        }
        else
        {
            for(int i = 0; i < AllEnergies.size(); i++)
            {
                if(fabs(Energy + Input.Integrals["0 0 0 0"] - AllEnergies[i]) < 10E-3) // Checks to see if new energy is equal to any previous energy.
                {
                    isUniqueSoln = false; // If it matches at least one, set this flag to false so the SCF procedure can repeat for this solution.
                }
            }
        }

        if(!isUniqueSoln) // If the flag is still false, we modify the bias and hope that this gives a better result.
        {
            std::cout << "SCF MetaD: Solution is not unique. Retrying solution " << SolnNum << "." << std::endl;
            ModifyBias(Bias); // Changes bias, usually means increase value of parameters.

            /* We should also change the density matrix to converge to different solution, but it is not
               obvious how we should do that. We could rotate two orbitals, but this may not be enough to
               find a different solution. We could randomize the density matrix, but then we get 
               unphysical results. */
            // NewDensityMatrix(DensityMatrix, CoeffMatrix, OccupiedOrbitals, VirtualOrbitals);
            DensityMatrix = Eigen::MatrixXd::Random(DensityMatrix.rows(), DensityMatrix.cols());
            // GenerateRandomDensity(DensityMatrix);
        }
    }

    AllEnergies.push_back(Energy + Input.Integrals["0 0 0 0"]);

	std::cout << "SCF MetaD: Solution " << SolnNum << " has converged with energy " << Energy + Input.Integrals["0 0 0 0"] << std::endl;
	std::cout << "SCF MetaD: This solution took " << (clock() - ClockStart) / CLOCKS_PER_SEC << " seconds." << std::endl;
	Output << "Solution " << SolnNum << " has converged with energy " << Energy + Input.Integrals["0 0 0 0"] << std::endl;
	Output << "This solution took " << (clock() - ClockStart) / CLOCKS_PER_SEC << " seconds." << std::endl;

    return Energy + Input.Integrals["0 0 0 0"];
}
