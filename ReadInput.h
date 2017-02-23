#include <Eigen/Dense>

class InputObj
{
    public:
        void GetInputName();
        void SetNames(char*, char*, char*);
        void Set();
        std::map< std::string, double > Integrals;
        Eigen::MatrixXd OverlapMatrix;
        std::string IntegralsInput;
        std::string OverlapInput;
        std::string OutputName;
		unsigned short int NumOcc;
		unsigned short int NumSoln;
        unsigned short int NumElectrons;
        unsigned short int NumAO;
};