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
    private:
        unsigned short int NumAO;
};