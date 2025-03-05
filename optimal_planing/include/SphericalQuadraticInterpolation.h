#include <iostream>
#include <Eigen/Core> 
#include <Eigen/SVD>
#include <Eigen/Dense>   
#include <vector>
#include <cmath>

class Squad {
    public:
        Squad();
        ~Squad();
        Eigen::Quaterniond QuadSquad(std::vector<Eigen::Quaterniond> &q, double t);
        double ComputeAngularVelocityAndAcceleration(std::vector<Eigen::Quaterniond> &q, Eigen::MatrixXd &dt);
    private:
        double sCurve(double t);
        bool AlmostEqual(double data1, double data2);
        Eigen::Quaterniond quatExp(Eigen::Quaterniond &v);
        Eigen::Quaterniond quatMultiply(Eigen::Quaterniond &q1, Eigen::Quaterniond &q2);
        Eigen::Quaterniond quatInverse(Eigen::Quaterniond &q);
        Eigen::Quaterniond quatNormalize(Eigen::Quaterniond &q);
        Eigen::Quaterniond quatlog(Eigen::Quaterniond &q);
        Eigen::Quaterniond Slerp(Eigen::Quaterniond &q1, Eigen::Quaterniond &q2, double t);
        Eigen::Quaterniond GetIntermediateControlPoint(int j, std::vector<Eigen::Quaterniond> &q);
        double EvalAlpha(double s, int i, int L);
        double eps = 1e-16;
        double FindMax(Eigen::MatrixXd &omegaNormal);
};