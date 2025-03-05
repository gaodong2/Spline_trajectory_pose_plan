#include <iostream>
#include <iostream>
#include <Eigen/Core> 
#include <Eigen/SVD>
#include <Eigen/Dense>   
#include <vector>
#include <cmath>


class TrapezoidalVelocity {


    public:
        TrapezoidalVelocity();
        ~TrapezoidalVelocity();
        Eigen::MatrixXd GetTotalS();
        Eigen::MatrixXd GetTotalT();
        Eigen::MatrixXd GetTotalV();
        Eigen::MatrixXd GetTotalA();
        Eigen::MatrixXd TrapezoidalVelPlan(Eigen::MatrixXd &L, Eigen::MatrixXd &kappa, double period, double Vx, double Ax, double threshold, double scale);

    private:
        Eigen::MatrixXd VelocityPlanning(Eigen::MatrixXd &L, Eigen::MatrixXd &v_max, double a_sup, double v0, double vt);
        Eigen::MatrixXd MaxVel(Eigen::MatrixXd &kappa, double threshold, double Vs);
        void TrapezoidalPlan(Eigen::MatrixXd &V_reachability, Eigen::MatrixXd &s, double a_max);
        Eigen::MatrixXd CalcSinPeriod(double period, const std::string &CurveType="B-Spline");
        Eigen::MatrixXd DelS;
        Eigen::MatrixXd sumS;
        Eigen::MatrixXd sumT;
        Eigen::MatrixXd sumV;
        Eigen::MatrixXd sumA;
};