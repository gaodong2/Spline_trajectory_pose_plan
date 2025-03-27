#include <iostream>
#include <iostream>
#include <Eigen/Core> 
#include <Eigen/SVD>
#include <Eigen/Dense>   
#include <vector>
#include <cmath>


class S_Type_Velocity {
    
    public:
        S_Type_Velocity();
        ~S_Type_Velocity();
        Eigen::MatrixXd STypeVelPlan(Eigen::MatrixXd &L, Eigen::MatrixXd &kappa, double threshold, double Vx, double Ax, double j_max, double scale, double period);

    private:
        void TrapezoidalPlan(Eigen::MatrixXd &V_reachability, Eigen::MatrixXd &s, double a_max);
        Eigen::MatrixXd P2PMultiAxisDoubleSTrajectory(double sampleTime, Eigen::MatrixXd input);
        Eigen::MatrixXd VelocityPlanning(Eigen::MatrixXd &L, Eigen::MatrixXd &v_max, double a_sup, double v0, double vt);
        Eigen::MatrixXd MaxVel(Eigen::MatrixXd &kappa, double threshold, double Vs);
        Eigen::MatrixXd STypeSpeedUp(Eigen::MatrixXd &TJAVS);
        Eigen::MatrixXd STypeReach(Eigen::MatrixXd &TJAVSf, Eigen::MatrixXd &newTJAVS);
        Eigen::MatrixXd FindSpeedUpStage(Eigen::MatrixXd &V_reachability, Eigen::MatrixXd &ssss, double Vmax, double a_sup, double j_max);
        Eigen::MatrixXd FindShutDownStage(Eigen::MatrixXd &V_reachability, Eigen::MatrixXd &ssss, double Vmax, double a_sup, double j_max);
        Eigen::MatrixXd Compose(double dt, Eigen::MatrixXd &newTJAVSf);  
        double Sum(Eigen::MatrixXd &TTTT, int ii, int jj);
        Eigen::MatrixXd Tp, Tpf;
        Eigen::MatrixXd Jp, Jpf;
        Eigen::MatrixXd Ap, Apf;
        Eigen::MatrixXd Vp, Vpf;
        Eigen::MatrixXd Sp, Spf;
        Eigen::MatrixXd DelS, sumS, sumT, sumV, sumA;
        Eigen::MatrixXd TJAVS, TJAVSf;
};