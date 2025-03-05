#include "trapezoidal_velocity.h"

TrapezoidalVelocity::TrapezoidalVelocity() {}

TrapezoidalVelocity::~TrapezoidalVelocity() {}

Eigen::MatrixXd TrapezoidalVelocity::TrapezoidalVelPlan(Eigen::MatrixXd &L, Eigen::MatrixXd &kappa, double period, double Vx, double Ax, double threshold, double scale) {
    double v0 = 0;
    double vt = 0;
    Eigen::MatrixXd Vmax =  MaxVel(kappa, threshold, Vx);
    Eigen::MatrixXd V_reachability = VelocityPlanning(L, Vmax, Ax, v0, vt);
    for(int i=0;i<V_reachability.rows();i++){
        V_reachability(i,0) = V_reachability(i,0)*scale;
    }
    TrapezoidalPlan(V_reachability, L, Ax);
    Eigen::MatrixXd uus = CalcSinPeriod(period);
    return uus;
}

Eigen::MatrixXd TrapezoidalVelocity::VelocityPlanning(Eigen::MatrixXd &L, Eigen::MatrixXd &v_max, double a_sup, double v0, double vt) {
    // 输入参数
    int rows = L.rows();
    Eigen::MatrixXd a_max = Eigen::MatrixXd::Ones(rows, 1) * a_sup;  // 各点的最大加速度 (m/s^2)
    Eigen::MatrixXd v = Eigen::MatrixXd::Ones(rows, 1);
    double delta_s;
    // 前向传播
    Eigen::MatrixXd v_forward = Eigen::MatrixXd::Ones(rows, 1);
    v_forward(0,0) = v0;
    for(int i=1;i<rows;i++){
        delta_s = L(i,0);
        v_forward(i,0) = std::min(v_max(i,0), sqrt(v_forward(i-1,0)*v_forward(i-1,0) + 2 * a_max(i,0) * delta_s));
    }

    // 反向传播
    Eigen::MatrixXd v_backward = Eigen::MatrixXd::Ones(rows, 1);
    v_backward(rows-1,0) = vt;
    for (int i=rows-2;i>0;i--){
        delta_s = L(i+1,0);
        v_backward(i,0) = std::min(v_max(i,0), sqrt(v_backward(i+1,0)*v_backward(i+1,0) + 2 * a_max(i,0) * delta_s));
        
    }

    // 可达速度
    for(int i=0;i<rows;i++) {
        v(i,0) = std::min(v_forward(i,0), v_backward(i,0));   
    }
    
    return v;  
}

Eigen::MatrixXd TrapezoidalVelocity::MaxVel(Eigen::MatrixXd &kappa, double threshold, double Vs) {
    int rows = kappa.rows();
    int cols = kappa.cols();
    Eigen::MatrixXd v_max = Eigen::MatrixXd::Zero(rows, cols);
    for(int i=0;i<rows;i++) {
        v_max(i,0) = Vs;
        if(kappa(i,0) < threshold)
            v_max(i,0) = Vs * kappa(i,0) / threshold;
    }
    return v_max;
}

void TrapezoidalVelocity::TrapezoidalPlan(Eigen::MatrixXd &V_reachability, Eigen::MatrixXd &s, double a_max) {
    int rows = V_reachability.rows();
    int cols = V_reachability.cols();
    DelS = Eigen::MatrixXd::Zero(rows, cols);
    sumS = Eigen::MatrixXd::Zero(rows, cols);
    sumT = Eigen::MatrixXd::Zero(rows, cols);
    sumV = Eigen::MatrixXd::Zero(rows, cols);
    sumA = Eigen::MatrixXd::Zero(rows, cols);
    double sum_T = 0, sum_S = 0;
    for(int i=0;i<rows-1;i++){
        double v0 = V_reachability(i,0);
        double vt = V_reachability(i+1,0); 
        double v_max = std::max(v0, vt); 
        double TTT = 2 * s(i+1,0)  / (v0 + vt);
        double a = (vt - v0) / TTT;
        
        double s = v0 * TTT + 0.5 * a * TTT*TTT;
        sum_T += TTT;
        sum_S += s;
        DelS(i+1,0) = s;
        sumS(i+1,0) = sum_S;
        sumT(i+1,0) = sum_T;
        sumV(i+1,0) = vt;
        sumA(i+1,0) = a;
    }
}

Eigen::MatrixXd TrapezoidalVelocity::GetTotalS() {
    return sumS;
}

Eigen::MatrixXd TrapezoidalVelocity::GetTotalT() {
    return sumT;
}

Eigen::MatrixXd TrapezoidalVelocity::GetTotalV() {
    return sumV;
}

Eigen::MatrixXd TrapezoidalVelocity::GetTotalA() {
    return sumA;
}

Eigen::MatrixXd TrapezoidalVelocity::CalcSinPeriod(double period, const std::string &CurveType) {
    int rows = sumT.rows();
    int urows = (sumT(rows-1,0) + period) / period;
    double ss = 0;
    Eigen::MatrixXd uus(urows, 1);
    if (CurveType == "B-Spline") {
        for(int i=0;i<urows-1;i++) {
            double toltalT = period * i;
            for(int j=0;j<rows-1;j++){
                if((sumT(j,0)<=toltalT) && (sumT(j+1,0)>toltalT)) {
                    double deltat = toltalT - sumT(j,0);
                    ss = sumV(j,0) * deltat + 0.5 * sumA(j+1,0) * deltat * deltat;
                    double uu = ss / DelS(j+1,0) / (rows-1) + (double) (j) / (double)(rows-1);
                    uus(i,0) = uu;                                             
                }
            }     
        }
        uus(urows-1, 0) = 1.0;   
    }
    return uus;
}