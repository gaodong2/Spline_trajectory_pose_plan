#include <iostream>
#include <Eigen/Core> 
#include <Eigen/SVD>
#include <Eigen/Dense>   
#include <vector>
#include <cmath>

class BSplines {
    public:
        BSplines();
        ~BSplines();
        void UniformCubicBSplines(Eigen::MatrixXd &points, int seg=500);
        void UniformCubicBSplines(Eigen::MatrixXd &points, Eigen::MatrixXd &u);
        double GetTotalL();
        Eigen::MatrixXd GetCurvePos();
        Eigen::MatrixXd GetCurveVel();
        Eigen::MatrixXd GetCurveAcc();
        Eigen::MatrixXd GetCurveJerk();
        Eigen::MatrixXd GetCurvature();
        Eigen::MatrixXd GetDeltaL();

    private:
        Eigen::MatrixXd KnotsUniform(int len, int dim);
        Eigen::MatrixXd LinSpace(int seg);
        Eigen::MatrixXd BSplineBasis(int p, Eigen::MatrixXd &knots, double u);
        double BSsplineBasisDerivative(int p, Eigen::MatrixXd &knots, int i, double u, int order);
        double BSplineBasisValue(int p, Eigen::MatrixXd &knots, int i, double u);
        bool AlmostEqual(double data1, double data2, double eps=1e-6);
        double CalTotalL();
        double TrapezoidalIntegral(double deltaX_1, double deltaX, double Y);
        
        Eigen::MatrixXd uniformBasisMatrix;
        Eigen::MatrixXd C;  // 位置
        Eigen::MatrixXd V;  // 速度
        Eigen::MatrixXd A;  // 加速度
        Eigen::MatrixXd J;  // 加加速度
        Eigen::MatrixXd kappa;  // 曲率半径
        Eigen::MatrixXd Ls;
};