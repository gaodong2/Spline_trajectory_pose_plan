#include <iostream>
#include <Eigen/Core> 
#include <Eigen/SVD>
#include <Eigen/Dense>   
#include <vector>
#include "kinematics_6axis_ur.h"
#include "B_Spline.h"
#include "trapezoidal_velocity.h"
#include "rotational_interpolation.h"
#include "SphericalQuadraticInterpolation.h"
#include "S_Type_Velocity.h"


class MultiPoints {

    public:
        MultiPoints();
        ~MultiPoints();
        std::vector<Eigen::Isometry3d> GenerateTrajectory(Eigen::MatrixXd &points, double VxPrm, double AxPrm, double thresholdPrm, double dtPram, const std::string &velType="T");

    private:
        std::vector<Eigen::Isometry3d> TTypeVelPlan(Eigen::MatrixXd &points);
        std::vector<Eigen::Isometry3d> STypeVelPlan(Eigen::MatrixXd &points);
        std::vector<Eigen::Quaterniond> EulerToQuats(Eigen::MatrixXd &points);
        BSplines curve;
        TrapezoidalVelocity Tvel; 
        S_Type_Velocity Svel;
        Squad RI;  
        int seg;
        double scale;
        double period;
        double Vx, Ax;
        double threshold;
        double allow_omega;
        double j_max = 0.1;
};
