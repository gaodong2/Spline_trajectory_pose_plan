#include <toppra/algorithm/toppra.hpp>
#include <toppra/constraint/linear_joint_acceleration.hpp>
#include <toppra/constraint/linear_joint_velocity.hpp>
#include <toppra/geometric_path/piecewise_poly_path.hpp>
#include <toppra/parametrizer/const_accel.hpp>
#include <toppra/toppra.hpp>
#include <Eigen/Core> 
#include <Eigen/SVD>
#include <Eigen/Dense>   
#include "kinematics_6axis_ur.h"

namespace toppra {
class Planning

{
 public:
    Planning();
    ~Planning();
    void JointSpacePlanning(toppra::Vectors &positions, toppra::Vectors &velocities, 
                              toppra::Vector velLimitLower,
                              toppra::Vector velLimitUpper,
                              toppra::Vector accLimitLower,
                              toppra::Vector accLimitUpper,
                              double period, bool printInfo, Eigen::MatrixXd &path_pos2_);
    int CartesianSpacePlanning(toppra::Vector &currentQ, toppra::Vectors &positions, toppra::Vectors &velocities, 
                            toppra::Vector velLimitLower,
                            toppra::Vector velLimitUpper,
                            toppra::Vector accLimitLower,
                            toppra::Vector accLimitUpper,
                            double period, bool printInfo, Eigen::MatrixXd &path_pos2_);
  private:
    void formatVecToMat(const std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>>& vec, Eigen::MatrixXd& mat);
    int TransCartesianToJoint(toppra::Vector &currentQ, Eigen::MatrixXd &path_pos2_, Eigen::MatrixXd &path);
    bool CheckSingularity(toppra::Vector &currentQ, toppra::Vectors &positions);
    Eigen::Vector3d getPosEulerFromZYX(const Eigen::Matrix3d &R, double Rref[3]);
    Eigen::Vector3d getPosEulerFromXYZ(const Eigen::Matrix3d &R, double Rref[3]);
    Eigen::Matrix3d getRotationMatrixFromEuler(const Eigen::Vector3d &eulerAngle, int eulerType);
    // DH_6axis_ur dh;
    DH_6axis_ur dh;
    kinematics_6axis_ur_two aubo_i3 = kinematics_6axis_ur_two(dh);
    double threshold = 6.28;
    double Rref[3] = {0, 0, 0};
    int eulerType = 210;;
};

}
