#include <toppra/algorithm/toppra.hpp>
#include <toppra/constraint/linear_joint_acceleration.hpp>
#include <toppra/constraint/linear_joint_velocity.hpp>
#include <toppra/geometric_path/piecewise_poly_path.hpp>
#include <toppra/parametrizer/const_accel.hpp>
#include <toppra/toppra.hpp>
#include "interface_toppra.h"

namespace toppra {
  
  Planning::Planning(){};
  Planning::~Planning(){};

  void Planning::formatVecToMat(const std::vector<Eigen::VectorXd,
                                        Eigen::aligned_allocator<Eigen::VectorXd>>& vec,
                      Eigen::MatrixXd& mat) {
    mat.resize(vec.at(0).rows(), vec.size());
    for (size_t i = 0; i < vec.size(); i++) {
      mat.col(i) = vec.at(i);
    }
  };


  void Planning::JointSpacePlanning(toppra::Vectors &positions, toppra::Vectors &velocities, 
                            toppra::Vector velLimitLower,
                            toppra::Vector velLimitUpper,
                            toppra::Vector accLimitLower,
                            toppra::Vector accLimitUpper,
                            double period, bool printInfo, Eigen::MatrixXd &path_pos2_) {
    // const bool printInfo = true;
    int numJoints = positions[0].size();
    std::cout << numJoints << std::endl;
    // //#### create piecewise polynomial geometric path ####
    // const int numJoints = 2;
    // toppra::Vector position0{numJoints}, position1{numJoints}, position2{numJoints};
    // position0 << 0.0, 0.0;
    // position1 << 1.0, 2.0;
    // position2 << 2.0, 0.0;
    // toppra::Vector velocity0{numJoints}, velocity1{numJoints}, velocity2{numJoints};
    // velocity0 << 0.0, 0.0;
    // velocity1 << 1.0, 1.0;
    // velocity2 << 0.0, 0.0;
    // toppra::Vectors positions = {position0, position1,
    //                             position2};  //[(0, 0), (1, 2), (2, 0)]
    // toppra::Vectors velocities = {velocity0, velocity1,
    //                               velocity2};  //[(0, 0), (1, 1), (0, 0)]
    int length = positions.size();
    std::cout << length << std::endl;
    std::vector<toppra::value_type> steps;
    // for(int i=0;i<length;i++) {
    //   steps.push_back(i);
    // }
    std::cout << length << std::endl;
    steps = std::vector<toppra::value_type>{0, 1, 2, 3};

    std::cout << positions.size() << std::endl;
    std::cout << velocities.size() << std::endl;

    toppra::PiecewisePolyPath hermite =
        toppra::PiecewisePolyPath::CubicHermiteSpline(positions, velocities, steps);
    toppra::GeometricPathPtr path;
    path = std::make_shared<toppra::PiecewisePolyPath>(hermite);
    std::cout << length << std::endl;
    //#### create linear joint-space constraints ####
    // toppra::Vector velLimitLower = Eigen::VectorXd::Zero(numJoints);
    // toppra::Vector velLimitUpper = Eigen::VectorXd::Zero(numJoints);
    // toppra::Vector accLimitLower = Eigen::VectorXd::Zero(numJoints);
    // toppra::Vector accLimitUpper = Eigen::VectorXd::Zero(numJoints);
    // velLimitLower << -1, -1;
    // velLimitUpper << 1, 1;
    // accLimitLower << -0.2, -0.3;
    // accLimitUpper << 0.2, 0.3;

    toppra::LinearConstraintPtr ljv, lja;
    ljv = std::make_shared<toppra::constraint::LinearJointVelocity>(
        velLimitLower, velLimitUpper);  //[[-1, 1], [-0.5, 0.5]]
    lja = std::make_shared<toppra::constraint::LinearJointAcceleration>(
        accLimitLower, accLimitUpper);  //[[-0.05, 0.2], [-0.1, 0.3]]
    lja->discretizationType(toppra::DiscretizationType::Interpolation);
    toppra::LinearConstraintPtrs constraints{ljv, lja};
    std::cout << length << std::endl;
    //#### create toppra ####
    toppra::PathParametrizationAlgorithmPtr algo =
        std::make_shared<toppra::algorithm::TOPPRA>(constraints, path);
    toppra::ReturnCode rc1 = algo->computePathParametrization(0, 0);
    if (printInfo) std::cout << "rc1 = " << int(rc1) << std::endl;


    toppra::ParametrizationData pd = algo->getParameterizationData();
    if (printInfo)
      std::cout << "pd.gridpoints \n " << pd.gridpoints.transpose() << std::endl;
    if (printInfo)
      std::cout << "pd.parametrization \n " << pd.parametrization.transpose()
                << std::endl;
    if (printInfo)
      std::cout << "pd.controllable_sets \n " << pd.controllable_sets.transpose()
                << std::endl;
    if (printInfo)
      std::cout << "pd.feasible_sets \n " << pd.feasible_sets.transpose() << std::endl;
    if (printInfo) std::cout << "pd.ret_code = " << int(pd.ret_code) << std::endl;

    //#### create constant accelaration parametrizer ####
    toppra::Vector gridpoints =
        pd.gridpoints;  // Grid-points used for solving the discretized problem.
    toppra::Vector vsquared =
        pd.parametrization;  // Output parametrization (squared path velocity)
    std::shared_ptr<toppra::parametrizer::ConstAccel> ca =
        std::make_shared<toppra::parametrizer::ConstAccel>(path, gridpoints, vsquared);
    if (printInfo) std::cout << "ca->validate() = " << ca->validate() << std::endl;

    Eigen::Matrix<toppra::value_type, 1, 2> interval2;
    interval2 = ca->pathInterval();
    if (printInfo) std::cout << "interval2 = " << interval2 << std::endl;

    const int length2 = int((interval2(1) - interval2(0)) / period);
    toppra::Vector times2 =
        toppra::Vector::LinSpaced(length2, interval2(0), interval2(1));
    toppra::Vectors path_pos2;
    path_pos2 = ca->eval(times2, 0);  // TODO this function call fails
    toppra::Vectors path_vel2;
    path_vel2 = ca->eval(times2, 1);  // TODO this function call fails
    toppra::Vectors path_acc2;
    path_acc2 = ca->eval(times2, 2);  // TODO this function call fails
    // Eigen::MatrixXd path_pos2_ = Eigen::MatrixXd::Zero(numJoints, length2);
    Eigen::MatrixXd path_vel2_ = Eigen::MatrixXd::Zero(numJoints, length2);
    Eigen::MatrixXd path_acc2_ = Eigen::MatrixXd::Zero(numJoints, length2);
    this->formatVecToMat(path_pos2, path_pos2_);
    this->formatVecToMat(path_vel2, path_vel2_);
    this->formatVecToMat(path_acc2, path_acc2_);

    // std::cout << "path_vel2::" << path_vel2_ << std::endl;
    // std::cout << "path_acc2::" << path_acc2_ << std::endl;

    if (printInfo) std::cout << "path_pos2_\n " << path_pos2_ << std::endl;
    if (printInfo) std::cout << "path_vel2_\n " << path_vel2_ << std::endl;
    if (printInfo) std::cout << "path_acc2_\n " << path_acc2_ << std::endl;
    if (printInfo) std::cout << "times2 \n " << times2.transpose() << std::endl;
  }

  int Planning::CartesianSpacePlanning(toppra::Vector &currentQ, toppra::Vectors &positions, toppra::Vectors &velocities, 
                            toppra::Vector velLimitLower,
                            toppra::Vector velLimitUpper,
                            toppra::Vector accLimitLower,
                            toppra::Vector accLimitUpper,
                            double period, bool printInfo, Eigen::MatrixXd &path) {
    Eigen::Matrix<double, 1, 6> cur_Q, targetQ;
    cur_Q << currentQ(0,0), currentQ(1,0), currentQ(2,0), currentQ(3,0), currentQ(4,0), currentQ(5,0);

    Eigen::Isometry3d Obj;
    Obj = aubo_i3.baseTend(currentQ);
    aubo_i3.getQ(Obj, cur_Q, targetQ, threshold);           
    if(CheckSingularity(currentQ, positions)) {
      if (printInfo) std::cout << "Planning failed, waypoints singularity\n "<< std::endl;
      return -1;
    }

    Eigen::MatrixXd path_pos2_;
    std::cout << path_pos2_ << std::endl;
    JointSpacePlanning(positions, velocities, velLimitLower, velLimitUpper, accLimitLower, accLimitUpper, period, printInfo, path_pos2_);
    std::cout << path_pos2_ << std::endl;
    if(TransCartesianToJoint(currentQ, path_pos2_, path) < 0){
      if (printInfo) std::cout << "Planning failed, waypoints singularity\n "<< std::endl;
      return -1;
    }
    return 0;
  }

  int Planning::TransCartesianToJoint(toppra::Vector &currentQ, Eigen::MatrixXd &path_pos2_, Eigen::MatrixXd &path) {
    int rows = path_pos2_.rows();
    int cols = path_pos2_.cols();
    path.resize(rows, cols);
    Eigen::Isometry3d Target = Eigen::Isometry3d::Identity();
    Eigen::Matrix<double, 1, 6> cur_Q, targetQ;
    cur_Q << currentQ(0,0), currentQ(1,0), currentQ(2,0), currentQ(3,0), currentQ(4,0), currentQ(5,0);
    for(int i=0;i<cols;i++) {
      Eigen::Vector3d Euler = {path_pos2_(3,i), path_pos2_(4,i), path_pos2_(5,i)};
      Eigen::Vector3d trans = {path_pos2_(0,i), path_pos2_(1,i), path_pos2_(2,i)};
      Eigen::Matrix3d rotm = getRotationMatrixFromEuler(Euler, eulerType);
      Target.translate(trans);
      Target.rotate(rotm);
      if(aubo_i3.getQ(Target, cur_Q, targetQ, threshold) < 0) {
        return -1;
      }
      for(int j=0;j<6;j++) {
        
        path(j,i) = targetQ(0, j);
      }
      Target = Eigen::Isometry3d::Identity();
      cur_Q << targetQ(0,0), targetQ(0,1), targetQ(0,2), targetQ(0,3), targetQ(0,4), targetQ(0,5);
    }
    return 0;
  }

  bool Planning::CheckSingularity(toppra::Vector &currentQ, toppra::Vectors &positions){
    int cols = positions.size();
    std::cout << cols << std::endl;
    Eigen::Isometry3d Target = Eigen::Isometry3d::Identity();
    Eigen::Matrix<double, 1, 6> cur_Q, targetQ;
    cur_Q << currentQ(0,0), currentQ(1,0), currentQ(2,0), currentQ(3,0), currentQ(4,0), currentQ(5,0);
    std::cout << cur_Q << std::endl;
    for(int i=0;i<cols;i++) {
      Eigen::Vector3d Euler = {positions[i][3], positions[i][4], positions[i][5]};
      std::cout << Euler << std::endl;
      Eigen::Vector3d trans = {positions[i][0], positions[i][1], positions[i][2]};
      Eigen::Matrix3d rotm = getRotationMatrixFromEuler(Euler, eulerType);
      Target.translate(trans);
      Target.rotate(rotm);
      
      std::cout << Target.rotation() << std::endl;
      std::cout << Target.translation() << std::endl;
      if(aubo_i3.getQ(Target, cur_Q, targetQ, threshold) < 0) {
        std::cout << "CheckSingularity" << std::endl;
        return true;
      }
      std::cout << "targetQ: "<< targetQ*180/M_PI << std::endl;
      Target = Eigen::Isometry3d::Identity();
    }
    return false;
  }

  Eigen::Matrix3d Planning::getRotationMatrixFromEuler(const Eigen::Vector3d &eulerAngle, int eulerType) {
    Eigen::Matrix3d m = Eigen::Matrix3d::Identity();

    if (eulerType == 210)
    {
        m = Eigen::AngleAxisd(eulerAngle[0], Eigen::Vector3d::UnitZ()) * Eigen::AngleAxisd(eulerAngle[1], Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(eulerAngle[2], Eigen::Vector3d::UnitX());
    }
    else if (eulerType == 012)
    {
        m = Eigen::AngleAxisd(eulerAngle[0], Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(eulerAngle[1], Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(eulerAngle[2], Eigen::Vector3d::UnitZ());
    }
    else
    {
        cout << "Wrong euler transform number: " << eulerType << endl;
    }

    return m;
}

Eigen::Vector3d Planning::getPosEulerFromXYZ(const Eigen::Matrix3d &R, double Rref[3])
{ // Rx*Ry*Rz out: xyz

    Eigen::Vector3d pr;
    double x[2], y[2], z[2];

    //-2Pi~2Pi
    if (R(0, 2) > 0) // pr[4] 0~Pi, -2*Pi~-Pi
    {                // total 4 result
        // 1:0~Pi/2
        x[0] = atan2(-R(1, 2), R(2, 2));
        y[0] = atan2(R(0, 2), sqrt(R(0, 0) * R(0, 0) + R(0, 1) * R(0, 1)));
        z[0] = atan2(-R(0, 1), R(0, 0));
        // 2:Pi/2~Pi
        x[1] = atan2(R(1, 2), -R(2, 2));
        y[1] = M_PI - y[0];
        z[1] = atan2(R(0, 1), -R(0, 0));
    }
    else // pr[4] 0~-Pi, Pi~2*Pi
    {
        // 1:0~-Pi/2
        x[0] = atan2(-R(1, 2), R(2, 2));
        y[0] = atan2(R(0, 2), sqrt(R(0, 0) * R(0, 0) + R(0, 1) * R(0, 1)));
        z[0] = atan2(-R(0, 1), R(0, 0));
        // 2:-Pi/2~-Pi
        x[1] = atan2(R(1, 2), -R(2, 2));
        y[1] = -M_PI - y[0];
        z[1] = atan2(R(0, 1), -R(0, 0));
    }

    int sel = 0;
    if (fabs(Rref[1] - y[0]) < 1.57)
    {
        pr[1] = y[0];
    }
    else if (fabs(Rref[1] - y[1]) < 1.57)
    {
        pr[1] = y[1];
        sel = 1;
    }
    else if (Rref[1] - y[1] > 1.5 * M_PI)
    {
        pr[1] = 2 * M_PI + y[1];
        sel = 1;
    }
    else if (Rref[1] - y[1] < -1.5 * M_PI)
    {
        pr[1] = -2 * M_PI + y[1];
        sel = 1;
    }
    else if (Rref[1] - y[0] > 1.5 * M_PI)
    {
        pr[1] = 2 * M_PI + y[0];
    }
    else if (Rref[1] - y[0] < -1.5 * M_PI)
    {
        pr[1] = -2 * M_PI + y[0];
    }

    if (Rref[0] - x[sel] > 1.5 * M_PI)
    {
        pr[0] = 2 * M_PI + x[sel];
    }
    else if (Rref[0] - x[sel] < -1.5 * M_PI)
    {
        pr[0] = -2 * M_PI + x[sel];
    }
    else
        pr[0] = x[sel];
    if (Rref[2] - z[sel] > 1.5 * M_PI)
    {
        pr[2] = 2 * M_PI + z[sel];
    }
    else if (Rref[2] - z[sel] < -1.5 * M_PI)
    {
        pr[2] = -2 * M_PI + z[sel];
    }
    else
        pr[2] = z[sel];
    return pr;
}

Eigen::Vector3d Planning::getPosEulerFromZYX(const Eigen::Matrix3d &R, double Rref[3]) {
    // Rz*Ry*Rx out: xyz
    Eigen::Vector3d pr;
    Eigen::Matrix3d Rt = R.transpose();
    pr = getPosEulerFromXYZ(Rt, Rref);
    return -pr;
}

}