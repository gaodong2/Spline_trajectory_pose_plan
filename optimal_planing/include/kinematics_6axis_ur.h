#ifndef KINEMITICS_6AXIS_UR_ONE_H
#define KINEMITICS_6AXIS_UR_ONE_H

#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <iostream>
using namespace std;
//
//Kinematics data struct.
//
typedef struct DH_6axis_ur
{
    double d2 = 157;
    double a3 = 266;
    double a4 = 256.5;
    double d5 = 119;
    double d6 = 102.5;
    double d7 = 94;

} DH_6axis_ur;


//
//Kinematics for 6 axises ur type. The forward rotation of axis z is from cap to axis.
//
class kinematics_6axis_ur_two {

public:
    //Type
    using Vector6 = Eigen::Matrix<double, 1, 6>;
    using Tranform44 = Eigen::Isometry3d;
    const std::string prefix = "IK error: ";

    //Construct function
    kinematics_6axis_ur_two(const DH_6axis_ur &dh);

    //Config DH parameters
    void change_DH(const DH_6axis_ur &dh);

    //Forward kinematics
    Eigen::Isometry3d baseTend(const Vector6 &q);
    Eigen::Isometry3d Transform3D_DH(const double &a, const double &alpha, const double &d, const double &theta);
    void Transfer2Zero(const double &x, const double &y, const double &z, const double &Rx, const double &Ry, const double &Rz, double res[6]);
    void Transfer2Aubo(const double &x, const double &y, const double &z, const double &Rx, const double &Ry, const double &Rz, double res[6]);
    void DeltaTransform(const double base[6], const double delta[6], double res[6]);
    //Inverse kinematics
    int getQ(const Eigen::Isometry3d& transform3D, const Vector6 &q, Vector6 &targetQ, double threahold = 6.28);
    void SixDMoseTransfer(const Eigen::Vector3d &origin, Eigen::Vector3d &outPut);
private:
    int solve(const Eigen::Isometry3d& transform3D, std::vector<Vector6>& joints, double eps = 1e-10);
    double norm2(Vector6 q1, Vector6 q2);
    int filter_joints(const Vector6 &cQ, const Vector6 &tQ, double threahold = 6.28);
    double solveQ0(bool positive);
    double solveQ4(double q0,bool positive);
    double solveQ5(double q0,double q4);
    double solveQ123(double q0,double q4);
    double solveQ2(double q0,double q4,double q123,bool positive);
    double solveQ1(double q0,double q4,double q123,double q2);
    double solveQ3(double q123,double q2,double q1);
    double M(double q0, double q4,double q123);
    double N(double q0, double q4,double q123);
    Eigen::Isometry3d AuboToZero(const Eigen::Vector3d &pos, const Eigen::Matrix3d &ori);
    Eigen::Isometry3d ZeroToAubo(const Eigen::Vector3d &pos, const Eigen::Matrix3d &ori);
    Eigen::Isometry3d VecToMat(const Eigen::Vector3d &pos, const Eigen::Matrix3d &ori);
    Eigen::Isometry3d RotaMpi();
    Eigen::Vector3d rotationMatrixToEulerAngles(Eigen::Matrix3d& R);
    

private:
    //DH
    Eigen::Matrix<double,7,1> _a;
    Eigen::Matrix<double,7,1> _alpha;
    Eigen::Matrix<double,7,1> _d;
    Eigen::Matrix<double,7,1> _theta;
    double d2;
    double a3;
    double a4;
    double d5;
    double d6;
    double d7;

    //tem
    Eigen::Matrix<double, 3, 3> r;
    Eigen::Matrix<double, 3, 1> p;

    //eps
    double _eps = 1e-10;



};

inline double kinematics_6axis_ur_two::norm2(Vector6 q1, Vector6 q2)
{
    Vector6 q = q2 - q1;
    return sqrt(q * q.transpose());
}

#endif // KINEMITICS_6AXIS_UR_ONE_H
