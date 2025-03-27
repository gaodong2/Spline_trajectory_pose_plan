#include <iostream>
#include <math.h>
#include "kinematics_6axis_ur.h"

using namespace std;

kinematics_6axis_ur_two::kinematics_6axis_ur_two(const DH_6axis_ur &dh) :
    d2(dh.d2),
    a3(dh.a3),
    a4(dh.a4),
    d5(dh.d5),
    d6(dh.d6),
    d7(dh.d7)
{
    _a     << 0,    0,      a3,      a4,   0,       0,      0;
    _alpha << M_PI, M_PI_2, 0,       M_PI, -M_PI_2, M_PI_2, 0;
    _d     << 0,    -d2,    0,       0,    d5,      d6,     d7;
    _theta << 0,    M_PI_2, -M_PI_2, 0,    -M_PI_2, 0,      0;
}

void kinematics_6axis_ur_two::change_DH(const DH_6axis_ur &dh)
{
    d2 = dh.d2;
    a3 = dh.a3;
    a4 = dh.a4;
    d5 = dh.d5;
    d6 = dh.d6;
    d7 = dh.d7;

    _a     << 0,    0,      a3,      a4,   0,       0,      0;
    _alpha << M_PI, M_PI_2, 0,       M_PI, -M_PI_2, M_PI_2, 0;
    _d     << 0,    -d2,    0,       0,    d5,      d6,     d7;
    _theta << 0,    M_PI_2, -M_PI_2, 0,    -M_PI_2, 0,      0;

}

Eigen::Vector3d kinematics_6axis_ur_two::rotationMatrixToEulerAngles(Eigen::Matrix3d& R)
{
	double sy = sqrt(R(0, 0) * R(0, 0) + R(1, 0) * R(1, 0));
	bool singular = sy < 1e-6;
	double x, y, z;
	if (!singular) {
		x = atan2(R(2, 1), R(2, 2));
		y = atan2(-R(2, 0), sy);
		z = atan2(R(1, 0), R(0, 0));
	}
	else {
		x = atan2(-R(1, 2), R(1, 1));
		y = atan2(-R(2, 0), sy);
		z = 0;
	}
	return { x,y,z };
}

void kinematics_6axis_ur_two::Transfer2Zero(const double &x, const double &y, const double &z, const double &Rx, const double &Ry, const double &Rz, double res[6]) {
    // Get pos and ori 
    Eigen::Vector3d pos = {x, y, z};
    Eigen::Matrix3d ori;
    ori = Eigen::AngleAxisd(Rz, Eigen::Vector3d::UnitZ())*
		  Eigen::AngleAxisd(Ry, Eigen::Vector3d::UnitY())*
		  Eigen::AngleAxisd(Rx, Eigen::Vector3d::UnitX());

    Eigen::Isometry3d Zero;
    Zero = AuboToZero(pos, ori);
    Eigen::Matrix3d Mat = Zero.rotation();
    Eigen::Vector3d rpy = rotationMatrixToEulerAngles(Mat);
    // return nums [x, y, x, rx, ry, ry]
    for(int i=0;i<3;i++) {
       res[i] = Zero.translation()[i];
       res[i+3] = rpy[i]; 
    }
}

void kinematics_6axis_ur_two::Transfer2Aubo(const double &x, const double &y, const double &z, const double &Rx, const double &Ry, const double &Rz, double res[6]) {
    // Get pos and ori 
    Eigen::Vector3d pos = {x, y, z};
    Eigen::Matrix3d ori =
	ori = Eigen::AngleAxisd(Rz, Eigen::Vector3d::UnitZ())*
		  Eigen::AngleAxisd(Ry, Eigen::Vector3d::UnitY())*
		  Eigen::AngleAxisd(Rx, Eigen::Vector3d::UnitX());

    Eigen::Isometry3d Aubo;
    Aubo = ZeroToAubo(pos, ori);
    Eigen::Matrix3d Mat = Aubo.rotation();
    Eigen::Vector3d rpy = rotationMatrixToEulerAngles(Mat);
    // return nums [x, y, x, rx, ry, rz]
    for(int i=0;i<3;i++) {
       res[i] = Aubo.translation()[i];
       res[i+3] = rpy[i]; 
    }
}

void kinematics_6axis_ur_two::DeltaTransform(const double base[6], const double delta[6], double res[6]) {
    // data form [x, y, x, rx, ry, rz]
    Eigen::Vector3d basePos = {base[0], base[1], base[2]};
    Eigen::Matrix3d baseOri =
	baseOri = Eigen::AngleAxisd(base[5], Eigen::Vector3d::UnitZ())*
		  Eigen::AngleAxisd(base[4], Eigen::Vector3d::UnitY())*
		  Eigen::AngleAxisd(base[3], Eigen::Vector3d::UnitX());

    Eigen::Vector3d deltaPos = {delta[0], delta[1], delta[2]};
    Eigen::Matrix3d deltaOri =
	deltaOri = Eigen::AngleAxisd(delta[5], Eigen::Vector3d::UnitZ())*
		  Eigen::AngleAxisd(delta[4], Eigen::Vector3d::UnitY())*
		  Eigen::AngleAxisd(delta[3], Eigen::Vector3d::UnitX());

    Eigen::Isometry3d baseMat, delteMat, outPut;
    baseMat = VecToMat(basePos, baseOri);
    delteMat = VecToMat(deltaPos, deltaOri);
    outPut = baseMat * delteMat;

    Eigen::Matrix3d Mat = outPut.rotation();
    Eigen::Vector3d rpy = rotationMatrixToEulerAngles(Mat);
    // return nums [x, y, x, rx, ry, rz]
    for(int i=0;i<3;i++) {
       res[i] = outPut.translation()[i];
       res[i+3] = rpy[i]; 
    }
}

// after forward kinematics 
Eigen::Isometry3d kinematics_6axis_ur_two::AuboToZero(const Eigen::Vector3d &pos, const Eigen::Matrix3d &ori) {
    Eigen::Isometry3d Zero2Base = Eigen::Isometry3d::Identity();
    Eigen::Isometry3d Zero2Aubo  = RotaMpi();
    Eigen::Isometry3d Aubo2Base = VecToMat(pos, ori);

    Zero2Base = Zero2Aubo * Aubo2Base;
    return Zero2Base;
}
// before inverse kinematics 
Eigen::Isometry3d kinematics_6axis_ur_two::ZeroToAubo(const Eigen::Vector3d &pos, const Eigen::Matrix3d &ori) {
    Eigen::Isometry3d Aubo2Base = Eigen::Isometry3d::Identity();
    Eigen::Isometry3d Zero2Aubo = RotaMpi();
    Eigen::Isometry3d Zero2Base = VecToMat(pos, ori);

    Aubo2Base = Zero2Aubo.inverse() * Zero2Base;
    return Aubo2Base;
}
void kinematics_6axis_ur_two::SixDMoseTransfer(const Eigen::Vector3d &origin, Eigen::Vector3d &outPut) {
    // origin [Rx, Ry, Rx]
    // outPut [Rx, Ry, Rx]
    Eigen::Matrix3d Transform, output;
    Transform(0,0) = 1;
    Transform(1,0) = 0;         
    Transform(2,0) = 0;

    Transform(0,1) = 0;
    Transform(1,1) = cos(-M_PI/4);
    Transform(2,1) = -sin(-M_PI/4);

    Transform(0,2) = 0;
    Transform(1,2) = sin(-M_PI/4);
    Transform(2,2) = cos(-M_PI/4);

    Eigen::Matrix3d ori =
	ori = Eigen::AngleAxisd(origin[2], Eigen::Vector3d::UnitZ())*
		  Eigen::AngleAxisd(origin[1], Eigen::Vector3d::UnitY())*
		  Eigen::AngleAxisd(origin[0], Eigen::Vector3d::UnitX());
    
    output = ori * Transform;
    outPut = output.eulerAngles(0, 1, 2);
    std::cout << outPut[0] << " " << outPut[1] << " " << outPut[2] << std::endl;
} 

Eigen::Isometry3d kinematics_6axis_ur_two::RotaMpi() {
    Eigen::Isometry3d tranDH;
    tranDH(0, 0) = 0;
    tranDH(1, 0) = -1;
    tranDH(2, 0) = 0;
    tranDH(3, 0) = 0;

    tranDH(0, 1) = 1;
    tranDH(1, 1) = 0;
    tranDH(2, 1) = 0;
    tranDH(3, 1) = 0;

    tranDH(0, 2) = 0;
    tranDH(1, 2) = 0;
    tranDH(2, 2) = 1;
    tranDH(3, 2) = 0;

    tranDH(0, 3) = 0;
    tranDH(1, 3) = 0;
    tranDH(2, 3) = 0;
    tranDH(3, 3) = 1;

    return tranDH;
}

Eigen::Isometry3d kinematics_6axis_ur_two::VecToMat(const Eigen::Vector3d &pos, const Eigen::Matrix3d &ori) {
    Eigen::Isometry3d TMat = Eigen::Isometry3d::Identity();
    TMat.rotate(ori);
    TMat.pretranslate(pos);
    return TMat;
}

Eigen::Isometry3d kinematics_6axis_ur_two::Transform3D_DH(const double &a, const double &alpha, const double &d, const double &theta)
{
    Eigen::Isometry3d tranDH;
    tranDH(0, 0) = cos(theta);
    tranDH(1, 0) = sin(theta);
    tranDH(2, 0) = 0;
    tranDH(3, 0) = 0;

    tranDH(0, 1) = -sin(theta)*cos(alpha);
    tranDH(1, 1) = cos(theta)*cos(alpha);
    tranDH(2, 1) = sin(alpha);
    tranDH(3, 1) = 0;

    tranDH(0, 2) = sin(theta)*sin(alpha);
    tranDH(1, 2) = -cos(theta)*sin(alpha);
    tranDH(2, 2) = cos(alpha);
    tranDH(3, 2) = 0;

    tranDH(0, 3) = a*cos(theta);
    tranDH(1, 3) = a*sin(theta);
    tranDH(2, 3) = d;
    tranDH(3, 3) = 1;

    return tranDH;
}
Eigen::Isometry3d kinematics_6axis_ur_two::baseTend(const Vector6 &q)
{
    Eigen::Isometry3d transTCPInBase = Transform3D_DH(_a(0), _alpha(0), _d(0), _theta(0));
    for(int ii = 1; ii < 7; ii++)
        transTCPInBase = transTCPInBase * Transform3D_DH(_a(ii), _alpha(ii), _d(ii), _theta(ii) + q(ii-1));

    return transTCPInBase;
}



int kinematics_6axis_ur_two::solve(const Eigen::Isometry3d& transform3D, std::vector<Vector6> &joints, double eps){
    _eps =eps;
    Eigen::Matrix<double, 8, 6> solutions, solutionsPlus, solutionsMinus;
    r=transform3D.rotation();
    p=transform3D.translation();

    double q0_p=solveQ0(true);
    double q0_n=solveQ0(false);

    double q4_pp=solveQ4(q0_p,true);
    double q4_pn=solveQ4(q0_p,false);
    double q4_np=solveQ4(q0_n,true);
    double q4_nn=solveQ4(q0_n,false);

    if(fabs(sin(q4_pp))<eps){
        std::cout << "Robot is in singular state q4_pp, joint0 = "
                  << q0_p/M_PI*180.0
                  << " or " << q0_n/M_PI*180.0
                  << ",joint4 = 0, 关节2,3,4,6平行"<<"\n";
        return -1;
    }
    if(fabs(sin(q4_pn))<eps){
        std::cout << "Robot is in singular state q4_pn, joint0 = "
                  << q0_p/M_PI*180.0
                  << " or " << q0_n/M_PI*180.0
                  << ",joint4 = 0, 关节2,3,4,6平行"<<"\n";
        return -1;
    }
    if(fabs(sin(q4_np))<eps){
        std::cout << "Robot is in singular state q4_np, joint0 = "
                  << q0_p/M_PI*180.0
                  << " or " << q0_n/M_PI*180.0
                  << ",joint4 = 0, 关节2,3,4,6平行"<<"\n";
        return -1;
    }
    if(fabs(sin(q4_nn))<eps){
        std::cout << "Robot is in singular state q4_nn, joint0 = "
                  << q0_p/M_PI*180.0
                  << " or " << q0_n/M_PI*180.0
                  << ",joint4 = 0, 关节2,3,4,6平行"<<"\n";
        return -1;
    }

    double q5_pp=solveQ5(q0_p,q4_pp);
    double q5_pn=solveQ5(q0_p,q4_pn);
    double q5_np=solveQ5(q0_n,q4_np);
    double q5_nn=solveQ5(q0_n,q4_nn);

    double q123_pp=solveQ123(q0_p,q4_pp);
    double q123_pn=solveQ123(q0_p,q4_pn);
    double q123_np=solveQ123(q0_n,q4_np);
    double q123_nn=solveQ123(q0_n,q4_nn);

    double q2_ppp=solveQ2(q0_p,q4_pp,q123_pp, true);
    double q2_ppn=solveQ2(q0_p,q4_pp,q123_pp, false);
    double q2_pnp=solveQ2(q0_p,q4_pn,q123_pn, true);
    double q2_pnn=solveQ2(q0_p,q4_pn,q123_pn, false);
    double q2_npp=solveQ2(q0_n,q4_np,q123_np, true);
    double q2_npn=solveQ2(q0_n,q4_np,q123_np, false);
    double q2_nnp=solveQ2(q0_n,q4_nn,q123_nn, true);
    double q2_nnn=solveQ2(q0_n,q4_nn,q123_nn, false);

    double q1_ppp=solveQ1(q0_p,q4_pp,q123_pp,q2_ppp);
    double q1_ppn=solveQ1(q0_p,q4_pp,q123_pp,q2_ppn);
    double q1_pnp=solveQ1(q0_p,q4_pn,q123_pn,q2_pnp);
    double q1_pnn=solveQ1(q0_p,q4_pn,q123_pn,q2_pnn);
    double q1_npp=solveQ1(q0_n,q4_np,q123_np,q2_npp);
    double q1_npn=solveQ1(q0_n,q4_np,q123_np,q2_npn);
    double q1_nnp=solveQ1(q0_n,q4_nn,q123_nn,q2_nnp);
    double q1_nnn=solveQ1(q0_n,q4_nn,q123_nn,q2_nnn);

    double q3_ppp=solveQ3(q123_pp,q2_ppp,q1_ppp);
    double q3_ppn=solveQ3(q123_pp,q2_ppn,q1_ppn);
    double q3_pnp=solveQ3(q123_pn,q2_pnp,q1_pnp);
    double q3_pnn=solveQ3(q123_pn,q2_pnn,q1_pnn);
    double q3_npp=solveQ3(q123_np,q2_npp,q1_npp);
    double q3_npn=solveQ3(q123_np,q2_npn,q1_npn);
    double q3_nnp=solveQ3(q123_nn,q2_nnp,q1_nnp);
    double q3_nnn=solveQ3(q123_nn,q2_nnn,q1_nnn);

    solutions <<
        q0_p,q1_ppp,q2_ppp,q3_ppp,q4_pp,q5_pp,
        q0_p,q1_ppn,q2_ppn,q3_ppn,q4_pp,q5_pp,
        q0_p,q1_pnp,q2_pnp,q3_pnp,q4_pn,q5_pn,
        q0_p,q1_pnn,q2_pnn,q3_pnn,q4_pn,q5_pn,
        q0_n,q1_npp,q2_npp,q3_npp,q4_np,q5_np,
        q0_n,q1_npn,q2_npn,q3_npn,q4_np,q5_np,
        q0_n,q1_nnp,q2_nnp,q3_nnp,q4_nn,q5_nn,
        q0_n,q1_nnn,q2_nnn,q3_nnn,q4_nn,q5_nn;
    solutionsPlus <<
        q0_p,q1_ppp,q2_ppp,q3_ppp,q4_pp,q5_pp + 2*M_PI,
        q0_p,q1_ppn,q2_ppn,q3_ppn,q4_pp,q5_pp + 2*M_PI,
        q0_p,q1_pnp,q2_pnp,q3_pnp,q4_pn,q5_pn + 2*M_PI,
        q0_p,q1_pnn,q2_pnn,q3_pnn,q4_pn,q5_pn + 2*M_PI,
        q0_n,q1_npp,q2_npp,q3_npp,q4_np,q5_np + 2*M_PI,
        q0_n,q1_npn,q2_npn,q3_npn,q4_np,q5_np + 2*M_PI,
        q0_n,q1_nnp,q2_nnp,q3_nnp,q4_nn,q5_nn + 2*M_PI,
        q0_n,q1_nnn,q2_nnn,q3_nnn,q4_nn,q5_nn + 2*M_PI;
    solutionsMinus <<
        q0_p,q1_ppp,q2_ppp,q3_ppp,q4_pp,q5_pp - 2*M_PI,
        q0_p,q1_ppn,q2_ppn,q3_ppn,q4_pp,q5_pp - 2*M_PI,
        q0_p,q1_pnp,q2_pnp,q3_pnp,q4_pn,q5_pn - 2*M_PI,
        q0_p,q1_pnn,q2_pnn,q3_pnn,q4_pn,q5_pn - 2*M_PI,
        q0_n,q1_npp,q2_npp,q3_npp,q4_np,q5_np - 2*M_PI,
        q0_n,q1_npn,q2_npn,q3_npn,q4_np,q5_np - 2*M_PI,
        q0_n,q1_nnp,q2_nnp,q3_nnp,q4_nn,q5_nn - 2*M_PI,
        q0_n,q1_nnn,q2_nnn,q3_nnn,q4_nn,q5_nn - 2*M_PI;

    //LOG(INFO) <<"\nSolutions:\n\n"<<solutions/M_PI*180.0;

     for(int ii=0;ii<8;ii++){
         if(!std::isnan(solutions(ii,0)) &&
            !std::isnan(solutions(ii,1)) &&
            !std::isnan(solutions(ii,2)) &&
            !std::isnan(solutions(ii,3)) &&
            !std::isnan(solutions(ii,4)) &&
            !std::isnan(solutions(ii,5)) )
         {
            joints.push_back(solutions.block(ii, 0, 1, 6));
            joints.push_back(solutionsPlus.block(ii, 0, 1, 6));
            joints.push_back(solutionsMinus.block(ii, 0, 1, 6));
             //joints.push_back(solutions.block(ii,0,1,6));
//             joints.push_back(solutionsPlus.block(ii,0,1,6));
//             joints.push_back(solutionsMinus.block(ii,0,1,6));
//             auto q = solutions.block(ii,0,1,6);
//             LOG_INFO << "Index "<<ii <<": "<<q/M_PI*180.0<<"\n";
//             baseTend(q);
         }
    }

    //Solutions are empty.
    if(joints.empty()){
        return -2;
    }

    return 0;
}

double kinematics_6axis_ur_two::solveQ0(bool positive){
    double A=r(0,2)*d7-p(0);
    double B=p(1)-r(1,2)*d7;
    double C=d5;
    double q;
    double sq = A*A+B*B-C*C;
    if(fabs(sq) < _eps){
        sq = 0.0;
    }
    if(positive){
        q=atan2(C,sqrt(sq))-atan2(A,B);
    }else{
        q=atan2(C,-sqrt(sq))-atan2(A,B);
    }
    if(q>M_PI){
        q-=M_PI*2;
    }
    if(q<-M_PI){
        q+=M_PI*2;
    }
    return q;
}

double kinematics_6axis_ur_two::solveQ4(double q0,bool positive){
    double b6=-cos(q0)*r(0,2)+sin(q0)*r(1,2);
    double q;
    double sq = 1-b6*b6;
    if(fabs(sq) < _eps){
        sq = 0.0;
    }
    if(positive){
        q=atan2(sqrt(sq),b6);
    }else{
        q=atan2(-sqrt(sq),b6);
    }
    if(q>M_PI){
        q-=M_PI*2;
    }
    if(q<-M_PI){
        q+=M_PI*2;
    }
    return q;
}

double kinematics_6axis_ur_two::solveQ5(double q0,double q4){
    double A6=-(cos(q0)*r(0,1)-sin(q0)*r(1,1))/sin(q4);
    double B6=(cos(q0)*r(0,0)-sin(q0)*r(1,0))/sin(q4);
    double q=atan2(A6,B6);
    if(q>M_PI){
        q-=M_PI*2;
    }
    if(q<-M_PI){
        q+=M_PI*2;
    }
    return q;
}

double kinematics_6axis_ur_two::solveQ123(double q0,double q4){
    double A345=-r(2,2)/sin(q4);
    double B345=(-sin(q0)*r(0,2)-cos(q0)*r(1,2))/sin(q4);
    double q=atan2(A345,B345);
    if(q>M_PI){
        q-=M_PI*2;
    }
    if(q<-M_PI){
        q+=M_PI*2;
    }
    return q;
}

double kinematics_6axis_ur_two::solveQ2(double q0,double q4,double q123,bool positive){
    double m = M(q0,q4,q123);
    double n = N(q0,q4,q123);
    double b=(m*m+n*n-a3*a3-a4*a4)/(2.0*a3*a4);
    double sq = 1 - b*b;
    if(fabs(sq) < _eps){
        sq = 0.0;
    }
    double q;
    if(positive){
        q=atan2(sqrt(sq),b);
    }else{
        q=atan2(-sqrt(sq),b);
    }
    if(q>M_PI){
        q-=M_PI*2;
    }
    if(q<-M_PI){
        q+=M_PI*2;
    }
    return q;
}

double kinematics_6axis_ur_two::solveQ1(double q0,double q4,double q123,double q2){
    double m = M(q0,q4,q123);
    double A = a4*sin(q2);
    double B = a4*cos(q2)+a3;
    double C = m;
    double sq = A*A+B*B-C*C;
    if(fabs(sq) < _eps){
        sq = 0.0;
    }
    //TODO 我不知道为什么这里不是同时存在两种情况都有解，而是只有一种情况是有解的。
    double qp=atan2(C,sqrt(sq))-atan2(A,B);
    double qn=atan2(C,-sqrt(sq))-atan2(A,B);
    double q;
    double pzp=a3*cos(qp)+a4*cos(q2+qp)+d2+d6*cos(q123)-d7*sin(q4)*sin(q123);
    double pzn=a3*cos(qn)+a4*cos(q2+qn)+d2+d6*cos(q123)-d7*sin(q4)*sin(q123);
    if(fabs(pzp-p(2))<_eps){
        q = qp;
    } else if(fabs(pzn-p(2))<_eps){
        q = qn;
    } else {
        q = NAN;
        return q;
    }

    if(q>M_PI){
        q-=M_PI*2;
    }
    if(q<-M_PI){
        q+=M_PI*2;
    }
    return q;
}

double kinematics_6axis_ur_two::solveQ3(double q123,double q2,double q1){
    double q = q1+q2-q123;
    if(q>M_PI){
        q-=M_PI*2;
    }
    if(q<-M_PI){
        q+=M_PI*2;
    }
    return q;
}

double kinematics_6axis_ur_two::M(double q0, double q4,double q123){
    double m=-p(0)*sin(q0)-p(1)*cos(q0)-d6*sin(q123)-d7*sin(q4)*cos(q123);
    return m;
}

double kinematics_6axis_ur_two::N(double q0, double q4,double q123){
    double n=-d6*cos(q123)+d7*sin(q4)*sin(q123)-d2+p(2);
    return n;
}

int kinematics_6axis_ur_two::getQ(const Eigen::Isometry3d& transform3D, const Vector6 &cur_Q, Vector6 &targetQ, double threshold)
{
    double norm = 200 * M_PI, _temnorm = 0;
    int choseNum = 0;

    //Get all solutions.
    std::vector<Vector6> solutions;
    if(solve(transform3D, solutions) < 0){
        std::cout << prefix << "NO IK solutions." << std::endl;
        return -1;
    }

    //Find bast solution.
    for(int i = 0, solNum = solutions.size(); i < solNum; ++i){
        if(norm > (_temnorm = norm2(solutions[i], cur_Q))){
            norm = _temnorm;
            choseNum = i;
        }
    }
    targetQ = solutions[choseNum];

    //Deal with the case when got a discontinuous solution.
    if(filter_joints(cur_Q, targetQ, threshold)){
        std::cout << prefix << "Get discontinuous solution(threshold " << threshold << ")." <<  std::endl;
        return -2;
    }

    return 0;
}

int kinematics_6axis_ur_two::filter_joints(const Vector6 &cQ, const Vector6 &tQ, double threshold)
{
    Vector6 offset = cQ - tQ;
    for(int i = 0; i < 6; ++i)
        if(fabs(offset[i]) >= threshold)
            return -1;

    return 0;
}
