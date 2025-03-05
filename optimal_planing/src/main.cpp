#include <toppra/algorithm/toppra.hpp>
#include <toppra/constraint/linear_joint_acceleration.hpp>
#include <toppra/constraint/linear_joint_velocity.hpp>
#include <toppra/geometric_path/piecewise_poly_path.hpp>
#include <toppra/parametrizer/const_accel.hpp>
#include <toppra/toppra.hpp>
#include "interface_toppra.h"
#include "kinematics_6axis_ur.h"
// #include "B_Spline.h"
// #include "trapezoidal_velocity.h"
#include "rotational_interpolation.h"
// #include "SphericalQuadraticInterpolation.h"
#include "MultiPointsTrajectory.h"

// #include <Eigen/Dense>
// #include <iostream>
// #include <cmath>
// #include <gl/glut.h>
// using namespace Eigen;
// using namespace std;

// void Display();
// void Caculate_R(int i);
// Vector3d Caculate_Curve(int i, double u);
// double delta_u(int i, double u_j);

// MatrixXd P(10, 3); //10个型值点，n+1=10
// MatrixXd B(12, 3); //12个控制点，10+2=12
// Vector3d R3, R2, R1, R0; //B样条曲线的参数矩阵，Pi = R0 + R1*u + R2*u*u + R3*u*u*u

// void Caculate_R(int i)
// {
// 	R3 = (-B.row(i) + 3 * B.row(i + 1) - 3 * B.row(i + 2) + B.row(i + 3)) / 6;
// 	R2 = (3 * B.row(i) - 6 * B.row(i + 1) + 3 * B.row(i + 2)) / 6;
// 	R1 = (-3 * B.row(i) + 3 * B.row(i + 2)) / 6;
// 	R0 = (B.row(i) + 4 * B.row(i + 1) + B.row(i + 2)) / 6;
// }
// Vector3d Caculate_Curve(int i,double u)
// {
// 	Caculate_R(i);
// 	return (R0 + R1*u + R2*u*u + R3*u*u*u);
// }
// double delta_u(int i, double u_j)
// {
// 	Caculate_R(i);
// 	double a4, a3, a2, a1, a0;
// 	a4 = 9 * R3.transpose()*R3;
// 	a3 = 12 * R3.transpose()*R2;
// 	a2 = 4 * R2.transpose()*R2;
// 	a2= a2 + (6 * R3.transpose()*R1);
// 	a1 = 4 * R2.transpose()*R1;
// 	a0 = R1.transpose()*R1;
// 	double delta;
// 	delta = sqrt(a4*pow(u_j,4)+a3*pow(u_j,3)+a2*pow(u_j,2)+a1*u_j+a0);
// 	return delta;
// }

// void Display(void)//OpenGL来绘制曲线
// {
// 	glClear(GL_COLOR_BUFFER_BIT);
// 	glEnableClientState(GL_VERTEX_ARRAY);
// 	float g_CP[24];
// 	int k = 0;
// 	for (int i = 0; i < 12; i++){
// 		for (int j = 0; j < 2; j++){
// 			g_CP[k] = B(i, j);
// 			cout << g_CP[k] << endl;
// 			k++;
// 		}
// 	}
// 	glPointSize(10);
// 	glColor3d(255, 0, 0);
// 	glBegin(GL_POINTS);
// 	k = 0;
// 	while (k < 22){
// 		glVertex2f(g_CP[k], g_CP[k + 1]);
// 		k = k + 2;
// 	}
// 	glEnd();

// 	glVertexPointer(2, GL_FLOAT, 0, g_CP);
// 	glPointSize(2);
// 	glColor3d(255, 0, 0);
// 	glDrawArrays(GL_LINE_STRIP, 0, 12); //不闭合折线绘制

// 	float u;
// 	Vector3d point_3d;
// 	GLfloat point[2];
// 	glColor3d(0, 255, 0);
// 	glBegin(GL_LINE_STRIP);
// 	glPointSize(4);
// 	for (int i = 0; i < 9; i++) {
// 		u = 0;
// 		Caculate_R(i);
// 		for (int j = 0; j < 20; j++){
// 			//u = u + delta_u(i, u);
// 			point_3d = Caculate_Curve(i,u);

// 			u = u + 0.05;
// 			point[0] = point_3d(0), point[1] = point_3d(1);
// 			glVertex2f(point[0], point[1]);
// 		}
		
// 	}
// 	glEnd();
// 	glPointSize(10);
// 	glBegin(GL_POINTS);
// 	for (int i = 0; i<10; i++){
// 		glVertex2f(P(i, 0), P(i, 1));
// 	}
// 	glEnd();
// 	glFlush();
// }

Eigen::Vector3d getPosEulerFromXYZ(const Eigen::Matrix3d &R, double Rref[3])
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

    // if (fabs(pr[0] - Rref[0]) + fabs(pr[1] - Rref[1]) + fabs(pr[2] - Rref[2]) > 1.57)
    // {
    //     cout << "xyz1=" << x[0] * 180 / M_PI << "," << y[0] * 180 / M_PI << "," << z[0] * 180 / M_PI << endl;
    //     cout << "xyz2=" << x[1] * 180 / M_PI << "," << y[1] * 180 / M_PI << "," << z[1] * 180 / M_PI << endl;
    //     cout << "pr=" << pr[0] * 180 / M_PI << "," << pr[1] * 180 / M_PI << "," << pr[2] * 180 / M_PI << endl;
    //     cout << "Rref=" << Rref[0] * 180 / M_PI << "," << Rref[1] * 180 / M_PI << "," << Rref[2] * 180 / M_PI << endl;
    // }
    return pr;
}

Eigen::Matrix3d QuatToRotm(Eigen::Quaterniond &IMUquat){
    return IMUquat.toRotationMatrix();
}

bool IsIMUCollision(Eigen::Matrix<double, 1,7> &ref, Eigen::Matrix<double, 1,7> &current, double threshold) {

        /* Input
        ref (w,x,y,z,RZ,RY,RX) degree
        cur (w,x,y,z,RZ,RY,RX) degree
        threshold degree */
        double Rref[3];
        Rref[0] = 0;
        Rref[1] = 0;
        Rref[2] = 0;
        Eigen::Matrix3d NTRom, NRRom, TR, RR, deltaTR, deltaRR;
        Eigen::Vector3d bNTR, NTR, NRR;
        Eigen::Matrix<double, 3,3> Pin;
        Eigen::Matrix<double, 3,1> CTR;
        Pin <<    0.985519,   0.120665, -0.0871064,
                 -0.188759,   0.970245, -0.0603612,
                  0.153308,  0.0639084,   0.992154;

        // Reference IMU
        double w = ref(0,0);
        double x = ref(0,1);
        double y = ref(0,2);
        double z = ref(0,3);
        Eigen::Quaterniond refIMU(w,x,y,z);
        Eigen::AngleAxisd Rz(M_PI/2, Eigen::Vector3d::UnitZ()); 
        Eigen::AngleAxisd Ry(0, Eigen::Vector3d::UnitY()); 
        Eigen::AngleAxisd Rx(M_PI/2, Eigen::Vector3d::UnitX()); 
        TR = QuatToRotm(refIMU) * Rz.toRotationMatrix()*Ry.toRotationMatrix()*Rx.toRotationMatrix();
        std::cout << "TR" << std::endl;
        std::cout << TR << std::endl;
        // Reference Robot
        double rz = ref(0,4) /180 *M_PI;
        double ry = ref(0,5) /180 *M_PI;
        double rx = ref(0,6) /180 *M_PI;
        Eigen::AngleAxisd roll(Eigen::AngleAxisd(rx,Eigen::Vector3d::UnitX()));
        Eigen::AngleAxisd pitch(Eigen::AngleAxisd(ry,Eigen::Vector3d::UnitY()));
        Eigen::AngleAxisd yaw(Eigen::AngleAxisd(rz,Eigen::Vector3d::UnitZ()));
        RR = yaw*pitch*roll;
        std::cout << "RR" << std::endl;
        std::cout << RR << std::endl;
		std::cout << rz << " " << ry << " " << rx << std::endl;
        //Current IMU
        w = current(0,0);
        x = current(0,1);
        y = current(0,2);
        z = current(0,3);
        Eigen::Quaterniond curIMU(w,x,y,z);
        // NTRom = QuatToRotm(curIMU) *Rz.toRotationMatrix()*Ry.toRotationMatrix()*Rx.toRotationMatrix(); 
        NTRom = QuatToRotm(curIMU);
        deltaTR = TR.transpose() * NTRom;
        std::cout << "deltaTR" << std::endl;
        std::cout << deltaTR << std::endl;
        std::cout << "NTRom" << std::endl;
        std::cout << NTRom << std::endl;
        NTR = getPosEulerFromXYZ(deltaTR, Rref);
        CTR << NTR[0], NTR[1], NTR[2];
        bNTR = Pin * CTR;
        std::cout << NTR[0]/M_PI*180 << " " << NTR[1]/M_PI*180 << " " << NTR[2]/M_PI*180 << std::endl;
        std::cout << bNTR[0]/M_PI*180 << " " << bNTR[1]/M_PI*180 << " " << bNTR[2]/M_PI*180 << std::endl;
        //Current Robot
        rz = current(0,4) /180 *M_PI;
        ry = current(0,5) /180 *M_PI;
        rx = current(0,6) /180 *M_PI;
        Eigen::AngleAxisd Croll(Eigen::AngleAxisd(rx,Eigen::Vector3d::UnitX()));
        Eigen::AngleAxisd Cpitch(Eigen::AngleAxisd(ry,Eigen::Vector3d::UnitY()));
        Eigen::AngleAxisd Cyaw(Eigen::AngleAxisd(rz,Eigen::Vector3d::UnitZ()));
        NRRom = Cyaw*Cpitch*Croll;        
        deltaRR = RR.transpose() * NRRom;
        std::cout << "deltaRR" << std::endl;
		std::cout << deltaRR << std::endl; 

        NRR = getPosEulerFromXYZ(deltaRR, Rref);
        std::cout << NRR[0]/M_PI*180 << " " << NRR[1]/M_PI*180 << " " << NRR[2]/M_PI*180 << std::endl; 
        for(int i=0;i<3;i++) {
            if(fabs(bNTR[i] - NRR[i])/M_PI*180 > threshold) {
                return true;
            }
        }
        return false;
}

int main() {
    // const int joint_num = 6;
    // toppra::Vectors positions, velocities;
    // toppra::Vector velLimitLower(joint_num), velLimitUpper(joint_num), accLimitLower(joint_num), accLimitUpper(joint_num);
    // toppra::Vector position0(joint_num), position1(joint_num), position2(joint_num), position3(joint_num), position4(joint_num);
    // toppra::Vector velocitie0(joint_num), velocitie1(joint_num), velocitie2(joint_num), velocitie3(joint_num), velocitie4(joint_num);
    // double period = 0.005;
    // bool printInfo = false;
    // Eigen::MatrixXd path_pos2_;
    // position0 << 10.471,  2.502,  113.411,  -50.126,  -8.889, -6.894;
    // position1 <<  9.890, -11.586, 127.216,  -59.438,  -36.698, -1.409;
    // position2 << -52.524, -0.068, 114.353,  -59.520,  -31.463,  3.288;
    // position3 << -52.525, 14.059, 120.009,  -46.632,  -31.479,  3.324;

    // positions = {position0/180*M_PI, position1/180*M_PI, position2/180*M_PI, position3/180*M_PI, position4/180*M_PI};

    // velocitie0 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    // velocitie1 << 0.00512708, -0.12022767, -0.18698188,  0.07181939, -0.09782292,  0.07245699;
    // velocitie2 << 0.09595612, -0.1151688,   0.14235637, -0.12214979,  0.01251548, -0.01386915;
    // velocitie3 <<  0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    
    // velocities = {velocitie0, velocitie1, velocitie2, velocitie3};


    // velLimitLower << -0.2, -0.2, -0.2, -0.1, -0.1, -0.1;
    // velLimitUpper <<  0.2,  0.2,  0.2,  0.1,  0.1,  0.1;
    // accLimitLower << -0.3, -0.3, -0.3, -0.1, -0.1, -0.1;
    // accLimitUpper <<  0.3,  0.3,  0.3,  0.1,  0.1,  0.1;

    // toppra::Planning optimizer; 
    // // optimizer.JointSpacePlanning(positions, velocities, 
    // //                           velLimitLower, velLimitUpper,
    // //                           accLimitLower, accLimitUpper,
    // //                           period, printInfo, path_pos2_);
    // // std::cout << path_pos2_ << std::endl;
    // Eigen::Matrix<double, 1, 6> cur_Q, targetQ;
    
    // DH_6axis_ur dh;
    // kinematics_6axis_ur_two aubo_i3 = kinematics_6axis_ur_two(dh);
    // double threshold = 6.28;
    // toppra::Vector currrentQ(joint_num), currrentQ2(joint_num), currrentQ3(joint_num), currrentQ4(joint_num);
    // currrentQ << -35.027,  75.625,  66.433, 142.257, 54.974, 0.001;
    // currrentQ2 <<   5.376,  74.235,  66.423, 142.259, 54.974, 0.001;
    // currrentQ3 <<   5.376,  65.282,  66.410, 131.430, 54.975, 0.001;
    // currrentQ4 << -17.249,  53.112,  66.425, 122.126, 68.203, 0.001;
    // currrentQ = currrentQ/180*M_PI;
    // currrentQ2 = currrentQ2/180*M_PI;
    // currrentQ3 = currrentQ3/180*M_PI;
    // currrentQ4 = currrentQ4/180*M_PI;
    // cur_Q << currrentQ2(0,0), currrentQ2(1,0), currrentQ2(2,0), currrentQ2(3,0), currrentQ2(4,0), currrentQ2(5,0);
    // Eigen::Isometry3d Obj;
    // Obj = aubo_i3.baseTend(currrentQ2);
    // std::cout << Obj.rotation() << std::endl;
    // std::cout << Obj.translation() << std::endl;
    // aubo_i3.getQ(Obj, cur_Q, targetQ, threshold);
    // std::cout << "targetQ: "<<targetQ*180/M_PI << std::endl;    


    // position0 << 140.765,-502.160, 123.521,  0.001/180*M_PI, -0.001/180*M_PI, 90.000/180*M_PI;
    // position1 <<-218.351,-474.313, 135.510,-40.391/180*M_PI, -0.920/180*M_PI, 88.689/180*M_PI;
    // position2 <<-220.025,-492.161, 199.772,-40.401/180*M_PI,  0.149/180*M_PI, 90.215/180*M_PI;
    // position3 <<   6.761,-540.792, 296.533, -4.528/180*M_PI, -0.962/180*M_PI, 87.596/180*M_PI;

    // positions = {position0, position1, position2, position3};
    // std::cout << position0 << std::endl;
    // velocitie0 << 0.0, 0.0, 0.0, 0.0,  0.0, 0.0;
    // velocitie1 << 80, 30, 30, 0.03, 0.03, 0.1;
    // velocitie2 << 30, 30, 80, 0.03, 0.03, 0.1;
    // velocitie3 <<  0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    
    // velocities = {velocitie0, velocitie1, velocitie2, velocitie3};


    // velLimitLower << -300, -300, -300, -0.2, -0.2, -0.2;
    // velLimitUpper <<  300,  300,  300,  0.2,  0.2,  0.2;
    // accLimitLower << -80, -80, -80, -0.2, -0.2, -0.2;
    // accLimitUpper <<  80,  80,  80,  0.2,  0.2,  0.2;
    // std::cout << accLimitUpper << std::endl;
    // optimizer.CartesianSpacePlanning(currrentQ, positions, velocities, 
    //                           velLimitLower, velLimitUpper,
    //                           accLimitLower, accLimitUpper,
    //                           period, printInfo, path_pos2_);
    // std::cout << path_pos2_ << std::endl;


    // 	P << -0.5, -0.5, 0.0,
	// 	-0.8, -0.1, 0.0,
	// 	-0.4, 0.5, 0.0,
	// 	-0.1, 0.6, 0.0,
	// 	0.3, 0.4, 0.0,
	// 	0.5, 0.25, 0.0,
	// 	0.7, 0.0, 0.0,
	// 	0.8, -0.2, 0.0,
	// 	0.3, -0.4, 0.0,
	// 	0.4, -0.5, 0.0; 
	// MatrixXd A(12, 12);
	// A = MatrixXd::Zero(12, 12);
	// cout << A << endl;
	// for (int i = 0; i < 12; i++){
	// 	if (i == 0){
	// 		A(0, 0) = -1;
	// 		A(0, 1) = 1;
	// 	}
	// 	else if (i == 11){
	// 		A(11, 10) = 1;
	// 		A(11, 11) = -1;
	// 	}
	// 	else{
	// 		A(i, i - 1) = 1;
	// 		A(i, i) = 4;
	// 		A(i, i + 1) = 1;
	// 	}
	// }
	// cout << A << endl;
	// MatrixXd PP(12, 3);
	// PP = MatrixXd::Zero(12, 3);
	// for (int i = 1; i < 11; i++){
	// 	for (int j = 0; j < 2; j++){
	// 		PP(i, j) = 6*P(i - 1, j);
	// 	}
	// }
	// cout << PP << endl;
	// B = A.inverse()*PP;
	// cout << "控制点" << endl;
	// cout << B << endl;

	// glutInit(&argc, argv);
	// glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	// glutInitWindowPosition(400, 200);
	// glutInitWindowSize(600, 600);
	// glutCreateWindow("B样条曲线");
	// glutDisplayFunc(&Display);
// 	glutMainLoop()
//     double p = -38.26149368286133 / 180 * M_PI;
//     double r = -3.4206910133361816/ 180 * M_PI;
//     double y = 0.1475760042667389 / 180 * M_PI;
// -38.26149368286133, -3.4206910133361816, 0.1475760042667389,
// -33.572227478027344, -2.2242000102996826, 346.8136901855469
//     Eigen::Quaterniond qq;
//         Eigen::AngleAxisd roll(Eigen::AngleAxisd(r,Eigen::Vector3d::UnitX()));
//         Eigen::AngleAxisd pitch(Eigen::AngleAxisd(p,Eigen::Vector3d::UnitY()));
//         Eigen::AngleAxisd yaw(Eigen::AngleAxisd(y,Eigen::Vector3d::UnitZ()));
//         qq = yaw*pitch*roll;
//     std::cout << "qq" << qq << std::endl;


// Eigen::AngleAxisd roll(Eigen::AngleAxisd(r, Eigen::Vector3d::UnitX()));
// Eigen::AngleAxisd pitch(Eigen::AngleAxisd(p, Eigen::Vector3d::UnitY()));
// Eigen::AngleAxisd yaw(Eigen::AngleAxisd(y, Eigen::Vector3d::UnitZ()));
// Eigen::Quaterniond q = yaw * pitch * roll;

// 	Eigen::Matrix<double, 1,7> ref;
//     Eigen::Matrix<double, 1,7> current;
//     ref << 0.944364, -0.0277765, -0.327615, -0.00856533, 99.46201721134291, -1.136849525403726, 123.87779164364814;
// current << -0.950235, 0.051612, 0.284702, 0.115473, 112.69802497437217, -1.7031973099832007, 119.00630143602517;
//    double threshold = 11;
//     bool isCol = IsIMUCollision(ref, current, threshold);
//     std::cout << "isCol" << isCol << std::endl;


	Eigen::MatrixXd points(6, 6);
    std::vector<Eigen::Quaterniond> quats;
	//points << 1, 1, 1, 2, 3, 2, 4, 5, 5, 2, 3, 3, 5, 4, 3, 6, 7, 1, 9, 9, 8, 12, 15, 11;
// 	points <<1.62338281300000,	0.704965027000000,	-0.0125106200000000,
// 1.28713110400000,	0.689627258000000,	0.251928284000000,
// 1.20999658200000,	0.694901245000000,	0.217869019000000,
// 1.06577929700000,	0.700669373000000,	0.193214111000000,
// 0.983905090000000,	0.704640991000000,	0.111798584000000,
// 0.892094543000000,	0.707915222000000,	0.0112802730000000,
// 0.835077087000000,	0.709012573000000,	0.00817981000000000,
// 0.796471741000000,	0.703724487000000,	0.123177490000000,
// 0.751769348000000,	0.698676758000000,	0.215058594000000,
// 0.713332092000000,	0.692176392000000,	0.273975037000000,
// 0.669832581000000,	0.690623779000000,	0.304385864000000,
// 0.601487244000000,	0.690779785000000,	0.324083374000000,
// 0.556960999000000,	0.690075928000000,	0.296066772000000,
// 0.475868347000000,	0.691352539000000,	0.270943787000000,
// 0.418927307000000,	0.695278809000000,	0.198506470000000,
// 0.344341064000000,	0.699948914000000,	0.117120605000000,
// 0.300111877000000,	0.702718506000000,	0.0414648440000000,
// 0.253929565000000,	0.688662476000000,	0.157157593000000,
// 0.200157654000000,	0.677834961000000,	0.251385193000000,
// 0.161352356000000,	0.669927490000000,	0.294979492000000,
// 0.106556396000000,	0.661460449000000,	0.327523560000000,
// 0.0496652830000000,	0.654651489000000,	0.331231323000000,
// -0.00696710200000000,	0.651501038000000,	0.307992493000000;
// std::cout << points << std::endl;

    points << -0.004033,0.3055,0.3909, 1.1784,1.5462,-0.83048,
-0.002464,0.30846,0.28987, 1.2859,1.5521,-0.72293,
-0.084906,0.47651,0.29636, 1.2525,1.5337,-0.76819,
-0.020277,0.50181,0.4644, 1.2534,1.5244,-0.76136,
0.030274,0.39687,0.46005, 1.2622,1.5269,-0.74918,
-0.006156,0.30654,0.59421, 1.2058,1.5164,-0.80915;

	// double period = 0.05; 
	double Vx = 0.3, allow_omega = 0.2, Ax = 0.05;
	double threshold = 0.3;
    // double scale = 1.0;
    MultiPoints mp;
    std::vector<Eigen::Isometry3d> output;
    output = mp.GenerateTrajectory(points, Vx, Ax, threshold);
    for(int i=0;i<output.size();i++) {
        std::cout << i << ": " << std::endl;
        std::cout << output[i].matrix() << std::endl;
    }
//  Eigen::Quaterniond R0;
 Squad RI;
//  for(int i=0;i<6;i++) {
//     Eigen::AngleAxisd roll(Eigen::AngleAxisd(points(i,5),Eigen::Vector3d::UnitX()));
//     Eigen::AngleAxisd pitch(Eigen::AngleAxisd(points(i,4),Eigen::Vector3d::UnitY()));
//     Eigen::AngleAxisd yaw(Eigen::AngleAxisd(points(i,3),Eigen::Vector3d::UnitZ()));
//     R0 = yaw*pitch*roll;
//     // std::cout << "R0" << std::endl;
//     // std::cout << R0 << std::endl;
//     quats.push_back(R0);
//  }



//  RI.SetStartEnd();
//  double angle = RI.Angle();


// std::cout << points << std::endl;
	BSplines curve;
	TrapezoidalVelocity Tvel;
	// // std::cout << points << std::endl;
	// curve.UniformCubicBSplines(points);
	// // Eigen::MatrixXd C = curve.GetCurvePos();
	// // std::cout << C << std::endl;

    // // double L =  curve.GetTotalL();
    // Eigen::MatrixXd LL = curve.GetDeltaL();
	// Eigen::MatrixXd kappa = curve.GetCurvature();
    // int seg;


    // for(int f=0;f<10;f++) {
    //     std::vector<Eigen::Quaterniond> dQuats;
    //     Eigen::MatrixXd uus = Tvel.TrapezoidalVelPlan(LL, kappa, period, Vx, Ax, threshold, scale);
    //     // Eigen::MatrixXd T = Tvel.GetTotalT();
    //     // Eigen::MatrixXd S = Tvel.GetTotalS();
    //     // Eigen::MatrixXd V = Tvel.GetTotalV();
    //     // Eigen::MatrixXd A = Tvel.GetTotalA();
    //     curve.UniformCubicBSplines(points, uus);
    //     Eigen::MatrixXd C = curve.GetCurvePos();
    //     // std::cout << T << std::endl;

    //     int seg = uus.rows();
    //     std::cout << seg << std::endl;

    //     dQuats.push_back(quats[0]);
    //     for(int i=1;i<seg-1;i++) {
    //         double s = i / (seg-1.0);
    //         R0 = RI.QuadSquad(quats, s);
    //         dQuats.push_back(R0);
    //     }
    //     dQuats.push_back(quats[5]);
    //     Eigen::MatrixXd dt(seg, 1);
        
    //     for(int i=0;i<seg;i++) {
    //         std::cout << i+1 << ": " << dQuats[i] << std::endl;
    //         dt(i,0) = period;
    //     }
    //     double maxD = RI.ComputeAngularVelocityAndAcceleration(dQuats, dt);
    //     std::cout << maxD << std::endl;
    //     scale = allow_omega / maxD;
    //     std::cout << scale << std::endl;
    //     if(scale >= 1) {   
    //         break;
    //     }
    // }


    // DH_6axis_ur dh;
    // kinematics_6axis_ur_two aubo_i3 = kinematics_6axis_ur_two(dh);
    // double threshold = 6.28;
    // Eigen::Isometry3d Obj;
    
    // for(int i=0;i<C.rows();i++) {
    //     obj = Eigen::Isometry3d::Identity();
    // }
    // Eigen::Isometry3d Obj;
    // Obj = aubo_i3.baseTend(currrentQ2);
    // std::cout << Obj.rotation() << std::endl;
    // std::cout << Obj.translation() << std::endl;
    // aubo_i3.getQ(Obj, cur_Q, targetQ, threshold);
    // std::cout << "targetQ: "<<targetQ*180/M_PI << std::endl;    
    return 0;

}
