#include "MultiPointsTrajectory.h"

MultiPoints::MultiPoints() : scale(1.0), period(0.005), Vx(1.0), Ax(0.1), threshold(0.3), allow_omega(0.1) {}

MultiPoints::~MultiPoints() {}

std::vector<Eigen::Isometry3d> MultiPoints::GenerateTrajectory(Eigen::MatrixXd &points, double VxPrm, double AxPrm, double thresholdPrm, const std::string &velType) {
    std::vector<Eigen::Isometry3d> diecretePoints;
    Vx = VxPrm;
    Ax = AxPrm;
    threshold = thresholdPrm;
    if(velType == "T") {
        diecretePoints = TTypeVelPlan(points);
    } else if (velType == "S") {
        diecretePoints = STypeVelPlan(points);
    }
    return diecretePoints;
}

std::vector<Eigen::Isometry3d> MultiPoints::TTypeVelPlan(Eigen::MatrixXd &points) {
    std::vector<Eigen::Isometry3d> output;
	std::vector<Eigen::Quaterniond> quats = EulerToQuats(points);
    int rows = points.rows();
    Eigen::MatrixXd pos = points.block(0,0,rows,3);
    curve.UniformCubicBSplines(pos);

    Eigen::MatrixXd LL = curve.GetDeltaL();
	Eigen::MatrixXd kappa = curve.GetCurvature();
    std::vector<Eigen::Quaterniond> finalQ;
    for(int f=0;f<10;f++) {
        std::vector<Eigen::Quaterniond> dQuats;
        Eigen::MatrixXd uus = Tvel.TrapezoidalVelPlan(LL, kappa, period, Vx, Ax, threshold, scale);

        curve.UniformCubicBSplines(pos, uus);
        Eigen::MatrixXd C = curve.GetCurvePos();

        seg = uus.rows();

        dQuats.push_back(quats[0]);
        for(int i=1;i<seg-1;i++) {
            double s = i / (seg-1.0);
            dQuats.push_back(RI.QuadSquad(quats, s));
        }
        dQuats.push_back(quats[5]);
        Eigen::MatrixXd dt(seg, 1);
        
        for(int i=0;i<seg;i++) {
            std::cout << i+1 << ": " << dQuats[i] << std::endl;
            dt(i,0) = period;
        }
        double maxD = RI.ComputeAngularVelocityAndAcceleration(dQuats, dt);

        scale = allow_omega / maxD;
        if(scale >= 1) {   
            finalQ = dQuats;
            break;
        }
    }    

    Eigen::MatrixXd C = curve.GetCurvePos();
    for(int i=0;i<seg;i++) {
        Eigen::Isometry3d single = Eigen::Isometry3d::Identity();
        Eigen::Vector3d vec = {C(i,0), C(i,1), C(i,2)};
        single.translate(vec);
        single.rotate(finalQ[i]);
        output.push_back(single);
    }
    return output;
}

std::vector<Eigen::Isometry3d> MultiPoints::STypeVelPlan(Eigen::MatrixXd &points) {
    std::vector<Eigen::Isometry3d> output;
    return output;
}

std::vector<Eigen::Quaterniond> MultiPoints::EulerToQuats(Eigen::MatrixXd &points) {
    std::vector<Eigen::Quaterniond> quats;
    int rows = points.rows();
    Eigen::Quaterniond R0;
    for(int i=0;i<rows;i++) {
        Eigen::AngleAxisd roll(Eigen::AngleAxisd(points(i,5),Eigen::Vector3d::UnitX()));
        Eigen::AngleAxisd pitch(Eigen::AngleAxisd(points(i,4),Eigen::Vector3d::UnitY()));
        Eigen::AngleAxisd yaw(Eigen::AngleAxisd(points(i,3),Eigen::Vector3d::UnitZ()));
        R0 = yaw*pitch*roll;
        quats.push_back(R0);
    }
    return quats;
}
