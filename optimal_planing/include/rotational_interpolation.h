

#ifndef KDL_ROTATIONALINTERPOLATION_H
#define KDL_ROTATIONALINTERPOLATION_H

#include <Eigen/Eigen>
#include <vector>
//#include "TrajDefine.h"
using namespace Eigen;
using Rotation33 = Eigen::Matrix<double, 3, 3>;
class RotationalInterpolation
{
	Rotation33 R_base_start;
	Rotation33 R_base_end;

	Vector3d rot_start_end;
	double angle;

public:
	RotationalInterpolation();
	virtual void SetStartEnd(Rotation33 start, Rotation33 end);
	virtual double Angle();
	virtual Rotation33 Pos(double th) const;
	virtual Vector3d Vel(double th, double thd) const;
	virtual Vector3d Acc(double th, double thd, double thdd) const;
	virtual void Write(std::ostream &os) const;
	virtual RotationalInterpolation *Clone() const;
	virtual ~RotationalInterpolation();

	double GetRotAngle(Rotation33 data, Vector3d &axis, double eps = 1e-6) const;
	static Rotation33 Rot2(const Vector3d &rotvec, double angle);
};

#endif
