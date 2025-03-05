

#include "rotational_interpolation.h"
#include <iostream>

#include <unistd.h>

RotationalInterpolation::RotationalInterpolation()
{
}

void RotationalInterpolation::SetStartEnd(Rotation33 start, Rotation33 end)
{
	R_base_start = start;
	R_base_end = end;
	Rotation33 R_start_end = R_base_start.inverse() * R_base_end;
	angle = GetRotAngle(R_start_end, rot_start_end);
}

Rotation33 RotationalInterpolation::Pos(double theta) const
{
	Rotation33 R_t = Rot2(rot_start_end, theta);
	return R_base_start * R_t;
}

Vector3d RotationalInterpolation::Vel(double theta, double thetad) const
{
	return R_base_start * (rot_start_end * thetad);
}

Vector3d RotationalInterpolation::Acc(double theta, double thetad, double thetadd) const
{
	return R_base_start * (rot_start_end * thetadd);
}

double RotationalInterpolation::Angle()
{
	return angle;
}

void RotationalInterpolation::Write(std::ostream &os) const
{
	os << "SingleAxis[] " << std::endl;
}

RotationalInterpolation::~RotationalInterpolation()
{
}

RotationalInterpolation *RotationalInterpolation::Clone() const
{
	return new RotationalInterpolation();
}

double RotationalInterpolation::GetRotAngle(Rotation33 data, Vector3d &axis, double eps) const
{
	double angle, x, y, z;		// variables for result
	double epsilon = eps;		// margin to allow for rounding errors
	double epsilon2 = eps * 10; // margin to distinguish between 0 and 180 degrees

	if (
		(std::fabs(data(0, 1) - data(1, 0)) < epsilon) && (std::fabs(data(0, 2) - data(2, 0)) < epsilon) && (std::fabs(data(1, 2) - data(2, 1)) < epsilon))
	{
		// singularity found
		// first check for identity matrix which must have +1 for all terms
		//  in leading diagonal and zero in other terms
		if ((std::fabs(data(0, 1) + data(1, 0)) < epsilon2) && (std::fabs(data(0, 2) + data(2, 0)) < epsilon2) && (std::fabs(data(1, 2) + data(2, 1)) < epsilon2) && (std::fabs(data(0, 0) + data(1, 1) + data(2, 2) - 3) < epsilon2))
		{
			// this singularity is identity matrix so angle = 0, axis is arbitrary
			// Choose 0, 0, 1 to pass orocos tests
			axis = {0, 0, 1};
			angle = 0.0;
			return angle;
		}

		// otherwise this singularity is angle = 180
		angle = 3.1415926;
		double xx = (data(0, 0) + 1) / 2;
		double yy = (data(1, 1) + 1) / 2;
		double zz = (data(2, 2) + 1) / 2;
		double xy = (data(0, 1) + data(1, 0)) / 4;
		double xz = (data(0, 2) + data(2, 0)) / 4;
		double yz = (data(1, 2) + data(2, 1)) / 4;

		if ((xx > yy) && (xx > zz))
		{
			// data[0] is the largest diagonal term
			x = sqrt(xx);
			y = xy / x;
			z = xz / x;
		}
		else if (yy > zz)
		{
			// data[4] is the largest diagonal term
			y = sqrt(yy);
			x = xy / y;
			z = yz / y;
		}
		else
		{
			// data[8] is the largest diagonal term so base result on this
			z = sqrt(zz);
			x = xz / z;
			y = yz / z;
		}
		axis = {x, y, z};
		return angle; // return 180 deg rotation
	}

	double f = (data(0, 0) + data(1, 1) + data(2, 2) - 1) / 2;

	x = (data(2, 1) - data(1, 2));
	y = (data(0, 2) - data(2, 0));
	z = (data(1, 0) - data(0, 1));
	axis = {x, y, z};

	angle = atan2(axis.norm() / 2, f);
	axis.normalize();
	return angle;
}

Rotation33 RotationalInterpolation::Rot2(const Vector3d &rotvec, double angle)
{
	// rotvec should be normalized !
	// The formula is
	// V.(V.tr) + st*[V x] + ct*(I-V.(V.tr))
	// can be found by multiplying it with an arbitrary vector p
	// and noting that this vector is rotated.
	double ct = cos(angle);
	double st = sin(angle);
	double vt = 1 - ct;
	double m_vt_0 = vt * rotvec(0);
	double m_vt_1 = vt * rotvec(1);
	double m_vt_2 = vt * rotvec(2);
	double m_st_0 = rotvec(0) * st;
	double m_st_1 = rotvec(1) * st;
	double m_st_2 = rotvec(2) * st;
	double m_vt_0_1 = m_vt_0 * rotvec(1);
	double m_vt_0_2 = m_vt_0 * rotvec(2);
	double m_vt_1_2 = m_vt_1 * rotvec(2);
	Rotation33 Rre;
	Rre << ct + m_vt_0 * rotvec(0),
		-m_st_2 + m_vt_0_1,
		m_st_1 + m_vt_0_2,
		m_st_2 + m_vt_0_1,
		ct + m_vt_1 * rotvec(1),
		-m_st_0 + m_vt_1_2,
		-m_st_1 + m_vt_0_2,
		m_st_0 + m_vt_1_2,
		ct + m_vt_2 * rotvec(2);

	return Rre;
}
