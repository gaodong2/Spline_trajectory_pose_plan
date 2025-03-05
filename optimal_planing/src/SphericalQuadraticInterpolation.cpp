#include "SphericalQuadraticInterpolation.h"

Squad::Squad(){};

Squad::~Squad(){};

double Squad::sCurve(double t) {
    return 3 * pow(t,2) - 2 * pow(t,3);
}

Eigen::Quaterniond Squad::quatMultiply(Eigen::Quaterniond &q1, Eigen::Quaterniond &q2) {
    Eigen::Quaterniond q(q1.w()*q2.w() - q1.x()*q2.x() - q1.y()*q2.y() - q1.z()*q2.z(),
         q1.w()*q2.x() + q1.x()*q2.w() + q1.y()*q2.z() - q1.z()*q2.y(),
         q1.w()*q2.y() - q1.x()*q2.z() + q1.y()*q2.w() + q1.z()*q2.x(),
         q1.w()*q2.z() + q1.x()*q2.y() - q1.y()*q2.x() + q1.z()*q2.w()); 
    return q;
}

Eigen::Quaterniond Squad::quatInverse(Eigen::Quaterniond &q) {
    Eigen::Quaterniond q_inv(q.w(), -q.x(), -q.y(), -q.z());
    return q_inv;    
}

Eigen::Quaterniond Squad::quatNormalize(Eigen::Quaterniond &q) {
    Eigen::Quaterniond q_normal = q.normalized();
    return q_normal;    
}

Eigen::Quaterniond Squad::quatlog(Eigen::Quaterniond &q) {
    Eigen::Quaterniond ln_q;
    
    double theta = acos(q.w());
    double t = 0.1;
    if(fabs(theta) < eps) {
        ln_q = Eigen::Quaterniond(0, 0, 0, 0);
    } else {
        double x = q.x()* (theta / sin(theta));
        double y = q.y()* (theta / sin(theta));
        double z = q.z()* (theta / sin(theta));
        ln_q = Eigen::Quaterniond(0, x, y, z);        
    }
    return ln_q;
}

Eigen::Quaterniond Squad::quatExp(Eigen::Quaterniond &v) {
    // 四元数的指数映射
    double theta = sqrt(v.x()*v.x()+v.y()*v.y()+v.z()*v.z());
    Eigen::Quaterniond q;
    if(theta < eps) {
        q = Eigen::Quaterniond(1, 0, 0, 0);    
    } else {
        q = Eigen::Quaterniond(cos(theta), v.x()/theta*sin(theta), v.y()/theta*sin(theta), v.z()/theta*sin(theta));    
    }
    return q;
}

Eigen::Quaterniond Squad::Slerp(Eigen::Quaterniond &q1, Eigen::Quaterniond &q2, double t) {
    return q1.slerp(t, q2);
}


double Squad::EvalAlpha(double s, int i, int L) {
    double k = s*(L-1)+1;
    double alpha;
    if((i>=k-1)&&(i<k)) {
        alpha= (double)k-((double)i);        
    } else {
        alpha=0.0;        
    }
    return alpha;
}

Eigen::Quaterniond Squad::QuadSquad(std::vector<Eigen::Quaterniond> &q, double t) {
    int L = q.size();
    Eigen::Quaterniond val = q[0];
    double s = sCurve(t);
    for(int j=1;j<L;j++){
        double C = q[j-1].dot(q[j]);
        if(C<0) {
            Eigen::Quaterniond new_q(-q[j].w(), -q[j].x(), -q[j].y(), -q[j].z());
            q[j] = new_q;           
        }
    }

    if(AlmostEqual(s,0.0)) {
        return q[0];
    } else if (AlmostEqual(s,1.0)) {
        return q[L-1];   
    }

    for(int j=1;j<L;j++) {
        // 全局细分映射到局部细分
        double alpha=EvalAlpha(s,j,L);
        double t= alpha;
        std::cout << "alpha: " << alpha << std::endl;
        if(alpha>0) {
            double EPS = 1e-9;

            // 计算两个姿态的点积，结果范围[-1,1]
            double C = q[j-1].dot(q[j]);
            if ((1.0 - C) <= EPS) {
                // 姿态之间过于接近采用线性插补
                val=Eigen::Quaterniond(q[j-1].w()*(1-s)+q[j].w()*s, q[j-1].x()*(1-s)+q[j].x()*s, q[j-1].y()*(1-s)+q[j].y()*s, q[j-1].z()*(1-s)+q[j].z()*s); 
                val = quatNormalize(val);
            return val;
            } 
            if((1.0 + C) <= EPS) {
                // 当姿态夹角接近180，无最短路径，结果不确定
                // 将角度旋转90
                Eigen::Quaterniond qtemp(q[j].z(), -q[j].y(), q[j].x(), -q[j].w());
                q[j] = qtemp;                
            }
            // Calculate interplations
            Eigen::Quaterniond qa = GetIntermediateControlPoint(j-1,q);
            Eigen::Quaterniond qap1 = GetIntermediateControlPoint(j,q);

            // 插补
            Eigen::Quaterniond qtemp1 = Slerp(q[j-1], q[j], t);
            Eigen::Quaterniond qtemp2 = Slerp(qa, qap1, t);
            Eigen::Quaterniond squad = Slerp(qtemp1, qtemp2, 2*t*(1-t));
            val = squad;val = quatNormalize(val);
            return val;
        }
    }
    return val;
}


bool Squad::AlmostEqual(double data1, double data2) {
    bool isEqual = std::fabs(data1-data2) <= eps ? true : false;
    return isEqual;
}

Eigen::Quaterniond Squad::GetIntermediateControlPoint(int j, std::vector<Eigen::Quaterniond> &q) {
// 插补中间点
// 若点为起点和终点，则直接返回当前值
    Eigen::Quaterniond qa;
    int L = q.size();
    if(j==0) {
        return q[0];    
    } else if (j==L-1) {
        return q[L-1];
    } else {
        Eigen::Quaterniond qji = quatInverse(q[j]);
        Eigen::Quaterniond qiqm1 = quatMultiply(qji,q[j-1]);
        Eigen::Quaterniond qiqp1 = quatMultiply(qji,q[j+1]);
        Eigen::Quaterniond ang_vel((-quatlog(qiqp1).w()-quatlog(qiqm1).w())/4, (-quatlog(qiqp1).x()-quatlog(qiqm1).x())/4, (-quatlog(qiqp1).y()-quatlog(qiqm1).y())/4, (-quatlog(qiqp1).z()-quatlog(qiqm1).z())/4); 
        Eigen::Quaterniond q2 = quatExp(ang_vel);
        qa = quatMultiply(q[j],q2);
        qa = quatNormalize(qa);    
    }
    return qa;
}

double Squad::ComputeAngularVelocityAndAcceleration(std::vector<Eigen::Quaterniond> &q, Eigen::MatrixXd &dt) {
    // 输入：
    // q: Nx4 矩阵，N 个时间步的四元数（每行是一个四元数 [w, x, y, z]）
    // dt: 时间步长
    // 输出：
    // omega: Nx3 矩阵，角速度（每行是 [ωx, ωy, ωz]）
    // alpha: Nx3 矩阵，角加速度（每行是 [αx, αy, αz]）

    int N = q.size(); // 时间步数
    Eigen::MatrixXd omega(N, 3); // 初始化角速度
    Eigen::MatrixXd alpha(N, 3); // 初始化角加速度
    Eigen::MatrixXd omegaNormal(N, 1); // 初始化角加速度

    // 计算四元数导数 dq/dt
    Eigen::MatrixXd dq(N, 4);
    dq(0, 0) = (q[1].w() - q[0].w()) / (dt(0));   
    dq(0, 1) = (q[1].x() - q[0].x()) / (dt(0)); 
    dq(0, 2) = (q[1].y() - q[0].y()) / (dt(0)); 
    dq(0, 3) = (q[1].z() - q[0].z()) / (dt(0));      

    for(int i=1;i<N-1;i++) {
        dq(i, 0) = (q[i+1].w() - q[i-1].w()) / (2 * dt(i-1)); // 中心差分
        dq(i, 1) = (q[i+1].x() - q[i-1].x()) / (2 * dt(i-1)); 
        dq(i, 2) = (q[i+1].y() - q[i-1].y()) / (2 * dt(i-1)); 
        dq(i, 3) = (q[i+1].z() - q[i-1].z()) / (2 * dt(i-1));     
    }

    dq(N-1, 0) = (q[N-1].w() - q[N-2].w()) / (dt(N-1));
    dq(N-1, 1) = (q[N-1].x() - q[N-2].x()) / (dt(N-1));
    dq(N-1, 2) = (q[N-1].y() - q[N-2].y()) / (dt(N-1));
    dq(N-1, 3) = (q[N-1].z() - q[N-2].z()) / (dt(N-1));

    // 计算角速度 omega
    for(int i=0;i<N;i++) {
        Eigen::Quaterniond q_inv = quatInverse(q[i]); // 四元数的逆
        Eigen::Quaterniond delatq(dq(i,0), dq(i,1), dq(i,2), dq(i,3));
        Eigen::Quaterniond omega_q = quatMultiply(q_inv, delatq); // 角速度（四元数形式）
        omega(i, 0) = 2.0 * omega_q.x(); // 提取虚部（角速度向量）
        omega(i, 1) = 2.0 * omega_q.y();
        omega(i, 2) = 2.0 * omega_q.z();
        omegaNormal(i,0) = sqrt(omega(i,0)*omega(i,0) + omega(i,1)*omega(i,1) + omega(i,2)*omega(i,2));
    }

    std::cout << omega << std::endl;
    // 计算角加速度 alpha
    for(int i=1;i<N-1;i++) {
        for(int j=0;j<3;j++) {
            alpha(i, j) = (omega(i+1, j) - omega(i-1, j)) / (2 * dt(i-1)); // 中心差分            
        }
    }
    double maxD = FindMax(omegaNormal);
    return maxD;
}

double Squad::FindMax(Eigen::MatrixXd &omegaNormal) {
    double maxData = 0.0;
    for(int i=0;i<omegaNormal.rows();i++) {
        if(maxData < omegaNormal(i,0)) {
            maxData = omegaNormal(i,0);
        }
    }
    return maxData;
}