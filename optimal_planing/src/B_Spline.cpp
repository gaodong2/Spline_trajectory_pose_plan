#include "B_Spline.h"

BSplines::BSplines(){}
BSplines::~BSplines(){}

void BSplines::UniformCubicBSplines(Eigen::MatrixXd &points, int seg) {
    int cols = points.cols();
    int rows = points.rows(); 
    // 离散化参数 u
    Eigen::MatrixXd u = LinSpace(seg);
   
    // 初始化位置、速度和加速度
    C = Eigen::MatrixXd::Zero(seg, 3);  // 位置
    V = Eigen::MatrixXd::Zero(seg, 3);  // 速度
    A = Eigen::MatrixXd::Zero(seg, 3);  // 加速度
    J = Eigen::MatrixXd::Zero(seg, 3);  // 加加速度
    Ls = Eigen::MatrixXd::Zero(seg, 1);  // 总长
    // 初始化速度向量的模
    Eigen::MatrixXd speed_norm(seg, 1);
    kappa = Eigen::MatrixXd::Zero(seg, 1);  // 曲率半径
    // std::cout << kappa << std::endl;
    Eigen::MatrixXd knots = KnotsUniform(rows-1, cols);

    // 计算位置、速度和加速度
    for(int i=0;i<seg;i++) {
        // 计算基函数及其导数
        Eigen::MatrixXd N = BSplineBasis(cols, knots, u(i,0));

        // 计算位置
        C.block(i,0,1,3) = N.block(0,0,1,rows) * points;
        
        // 计算速度
        V.block(i,0,1,3) = N.block(1,0,1,rows) * points;
        
        // 计算加速度
        A.block(i,0,1,3) = N.block(2,0,1,rows) * points;

        // 计算加速度
        J.block(i,0,1,3) = N.block(3,0,1,rows) * points;

        // 计算曲率半径
        Eigen::Vector3d dV = {V(i,0),V(i,1),V(i,2)};
        Eigen::Vector3d dA = {A(i,0),A(i,1),A(i,2)};
        Eigen::MatrixXd cross_prod = dV.cross(dA);  // 二维曲线需要补零

        kappa(i, 0) = 1 / (sqrt(cross_prod(2, 0) * cross_prod(2, 0) ) / (dV.norm() * dV.norm() * dV.norm()));
        // 计算速度向量的模
        speed_norm(i, 0) =  V.block(i,0,1,3).norm();
        if(i>0) {
            Ls(i,0) = TrapezoidalIntegral(speed_norm(i-1, 0), speed_norm(i, 0), u(1));   
        }
        
    }
}

 void BSplines::UniformCubicBSplines(Eigen::MatrixXd &points, Eigen::MatrixXd &u) {
    int cols = points.cols();
    int rows = points.rows(); 
    // 离散化参数 u
    // Eigen::MatrixXd u = LinSpace(seg);
    int urows = u.rows();

    // 初始化位置、速度和加速度
    C = Eigen::MatrixXd::Zero(urows, 3);  // 位置
    V = Eigen::MatrixXd::Zero(urows, 3);  // 速度
    A = Eigen::MatrixXd::Zero(urows, 3);  // 加速度
    J = Eigen::MatrixXd::Zero(urows, 3);  // 加加速度
    Ls = Eigen::MatrixXd::Zero(urows, 1);  // 总长
    // 初始化速度向量的模
    Eigen::MatrixXd speed_norm(urows, 1);
    kappa = Eigen::MatrixXd::Zero(urows, 1);  // 曲率半径
    Eigen::MatrixXd knots = KnotsUniform(rows-1, cols);

    // 计算位置、速度和加速度
    for(int i=0;i<urows;i++) {

        // 计算基函数及其导数
        Eigen::MatrixXd N = BSplineBasis(cols, knots, u(i,0));

        // 计算位置
        C.block(i,0,1,3) = N.block(0,0,1,rows) * points;
      
        // 计算速度
        V.block(i,0,1,3) = N.block(1,0,1,rows) * points;
        
        // 计算加速度
        A.block(i,0,1,3) = N.block(2,0,1,rows) * points;

        // 计算加速度
        J.block(i,0,1,3) = N.block(3,0,1,rows) * points;

        // 计算曲率半径
        Eigen::Vector3d dV = {V(i,0),V(i,1),V(i,2)};
        Eigen::Vector3d dA = {A(i,0),A(i,1),A(i,2)};
        Eigen::MatrixXd cross_prod = dV.cross(dA);  // 二维曲线需要补零

        kappa(i, 0) = 1 / (sqrt(cross_prod(2, 0) * cross_prod(2, 0) ) / (dV.norm() * dV.norm() * dV.norm()));
        // 计算速度向量的模
        speed_norm(i, 0) =  V.block(i,0,1,3).norm();
        if(i>0) {
            Ls(i,0) = TrapezoidalIntegral(speed_norm(i-1, 0), speed_norm(i, 0), u(1));   
        }
        
    }
 }

Eigen::MatrixXd BSplines::KnotsUniform(int len, int dim) {
    Eigen::MatrixXd knots(1, len+dim+2);
    for(int i=0; i<len+dim+2; i++) {
        if(i < dim + 1) {
            knots(0,i) = 0;            
        } else if(i >= len + 1) {
            knots(0,i) = 1;            
        } else {
            knots(0,i) = (i - dim) / (len - dim + 1.0);      
        }
    }
    return knots;
}

Eigen::MatrixXd BSplines::BSplineBasis(int p, Eigen::MatrixXd &knots, double u) {
    if(u == 1.0) {
        u = u - 1e-10;
    }
    // 计算基函数值
    int cols = knots.cols();
    Eigen::MatrixXd N(4, cols - p - 1);
    for(int i=0;i<N.cols();i++) {
        N(0,i) = BSplineBasisValue(p, knots, i, u);
        N(1,i) = BSsplineBasisDerivative(p, knots, i, u, 1);
        N(2,i) = BSsplineBasisDerivative(p, knots, i, u, 2);
        N(3,i) = BSsplineBasisDerivative(p, knots, i, u, 3);
    }
    
    return N;
}

// 定义 B 样条基函数及其导数
double  BSplines::BSplineBasisValue(int p, Eigen::MatrixXd &knots, int i, double u){
    double val;
    if(p == 0) {
        val = (knots(0, i) <= u) && (u < knots(0, i + 1));        
    } else {
        double denom1 = knots(0, i + p) - knots(0, i);
        double denom2 = knots(0, i + p + 1) - knots(0, i + 1);
        double term1 = 0, term2 = 0;
        // std::cout << denom1 << std::endl;
        // std::cout << denom2 << std::endl;
        if(!AlmostEqual(denom1,0)) {
            term1 = (u - knots(0, i)) / denom1 * BSplineBasisValue(p - 1, knots, i, u);            
        }
        if(!AlmostEqual(denom2,0)) {
            term2 = (knots(0, i + p + 1) - u) / denom2 * BSplineBasisValue(p - 1, knots, i + 1, u);            
        }
        // std::cout << "val: "<< val << std::endl;
        val = term1 + term2;
    }
    return val;

}

// 计算 B 样条基函数的导数
double BSplines::BSsplineBasisDerivative(int p, Eigen::MatrixXd &knots, int i, double u, int order) {
    double val;
    if(order == 0) {
        val = BSplineBasisValue(p, knots, i, u);
    } else {
        double denom1 = knots(0, i + p) - knots(0, i);
        double denom2 = knots(0, i + p + 1) - knots(0, i + 1);

        double term1 = 0, term2 = 0;
        
        if(!AlmostEqual(denom1,0)) {
            term1 = p / denom1 * BSsplineBasisDerivative(p - 1, knots, i, u, order - 1);            
        }
        if(!AlmostEqual(denom2,0)) {
            term2 = -p / denom2 * BSsplineBasisDerivative(p - 1, knots, i + 1, u, order - 1);            
        }
        val = term1 + term2;
    } 
    return val;    
}

double BSplines::TrapezoidalIntegral(double deltaX_1, double deltaX, double Y) {
    return (deltaX_1 + deltaX)*Y/2;
}


bool BSplines::AlmostEqual(double data1, double data2, double eps) {
    bool isEqual = std::fabs(data1-data2) <= eps ? true : false;
    return isEqual;
}

Eigen::MatrixXd BSplines::LinSpace(int seg) {
    Eigen::MatrixXd u(seg, 1);
    double delta = 1.0/(seg-1);
    for(int i=0;i<seg;i++) {
        u(i,0) = i * delta;  // 100 个离散点
    }
    return u;
}

Eigen::MatrixXd BSplines::GetCurvePos() {
    return C;
}

Eigen::MatrixXd BSplines::GetCurveVel() {
    return V;
}

Eigen::MatrixXd BSplines::GetCurveAcc() {
    return A;
}

Eigen::MatrixXd BSplines::GetCurveJerk() {
    return J;
}

Eigen::MatrixXd BSplines::GetCurvature() {
    return kappa;
}

Eigen::MatrixXd BSplines::GetDeltaL() {
    return Ls;    
}

double BSplines::GetTotalL() {
    return CalTotalL();    
}      
double BSplines::CalTotalL() {
    double res = 0.0;
    for(int i=0;i<Ls.rows();i++){
        res += Ls(i,0);
    }
    return res;
}