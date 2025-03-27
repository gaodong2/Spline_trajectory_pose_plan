#include "S_Type_Velocity.h"

S_Type_Velocity::S_Type_Velocity() {};

S_Type_Velocity::~S_Type_Velocity() {};

Eigen::MatrixXd S_Type_Velocity::STypeVelPlan(Eigen::MatrixXd &L, Eigen::MatrixXd &kappa, double threshold, double Vx, double Ax, double j_max, double scale, double period) {
    double v0 = 0;
    double vt = 0;
    Eigen::MatrixXd TJAV, TJAVf;
    Eigen::MatrixXd Vmax =  MaxVel(kappa, threshold, Vx);
    Eigen::MatrixXd V_reachability = VelocityPlanning(L, Vmax, Ax, v0, vt);
    for(int i=0;i<V_reachability.rows();i++){
        V_reachability(i,0) = V_reachability(i,0)*scale;
    }
    Vx = Vx * scale;
    TrapezoidalPlan(V_reachability, L, Ax);
    TJAV = FindSpeedUpStage(V_reachability, L, Vx, Ax, j_max);
    TJAVf = FindShutDownStage(V_reachability, L, Vx, Ax, j_max);      
    
    Eigen::MatrixXd newTAVS = STypeSpeedUp(TJAV);
    Eigen::MatrixXd newTAVSf = STypeReach(TJAVf, newTAVS);
    Eigen::MatrixXd uus = Compose(period, newTAVSf);
    return uus;
}

void S_Type_Velocity::TrapezoidalPlan(Eigen::MatrixXd &V_reachability, Eigen::MatrixXd &s, double a_max) {
    int rows = V_reachability.rows();
    int cols = V_reachability.cols();
    DelS = Eigen::MatrixXd::Zero(rows, cols);
    sumS = Eigen::MatrixXd::Zero(rows, cols);
    sumT = Eigen::MatrixXd::Zero(rows, cols);
    sumV = Eigen::MatrixXd::Zero(rows, cols);
    sumA = Eigen::MatrixXd::Zero(rows, cols);
    double sum_T = 0, sum_S = 0;
    for(int i=0;i<rows-1;i++){
        double v0 = V_reachability(i,0);
        double vt = V_reachability(i+1,0); 
        double v_max = std::max(v0, vt); 
        double TTT = 2 * s(i+1,0)  / (v0 + vt);
        double a = (vt - v0) / TTT;
        
        double s = v0 * TTT + 0.5 * a * TTT*TTT;
        sum_T += TTT;
        sum_S += s;
        DelS(i+1,0) = s;
        sumS(i+1,0) = sum_S;
        sumT(i+1,0) = sum_T;
        sumV(i+1,0) = vt;
        sumA(i+1,0) = a;
    }
    // std::cout << "sumS" << std::endl;
    // std::cout << sumS << std::endl;
}


Eigen::MatrixXd S_Type_Velocity::MaxVel(Eigen::MatrixXd &kappa, double threshold, double Vs) {
    int rows = kappa.rows();
    int cols = kappa.cols();

    Eigen::MatrixXd v_max = Eigen::MatrixXd::Zero(rows, cols);
    for(int i=0;i<rows;i++) {
        v_max(i,0) = Vs;      
        if(kappa(i,0) < threshold)
            v_max(i,0) = Vs * pow((kappa(i,0) / threshold),0.5);
    }
    return v_max;
}

Eigen::MatrixXd S_Type_Velocity::P2PMultiAxisDoubleSTrajectory(double sampleTime, Eigen::MatrixXd input) {
    int rows = input.rows();
    int cols = input.cols();
    

    double T, T1, T2, T3, T4, T5, T6, T7, Tj;
    double Da, Dup;
    
    double Sr = input(1) - input(0);
    double Dr = std::fabs(Sr);
    double Dir = Sr>0 ? 1.0 : -1.0;
    double Vm = input(2);
    double Am = input(3);
    double Jm = input(4);
    double Vj = Am*Am / (2 * Jm);  
    if (Vm - 2 * Vj >= 0) {
        Tj = Am / Jm;
        Da = 2 * Am*Am / Jm*Jm;
        Dup = Vm*Vm / (2 * Am) + Am * Vm / (2 * Jm);
        if (Dr >= 2 * Dup) {
            T1 = Tj;  T3 = Tj;  T5 = Tj;  T7 = Tj;
            T2 = (Vm - 2 * Vj) / Am;    T6 = T2;
            T4 = (Dr - 2 * Dup) / Vm;            
        } else {
            if (Da <= Dr) {
                T1 = Tj;  T3 = Tj;  T5 = Tj;  T7 = Tj;
                T2 = (sqrt(pow(Am, 4) + 4 * Dr * Am * pow(Jm, 2)) - 3 * pow(Am, 2)) / (2 *Jm * Am);
                T6 = T2;
                T4 = 0;                
            } else {
                T1 = pow(Dr / (2 * Jm), 1/3);  T3 = T1;  T5 = T1;  T7 = T1;
                T2 = 0; T4 = 0; T6 = 0;                 
            }     
        }  
    } else {
        Am = sqrt(Vm * Jm);
        Tj = Am / Jm;
        Da = 2 * pow(Am, 3) / pow(Jm, 2);
        Vj = pow(Am, 2) / (2 * Jm);
        T2 = 0; T6 = 0;
        if (Da <= Dr) {
            T1 = Tj;  T3 = Tj;  T5 = Tj;  T7 = Tj;
            T4 = (Dr - Da) / Vm;            
        } else {
            T1 = pow(Dr / (2 * Jm), 1/3);  T3 = T1;  T5 = T1;  T7 = T1;
            T4 = 0;            
        }
    }
    T = 4 * T1 + 2* T2 + T4;
    int N = T / sampleTime;
    Tp = Eigen::MatrixXd::Zero(1,N);
    Jp = Eigen::MatrixXd::Zero(1,N);
    Ap = Eigen::MatrixXd::Zero(1,N);
    Vp = Eigen::MatrixXd::Zero(1,N);
    Sp = Eigen::MatrixXd::Zero(1,N);
    Eigen::MatrixXd TJAVS(5,N);
    for(int j=0;j<N;j++) {
        Tp(0,j) = j * sampleTime;
        if (j * sampleTime < T1) {
            Jp(0,j) = Dir * Jm;            
        } else if (j * sampleTime >= T1 && j * sampleTime < T1 + T2) {
            Jp(0,j) = 0;            
        } else if (j * sampleTime >= T1 + T2 && j * sampleTime < T1 + T2 + T3) {
            Jp(0,j) = -1 * Dir * Jm;            
        } else if (j * sampleTime >= T1 + T2 + T3 && j * sampleTime < T1 + T2 + T3 + T4) {
            Jp(0,j) = 0;            
        } else if (j * sampleTime >= T1 + T2 + T3 + T4 && j * sampleTime < T1 + T2 + T3 + T4 + T5) {
            Jp(0,j) = -1 * Dir * Jm;            
        } else if (j * sampleTime >= T1 + T2 + T3 + T4 + T5 && j * sampleTime < T1 + T2 + T3 + T4 + T5 + T6) {
            Jp(0,j) = 0;            
        } else if (j * sampleTime >= T1 + T2 + T3 + T4 + T6) {
            Jp(0,j) = Dir * Jm;                 
        }
    }

    int pp = 0;
    int pm = 0;
    int mp = 0;
    int mm = 0;
    std::vector<int> pmj, mpj, mmj;
    for(int j=0;j<N-1;j++) {
        if(j < N/2) {
            if (Jp(0,j) > 0) {
                pp ++;                
            } else if(Jp(0,j) < 0) {
                pm ++;
                pmj.push_back(j);                
            }            
        } else {
            if(Jp(0,j) < 0) {
                mp ++;
                mpj.push_back(j);                
            } else if (Jp(0,j) > 0) {
                mm ++;
                mmj.push_back(j);                
            }            
        }
    }

    
    if(pp!=pm) {
        Jp(pmj.back()) = 0;         
    }
    if(mp!=mm) {
        Jp(0,mpj[0]) = 0;
        Jp(0,mmj[0]) = Jp(0, mmj[0]-1);            
    }
    for(int j=0;j<N;j++) {
        if(j == 0) {
            Ap(0,j) = Jp(0,j) * sampleTime;
            Vp(0,j) = 1 / 2 * Jp(0,j) * pow(sampleTime, 2);
            Sp(0,j) = 1 / 6 * Jp(0,j) * pow(sampleTime, 3);                
        } else {
            Ap(0,j) = Ap(0,j-1) + Jp(0,j) * sampleTime;
            Vp(0,j) = Vp(0,j-1) + Ap(0,j-1) * sampleTime + 1 / 2 * Jp(0,j) * pow(sampleTime, 2);
            Sp(0,j) = Sp(0,j-1) + Vp(0,j-1) * sampleTime + 1 / 2 * Ap(0,j-1) * pow(sampleTime, 2) + 1 / 6 * Jp(0,j) * pow(sampleTime, 3);                   
        }
    }

    // Jp(0,N) = 0;
    // Ap(0,N) = 0;
    // Vp(0,N) = 0;
    // Sp(0,N) = Sr;


    TJAVS.block(0,0,1,N) = Tp; 
    TJAVS.block(1,0,1,N) = Jp; 
    TJAVS.block(2,0,1,N) = Ap; 
    TJAVS.block(3,0,1,N) = Vp; 
    TJAVS.block(4,0,1,N) = Sp; 
    // std::cout << "TJAVS" << std::endl;
    // std::cout << TJAVS << std::endl;    
    return TJAVS.transpose();
}

double S_Type_Velocity::Sum(Eigen::MatrixXd &deltaT, int ii, int jj) {
    double res = 0.0;
    for(int i=ii;i<=jj;i++) {
        res = res + deltaT(i,0);
    }
    return res;
}

Eigen::MatrixXd S_Type_Velocity::STypeSpeedUp(Eigen::MatrixXd &TJAVS) {
    // S型速度加速
    Eigen::MatrixXd TTp, SSp, VVp, AAp;
    Eigen::MatrixXd TTTT = sumT;
    Eigen::MatrixXd AAAA = sumA;
    Eigen::MatrixXd VVVV = sumV;
    Eigen::MatrixXd SSSS = sumS;
    int steps = SSSS.rows();
    int w = TJAVS.rows();
    
    for(int jj=1;jj<=w/2;jj++) {
        if((fabs(TJAVS(jj,1)) < 1e-6 && TJAVS(jj-1,1) < 0) || jj==w/2) {
            TTp = TJAVS.block(0,0,jj-1,1);
            SSp = TJAVS.block(0,4,jj-1,1);
            VVp = TJAVS.block(0,3,jj-1,1);
            AAp = TJAVS.block(0,2,jj-1,1);
            break;              
        }
    }
    int Tros = SSp.rows();
    Eigen::MatrixXd NTT, deltaT(steps, 1);
    for(int jj=1;jj<steps;jj++) {
        deltaT(jj, 0) = TTTT(jj,0)-TTTT(jj-1,0);
    }
    Eigen::MatrixXd newTAVS;
    for(int jj=1;jj<steps;jj++) {
        // std::cout << SSSS(jj,0) << "  " << SSp(Tros-1,0) << std::endl;
        if(SSSS(jj,0)>SSp(Tros-1,0) && SSSS(jj-1,0)<SSp(Tros-1,0)) {
            double deltaS = SSSS(jj,0) - SSp(Tros-1,0);
            // std::cout << "deltaS" << std::endl;
            // std::cout << deltaS << std::endl; 
            double Tjj = 2 * deltaS / (VVVV(jj,0) + VVp(Tros-1,0));
            // std::cout << "Tjj" << std::endl;
            // std::cout << Tjj << std::endl;  
            NTT = Eigen::MatrixXd(steps-jj,1);
            for(int ii=jj;ii<steps;ii++) {
                deltaT(jj,0) = Tjj;
                NTT(ii-jj,0)  = TTp(Tros-1,0) + Sum(deltaT, jj, ii);                
            } 
            newTAVS = Eigen::MatrixXd(steps-jj+Tros, 4);
  
            newTAVS.block(0,3,Tros,1) = SSp;
            newTAVS.block(Tros,3,steps-jj,1) = SSSS.block(jj,0,steps-jj,1);

            newTAVS.block(0,2,Tros,1) = VVp;
            newTAVS.block(Tros,2,steps-jj,1) = VVVV.block(jj,0,steps-jj,1);

            newTAVS.block(0,1,Tros,1) = AAp;
            newTAVS.block(Tros,1,steps-jj,1) = AAAA.block(jj,0,steps-jj,1);

            newTAVS.block(0,0,Tros,1) = TTp;
            newTAVS.block(Tros,0,steps-jj,1) = NTT;
            int begins = jj;
            break;               
        }
      
    }
    return newTAVS;
}

Eigen::MatrixXd S_Type_Velocity::STypeReach(Eigen::MatrixXd &TJAVSf, Eigen::MatrixXd &newTJAVS) {
    Eigen::MatrixXd TTpf, SSpf, VVpf, AApf, NTT, NSS;
    int Srows = newTJAVS.rows();
    Eigen::MatrixXd TTTT = newTJAVS.block(0,0,Srows,1);
    Eigen::MatrixXd AAAA = newTJAVS.block(0,1,Srows,1);
    Eigen::MatrixXd VVVV = newTJAVS.block(0,2,Srows,1);
    Eigen::MatrixXd SSSS = newTJAVS.block(0,3,Srows,1);
    
    int rows = TJAVSf.rows();
    // S型速度减速 
    for(int jj=rows-1;jj>1;jj--) {
        if((fabs(TJAVSf(jj,1)) < 1e-6 && TJAVSf(jj+1,1) < 0) || jj==rows/2) {
            // std::cout << "jj" << std::endl;
            // std::cout << jj << std::endl;           
            TTpf = TJAVSf.block(jj-1,0,rows-jj,1);
            SSpf = TJAVSf.block(jj-1,4,rows-jj,1);
            VVpf = TJAVSf.block(jj-1,3,rows-jj,1);
            AApf = TJAVSf.block(jj-1,2,rows-jj,1);
            break;        
        }        
    }


    double deltas = TJAVSf(rows-1,4) - SSpf(0,0);
    std::cout << "VVpf" << std::endl;
    std::cout << VVpf.rows() << std::endl;    
    std::cout << "deltas" << std::endl;
    std::cout << deltas << std::endl;
    std::cout << TJAVSf(rows-1,4) << std::endl;
    std::cout << SSpf(0,0) << std::endl;
    
    int lenpf = AApf.rows();
    Eigen::MatrixXd newTAVSf;
    for(int jj=Srows-1;jj>1;jj--) {
        // std::cout << "deltas" << std::endl;
        // std::cout << deltas << std::endl;
        // std::cout << (SSSS(Srows-1,0)-SSSS(jj,0)) << std::endl;
        
        if ((SSSS(Srows-1,0)-SSSS(jj,0))<deltas && (SSSS(Srows-1,0)-SSSS(jj-1,0))>deltas) {
            double deltaS = SSSS(Srows-1,0) - SSSS(jj-1,0) - deltas;
            double deltaT = deltaS / VVpf(0,0);
            double importS = SSSS(jj-1,0) + deltaS;
            double importT = TTTT(jj-1,0) + deltaT;

            // std::cout << "deltas" << std::endl;
            // std::cout << SSSS(Srows-1,0) << std::endl;
            // std::cout << SSSS(jj-1,0) << std::endl;
            // std::cout << deltas << std::endl;
            // std::cout << deltaS << std::endl;
            // std::cout << deltaT << std::endl;
            // std::cout << importS << std::endl;
            // std::cout << importT << std::endl;
            NSS = Eigen::MatrixXd(lenpf,1);
            NTT = Eigen::MatrixXd(lenpf,1);
            for(int ii=0;ii<lenpf;ii++) {              
                NSS(ii,0) = importS + SSpf(ii,0) - SSpf(0,0);
                NTT(ii,0) = importT + TTpf(ii,0) - TTpf(0,0);                
            }

            std::cout << "NSS" << std::endl;
            std::cout << NSS.rows() << std::endl;   
            std::cout << NTT.rows() << std::endl;  
            newTAVSf = Eigen::MatrixXd(jj-1+lenpf, 4);
            std::cout << "SSSS.block" << std::endl;
            std::cout << jj << std::endl; 
            std::cout << SSSS.block(0,0,jj-1,0) << std::endl; 
            newTAVSf.block(0,3,jj-1,1) = SSSS.block(0,0,jj-1,1);
            newTAVSf(jj-1,3) = importS;
            newTAVSf.block(jj,3,lenpf-1,1) = NSS.block(1,0,lenpf-1,1);

            newTAVSf.block(0,2,jj-1,1) = VVVV.block(0,0,jj-1,1);
            newTAVSf(jj-1,2) = VVpf(0,0);
            newTAVSf.block(jj,2,lenpf-1,1) = VVpf.block(1,0,lenpf-1,1);

            newTAVSf.block(0,1,jj-1,1) = AAAA.block(0,0,jj-1,1);
            newTAVSf(jj-1,1) = AApf(0,0);
            newTAVSf.block(jj,1,lenpf-1,1) = AApf.block(1,0,lenpf-1,1);

            newTAVSf.block(0,0,jj-1,1) = TTTT.block(0,0,jj-1,1);
            newTAVSf(jj-1,0) = importT;
            newTAVSf.block(jj,0,lenpf-1,1) = NTT.block(1,0,lenpf-1,1);
            int ends = Srows - jj;
            break;            
        }        
    }
    // std::cout << "newTAVSf" << std::endl;
    // std::cout << newTAVSf << std::endl;   
    return newTAVSf;
}

Eigen::MatrixXd S_Type_Velocity::FindSpeedUpStage(Eigen::MatrixXd &V_reachability, Eigen::MatrixXd &ssss, double Vmax, double a_sup, double j_max) {
    int steps = sumS.rows();
    Eigen::MatrixXd List(5,1);
    Eigen::MatrixXd TJAVS;
    for(int i=1;i<steps-2;i++) {
        if((V_reachability(i+1)-V_reachability(i)<0) && (V_reachability(i)-V_reachability(i-1)>=0)) {

            double suml = sumS(i, 0);
            List << 0, 2*suml, Vmax, a_sup, j_max;
            std::cout << List << std::endl;
            TJAVS = P2PMultiAxisDoubleSTrajectory(0.001, List);
            break;            
        }
    }    
    // std::cout << "TJAVS" << std::endl;
    // std::cout << TJAVS.transpose() << std::endl; 
    return TJAVS;
}

Eigen::MatrixXd S_Type_Velocity::FindShutDownStage(Eigen::MatrixXd &V_reachability, Eigen::MatrixXd &ssss, double Vmax, double a_sup, double j_max) {
    int steps = sumS.rows();
    Eigen::MatrixXd List(5,1);
    Eigen::MatrixXd TJAVSf;
    for(int i=steps-1;i>2;i--){
        if((V_reachability(i)-V_reachability(i-1)>0) && (V_reachability(i+1)-V_reachability(i)<=0)) {
            double suml = sumS(steps-1) - sumS(i);
            List << 0, 2*suml, Vmax, a_sup, j_max;
            std::cout << List << std::endl;
            TJAVSf = P2PMultiAxisDoubleSTrajectory(0.001, List);
            break;   
        }
    }   
    // std::cout << "TJAVSf" << std::endl;
    // std::cout << TJAVSf.transpose() << std::endl; 
    return TJAVSf;
}


Eigen::MatrixXd S_Type_Velocity::Compose(double dt, Eigen::MatrixXd &newTJAVSf) {
    Eigen::MatrixXd uus;
    int mm = newTJAVSf.rows();
    Eigen::MatrixXd TTTT = newTJAVSf.block(0,0,mm,1);
    Eigen::MatrixXd SSSS = newTJAVSf.block(0,3,mm,1);
    // std::cout << "TTTT" <<std::endl;
    // std::cout << TTTT <<std::endl;
    // std::cout << "SSSS" <<std::endl;
    // std::cout << SSSS <<std::endl;
    // 不同规划方式融合
    int rows = int(TTTT(mm-1,0) / dt)+1;
    uus = Eigen::MatrixXd(rows,1);
    double delta, percent, lucheng, uu;
    int j = 0;
    for(double i=0;i<TTTT(mm-1,0);i=i+dt) {
        // std::cout << "i" <<std::endl;
        // std::cout << i <<std::endl; 
        for(int j=0;j<mm-1;j++) {
            if ((i >= TTTT(j,0)) && (i < TTTT(j+1,0))) {
                delta = i - TTTT(j,0);
                // std::cout << "delta" <<std::endl;
                // std::cout << delta <<std::endl;
                percent = delta/(TTTT(j+1,0) - TTTT(j,0));
                // std::cout << percent <<std::endl;
                lucheng = percent * (SSSS(j+1,0) - SSSS(j,0)) + SSSS(j,0);
                // std::cout << lucheng <<std::endl;
                uu = lucheng / SSSS(mm-1,0);
                // std::cout << uu <<std::endl;
                if(uu > 1) {
                uu = 1;                    
                }
                break;               
            }      
        }
        uus(j,0) = uu;   
        j++; 
        // std::cout << j <<std::endl;
    }
    return uus;
}

Eigen::MatrixXd S_Type_Velocity::VelocityPlanning(Eigen::MatrixXd &L, Eigen::MatrixXd &v_max, double a_sup, double v0, double vt) {
    // 输入参数
    int rows = L.rows();
    Eigen::MatrixXd a_max = Eigen::MatrixXd::Ones(rows, 1) * a_sup;  // 各点的最大加速度 (m/s^2)
    Eigen::MatrixXd v = Eigen::MatrixXd::Ones(rows, 1);
    double delta_s;
    // 前向传播
    Eigen::MatrixXd v_forward = Eigen::MatrixXd::Ones(rows, 1);
    v_forward(0,0) = v0;
    for(int i=1;i<rows;i++){
        delta_s = L(i,0);
        v_forward(i,0) = std::min(v_max(i,0), sqrt(v_forward(i-1,0)*v_forward(i-1,0) + 2 * a_max(i,0) * delta_s));
    }

    // 反向传播
    Eigen::MatrixXd v_backward = Eigen::MatrixXd::Ones(rows, 1);
    v_backward(rows-1,0) = vt;
    for (int i=rows-2;i>0;i--){
        delta_s = L(i+1,0);
        v_backward(i,0) = std::min(v_max(i,0), sqrt(v_backward(i+1,0)*v_backward(i+1,0) + 2 * a_max(i,0) * delta_s));
        
    }

    // 可达速度
    for(int i=0;i<rows;i++) {
        v(i,0) = std::min(v_forward(i,0), v_backward(i,0));   
    }
    
    return v;  
}