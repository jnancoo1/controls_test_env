#ifndef NORMS_HPP
#define NORMS_HPP

#include <Eigen/Dense>
#include "discrete_state_space.hpp"
#include <cmath>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include <complex>
#include "analysis.hpp"
class Norms {
public:
    // Compute H2 norm for stable continuous-time systems: sqrt(trace(C * P * C^T))
    static double H2_norm_continuous(const Discrete_StateSpace_System& System){

        Eigen::MatrixXd Q=System.B*System.B.transpose();
        Eigen::MatrixXd X=(solve_lyapunov(System.A,Q));

        double norm=sqrt((System.C*X*System.C.transpose()).diagonal().sum());

        return norm;
    }
    // Compute H2 norm for stable discrete-time systems
    static double H2_norm_discrete(const Discrete_StateSpace_System& System) {
        Eigen::MatrixXd Q = System.B * System.B.transpose();
        Eigen::MatrixXd X = solve_discrete_lyapunov(System.A, Q);

        double norm = std::sqrt((System.C * X * System.C.transpose()).trace());

        return norm;
    }

    // Solve continuous Lyapunov equation A'P + PA + Q = 0
    static Eigen::MatrixXd solve_lyapunov(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Q){

        int n=A.rows();
        Eigen::VectorXd vecQ = Eigen::Map<const Eigen::VectorXd>(Q.data(), Q.size());
         Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);

        Eigen::MatrixXd kron1 = Eigen::kroneckerProduct(A, I);
        Eigen::MatrixXd kron2 = Eigen::kroneckerProduct(I,A);

        Eigen::VectorXd w = (kron1 + kron2).fullPivLu().solve(-vecQ);

        Eigen::MatrixXd P = Eigen::Map<Eigen::MatrixXd>(w.data(), n, n);

        return P;
    };

    // Solve discrete Lyapunov equation P = A*P*A' + Q
    static Eigen::MatrixXd solve_discrete_lyapunov(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Q){

        int n=A.rows();

        Eigen::VectorXd vecQ=Eigen::Map<const Eigen::VectorXd>(Q.data(),Q.size());
        Eigen::MatrixXd kron1 = Eigen::kroneckerProduct(A,A);

        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n * n, n * n);
        Eigen::MatrixXd K = kron1 - I;


        Eigen::VectorXd p=K.fullPivLu().solve(-vecQ);
        Eigen::MatrixXd P = Eigen::Map<Eigen::MatrixXd>(p.data(), n, n);

        return P;
    };

    // Compute Frobenius norm of a matrix
    static double frobenius_norm(const Eigen::MatrixXd& M){double sum = 0.0;
    for (int i = 0; i < M.rows(); ++i) {
        for (int j = 0; j < M.cols(); ++j) {
            sum += M(i, j) * M(i, j);
        }
    }
    return std::sqrt(sum);
};

    //generate range of frequencies
    static std::vector<std::complex<double>> frequency_range(double low,double high,int n){

        std::vector<std::complex<double>> out;
        std::complex<double> temp;
        std::complex<double> onei(0.0,1.0);
        double w;
        for(int i=0;i<n;i++){

            w = low * pow(high / low, static_cast<double>(i) / (n - 1));
            std::complex<double> s = std::complex<double>(0, w); // jω
            out.push_back(w);
        }
        return out;

    }

    // Compute Frobenius norm of a State space
    static double frobenius_norm_state_space(const Discrete_StateSpace_System& System,int n){
        Eigen::MatrixXcd Identity=Eigen::MatrixXcd::Identity(System.A.rows(),System.A.cols());
        std::vector<std::complex<double>> jw=frequency_range(0.001,1000,1000);
        std::vector<std::complex<double>> Gjw_all;
        Eigen::MatrixXcd Gjwi;
        double total=0;
        for(int i=0;i<n;i++){

            Gjwi=System.C*((jw[i]*Identity-System.A).fullPivLu().inverse())*System.B;
            total+=(Gjwi.squaredNorm());
        }
        return std::sqrt(total / n);
        
    }
    // Compute spectral norm (2-norm)
    static double spectral_norm(const Discrete_StateSpace_System& System,double freq){

        Eigen::MatrixXcd Identity=Eigen::MatrixXcd::Identity(System.A.rows(),System.A.cols());
        double w=2*M_PI*freq;
        std::complex<double> jw=std::complex<double>(0, w);
        Eigen::MatrixXcd Gjw;
        

        Gjw=System.C*((jw*Identity-System.A).fullPivLu().inverse())*System.B;
       double val=(Gjw.jacobiSvd().singularValues()(0));
        
        return val;
        

    };

    // Compute induced 1-norm of a matrix
    static double induced_one_norm(const Discrete_StateSpace_System& System,double freq){
    
        Eigen::MatrixXcd Identity=Eigen::MatrixXcd::Identity(System.A.rows(),System.A.cols());
        double w=2*M_PI*freq;
        std::complex<double> jw=std::complex<double>(0, w);
        Eigen::MatrixXcd Gjw;
        Gjw=System.C*((jw*Identity-System.A).fullPivLu().inverse())*System.B;
        int cols=Gjw.cols();
        int rows=Gjw.rows();
        double max_col_sum = 0.0;


        for (int j = 0; j < cols; ++j) {
            double col_sum = 0.0;
            for (int i = 0; i < rows; ++i) {
                col_sum += std::abs(Gjw(i, j)); // abs for complex or real
            }
            if (col_sum > max_col_sum) {
                max_col_sum = col_sum;
            }
    }
    return max_col_sum;
    } 

    // Compute induced infinity norm of a matrix
    static double induced_infinity_norm(const Discrete_StateSpace_System& System,const double& frequency){

        Eigen::MatrixXcd Identity=Eigen::MatrixXcd::Identity(System.A.rows(),System.A.cols());
        double w=2*M_PI*frequency;
        std::complex<double> jw=std::complex<double>(0, w);
        Eigen::MatrixXcd Gjw;
        Gjw=System.C*((jw*Identity-System.A).fullPivLu().inverse())*System.B;
        Eigen::MatrixXd Gjwabs=Gjw.cwiseAbs();

        double Max_R_sum=0;
        double curr_R_sum=0;
        for(int i=0;i<Gjw.rows();i++){

             curr_R_sum=Gjwabs.row(i).sum();

             if(Max_R_sum<curr_R_sum){
                Max_R_sum=curr_R_sum;
             }
        }
        return Max_R_sum;

    }

    static double H_inf_norm(const Discrete_StateSpace_System& System,
                                double low_guess,
                                double high_guess){

    double gamma_min=low_guess;
    double gamma_max=high_guess;
    double tol=10e-6;
    double iter =1000;
    double gamma_test = (gamma_min + gamma_max) / 2.0;
    double k=0;
    
while ((gamma_max - gamma_min) > tol && k < iter){
        gamma_test = (gamma_min + gamma_max) / 2.0;


            int n = System.A.rows();  // state dimension

    Eigen::MatrixXd Q = System.C.transpose() * System.C;
    Eigen::MatrixXd R = System.D.transpose() * System.D - gamma_test * gamma_test * Eigen::MatrixXd::Identity(System.D.cols(), System.D.cols());

    Eigen::MatrixXd R_inv = R.inverse();

    Eigen::MatrixXd Acl = System.A - System.B * R_inv * System.D.transpose() * System.C;
    Eigen::MatrixXd S = System.B * R_inv * System.B.transpose();
    Eigen::MatrixXd T = Q - System.C.transpose() * System.D * R_inv * System.D.transpose() * System.C;

    Eigen::MatrixXd H(2 * n, 2 * n);
    H.topLeftCorner(n, n)     = Acl;
    H.topRightCorner(n, n)    = -S;
    H.bottomLeftCorner(n, n)  = -T;
    H.bottomRightCorner(n, n) = -Acl.transpose();
    bool has_unit_circle_eigenvalue=false;
   
    Eigen::EigenSolver<Eigen::MatrixXd> es(H);
            
    for (int i = 0; i < es.eigenvalues().size(); ++i) {
        std::complex<double> lambda = es.eigenvalues()(i);
        if (std::abs(std::abs(lambda) - 1.0) < tol) {
            has_unit_circle_eigenvalue = true;
            break;
            }

        if (has_unit_circle_eigenvalue) {
        gamma_min = gamma_test;                } 
                else {
                    gamma_max = gamma_test;
                }
        
    }
    k++;
    if(k==iter){
            break;
        }
}
       return gamma_max;
 }


 static double H_inf_norm_discrete(const Discrete_StateSpace_System& System,
                                   double gamma_min_init,
                                   double gamma_max_init,
                                   double tol = 1e-6,
                                   int max_iter = 100) {
    using namespace Eigen;

    int n = System.A.rows();
    int m = System.B.cols();
    int p = System.C.rows();

    double gamma_min = gamma_min_init;
    double gamma_max = gamma_max_init;
    double gamma_test = 0.0;

    for (int iter = 0; iter < max_iter && (gamma_max - gamma_min) > tol; ++iter) {
        gamma_test = 0.5 * (gamma_min + gamma_max);

        MatrixXd R = System.D.transpose() * System.D - gamma_test * gamma_test * MatrixXd::Identity(m, m);

        FullPivLU<MatrixXd> lu(R);
        if (!lu.isInvertible()) {
            gamma_min = gamma_test;
            continue;
        }

        MatrixXd R_inv = R.inverse();
        MatrixXd Q = System.C.transpose() * System.C;

        MatrixXd Acl = System.A - System.B * R_inv * System.D.transpose() * System.C;
        MatrixXd S = System.B * R_inv * System.B.transpose();
        MatrixXd T = Q - System.C.transpose() * System.D * R_inv * System.D.transpose() * System.C;

        MatrixXd H(2 * n, 2 * n);
        H.topLeftCorner(n, n) = Acl;
        H.topRightCorner(n, n) = -S;
        H.bottomLeftCorner(n, n) = -T;
        H.bottomRightCorner(n, n) = Acl.transpose();

        EigenSolver<MatrixXd> es(H);
        bool has_unit_circle_eigenvalue = false;

        for (int i = 0; i < es.eigenvalues().size(); ++i) {
            double abs_val = std::abs(es.eigenvalues()[i]);
            if (std::abs(abs_val - 1.0) < tol) {
                has_unit_circle_eigenvalue = true;
                break;
            }
        }

        if (has_unit_circle_eigenvalue) {
            gamma_min = gamma_test;
        } else {
            gamma_max = gamma_test;
        }
    }

    return gamma_max;
}



    // Compute L2 norm of impulse response (energy)
    static double impulse_response_L2_norm(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& C) {
        // Example implementation: computes the L2 norm of the impulse response
        // For a simple system, this could be approximated as sqrt(trace(C * P * C^T)), where P solves the Lyapunov equation
        Eigen::MatrixXd Q = B * B.transpose();
        Eigen::MatrixXd P = Norms::solve_lyapunov(A, Q);
        return std::sqrt((C * P * C.transpose()).trace());
    }

    // Compute condition number w.r.t. a given norm (default spectral norm)
    static double condition_number(const Eigen::MatrixXd& M);

    // Compute Hankel singular values for balanced truncation
    static Eigen::VectorXd hankel_singular_values(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& C);

    // Compute system norm bounds (returns lower and upper bounds)
    static std::pair<double, double> system_norm_bounds(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& C);

    // Compute the Hankel norm of a stable continuous-time system
    static double hankel_norm(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& C);

};

#endif // NORMS_HPP
