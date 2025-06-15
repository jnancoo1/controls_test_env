#ifndef NORMS_HPP
#define NORMS_HPP

#include <Eigen/Dense>
#include "discrete_state_space.hpp"
#include <cmath>
#include <eigen3/unsupported/Eigen/KroneckerProduct>


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
    static double frobenius_norm(const Eigen::MatrixXd& M)}{double sum = 0.0;
    for (int i = 0; i < M.rows(); ++i) {
        for (int j = 0; j < M.cols(); ++j) {
            sum += M(i, j) * M(i, j);
        }
    }
    return std::sqrt(sum);
};

    // Compute spectral norm (2-norm) of a matrix
    static double spectral_norm(const Eigen::MatrixXd& M);

    // Compute induced 1-norm of a matrix
    static double induced_one_norm(const Eigen::MatrixXd& M);

    // Compute induced infinity norm of a matrix
    static double induced_infinity_norm(const Eigen::MatrixXd& M);

    // Approximate H-infinity norm (upper bound estimate)
    static double Hinf_norm_approximate(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& C);

    // Compute L2 norm of impulse response (energy)
    static double impulse_response_L2_norm(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& C);

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
