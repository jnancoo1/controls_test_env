#ifndef NORMS_HPP
#define NORMS_HPP

#include <Eigen/Dense>

class Norms {
public:
    // Compute H2 norm for stable continuous-time systems: sqrt(trace(C * P * C^T))
    static double H2_norm_continuous(const Eigen::MatrixXd& A, const Eigen::MatrixXd& C);

    // Compute H2 norm for stable discrete-time systems
    static double H2_norm_discrete(const Eigen::MatrixXd& A, const Eigen::MatrixXd& C);

    // Solve continuous Lyapunov equation A'P + PA + Q = 0
    static Eigen::MatrixXd solve_lyapunov(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Q);

    // Solve discrete Lyapunov equation P = A*P*A' + Q
    static Eigen::MatrixXd solve_discrete_lyapunov(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Q);

    // Check if continuous-time system matrix A is stable (all eigenvalues real part < 0)
    static bool is_stable_continuous(const Eigen::MatrixXd& A);

    // Check if discrete-time system matrix A is stable (all eigenvalues magnitude < 1)
    static bool is_stable_discrete(const Eigen::MatrixXd& A);

    // Compute Frobenius norm of a matrix
    static double frobenius_norm(const Eigen::MatrixXd& M);

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
