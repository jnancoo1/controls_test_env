#ifndef MODES_HPP
#define MODES_HPP

#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <utility>

class Modes {
public:
    // Perform modal decomposition: eigenvalues and right eigenvectors
    static std::pair<Eigen::VectorXcd, Eigen::MatrixXcd> modal_decomposition(const Eigen::MatrixXd& A);

    // Compute left eigenvectors of matrix A
    static Eigen::MatrixXcd left_eigenvectors(const Eigen::MatrixXd& A);

    // Compute damping ratios from eigenvalues
    static Eigen::VectorXd damping_ratios(const Eigen::VectorXcd& eigvals);

    // Compute natural frequencies from eigenvalues
    static Eigen::VectorXd natural_frequencies(const Eigen::VectorXcd& eigvals);

    // Compute participation factors given right and left eigenvectors
    static Eigen::MatrixXd participation_factors(const Eigen::MatrixXcd& right_eigvecs, const Eigen::MatrixXcd& left_eigvecs);

    // Modal controllability indices for each mode
    static Eigen::VectorXd modal_controllability_indices(const Eigen::MatrixXd& B, const Eigen::MatrixXcd& right_eigvecs);

    // Modal observability indices for each mode
    static Eigen::VectorXd modal_observability_indices(const Eigen::MatrixXd& C, const Eigen::MatrixXcd& left_eigvecs);

    // Sort modes by criteria: frequency, damping ratio, participation factor, etc.
    static void sort_modes(
        Eigen::VectorXcd& eigvals,
        Eigen::MatrixXcd& right_eigvecs,
        const Eigen::VectorXd* criteria = nullptr);

    // Normalize mode shapes (eigenvectors), e.g. max component = 1
    static void normalize_modes(Eigen::MatrixXcd& eigvecs);

    // Compute modal residues for frequency response analysis
    static std::vector<Eigen::MatrixXcd> modal_residues(
        const Eigen::MatrixXd& A,
        const Eigen::MatrixXd& B,
        const Eigen::MatrixXd& C);

    // Time-domain contribution of modes given initial state vector
    static Eigen::MatrixXd mode_time_response(
        const Eigen::VectorXcd& eigvals,
        const Eigen::MatrixXcd& right_eigvecs,
        const Eigen::VectorXcd& left_eigvecs,
        const Eigen::VectorXd& x0,
        const std::vector<double>& time_points);

    // Generate reduced order model by selecting dominant modes
    static void mode_truncation(
        const Eigen::VectorXcd& eigvals,
        const Eigen::MatrixXcd& right_eigvecs,
        const Eigen::MatrixXcd& left_eigvecs,
        int num_modes,
        Eigen::VectorXcd& truncated_eigvals,
        Eigen::MatrixXcd& truncated_right_eigvecs,
        Eigen::MatrixXcd& truncated_left_eigvecs);

    // Convert complex conjugate modes into real modal decomposition
    static void real_modal_decomposition(
        const Eigen::VectorXcd& eigvals,
        const Eigen::MatrixXcd& right_eigvecs,
        Eigen::MatrixXd& real_modes_out);
};

#endif // MODES_HPP
