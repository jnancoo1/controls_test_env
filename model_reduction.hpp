#ifndef MODEL_REDUCTION_HPP
#define MODEL_REDUCTION_HPP

#include <Eigen/Dense>
#include "discrete_state_space.hpp"
#include <vector>
#include <complex>
#include "analysis.hpp"

class ModelReduction {
public:

    /**
     * Balanced Truncation
     * - Uses Gramians to identify energy-dominant states.
     * - Preserves stability and provides error bounds.
     */
    static Discrete_StateSpace_System balanced_truncation(const Discrete_StateSpace_System& System, int reduced_order);

    /**
     * Hankel Norm Approximation
     * - Minimizes the Hankel norm of the error system.
     */
    static Discrete_StateSpace_System hankel_norm_approximation(const Discrete_StateSpace_System& System, int reduced_order);

    /**
     * Modal Truncation
     * - Removes modes with fast decay or low observability.
     */
    static Discrete_StateSpace_System modal_truncation(const Discrete_StateSpace_System& System, int reduced_order);

    /**
     * Singular Perturbation Approximation
     * - Truncates fast modes based on a time-scale threshold.
     */
    static Discrete_StateSpace_System singular_perturbation_approx(const Discrete_StateSpace_System& System, double time_scale_threshold);

    /**
     * Proper Orthogonal Decomposition (POD)
     * - Data-driven method based on SVD of snapshots.
     */
    static Discrete_StateSpace_System pod(const Eigen::MatrixXd& snapshot_matrix,
                                          const Discrete_StateSpace_System& full_model,
                                          int reduced_order);

    /**
     * Krylov Subspace Reduction (Arnoldi / Lanczos)
     * - Good for large-scale, sparse systems.
     */
    static Discrete_StateSpace_System krylov_reduction(const Discrete_StateSpace_System& System, int reduced_order);

    /**
     * Balanced Stochastic Truncation
     * - Minimizes stochastic approximation error.
     */
    static Discrete_StateSpace_System balanced_stochastic_truncation(const Discrete_StateSpace_System& System, int reduced_order);

    /**
     * Loewner Framework Reduction
     * - Data-driven method for frequency-domain interpolation.
     */
    static Discrete_StateSpace_System loewner_reduction(const std::vector<std::complex<double>>& frequencies,
                                                        const std::vector<Eigen::MatrixXcd>& responses,
                                                        int reduced_order);

    /**
     * Empirical Gramian-based Reduction
     * - Model-free approach using snapshot data.
     */
    static Discrete_StateSpace_System empirical_gramian_reduction(const std::vector<Eigen::VectorXd>& snapshots,
                                                                  const Discrete_StateSpace_System& full_model,
                                                                  int reduced_order);

    /**
     * Compute Hankel Singular Values
     * - Indicates state energy contributions; aids in selecting reduced order.
     */
    static Eigen::VectorXd hankel_singular_values(const Discrete_StateSpace_System& System);

    /**
     * Select reduced order based on cumulative energy threshold (e.g., 0.99).
     */
    static int select_order_from_energy(const Eigen::VectorXd& hsv, double energy_threshold);

    /**
     * Compute balancing transformation matrices (T, T_inv) from Gramians.
     */
    static std::pair<Eigen::MatrixXd, Eigen::MatrixXd> compute_balancing_transforms(const Eigen::MatrixXd& Wc,
                                                                                     const Eigen::MatrixXd& Wo);




        // Interpolatory methods (e.g., IRKA)
    static std::pair<Eigen::MatrixXd, Eigen::MatrixXd> irka(const Discrete_StateSpace_System& sys, int max_iter, double tol);

    // H-infinity norm reduction (advanced, optimization-based)
    static std::pair<Eigen::MatrixXd, Eigen::MatrixXd> hinf_norm_reduction(const Discrete_StateSpace_System& sys, double tolerance);

private:
    ModelReduction() = default;
};

#endif // MODEL_REDUCTION_HPP
