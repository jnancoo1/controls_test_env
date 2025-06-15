/**
 * @file discrete_state_space.hpp
 * @brief Discrete-time state-space system implementation
 * @details Provides a class for simulating discrete-time state-space systems,
 *          including methods for state updates, output computation, and sequence simulation.
 */

#ifndef DISCRETE_STATE_SPACE
#define DISCRETE_STATE_SPACE

#include <iostream>
#include <cmath>
#include <Eigen/Dense>

/**
 * @class Discrete_StateSpace_System
 * @brief Represents a discrete-time state-space system
 * @details This class models a system in the form:
 *          x[k+1] = Ax[k] + Bu[k]
 *          y[k] = Cx[k] + Du[k]
 */
class Discrete_StateSpace_System {
public:
    Eigen::MatrixXd A, B, C, D; ///< State-space matrices
    int n_states;              ///< Number of states
    int n_inputs;              ///< Number of inputs
    int n_outputs;             ///< Number of outputs

    /**
     * @brief Constructor for the state-space system
     * @param n_states_ Number of states
     * @param n_inputs_ Number of inputs
     * @param n_outputs_ Number of outputs
     */
    Discrete_StateSpace_System(int n_states_, int n_inputs_, int n_outputs_)
        : n_states(n_states_), n_inputs(n_inputs_), n_outputs(n_outputs_) {
        A = Eigen::MatrixXd::Zero(n_states, n_states);
        B = Eigen::MatrixXd::Zero(n_states, n_inputs);
        C = Eigen::MatrixXd::Zero(n_outputs, n_states);
        D = Eigen::MatrixXd::Zero(n_outputs, n_inputs);
    }

    /**
     * @brief Simulates one time step of the system
     * @param x Current state vector
     * @param u Input vector
     * @return Next state vector
     * @throw std::invalid_argument if input dimensions are invalid
     */
    Eigen::VectorXd simulate_step(const Eigen::VectorXd& x, const Eigen::VectorXd& u) {
        if (x.size() == A.rows() && u.size() == B.cols()) {
            Eigen::VectorXd X_next = A * x + B * u;
            return X_next;
        }
        throw std::invalid_argument("Invalid inputs");
    }

    /**
     * @brief Computes the output of the system
     * @param x Current state vector
     * @param u Input vector
     * @return Output vector
     * @throw std::invalid_argument if input dimensions are invalid
     */
    Eigen::VectorXd get_output(const Eigen::VectorXd& x, const Eigen::VectorXd& u) {
        if (x.size() == A.rows() && u.size() == B.cols()) {
            Eigen::VectorXd output = C * x + D * u;
            return output;
        }
        throw std::invalid_argument("Invalid inputs");
    }

    /**
     * @brief Simulates a sequence of outputs over multiple time steps
     * @param x_init Initial state vector
     * @param no_of_steps Number of time steps to simulate
     * @param u_seq Input sequence matrix (columns represent inputs at each step)
     * @param n_outputs Number of outputs
     * @return Matrix of outputs (columns represent outputs at each step)
     * @throw std::invalid_argument if input dimensions are invalid
     */
    Eigen::MatrixXd sim_seq_outputs(const Eigen::VectorXd& x_init, const int& no_of_steps, const Eigen::MatrixXd& u_seq, const int& n_outputs) {
        Eigen::VectorXd X_next = x_init;
        Eigen::MatrixXd Out(n_outputs, no_of_steps);
        if (u_seq.cols() == no_of_steps && u_seq.rows() == B.cols()) {

            for (int k = 0; k < no_of_steps; k++) {

                Out.col(k) = get_output(X_next, u_seq.col(k));
                X_next = simulate_step(X_next, u_seq.col(k));
            }

            return Out;
        }
        throw std::invalid_argument("Invalid inputs");
    }

    /**
     * @brief Simulates a sequence of states over multiple time steps
     * @param x_init Initial state vector
     * @param no_of_steps Number of time steps to simulate
     * @param u_seq Input sequence matrix (columns represent inputs at each step)
     * @param n_states Number of states
     * @return Matrix of states (columns represent states at each step)
     * @throw std::invalid_argument if input dimensions are invalid
     */
    Eigen::MatrixXd sim_seq_states(const Eigen::VectorXd& x_init, const int& no_of_steps, const Eigen::MatrixXd& u_seq, const int& n_states) {
        int i = no_of_steps;
        int checker = 0;
        Eigen::VectorXd X_next = x_init;
        Eigen::MatrixXd Out(n_states, no_of_steps);
        if (u_seq.cols() == no_of_steps && u_seq.rows() == B.cols()) {

            for (int k = 0; k < no_of_steps; k++) {

                Out.col(k) = X_next;
                X_next = simulate_step(X_next, u_seq.col(k));
            }

            return Out;
        }
        throw std::invalid_argument("Invalid inputs");
    };

    /**
     * @brief Performs a single step of simulation using the RK4 method
     * @param x Current state vector
     * @param u Input vector
     * @param dt Time step size
     * @return Next state vector
     */
    Eigen::VectorXd simulate_rk4_step(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double& dt) {
        Eigen::VectorXd k1, k2, k3, k4, x_next;

        k1 = A * x + B * u;
        k2 = A * (x + 0.5 * dt * k1) + B * u;
        k3 = A * (x + 0.5 * dt * k2) + B * u;
        k4 = A * (x + dt * k3) + B * u;
        x_next = x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        return x_next;
    }

    /**
     * @brief Simulates a sequence of states using the RK4 method
     * @param x_init Initial state vector
     * @param u_seq Input sequence matrix (columns represent inputs at each step)
     * @param dt Time step size
     * @param n_steps Number of time steps to simulate
     * @return Matrix of states (columns represent states at each step)
     * @throw std::invalid_argument if input dimensions are invalid
     */
    Eigen::MatrixXd sim_seq_rk4_states(const Eigen::VectorXd& x_init,
        const Eigen::MatrixXd& u_seq,
        double dt,
        int n_steps) {
        int n_states = A.rows();
        Eigen::MatrixXd X_seq(n_states, n_steps);
        Eigen::VectorXd x = x_init;

        if (u_seq.cols() != n_steps || u_seq.rows() != B.cols()) {
            throw std::invalid_argument("Invalid input sequence size.");
        }

        for (int k = 0; k < n_steps; ++k) {
            X_seq.col(k) = x;
            x = simulate_rk4_step(x, u_seq.col(k), dt);
        }

        return X_seq;
    }

    /**
     * @brief Performs a single step of simulation using the Euler method
     * @param x Current state vector
     * @param u Input vector
     * @param dt Time step size
     * @return Next state vector
     */
    Eigen::VectorXd sim_euler_step(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double& dt) {
        Eigen::VectorXd fxk;
        fxk = A * x + B * u;
        Eigen::VectorXd x_kp = x + dt * fxk;
        return x_kp;
    }

    /**
     * @brief Simulates a sequence of states using the Euler method
     * @param x_init Initial state vector
     * @param u_seq Input sequence matrix (columns represent inputs at each step)
     * @param dt Time step size
     * @param n_steps Number of time steps to simulate
     * @return Matrix of states (columns represent states at each step)
     * @throw std::invalid_argument if input dimensions are invalid
     */
    Eigen::MatrixXd sim_seq_euler_states(const Eigen::VectorXd& x_init,
        const Eigen::MatrixXd& u_seq,
        double dt,
        int n_steps) {
        int n_states = A.rows();
        Eigen::MatrixXd X_seq(n_states, n_steps);
        Eigen::VectorXd x = x_init;

        if (u_seq.cols() != n_steps || u_seq.rows() != B.cols()) {
            throw std::invalid_argument("Invalid input sequence size.");
        }

        for (int k = 0; k < n_steps; ++k) {
            X_seq.col(k) = x;
            x = sim_euler_step(x, u_seq.col(k), dt);
        }

        return X_seq;
    }

    /**
     * @brief Performs a single step of simulation using the Crank-Nicolson method
     * @param x Current state vector
     * @param u_k Current input vector
     * @param u_k1 Next input vector
     * @param dt Time step size
     * @return Next state vector
     */
    Eigen::VectorXd sim_crank_nelson_step(const Eigen::VectorXd& x, const Eigen::VectorXd& u_k, const Eigen::VectorXd& u_k1, const double& dt) {

        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(A.rows(), A.cols());
        Eigen::MatrixXd M = (I - (dt / 2) * A);

        Eigen::VectorXd RHS = (I + (dt / 2) * A) * x + ((dt / 2) * B) * (u_k + u_k1);
        Eigen::VectorXd output = M.colPivHouseholderQr().solve(RHS);

        return output;
    }

    /**
     * @brief Performs a single step of simulation using the backward Euler method
     * @param x Current state vector
     * @param u Input vector
     * @param dt Time step size
     * @return Next state vector
     */
    Eigen::VectorXd sim_backward_euler_step(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double& dt) {
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(A.rows(), A.cols());
        Eigen::MatrixXd M = (I - dt * A);
        Eigen::VectorXd RHS = x + dt * B * u;
        return M.colPivHouseholderQr().solve(RHS);
    }

    /**
     * @brief Performs a single step of simulation using the Heun's method
     * @param x Current state vector
     * @param u Input vector
     * @param dt Time step size
     * @return Next state vector
     */
    Eigen::VectorXd sim_heun_step(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const double& dt) {
        Eigen::VectorXd k1 = A * x + B * u;
        Eigen::VectorXd x_pred = x + dt * k1;
        Eigen::VectorXd k2 = A * x_pred + B * u;
        return x + dt * 0.5 * (k1 + k2);
    }

private:

};

#endif
