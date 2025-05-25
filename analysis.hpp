#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "discrete_state_space.hpp"

class Analysis {
public:
    Eigen::MatrixXd compute_controllability_matrix(const Discrete_StateSpace_System& System)
    {
        int n = System.A.rows();
        int m = System.B.cols();
        Eigen::MatrixXd controllability_mat(n, n * m);

        for (int i = 0; i < n; ++i) {
            controllability_mat.block(0, i * m, n, m) = System.A.pow(i) * System.B;
        }

        return controllability_mat;
    }

    bool is_controllable(const Discrete_StateSpace_System& System) 
    {
        Eigen::MatrixXd controllability_mat = compute_controllability_matrix(System);
        Eigen::FullPivLU<Eigen::MatrixXd> lu(controllability_mat);
        return lu.rank() == System.n_states;
    }

    Eigen::MatrixXd compute_observability_matrix(const Discrete_StateSpace_System& System)
    {
        int n = System.A.rows();
        int p = System.C.rows();
        Eigen::MatrixXd observability_mat(p * n, System.A.cols());

        for (int i = 0; i < n; ++i) {
            observability_mat.block(i * p, 0, p, System.A.cols()) = System.C * System.A.pow(i);
        }

        return observability_mat;
    }

    bool is_observable(const Discrete_StateSpace_System& System) 
    {
        Eigen::MatrixXd observability_mat = compute_observability_matrix(System);
        Eigen::FullPivLU<Eigen::MatrixXd> lu(observability_mat);
        return lu.rank() == System.n_states;
    }

private:

};

#endif
