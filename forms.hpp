#ifndef FORMS_HPP
#define FORMS_HPP

#include <iostream>
#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include "discrete_state_space.hpp"
#include "analysis.hpp"
#include <Eigen/SVD>
#include <eigen3/unsupported/Eigen/KroneckerProduct>

struct StateSpace_System {

    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    Eigen::MatrixXd C;
    Eigen::MatrixXd D;

};

/**
 * @brief A class for transforming state space systems into various canonical forms
 *
 * This class provides methods to transform discrete-time state space systems
 * into different canonical forms. Each transformation is achieved through similarity
 * transformations of the form:
 * - A_new = T^(-1)AT
 * - B_new = T^(-1)B
 * - C_new = CT
 * where T is the transformation matrix specific to each form.
 */
class Forms{

    public:
    /**
     * @brief Transforms a state space system into controllable canonical form
     * 
     * The controllable canonical form transforms the system into the form:
     * \f[
     * A = \begin{bmatrix} 
     * 0 & 1 & 0 & \cdots & 0 \\
     * 0 & 0 & 1 & \cdots & 0 \\
     * \vdots & \vdots & \vdots & \ddots & \vdots \\
     * 0 & 0 & 0 & \cdots & 1 \\
     * -a_0 & -a_1 & -a_2 & \cdots & -a_{n-1}
     * \end{bmatrix}
     * \f]
     * 
     * where \f$a_i\f$ are the coefficients of the characteristic polynomial:
     * \f[
     * p(s) = s^n + a_{n-1}s^{n-1} + \cdots + a_1s + a_0
     * \f]
     * 
     * The transformation uses the controllability matrix:
     * \f[
     * T = [B \quad AB \quad A^2B \quad \cdots \quad A^{n-1}B]
     * \f]
     * 
     * @param System The discrete state space system to transform
     * @return StateSpace_System The transformed system in controllable canonical form
     */
    StateSpace_System Cont_Cannonical_form(const StateSpace_System& System)
    {


        int n = System.A.rows();
        Eigen::MatrixXd B0 = Eigen::MatrixXd::Identity(n, n);
        Eigen::MatrixXd Bk;
        std::vector<double> coeffs;
        double ak_1 = 0;
        double trk = 0;

        for (int k = 1; k <= n; k++) {
            if (k == 1) {
                // First iteration: trace(A), a0 coefficient, initialize Bk
                trk = System.A.diagonal().sum();  // trace(A)
                ak_1 = -trk / k;
                coeffs.push_back(ak_1);
                Bk = B0;
            } else {
                // Subsequent iterations: compute Bk, trace, and coefficient
                Eigen::MatrixXd Bk_1 = System.A * Bk + ak_1 * B0;
                trk = (System.A * Bk_1).diagonal().sum();
                ak_1 = -trk / k;
                coeffs.push_back(ak_1);
                Bk = Bk_1;
            }
        };


        Eigen::MatrixXd Accf = Eigen::MatrixXd::Zero(n, n);
        for(int i=0;i<n-1;i++){

            Accf(i,i+1)=1;

        }
        for(int j=0;j<n;j++){

            Accf(n-1,j)=-coeffs[j];

        }
       
        //Compute Controlability Matrix
        Analysis A;
        Eigen::MatrixXd T = A.compute_controllability_matrix(System);
        Eigen::MatrixXd T_inv=T.inverse();

        // Ensure the constructor matches the expected signature, e.g. (A, B, C, D)
        StateSpace_System new_system = System;
        new_system.A = Accf;
        new_system.B = T_inv * System.B;
        new_system.C = System.C * T;
        new_system.D = System.D;




        return new_system;
    }


        /**
     * @brief Transforms a state space system into observable canonical form
     * 
     * The observable canonical form transforms the system into the form:
     * \f{equation*}{
     * A = \begin{pmatrix} 
     * -a_{n-1} & -a_{n-2} & \cdots & -a_1 & -a_0 \\
     * 1 & 0 & \cdots & 0 & 0 \\
     * 0 & 1 & \cdots & 0 & 0 \\
     * \vdots & \vdots & \ddots & \vdots & \vdots \\
     * 0 & 0 & \cdots & 1 & 0
     * \end{pmatrix}
     * \f}
     * 
     * The transformation uses the observability matrix:
     * \f[
     * T = \begin{bmatrix} C \\ CA \\ CA^2 \\ \vdots \\ CA^{n-1} \end{bmatrix}
     * \f]
     * 
     * @param System The discrete state space system to transform
     * @return StateSpace_System The transformed system in observable canonical form
     */
        StateSpace_System obs_Cannonical_form(const StateSpace_System& System){

        int n = System.A.rows();
        Eigen::MatrixXd B0 = Eigen::MatrixXd::Identity(n, n);
        Eigen::MatrixXd Bk;
        std::vector<double> coeffs;
        double ak_1 = 0;
        double trk = 0;

        for (int k = 1; k <= n; k++) {
            if (k == 1) {
                // First iteration: trace(A), a0 coefficient, initialize Bk
                trk = System.A.diagonal().sum();  // trace(A)
                ak_1 = -trk / k;
                coeffs.push_back(ak_1);
                Bk = B0;
            } else {
                // Subsequent iterations: compute Bk, trace, and coefficient
                Eigen::MatrixXd Bk_1 = System.A * Bk + ak_1 * B0;
                trk = (System.A * Bk_1).diagonal().sum();
                ak_1 = -trk / k;
                coeffs.push_back(ak_1);
                Bk = Bk_1;
            }
        };

        Eigen::MatrixXd Aocf = Eigen::MatrixXd::Zero(n, n);
        for(int i=0;i<n-1;i++){

            Aocf(i,i+1)=1;

        }
        for(int j=0;j<n;j++){

        Aocf(j,0) = -coeffs[n - j - 1];
        }

        //Compute Controlability Matrix
        Analysis A;
        Eigen::MatrixXd T = Analysis::compute_observability_matrix(System);
        Eigen::MatrixXd T_inv=T.inverse();

        // Ensure the constructor matches the expected signature, e.g. (A, B, C, D)
        StateSpace_System new_system = System;
        new_system.A = Aocf;
        new_system.B = T_inv * System.B;
        new_system.C = System.C * T;
        new_system.D = System.D;

        return new_system;


     }


    /**
     * @brief Transforms a state space system into phase variable form
     * 
     * The phase variable form represents the system in terms of a state vector
     * containing successive derivatives (or differences in discrete-time):
     * \f[
     * x = \begin{bmatrix} y \\ \Delta y \\ \Delta^2 y \\ \vdots \\ \Delta^{n-1} y \end{bmatrix}
     * \f]
     * 
     * This form requires the system to be controllable. The transformation matrix T
     * is constructed from the controllability matrix.
     * 
     * @param System The discrete state space system to transform
     * @return StateSpace_System The transformed system in phase variable form
     * @throw std::runtime_error if the system is not controllable
     */
    StateSpace_System Phase_Variable_Form(const StateSpace_System& System) {
    int n = System.A.rows();

    // Build controllability matrix and check invertibility
    Eigen::MatrixXd T = Analysis::compute_controllability_matrix(System);
    if (T.determinant() == 0) 
        throw std::runtime_error("System is not controllable. Cannot convert to phase variable form.");
    }
    Eigen::MatrixXd T_inv = T.inverse();

    // Transform system to phase variable form
    StateSpace_System new_system = System;
    new_system.A = T_inv * System.A * T;
    new_system.B = T_inv * System.B;
    new_system.C = System.C * T;
    // D remains unchanged

    return new_system;
}


    /**
     * @brief Transforms a state space system into Schur form
     * 
     * The Schur form uses an orthogonal similarity transformation Q to transform
     * the system matrix A into an upper triangular matrix T:
     * \f[
     * T = Q^TAQ
     * \f]
     * 
     * Properties of the Schur form:
     * - T is upper triangular
     * - The diagonal elements of T are the eigenvalues of A
     * - Q is orthogonal (Q^TQ = QQ^T = I)
     * 
     * This form is numerically stable and preserves the eigenvalues of the original system.
     * 
     * @param System The discrete state space system to transform
     * @return StateSpace_System The transformed system in Schur form
     */
    StateSpace_System Schur_Form(const StateSpace_System& System){


        Eigen::RealSchur<Eigen::MatrixXd> Schur(System.A)
        Eigen::MatrixXd Q = schur.matrixU();  
        Eigen::MatrixXd T = schur.matrixT();

        Eigen::MatrixXd Q_inv=Q.transpose();

        StateSpace_System new_system = System;
        new_system.A = T;
        new_system.B = Q_inv * System.B;
        new_system.C = System.C * Q;

        return new_system
}


    /**
     * @brief Transforms a state space system into diagonal form
     * 
     * The diagonal form represents the system using the eigenvector decomposition:
     * \f[
     * A = PDP^{-1}
     * \f]
     * where:
     * - P is the matrix of eigenvectors
     * - D is a diagonal matrix of eigenvalues
     * - P^{-1} is the inverse of the eigenvector matrix
     * 
     * The transformed system has the properties:
     * - A_new = D (diagonal matrix of eigenvalues)
     * - B_new = P^{-1}B
     * - C_new = CP
     * 
     * This form decouples the system into independent first-order subsystems,
     * making it useful for stability analysis and control design.
     * 
     * @param System The discrete state space system to transform
     * @return StateSpace_System The transformed system in diagonal form
     * @throw std::runtime_error if the matrix is not diagonalizable (eigenvectors are not linearly independent)
     */
        StateSpace_System Diagonalize(const StateSpace_System& System) {
            Eigen::EigenSolver<Eigen::MatrixXd> es(System.A);
            Eigen::MatrixXcd eigenvectors = es.eigenvectors();
            Eigen::VectorXcd eigenvalues = es.eigenvalues();

            // Check if matrix is diagonalizable (eigenvectors matrix invertible)
            if (eigenvectors.determinant() == 0) {
                throw std::runtime_error("Matrix is not diagonalizable");
            }

            Eigen::MatrixXcd P = eigenvectors;
            Eigen::MatrixXcd P_inv = P.inverse();

            Eigen::MatrixXcd D = eigenvalues.asDiagonal();

            StateSpace_System new_system;

            new_system.A = (P_inv * System.A.cast<std::complex<double>>() * P).real();
            new_system.B = (P_inv * System.B.cast<std::complex<double>>()).real();
            new_system.C = (System.C.cast<std::complex<double>>() * P).real();
            new_system.D = System.D;

            return new_system;
        }



    private:




};
#endif