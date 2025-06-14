#ifndef FORMS_HPP
#define FORMS_HPP

#include <iostream>
#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include "analysis.hpp"
#include "discrete_state_space.hpp"
#include <Eigen/SVD>
#include <eigen3/unsupported/Eigen/KroneckerProduct>



class Forms{

    public:

    Discrete_StateSpace_System Cont_Cannonical_form(const Discrete_StateSpace_System& System)
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
        Discrete_StateSpace_System new_system=System;
        new_system.A = Accf;
        new_system.B = T_inv * System.B;
        new_system.C = System.C * T;
        new_system.D = System.D;




        return new_system;
    }


        Discrete_StateSpace_System obs_Cannonical_form(const Discrete_StateSpace_System& System){

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
        Eigen::MatrixXd T = A.compute_observability_matrix(System);
        Eigen::MatrixXd T_inv=T.inverse();

        // Ensure the constructor matches the expected signature, e.g. (A, B, C, D)
        Discrete_StateSpace_System new_system=System;
        new_system.A = Aocf;
        new_system.B = T_inv * System.B;
        new_system.C = System.C * T;
        new_system.D = System.D;

        return new_system;


     }
    private:


};
