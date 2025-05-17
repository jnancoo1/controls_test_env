#include <iostream>
#include <cmath>
#include "matrix_math.hpp"

#ifndef LINEAR_SOLVER_HPP
#define LINEAR_SOLVER_HPP



class Linear_Solvers{

    static My_Vec SolveLU(const Matrix& A, const My_Vec& b);
    static My_Vec SolveQR(const Matrix& A, const My_Vec& b);
    static My_Vec Inverse(const Matrix& A);

    static My_Vec ForwardSubstitution(const Matrix& A, const My_Vec& b){

       My_Vec solution_VEC=My_Vec::ones(A.rows);
       solution_VEC.Scalar_Mul(0);
       LUResult decomp;
       decomp=A.L_U();

       Matrix L= decomp.L;
       std::vector<double> knowns;
       double new_unknown;
       for(int i=0;i<L.rows;i++){

        double known_sum=0;
        for(int j=0;j<i-1;j++){
            known_sum+=L.MyMAT[i][j]*solution_VEC.myvector[j];
        }
        new_unknown=(b.myvector[i]-known_sum)/L.MyMAT[i][i];


        solution_VEC.myvector[i]=new_unknown;


       }

       return solution_VEC;


    };
    static My_Vec BackwardSubstitution(const Matrix& A, const My_Vec& b){
        My_Vec solution_VEC=My_Vec::ones(A.rows);
        solution_VEC.Scalar_Mul(0);
        LUResult decomp;
        decomp=A.L_U();
 
        Matrix U= decomp.U;
        std::vector<double> knowns;
        double new_unknown;
        for(int i=U.rows;i>=0;i--){
 
         double known_sum=0;
         for (int j = i + 1; j < U.cols; j++){
             known_sum+=U.MyMAT[i][j]*solution_VEC.myvector[j];
         }
         new_unknown=(b.myvector[i]-known_sum)/U.MyMAT[i][i];
 
 
         solution_VEC.myvector[i]=new_unknown;
 
 
        }
 
        return solution_VEC;
 
 
     };
    };



};
#endif