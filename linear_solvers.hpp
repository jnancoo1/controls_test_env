
#ifndef LINEAR_SOLVER_HPP
#define LINEAR_SOLVER_HPP


#include <iostream>
#include <cmath>
#include "matrix_math.hpp"



class Linear_Solvers{
public:
    static My_Vec SolveLU(const Matrix& A, const My_Vec& b){

        LUResult LU=A.L_U();
        My_Vec fwd_sub=ForwardSubstitution(LU.L,b);
        My_Vec bck_sub=BackwardSubstitution(LU.U,fwd_sub);

        return bck_sub;

    }
    static My_Vec SolveQR(const Matrix& A, const My_Vec& b){

        QRresult decomp=A.QR_fact();
        My_Vec qTb=(decomp.Q.Transpose()).multiply(b);
        My_Vec x=BackwardSubstitution(decomp.R,qTb);
        return x;
    };
    static Matrix Inverse(const Matrix& A){

        Matrix I=Matrix::eye(A.rows);
        LUResult decomp=A.L_U();
        Matrix L=decomp.L;
        Matrix U=decomp.U;
        My_Vec e_i = My_Vec::ones(A.rows);
        Matrix inverse_matrix=Matrix::Zeros(A.rows,A.cols);        
        for(int j=0;j<I.cols;j++){
            e_i.Scalar_Mul(0);
            e_i.myvector[j] = 1;

            My_Vec F1=ForwardSubstitution(L,e_i);
            My_Vec B1=BackwardSubstitution(U,F1);
            for (int row = 0; row < A.rows; row++) {
                inverse_matrix.MyMAT[row][j] = B1.myvector[row];
            }
        }

        return inverse_matrix;

    };
    static My_Vec ForwardSubstitution(const Matrix& L_1, const My_Vec& b){

       My_Vec solution_VEC=My_Vec::ones(L_1.rows);
       solution_VEC.Scalar_Mul(0);

       Matrix L= L_1;
       std::vector<double> knowns;
       double new_unknown;
       for(int i=0;i<L.rows;i++){

        double known_sum=0;
        for(int j=0;j<i;j++){
            known_sum+=L.MyMAT[i][j]*solution_VEC.myvector[j];
        }
        new_unknown=(b.myvector[i]-known_sum)/L.MyMAT[i][i];


        solution_VEC.myvector[i]=new_unknown;


       }

       return solution_VEC;


    };
    static My_Vec BackwardSubstitution(const Matrix& U_1, const My_Vec& b){
        My_Vec solution_VEC=My_Vec::ones(U_1.rows);
        solution_VEC.Scalar_Mul(0);

        Matrix U= U_1;
        std::vector<double> knowns;
        double new_unknown;
        for(int i=U.rows-1;i>=0;i--){
 
         double known_sum=0;
         for (int j = i + 1; j < U.cols; j++){
             known_sum+=U.MyMAT[i][j]*solution_VEC.myvector[j];
         }
         new_unknown=(b.myvector[i]-known_sum)/U.MyMAT[i][i];
 
 
         solution_VEC.myvector[i]=new_unknown;
 
 
        }
 
        return solution_VEC;
 
 
     };
     static My_Vec ApplyPermutation(const std::vector<int>& P,const My_Vec& V){
            
        if (P.size() != V.myvector.size()) {
            throw std::invalid_argument("Permutation size must match vector size");
        }

        My_Vec result = My_Vec::zeros(v.myvector.size());
        
        for (size_t i = 0; i < P.size(); i++) {
            result.myvector[i] = v.myvector[P[i]];
        }
        
        return result;
    }
    };
#endif