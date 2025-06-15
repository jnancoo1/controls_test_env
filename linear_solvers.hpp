/**
 * @file linear_solvers.hpp
 * @brief Linear system solvers implementation
 * @details Contains implementations of various methods for solving linear systems,
 *          including LU decomposition, QR decomposition, and direct matrix inversion.
 */

#ifndef LINEAR_SOLVER_HPP
#define LINEAR_SOLVER_HPP

#include <iostream>
#include <cmath>
#include "matrix_math.hpp"

/**
 * @class Linear_Solvers
 * @brief Static class providing methods for solving linear systems
 * @details Implements various numerical methods for solving linear systems including:
 *          - LU decomposition based solver
 *          - QR decomposition based solver
 *          - Matrix inversion
 *          - Forward and backward substitution
 */
class Linear_Solvers {
public:
    /**
     * @brief Solves a linear system using LU decomposition
     * @param A The coefficient matrix
     * @param b The right-hand side vector
     * @return Solution vector x where Ax = b
     */
    static My_Vec SolveLU(const Matrix& A, const My_Vec& b) {

        LUResult LU_temp=A.L_U();
        LUResult_to_pass LU=conv_LU(LU_temp);
        My_Vec pb= ApplyPermutation(LU.P,b);
        My_Vec fwd_sub=ForwardSubstitution(LU.L,pb);
        My_Vec bck_sub=BackwardSubstitution(LU.U,fwd_sub);

        return bck_sub;

    }
    /**
     * @brief Solves a linear system using QR decomposition
     * @param A The coefficient matrix
     * @param b The right-hand side vector
     * @return Solution vector x where Ax = b
     */
    static My_Vec SolveQR(const Matrix& A, const My_Vec& b) {

        QRresult decomp_temp=A.QR_fact();
        QR_result_to_pass decomp=conv_QR(decomp_temp);
        My_Vec qTb=(decomp.Q.Transpose()).multiply(b);
        My_Vec x=BackwardSubstitution(decomp.R,qTb);
        return x;
    };
    /**
     * @brief Computes the inverse of a matrix using LU decomposition
     * @param A The matrix to invert
     * @return The inverse of matrix A
     * @throw std::runtime_error if matrix is singular
     */
    static Matrix Inverse(const Matrix& A) {

        Matrix I=Matrix::eye(A.rows);
        LUResult decomp_temp=A.L_U();
        LUResult_to_pass decomp=conv_LU(decomp_temp);
        Matrix L=decomp.L;
        Matrix U=decomp.U;
        My_Vec e_i = My_Vec::ones(A.rows);
        Matrix inverse_matrix=Matrix::Zeros(A.rows,A.cols);        
        for(int j=0;j<I.cols;j++){
            e_i.Scalar_Mul(0);
            e_i.myvector[j] = 1;
            My_Vec Pe_j = ApplyPermutation(decomp.P, e_i);

            My_Vec F1=ForwardSubstitution(L,Pe_j);
            My_Vec B1=BackwardSubstitution(U,F1);
            for (int row = 0; row < A.rows; row++) {
                inverse_matrix.MyMAT[row][j] = B1.myvector[row];
            }
        }

        return inverse_matrix;

    };
    /**
     * @brief Performs forward substitution to solve Lx = b
     * @param L_1 Lower triangular matrix
     * @param b Right-hand side vector
     * @return Solution vector x
     */
    static My_Vec ForwardSubstitution(const Matrix& L_1, const My_Vec& b) {

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
    /**
     * @brief Performs backward substitution to solve Ux = b
     * @param U_1 Upper triangular matrix
     * @param b Right-hand side vector
     * @return Solution vector x
     */
    static My_Vec BackwardSubstitution(const Matrix& U_1, const My_Vec& b) {
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
    /**
     * @brief Applies a permutation to a vector
     * @param P Permutation vector
     * @param V Vector to permute
     * @return Permuted vector
     * @throw std::invalid_argument if permutation size doesn't match vector size
     */
     static My_Vec ApplyPermutation(const std::vector<int>& P,const My_Vec& V){
        if (P.size() != V.myvector. size()) {
            throw std::invalid_argument("Permutation size must match vector size");
        }
        My_Vec result;

        My_Vec result = My_Vec::Zeros(V.myvector.size());
        
        for (size_t i = 0; i < P.size(); i++) {
            result.myvector[i] = V.myvector[P[i]];
        }
        
        return result;
    }

    /**
     * @brief Computes the determinant of a matrix using LU decomposition
     * @param A Input matrix
     * @return Determinant value
     */
    static double determinant(const Matrix& A) {
    LUResult lu_temp = A.L_U();
    LUResult_to_pass lu=conv_LU(lu_temp);
    double det = 1.0;
    for(int i = 0; i < A.rows; i++) {
        det *= lu.U.MyMAT[i][i];
    }
    return det;
}


    /**
     * @brief Alias for SolveLU for solving linear systems
     * @param A The coefficient matrix
     * @param b The right-hand side vector
     * @return Solution vector x where Ax = b
     */
    static My_Vec solve_linear_system_LU(const Matrix& A, const My_Vec& b) {
    return SolveLU(A, b);  
}

    /**
     * @brief Alias for SolveQR for solving linear systems
     * @param A The coefficient matrix
     * @param b The right-hand side vector
     * @return Solution vector x where Ax = b
     */
    static My_Vec solve_linear_system_QR(const Matrix& A, const My_Vec& b) {
    return SolveQR(A, b);  
}

    };
#endif