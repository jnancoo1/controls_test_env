/**
 * @file matrix_math.hpp
 * @brief Linear algebra and matrix mathematics library
 * @details This library provides implementation for basic linear algebra operations,
 *          including matrix operations, LU decomposition, and QR factorization.
 */

//gcc (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0
// Standard library includes for basic operations
#include <iostream>
#include <vector>
#include <cmath>

#ifndef MATRIX_MATH_HPP
#define MATRIX_MATH_HPP

/** 
 * @struct LUResult
 * @brief Result of LU decomposition of a matrix
 * @details Stores the lower triangular (L), upper triangular (U) matrices,
 *          and the permutation vector (P) from the decomposition
 */
struct LUResult {
    std::vector<std::vector<double>> L;    ///< Lower triangular matrix
    std::vector<std::vector<double>> U;    ///< Upper triangular matrix
    std::vector<int> P;                    ///< Permutation vector
};

/** 
 * @struct QRresult
 * @brief Result of QR decomposition of a matrix
 * @details Stores the orthogonal matrix Q and upper triangular matrix R
 */
struct QRresult {
     std::vector<std::vector<double>> Q;   ///< Orthogonal matrix
     std::vector<std::vector<double>> R;   ///< Upper triangular matrix
};

/**
 * @class My_Vec
 * @brief Vector class for mathematical operations
 * @details Implements a mathematical vector with common operations like
 *          addition, subtraction, dot product, and scalar multiplication
 */
class My_Vec {
public:
    int length;                           ///< Length of the vector
    std::vector<double> myvector;         ///< Vector data storage

    /**
     * @brief Construct a new vector
     * @param l Length of the vector (defaults to 1)
     */
    My_Vec(int l=1) {
        length = l;
        myvector.resize(length);
    }

    /**
     * @brief Add two vectors
     * @param other Vector to add to this one
     * @return New vector containing the sum
     * @throw std::invalid_argument if vectors have different lengths
     */
    My_Vec operator+ (const My_Vec& other) const {
        if(this->length != other.length) {
            throw std::invalid_argument("Invalid lengths");
        }
        My_Vec Vec(this->length);
        for (int i=0; i<this->length; i++) {
            Vec.myvector[i] = this->myvector[i] + other.myvector[i];
        }
        return Vec;  
    }

    /**
     * @brief Creates a vector of ones
     * @param a Length of the vector
     * @return Vector with all elements set to 1
     */
    static My_Vec ones(int a) {
        My_Vec ones_vec(a);
        for(int i=0; i<a; i++) {
            ones_vec.myvector[i] = 1;
        }
        return ones_vec;
    }

    /**
     * @brief Vector subtraction operator
     * @param other The vector to subtract
     * @return Result of vector subtraction
     * @throws std::invalid_argument if vectors have different lengths
     */
    My_Vec operator- (const My_Vec& other) const{

        if(this->length != other.length) {
            throw std::invalid_argument("Invalid lengths");
        }
        My_Vec Vec(this->length);

        for (int i=0; i<this->length; i++) {
            Vec.myvector[i] = this->myvector[i] - other.myvector[i];
        }
        return Vec;        
    }

    /**
     * @brief Computes the Euclidean norm (magnitude) of the vector
     * @return Norm of the vector
     */
    double Norm() const {
        double a=0;
        for (int i=0; i<this->length; i++) {
            a += pow(this->myvector[i], 2.0);
        }
        return sqrt(a);
    }

    /**
     * @brief Computes the dot product with another vector
     * @param other The other vector
     * @return Dot product result
     * @throws std::invalid_argument if vectors have different lengths
     */
    double dot (const My_Vec& other) const {
        if(this->length != other.length) {
            throw std::invalid_argument("Invalid lengths");
        }

        double a=0;
        for(int i=0; i<this->length; i++) {
            a += this->myvector[i] * other.myvector[i];
        }
        return a;
    }

    /**
     * @brief Scalar multiplication of the vector
     * @param k Scalar value
     * @return New vector resulting from scalar multiplication
     */
    My_Vec Scalar_Mul(double k) const {
        My_Vec New(this->length);
        for (int i=0; i<this->length; i++) {
            New.myvector[i] = k * this->myvector[i];
        }
        return New;
    }

    /**
     * @brief Copy constructor for My_Vec
     * @param other The vector to copy
     */
    My_Vec(const My_Vec& other)
    : length(other.length), myvector(other.myvector) {};

    /**
     * @brief Assignment operator for My_Vec
     * @param other The vector to assign
     * @return Reference to this vector
     */
    My_Vec& operator=(const My_Vec& other) {
        if (this != &other) {
            length = other.length;
            myvector = other.myvector;
        }
        return *this;
    }

    /**
     * @brief Creates a unit vector with a 1 at position i
     * @param i Position of the 1 in the unit vector
     * @param L Total length of the vector
     * @return Unit vector with 1 at position i
     */
    static My_Vec unit_vec(int i, int L) {
        My_Vec unit_vec(L);
        for (int j=0; j<L; j++) {
            unit_vec.myvector[j] = (j == i) ? 1 : 0;
        }
        return unit_vec;
    }

    /**
     * @brief Creates a zero vector of length i
     * @param i Length of the zero vector
     * @return Zero vector of length i
     */
    static My_Vec Zeros(const int& i) {
        My_Vec unit_vec(i);
        for (int j=0; j<i; j++) {
            unit_vec.myvector[j] = 0;
        }
        return unit_vec;
    }
};

/**
 * @class Matrix
 * @brief Matrix class implementation for mathematical operations
 * @details Provides basic matrix operations including addition, subtraction,
 *          multiplication, and various matrix factorizations (LU, QR)
 */
class Matrix {
public:
    int rows;                              ///< Number of rows in the matrix
    int cols;                              ///< Number of columns in the matrix
    std::vector<std::vector<double>> MyMAT;  ///< Internal matrix storage

    /**
     * @brief Constructor for Matrix
     * @param rs Number of rows (defaults to 1)
     * @param cs Number of columns (defaults to 1)
     */
    Matrix(int rs=1, int cs=1) {
        rows = rs;
        cols = cs;
        MyMAT.resize(rs, std::vector<double>(cs));
    }

    /**
     * @brief Copy constructor for Matrix
     * @param other The matrix to copy
     */
    Matrix(const Matrix& other)
    : rows(other.rows), cols(other.cols), MyMAT(other.MyMAT) {}

    /**
     * @brief Assignment operator for Matrix
     * @param other The matrix to assign
     * @return Reference to this matrix
     */
    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            rows = other.rows;
            cols = other.cols;
            MyMAT = other.MyMAT;
        }
        return *this;
    }

    /**
     * @brief Matrix addition operator
     * @param other The matrix to add
     * @return Resulting matrix after addition
     * @throws std::invalid_argument if matrices have different dimensions
     */
    Matrix operator+(const Matrix& other) const {
        Matrix New_Mat(other.rows, other.cols);
        if (this->cols != other.cols || other.rows != this->rows) {
            throw std::invalid_argument("Dimension mismatch");
        } else {
            for (int i=0; i<rows; i++) {
                for (int j=0; j<cols; j++) {
                    New_Mat.MyMAT[i][j] = other.MyMAT[i][j] + this->MyMAT[i][j];
                }
            }
        }
        return New_Mat;
    }

    /**
     * @brief Matrix subtraction operator
     * @param other The matrix to subtract
     * @return Resulting matrix after subtraction
     * @throws std::invalid_argument if matrices have different dimensions
     */
    Matrix operator- (const Matrix& other) const {
        Matrix New_Mat(other.rows, other.cols);
        if (this->cols != other.cols || other.rows != this->rows) {
            throw std::invalid_argument("Dimension mismatch");
        } else {
            for (int i=0; i<rows; i++) {
                for (int j=0; j<cols; j++) {
                    New_Mat.MyMAT[i][j] = this->MyMAT[i][j] - other.MyMAT[i][j];
                }
            }
        }
        return New_Mat;
    }

    /**
     * @brief Computes the LU Decomposition of the matrix
     * @return LUResult structure containing L, U matrices and P vector
     * @throws std::invalid_argument if the matrix is not square
     * @throws std::runtime_error if the matrix is singular or nearly singular
     */
    LUResult L_U() const {
        // Check for square matrix
        if(this->rows != this->cols) {
            throw std::invalid_argument("Not square");
        }
        
        Matrix L = eye(this->rows);
        Matrix U = Zeros(this->rows, this->cols);
        std::vector<int> P(this->rows);
        for(int i = 0; i < this->rows; i++) {
            P[i] = i;  
        }
        
        Matrix A_work(*this);
        
        for(int j = 0; j < this->cols; j++) {
            int pivot_row = j;
            double max_val = std::abs(A_work.MyMAT[j][j]);
            
            for(int i = j+1; i < this->rows; i++) {
                if(std::abs(A_work.MyMAT[i][j]) > max_val) {
                    pivot_row = i;
                    max_val = std::abs(A_work.MyMAT[i][j]);
                }
            }
            
            if(pivot_row != j) {
                for(int k = 0; k < this->cols; k++) {
                    std::swap(A_work.MyMAT[j][k], A_work.MyMAT[pivot_row][k]);
                }
                
                std::swap(P[j], P[pivot_row]);
                
                if(j > 0) {
                    for(int k = 0; k < j; k++) {
                        std::swap(L.MyMAT[j][k], L.MyMAT[pivot_row][k]);
                    }
                }
            }
            
            if(std::abs(A_work.MyMAT[j][j]) < 1e-10) {
                throw std::runtime_error("Matrix is singular or nearly singular");
            }
            
            for(int i = 0; i <= j; i++) {
                double sum = 0.0;
                for(int k = 0; k < i; k++) {
                    sum += L.MyMAT[i][k] * U.MyMAT[k][j];
                }
                U.MyMAT[i][j] = A_work.MyMAT[i][j] - sum;
            }
            
            for(int i = j+1; i < this->rows; i++) {
                double sum = 0.0;
                for(int k = 0; k < j; k++) {
                    sum += L.MyMAT[i][k] * U.MyMAT[k][j];
                }
                L.MyMAT[i][j] = (A_work.MyMAT[i][j] - sum) / U.MyMAT[j][j];
            }
        }
        
        LUResult result;
        result.L = L.MyMAT;
        result.U = U.MyMAT;
        result.P = P;    
        return result;
    }

    /**
     * @brief Creates a unit vector with a 1 at position i
     * @param i Position of the 1 in the unit vector
     * @param L Total length of the vector
     * @return Unit vector with 1 at position i
     */
    static My_Vec UV(int i, int L) {
        My_Vec unit_vec(L);
        for (int j=0; j<L; j++) {
            unit_vec.myvector[j] = (j == i) ? 1 : 0;
        }
        return unit_vec;
    }

    /**
     * @brief Embeds a Householder matrix into a larger identity matrix
     * @param Householder The Householder matrix to embed
     * @param A The original matrix
     * @return The embedded matrix
     */
    static Matrix Embed(const Matrix& Householder, const Matrix& A) {
        Matrix Hp = eye(A.rows);
        int i = A.rows - Householder.rows;
        for (int a = 0; a < Householder.rows; a++) {
            for (int b = 0; b < Householder.cols; b++) {
                Hp.MyMAT[i+a][i+b] = Householder.MyMAT[a][b];
            }
        }
        return Hp;
    }

    /**
     * @brief Computes the QR Decomposition of the matrix
     * @return QRresult structure containing Q, R matrices
     */
    QRresult QR_fact() const {
        Matrix Q = eye(this->rows);
        Matrix Aupdate = *this;
        
        for (int i = 0; i < std::min(this->rows-1, this->cols); i++) {
            My_Vec vecx(this->rows - i);
            for (int row = i; row < this->rows; row++) {
                vecx.myvector[row - i] = Aupdate.MyMAT[row][i];
            }
            
            double n = vecx.Norm();
            
            My_Vec unit = My_Vec::unit_vec(0, this->rows - i);  // Use consistent naming
            
            My_Vec reflec_vec;
            if (vecx.myvector[0] < 0) {
                reflec_vec = vecx + unit.Scalar_Mul(n);
            } else {
                reflec_vec = vecx - unit.Scalar_Mul(n);
            }
            
            double normalize_ref = reflec_vec.Norm();
            
            // Check for zero vector (avoid division by zero)
            if (normalize_ref < 1e-10) {
                continue;  // Skip this iteration if vector is too small
            }
            
            My_Vec V = reflec_vec.Scalar_Mul(1.0 / normalize_ref);
            
            // Create Householder matrix H = I - 2*vv^T
            Matrix I = Matrix::eye(V.length);
            Matrix vvT = Matrix::Outer_Product(V, V);
            Matrix Householder = I - vvT.Scalar_Mul(2.0);
            
            // Embed the Householder matrix into a larger identity matrix
            Matrix Hprime = Matrix::eye(this->rows);
            for (int r = 0; r < Householder.rows; r++) {
                for (int c = 0; c < Householder.cols; c++) {
                    Hprime.MyMAT[i + r][i + c] = Householder.MyMAT[r][c];
                }
            }
            
            Aupdate = Hprime * Aupdate;
            Q = Q * Hprime; 
        }
        
        QRresult QR;
        std::vector<std::vector<double>> Q_n = Q.Transpose().MyMAT;
        std::vector<std::vector<double>> R_n = Aupdate.MyMAT;

        QR.Q = Q_n; 
        QR.R = R_n;
        return QR;
    }  
    
    /**
     * @brief Matrix multiplication operator
     * @param other The matrix to multiply with
     * @return Resulting matrix after multiplication
     * @throws std::invalid_argument if dimensions are not compatible
     */
    Matrix operator*(const Matrix& other) const {
        if (this->cols != other.rows) {
            throw std::invalid_argument("Dimension mismatch");
        }
    
        Matrix Ans(this->rows, other.cols); // Ensure constructor handles allocation
    
        for (int i = 0; i < this->rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                double sum_of_multiples = 0.0;
                for (int k = 0; k < this->cols; ++k) {
                    sum_of_multiples += this->MyMAT[i][k] * other.MyMAT[k][j];
                }
                Ans.MyMAT[i][j] = sum_of_multiples;
            }
        }
    
        return Ans;
    }
    

    /**
     * @brief Multiplies the matrix with a vector
     * @param x The vector to multiply with
     * @return Resulting vector after multiplication
     * @throws std::invalid_argument if dimensions are not compatible
     */
    My_Vec multiply(const My_Vec& x) const {
        if(this->cols != x.length) {
            throw std::invalid_argument("Dimension mismatch");
        }

        My_Vec ans(this->rows);

        for (int i=0; i<this->rows; i++) {
            double dot_prod = 0;
            for(int j=0; j<this->cols; j++) {
                dot_prod += this->MyMAT[i][j] * x.myvector[j];
            }
            ans.myvector[i] = dot_prod;
        }

        return ans;
    }

    /**
     * @brief Static method to compute the outer product of two vectors
     * @param u First vector
     * @param v Second vector
     * @return Matrix resulting from the outer product
     */
    static Matrix Outer_Product(const My_Vec& u, const My_Vec& v) {
        Matrix Output(u.length, v.length);
        for(int i=0; i<u.length; i++) {
            for(int j=0; j<v.length; j++) {
                Output.MyMAT[i][j] = u.myvector[i] * v.myvector[j];
            }
        }
        return Output;
    }

    /**
     * @brief Transposes the matrix
     * @return New matrix that is the transpose of this matrix
     */
    Matrix Transpose() const {
        Matrix New_Mat(this->cols, this->rows);
        for (int i = 0; i<rows; i++) {
            for(int j=0; j<cols; j++) {
                New_Mat.MyMAT[j][i] = this->MyMAT[i][j];
            }
        }
        return New_Mat;
    }

    /**
     * @brief Scalar multiplication of the matrix
     * @param k Scalar value
     * @return New matrix resulting from scalar multiplication
     */
    Matrix Scalar_Mul(double k) const {
        Matrix new_mat(this->rows, this->cols);
        
        for (int i = 0; i<rows; i++) {
            for(int j=0; j<cols; j++) {
                new_mat.MyMAT[i][j] = this->MyMAT[i][j] * k;
            }
        }
        return new_mat;
    }

    /**
     * @brief Static method to create an identity matrix of size a
     * @param a Size of the matrix (number of rows and columns)
     * @return Identity matrix of size a
     */
    static Matrix eye(int a) {
        Matrix Identity(a, a);
        for(int i=0; i<a; i++) {
            Identity.MyMAT[i][i] = 1;
        }
        return Identity;
    }

    /**
     * @brief Static method to create a matrix of ones
     * @param a Number of rows
     * @param b Number of columns
     * @return Matrix of size a x b with all elements set to 1
     */
    static Matrix Ones(int a, int b) {
        Matrix Ones(a, b);
        for(int i=0; i<a; i++) {
            for(int j=0; j<b; j++) {
                Ones.MyMAT[i][j] = 1;
            }
        }
        return Ones;
    }

    /**
     * @brief Static method to create a matrix of zeros
     * @param a Number of rows
     * @param b Number of columns
     * @return Matrix of size a x b with all elements set to 0
     */
    static Matrix Zeros(int a, int b) { 
        Matrix Zeros(a, b);
        for(int i=0; i<a; i++) {
            for(int j=0; j<b; j++) {
                Zeros.MyMAT[i][j] = 0;
            }
        }
        return Zeros;
    }
};

/**
 * @class Linear_Solvers
 * @brief Class containing static methods for solving linear systems
 * @details using various methods including LU and QR decomposition
 */
class Linear_Solvers {
    static My_Vec SolveLU(const Matrix& A, const My_Vec& b);
    static My_Vec SolveQR(const Matrix& A, const My_Vec& b);
    static My_Vec Inverse(const Matrix& A);

    static My_Vec ForwardSubstitution(const Matrix& L, const My_Vec& b);
    static My_Vec BackwardSubstitution(const Matrix& U, const My_Vec& y);
};

struct LUResult_to_pass {
   Matrix L;
   Matrix U;
   std::vector<int> P;
};


struct QR_result_to_pass {
     Matrix Q;
     Matrix R;
};


LUResult_to_pass conv_LU(const LUResult& LU) {
    int a = LU.L.size();
    int b = LU.L[0].size();

    Matrix L(a, b);
    L.MyMAT = LU.L;


    int c = LU.U.size();
    int d = LU.U[0].size();
    
    Matrix U(c, d);
    U.MyMAT = LU.U;
    
    LUResult_to_pass LU_new;

    LU_new.L = L;
    LU_new.U = U;
    LU_new.P = LU.P;
    return LU_new;
}


QR_result_to_pass conv_QR(const QRresult& QR) {
    int a = QR.Q.size();
    int b = QR.Q[0].size();

    Matrix Q(a, b);
    Q.MyMAT = QR.Q;


    int c = QR.R.size();
    int d = QR.R[0].size();
    
    Matrix R(c, d);
    R.MyMAT = QR.R;
    
    QR_result_to_pass QR_new;

    QR_new.Q = Q;
    QR_new.R = R;
    return QR_new;
}


#endif
