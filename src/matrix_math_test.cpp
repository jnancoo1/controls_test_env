#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "matrix_math.hpp"

// Test utilities
const double EPSILON = 1e-10;

bool isEqual(double a, double b, double epsilon = EPSILON) {
    return std::abs(a - b) < epsilon;
}

bool vectorsEqual(const My_Vec& a, const My_Vec& b, double epsilon = EPSILON) {
    if (a.length != b.length) return false;
    for (int i = 0; i < a.length; i++) {
        if (!isEqual(a.myvector[i], b.myvector[i], epsilon)) {
            return false;
        }
    }
    return true;
}

bool matricesEqual(const Matrix& a, const Matrix& b, double epsilon = EPSILON) {
    if (a.rows != b.rows || a.cols != b.cols) return false;
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            if (!isEqual(a.MyMAT[i][j], b.MyMAT[i][j], epsilon)) {
                return false;
            }
        }
    }
    return true;
}

// Custom matchers for better error messages
::testing::AssertionResult VectorsAreEqual(const My_Vec& expected, const My_Vec& actual) {
    if (vectorsEqual(expected, actual)) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure() 
        << "Vectors are not equal. Expected length: " << expected.length 
        << ", Actual length: " << actual.length;
}

::testing::AssertionResult MatricesAreEqual(const Matrix& expected, const Matrix& actual) {
    if (matricesEqual(expected, actual)) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure() 
        << "Matrices are not equal. Expected: " << expected.rows << "x" << expected.cols
        << ", Actual: " << actual.rows << "x" << actual.cols;
}

// ==================== VECTOR TESTS ====================

class VectorTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for vector tests
    }
    
    void TearDown() override {
        // Cleanup
    }
};

TEST_F(VectorTest, DefaultConstructor) {
    My_Vec v;
    EXPECT_EQ(v.length, 1);
    EXPECT_EQ(v.myvector.size(), 1);
}

TEST_F(VectorTest, ParameterizedConstructor) {
    My_Vec v(5);
    EXPECT_EQ(v.length, 5);
    EXPECT_EQ(v.myvector.size(), 5);
}

TEST_F(VectorTest, CopyConstructor) {
    My_Vec v1(3);
    v1.myvector = {1.0, 2.0, 3.0};
    
    My_Vec v2(v1);
    EXPECT_EQ(v2.length, 3);
    EXPECT_TRUE(VectorsAreEqual(v1, v2));
}

TEST_F(VectorTest, AssignmentOperator) {
    My_Vec v1(3);
    v1.myvector = {1.0, 2.0, 3.0};
    
    My_Vec v2(2);
    v2 = v1;
    
    EXPECT_EQ(v2.length, 3);
    EXPECT_TRUE(VectorsAreEqual(v1, v2));
}

TEST_F(VectorTest, Addition) {
    My_Vec v1(3);
    My_Vec v2(3);
    v1.myvector = {1.0, 2.0, 3.0};
    v2.myvector = {4.0, 5.0, 6.0};
    
    My_Vec result = v1 + v2;
    My_Vec expected(3);
    expected.myvector = {5.0, 7.0, 9.0};
    
    EXPECT_TRUE(VectorsAreEqual(expected, result));
}

TEST_F(VectorTest, AdditionDimensionMismatch) {
    My_Vec v1(3);
    My_Vec v2(2);
    
    EXPECT_THROW(v1 + v2, std::invalid_argument);
}

TEST_F(VectorTest, Subtraction) {
    My_Vec v1(3);
    My_Vec v2(3);
    v1.myvector = {5.0, 7.0, 9.0};
    v2.myvector = {1.0, 2.0, 3.0};
    
    My_Vec result = v1 - v2;
    My_Vec expected(3);
    expected.myvector = {4.0, 5.0, 6.0};
    
    EXPECT_TRUE(VectorsAreEqual(expected, result));
}

TEST_F(VectorTest, SubtractionDimensionMismatch) {
    My_Vec v1(3);
    My_Vec v2(2);
    
    EXPECT_THROW(v1 - v2, std::invalid_argument);
}

TEST_F(VectorTest, ScalarMultiplication) {
    My_Vec v(3);
    v.myvector = {1.0, 2.0, 3.0};
    
    My_Vec result = v.Scalar_Mul(2.5);
    My_Vec expected(3);
    expected.myvector = {2.5, 5.0, 7.5};
    
    EXPECT_TRUE(VectorsAreEqual(expected, result));
}

TEST_F(VectorTest, DotProduct) {
    My_Vec v1(3);
    My_Vec v2(3);
    v1.myvector = {1.0, 2.0, 3.0};
    v2.myvector = {4.0, 5.0, 6.0};
    
    double result = v1.dot(v2);
    double expected = 1.0*4.0 + 2.0*5.0 + 3.0*6.0; // 32.0
    
    EXPECT_NEAR(result, expected, EPSILON);
}

TEST_F(VectorTest, DotProductDimensionMismatch) {
    My_Vec v1(3);
    My_Vec v2(2);
    
    EXPECT_THROW(v1.dot(v2), std::invalid_argument);
}

TEST_F(VectorTest, Norm) {
    My_Vec v(3);
    v.myvector = {3.0, 4.0, 0.0};
    
    double result = v.Norm();
    double expected = 5.0; // sqrt(3^2 + 4^2 + 0^2) = 5
    
    EXPECT_NEAR(result, expected, EPSILON);
}

TEST_F(VectorTest, OnesVector) {
    My_Vec ones = My_Vec::ones(4);
    My_Vec expected(4);
    expected.myvector = {1.0, 1.0, 1.0, 1.0};
    
    EXPECT_TRUE(VectorsAreEqual(expected, ones));
}

TEST_F(VectorTest, ZerosVector) {
    My_Vec zeros = My_Vec::Zeros(3);
    My_Vec expected(3);
    expected.myvector = {0.0, 0.0, 0.0};
    
    EXPECT_TRUE(VectorsAreEqual(expected, zeros));
}

TEST_F(VectorTest, UnitVector) {
    My_Vec unit = My_Vec::unit_vec(1, 4);
    My_Vec expected(4);
    expected.myvector = {0.0, 1.0, 0.0, 0.0};
    
    EXPECT_TRUE(VectorsAreEqual(expected, unit));
}

// ==================== MATRIX TESTS ====================

class MatrixTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for matrix tests
    }
};

TEST_F(MatrixTest, DefaultConstructor) {
    Matrix m;
    EXPECT_EQ(m.rows, 1);
    EXPECT_EQ(m.cols, 1);
    EXPECT_EQ(m.MyMAT.size(), 1);
    EXPECT_EQ(m.MyMAT[0].size(), 1);
}

TEST_F(MatrixTest, ParameterizedConstructor) {
    Matrix m(3, 4);
    EXPECT_EQ(m.rows, 3);
    EXPECT_EQ(m.cols, 4);
    EXPECT_EQ(m.MyMAT.size(), 3);
    EXPECT_EQ(m.MyMAT[0].size(), 4);
}

TEST_F(MatrixTest, CopyConstructor) {
    Matrix m1(2, 2);
    m1.MyMAT = {{1.0, 2.0}, {3.0, 4.0}};
    
    Matrix m2(m1);
    EXPECT_EQ(m2.rows, 2);
    EXPECT_EQ(m2.cols, 2);
    EXPECT_TRUE(MatricesAreEqual(m1, m2));
}

TEST_F(MatrixTest, AssignmentOperator) {
    Matrix m1(2, 2);
    m1.MyMAT = {{1.0, 2.0}, {3.0, 4.0}};
    
    Matrix m2(3, 3);
    m2 = m1;
    
    EXPECT_EQ(m2.rows, 2);
    EXPECT_EQ(m2.cols, 2);
    EXPECT_TRUE(MatricesAreEqual(m1, m2));
}

TEST_F(MatrixTest, Addition) {
    Matrix m1(2, 2);
    Matrix m2(2, 2);
    m1.MyMAT = {{1.0, 2.0}, {3.0, 4.0}};
    m2.MyMAT = {{5.0, 6.0}, {7.0, 8.0}};
    
    Matrix result = m1 + m2;
    Matrix expected(2, 2);
    expected.MyMAT = {{6.0, 8.0}, {10.0, 12.0}};
    
    EXPECT_TRUE(MatricesAreEqual(expected, result));
}

TEST_F(MatrixTest, AdditionDimensionMismatch) {
    Matrix m1(2, 2);
    Matrix m2(3, 3);
    
    EXPECT_THROW(m1 + m2, std::invalid_argument);
}

TEST_F(MatrixTest, Subtraction) {
    Matrix m1(2, 2);
    Matrix m2(2, 2);
    m1.MyMAT = {{5.0, 6.0}, {7.0, 8.0}};
    m2.MyMAT = {{1.0, 2.0}, {3.0, 4.0}};
    
    Matrix result = m1 - m2;
    Matrix expected(2, 2);
    expected.MyMAT = {{4.0, 4.0}, {4.0, 4.0}};
    
    EXPECT_TRUE(MatricesAreEqual(expected, result));
}

TEST_F(MatrixTest, Multiplication) {
    Matrix m1(2, 3);
    Matrix m2(3, 2);
    m1.MyMAT = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
    m2.MyMAT = {{7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0}};
    
    Matrix result = m1 * m2;
    Matrix expected(2, 2);
    expected.MyMAT = {{58.0, 64.0}, {139.0, 154.0}};
    
    EXPECT_TRUE(MatricesAreEqual(expected, result));
}

TEST_F(MatrixTest, MultiplicationDimensionMismatch) {
    Matrix m1(2, 3);
    Matrix m2(2, 2);
    
    EXPECT_THROW(m1 * m2, std::invalid_argument);
}

TEST_F(MatrixTest, ScalarMultiplication) {
    Matrix m(2, 2);
    m.MyMAT = {{1.0, 2.0}, {3.0, 4.0}};
    
    Matrix result = m.Scalar_Mul(2.5);
    Matrix expected(2, 2);
    expected.MyMAT = {{2.5, 5.0}, {7.5, 10.0}};
    
    EXPECT_TRUE(MatricesAreEqual(expected, result));
}

TEST_F(MatrixTest, Transpose) {
    Matrix m(2, 3);
    m.MyMAT = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
    
    Matrix result = m.Transpose();
    Matrix expected(3, 2);
    expected.MyMAT = {{1.0, 4.0}, {2.0, 5.0}, {3.0, 6.0}};
    
    EXPECT_TRUE(MatricesAreEqual(expected, result));
}

TEST_F(MatrixTest, IdentityMatrix) {
    Matrix identity = Matrix::eye(3);
    Matrix expected(3, 3);
    expected.MyMAT = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    
    EXPECT_TRUE(MatricesAreEqual(expected, identity));
}

TEST_F(MatrixTest, OnesMatrix) {
    Matrix ones = Matrix::Ones(2, 3);
    Matrix expected(2, 3);
    expected.MyMAT = {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};
    
    EXPECT_TRUE(MatricesAreEqual(expected, ones));
}

TEST_F(MatrixTest, ZerosMatrix) {
    Matrix zeros = Matrix::Zeros(2, 3);
    Matrix expected(2, 3);
    expected.MyMAT = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    
    EXPECT_TRUE(MatricesAreEqual(expected, zeros));
}

TEST_F(MatrixTest, OuterProduct) {
    My_Vec u(2);
    My_Vec v(3);
    u.myvector = {1.0, 2.0};
    v.myvector = {3.0, 4.0, 5.0};
    
    Matrix result = Matrix::Outer_Product(u, v);
    Matrix expected(2, 3);
    expected.MyMAT = {{3.0, 4.0, 5.0}, {6.0, 8.0, 10.0}};
    
    EXPECT_TRUE(MatricesAreEqual(expected, result));
}

TEST_F(MatrixTest, MatrixVectorMultiplication) {
    Matrix m(2, 3);
    My_Vec v(3);
    m.MyMAT = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
    v.myvector = {1.0, 2.0, 3.0};
    
    My_Vec result = m.multiply(v);
    My_Vec expected(2);
    expected.myvector = {14.0, 32.0}; // [1*1+2*2+3*3, 4*1+5*2+6*3]
    
    EXPECT_TRUE(VectorsAreEqual(expected, result));
}

TEST_F(MatrixTest, MatrixVectorMultiplicationDimensionMismatch) {
    Matrix m(2, 3);
    My_Vec v(2);
    
    EXPECT_THROW(m.multiply(v), std::invalid_argument);
}

//LU DECOMPOSITION TESTS
class LUDecompositionTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Setup common matrices for LU tests
        simple_matrix = Matrix(3, 3);
        simple_matrix.MyMAT = {{2.0, 1.0, 1.0}, 
                               {4.0, 3.0, 3.0}, 
                               {8.0, 7.0, 9.0}};
    }
    
    Matrix simple_matrix;
};

TEST_F(LUDecompositionTest, BasicLUDecomposition) {
    LUResult lu_result = simple_matrix.L_U();
    
    // Convert to Matrix objects for easier testing
    LUResult_to_pass lu_converted = conv_LU(lu_result);
    
    // Verify L is lower triangular
    for (int i = 0; i < lu_converted.L.rows; i++) {
        for (int j = i + 1; j < lu_converted.L.cols; j++) {
            EXPECT_NEAR(lu_converted.L.MyMAT[i][j], 0.0, EPSILON);
        }
    }
    
    // Verify U is upper triangular
    for (int i = 1; i < lu_converted.U.rows; i++) {
        for (int j = 0; j < i; j++) {
            EXPECT_NEAR(lu_converted.U.MyMAT[i][j], 0.0, EPSILON);
        }
    }
    
    // Verify L has 1s on diagonal
    for (int i = 0; i < lu_converted.L.rows; i++) {
        EXPECT_NEAR(lu_converted.L.MyMAT[i][i], 1.0, EPSILON);
    }
}

TEST_F(LUDecompositionTest, LUReconstructOriginal) {
    LUResult lu_result = simple_matrix.L_U();
    LUResult_to_pass lu_converted = conv_LU(lu_result);
    
    // Reconstruct the original matrix (with permutation)
    Matrix reconstructed = lu_converted.L * lu_converted.U;
    
    // Apply permutation to original matrix for comparison
    Matrix permuted_original(simple_matrix.rows, simple_matrix.cols);
    for (int i = 0; i < simple_matrix.rows; i++) {
        for (int j = 0; j < simple_matrix.cols; j++) {
            permuted_original.MyMAT[i][j] = simple_matrix.MyMAT[lu_converted.P[i]][j];
        }
    }
    
    EXPECT_TRUE(MatricesAreEqual(permuted_original, reconstructed));
}

TEST_F(LUDecompositionTest, NonSquareMatrixThrows) {
    Matrix non_square(2, 3);
    EXPECT_THROW(non_square.L_U(), std::invalid_argument);
}

TEST_F(LUDecompositionTest, SingularMatrixThrows) {
    Matrix singular(3, 3);
    singular.MyMAT = {{1.0, 2.0, 3.0}, 
                      {2.0, 4.0, 6.0}, 
                      {1.0, 1.0, 1.0}};
    
    EXPECT_THROW(singular.L_U(), std::runtime_error);
}

// QR DECOMPOSITION TESTS 

class QRDecompositionTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_matrix = Matrix(3, 3);
        test_matrix.MyMAT = {{1.0, 1.0, 0.0}, 
                             {1.0, 0.0, 1.0}, 
                             {0.0, 1.0, 1.0}};
    }
    
    Matrix test_matrix;
};

TEST_F(QRDecompositionTest, BasicQRDecomposition) {
    QRresult qr_result = test_matrix.QR_fact();
    QR_result_to_pass qr_converted = conv_QR(qr_result);
    
    // Verify R is upper triangular
    for (int i = 1; i < qr_converted.R.rows; i++) {
        for (int j = 0; j < i && j < qr_converted.R.cols; j++) {
            EXPECT_NEAR(qr_converted.R.MyMAT[i][j], 0.0, EPSILON);
        }
    }
}

TEST_F(QRDecompositionTest, QRReconstructOriginal) {
    QRresult qr_result = test_matrix.QR_fact();
    QR_result_to_pass qr_converted = conv_QR(qr_result);
    
    Matrix reconstructed = qr_converted.Q * qr_converted.R;
    EXPECT_TRUE(MatricesAreEqual(test_matrix, reconstructed));
}

TEST_F(QRDecompositionTest, QIsOrthogonal) {
    QRresult qr_result = test_matrix.QR_fact();
    QR_result_to_pass qr_converted = conv_QR(qr_result);
    
    Matrix Q_transpose = qr_converted.Q.Transpose();
    Matrix should_be_identity = Q_transpose * qr_converted.Q;
    Matrix identity = Matrix::eye(should_be_identity.rows);
    
    EXPECT_TRUE(MatricesAreEqual(identity, should_be_identity));
}

//CONVERSION FUNCTION TESTS

class ConversionTest : public ::testing::Test {};

TEST_F(ConversionTest, LUConversion) {
    Matrix test_matrix(2, 2);
    test_matrix.MyMAT = {{4.0, 3.0}, {6.0, 3.0}};
    
    LUResult lu_result = test_matrix.L_U();
    LUResult_to_pass converted = conv_LU(lu_result);
    
    EXPECT_EQ(converted.L.rows, 2);
    EXPECT_EQ(converted.L.cols, 2);
    EXPECT_EQ(converted.U.rows, 2);
    EXPECT_EQ(converted.U.cols, 2);
    EXPECT_EQ(converted.P.size(), 2);
}

TEST_F(ConversionTest, QRConversion) {
    Matrix test_matrix(2, 2);
    test_matrix.MyMAT = {{1.0, 1.0}, {0.0, 1.0}};
    
    QRresult qr_result = test_matrix.QR_fact();
    QR_result_to_pass converted = conv_QR(qr_result);
    
    EXPECT_EQ(converted.Q.rows, 2);
    EXPECT_EQ(converted.Q.cols, 2);
    EXPECT_EQ(converted.R.rows, 2);
    EXPECT_EQ(converted.R.cols, 2);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
