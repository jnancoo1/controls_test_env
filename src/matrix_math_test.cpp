#include <iostream>
#include <vector>
#include <cmath>
#include <gtest/gtest.h>
#include "matrix_math.hpp"

TEST(MatrixVectorTest, MultiplyByIdentity) {
    Matrix A(3, 3);
    Matrix I(3, 3);

    // Fill A with some values
    A.MyMAT = {
        {2, -1, 0},
        {4, 3, 1},
        {0, 5, 6}
    };

    // Create identity matrix I
    I.MyMAT = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    Matrix C = A * I;

    // Multiplying by identity matrix should return A unchanged
    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < A.cols; ++j) {
            EXPECT_NEAR(C.MyMAT[i][j], A.MyMAT[i][j], 1e-10);
        }
    }
}

TEST(MatrixVectorTest, MultiplyRectangularMatrices) {
    // 2x3 matrix
    Matrix A(2, 3);
    A.MyMAT = {
        {1, 2, 3},
        {4, 5, 6}
    };

    // 3x2 matrix
    Matrix B(3, 2);
    B.MyMAT = {
        {7, 8},
        {9, 10},
        {11, 12}
    };

    // Result should be 2x2 matrix
    Matrix C = A * B;

    EXPECT_EQ(C.rows, 2);
    EXPECT_EQ(C.cols, 2);

    EXPECT_NEAR(C.MyMAT[0][0], 1*7 + 2*9 + 3*11, 1e-10);  // 58
    EXPECT_NEAR(C.MyMAT[0][1], 1*8 + 2*10 + 3*12, 1e-10); // 64
    EXPECT_NEAR(C.MyMAT[1][0], 4*7 + 5*9 + 6*11, 1e-10);  // 139
    EXPECT_NEAR(C.MyMAT[1][1], 4*8 + 5*10 + 6*12, 1e-10); // 154
}



TEST(MatrixVectorTest, MatrixMultiplication) {
    Matrix A(3, 3);
    Matrix B(3, 3);

    A.MyMAT = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    B.MyMAT = {
        {4, 0, 0},
        {0, 13.5, 0},
        {0, 0, 24}
    };

    Matrix C = A * B;

    EXPECT_NEAR(C.MyMAT[0][0], 4.0, 1e-10);
    EXPECT_NEAR(C.MyMAT[1][1], 13.5, 1e-10);
    EXPECT_NEAR(C.MyMAT[2][2], 24.0, 1e-10);
}

TEST(MatrixVectorTest, AssociativityOfMultiplication) {
    // Define matrices with compatible sizes:
    Matrix A(2, 3);
    A.MyMAT = {
        {1, 2, 3},
        {4, 5, 6}
    };

    Matrix B(3, 2);
    B.MyMAT = {
        {7, 8},
        {9, 10},
        {11, 12}
    };

    Matrix C(2, 2);
    C.MyMAT = {
        {1, 0},
        {0, 1}
    };

    // Calculate (A * B) * C
    Matrix AB = A * B;      // 2x2 matrix
    Matrix left = AB * C;   // 2x2 * 2x2 = 2x2

    // Calculate A * (B * C)
    Matrix BC = B * C;      // 3x2 * 2x2 = 3x2
    Matrix right = A * BC;  // 2x3 * 3x2 = 2x2

    // Compare element-wise
    for (int i = 0; i < left.rows; ++i) {
        for (int j = 0; j < left.cols; ++j) {
            EXPECT_NEAR(left.MyMAT[i][j], right.MyMAT[i][j], 1e-10);
        }
    }
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
