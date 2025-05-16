#include <iostream>
#include <vector>
#include <cmath>

#ifndef MATRIX_MATH_HPP
#define MATRIX_MATH_HPP




 class My_Vec{
public:
    int length=1;
    std::vector<double> myvector;
    My_Vec(int l=1){
        length=l;
        myvector.resize(length);
    }
    My_Vec operator+ (const My_Vec& other) const{

        if(this->length!=other.length){
            throw std::invalid_argument("Invalid lengths");
        }
        My_Vec Vec(this->length);

        for (int i=0;i<this->length;i++){

            Vec.myvector[i]=this->myvector[i]+other.myvector[i];

        }
        return Vec;  
    }
    My_Vec operator- (const My_Vec& other) const{

        if(this->length!=other.length){
            throw std::invalid_argument("Invalid lengths");
        }
        My_Vec Vec(this->length);

        for (int i=0;i<this->length;i++){

            Vec.myvector[i]=this->myvector[i]-other.myvector[i];
        }
        return Vec;
        
}

 double Norm() const{
     double a=0;
    for (int i=0;i<this->length;i++){
        
        a+=pow(this->myvector[i],2.0);

    }
    double ans=sqrt(a);
    return ans;
}

double dot (const My_Vec& other) const{

        if(this->length!=other.length)
            throw std::invalid_argument("Invalid lengths");
        

        double a=0;
        for(int i=0;i<this->length;i++){
          a+=this->myvector[i]*other.myvector[i];

        }
        return a;
    }

    My_Vec Scalar_Mul(double k) const{

        My_Vec New(this->length);
        for (int i=0;i<this->length;i++){

            New.myvector[i]=k*this->myvector[i];
            
        }
        return New;
    }

    My_Vec(const My_Vec& other)
    : length(other.length), myvector(other.myvector){};

        My_Vec& operator=(const My_Vec& other) {
            if (this != &other) {
                length = other.length;
                myvector = other.myvector;
            }
            return *this;
        }

    
};

class Matrix{

    public:

        int rows=1;
        int cols=1;
    
        std::vector<std::vector<double>> MyMAT;
    
        Matrix(int rs=1,int cs=1){
    
            rows=rs;
            cols=cs;
            MyMAT.resize(rs,std::vector<double>(cs));
    
        }

        Matrix(const Matrix& other)
        : rows(other.rows), cols(other.cols), MyMAT(other.MyMAT) {}
        Matrix& operator=(const Matrix& other) {
            if (this != &other) {
                rows = other.rows;
                cols = other.cols;
                MyMAT = other.MyMAT;
            }
            return *this;
        }
        
        Matrix operator+(const Matrix& other)const{
    
            Matrix New_Mat(other.rows,other.cols);
            if (this->cols!=other.cols || other.rows!=this->rows){
    
                throw std::invalid_argument("Dimension mismatch");
            }
            else{
    
                for (int i=0;i<rows;i++){
    
                    for (int j=0;j<cols;j++){

                        New_Mat.MyMAT[i][j]=other.MyMAT[i][j] + this->MyMAT[i][j];
                    }
                }

        }
        return New_Mat;
    }

        Matrix operator- (const Matrix& other) const{
    
            Matrix New_Mat(other.rows,other.cols);
            if (this->cols!=other.cols || other.rows!=this->rows){
    
                throw std::invalid_argument("Dimension mismatch");
            }
            else{
    
                for (int i=0;i<rows;i++){
    
                    for (int j=0;j<cols;j++){

                        New_Mat.MyMAT[i][j]= this->MyMAT[i][j]-other.MyMAT[i][j];
                    }
                }

        }
            return New_Mat;


    }

    Matrix Transpose() const{
        Matrix New_Mat(this->cols,this->rows);
        for (int i =0;i<rows;i++){

            for(int j=0;j<cols;j++)

            New_Mat.MyMAT[j][i]=this->MyMAT[i][j];
        }

        return New_Mat;

    }


    Matrix Scalar_Mul(double k) const{
        Matrix new_mat(this->rows,this->cols);
        
        for (int i =0;i<rows;i++){

            for(int j=0;j<cols;j++)

            new_mat.MyMAT[i][j]=this->MyMAT[i][j]*k;
        }
        return new_mat;

    }

    static Matrix eye(int a) {
        Matrix Identity(a,a);
        for(int i=0;i<a;i++){

            Identity.MyMAT[i][i]=1;
        }
        return Identity;

    }

    static Matrix Ones(int a,int b) {
        Matrix Ones(a,b);

        for(int i=0;i<a;i++){
            for(int j=0;j<b;j++){

            Ones.MyMAT[i][j]=1;
        
        }
    }
        return Ones;

    }



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
    

    My_Vec multiply(const My_Vec& x) const{

        if(this->cols != x.length){
            throw std::invalid_argument("Dimension mismatch");
        }

        My_Vec ans(this->rows);

        for (int i=0;i<this->rows;i++){

            double dot_prod=0;
            for(int j=0;j<this->cols;j++){

                dot_prod+=this->MyMAT[i][j]*x.myvector[j];

            }

            ans.myvector[i]=dot_prod;
        }

        return ans;
    }
};


#endif