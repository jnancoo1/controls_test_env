#include <iostream>
#include <vector>
#include <cmath>

#ifndef MATRIX_MATH_HPP
#define MATRIX_MATH_HPP
struct LUResult {
    Matrix L;
    Matrix U;
};


struct QRresult {
    Matrix L;
    Matrix U;
};



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
    static My_Vec ones(int a){}
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


    static My_Vec unit_vec(int i,int L){

        My_Vec unit_vec(L);

        for (int j=0;j<L;j++){

            if(j==i){
                unit_vec.myvector[j]=1;
            }
            else{
                unit_vec.myvector[j]=0;
            }
        }

        return unit_vec;
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

    static Matrix Zeros(int a,int b){ 
        
        Matrix Zeros(a,b);
        for(int i=0;i<a;i++){

            for(int j=0;j<b;j++){

                Zeros.MyMAT[i][j]=0;
            }

        }
        return Zeros;
    }


    LUResult L_U() const{

        if(this->rows!=this->cols){

            throw std::invalid_argument("Not square");
        }

        Matrix L=eye(this->rows);
        Matrix U=Zeros(this->rows,this->cols);

    
        for(int i=0;i<this->rows;i++){
            for ( int j=0;j<this->cols;j++){
                if(i==0){

                    U.MyMAT[i][j]=this->MyMAT[i][j];

                }


                else{
                    double sumterm=0;
                    for(int k=0;k<i;k++){

                        sumterm+=(L.MyMAT[i][k]*U.MyMAT[k][j]);
                    }
                    U.MyMAT[i][j]=this->MyMAT[i][j]-sumterm;
                }
            }

        }


        for(int i=0;i<this->rows;i++){
            for(int j=0;j<this->cols;j++){
                
                if(j<=i){
                    double sumterm2=0;
                    for(int k=0;k<j;k++){
                        sumterm2+=(L.MyMAT[i][k]*U.MyMAT[k][j]);
                    }
                    L.MyMAT[i][j]=(this->MyMAT[i][j]-sumterm2)/U.MyMAT[j][j];
                }

                else{

                    L.MyMAT[i][j]=0;  
                }
            }
           
        }


        LUResult ans;
        ans.L=L;
        ans.U=U;
        return ans;
    }

    static My_Vec UV(int i,int L){

        My_Vec unit_vec(L);

        for (int j=0;j<L;j++){

            if(j==i){
                unit_vec.myvector[j]=1;
            }
            else{
                unit_vec.myvector[j]=0;
            }
        }

        return unit_vec;
    }

    QRresult QR_fact() const {

        Matrix Q(this->rows,this->rows);
        Matrix R(this->rows,this->cols);

        for (int i=0;i<this->cols;i++){
        My_Vec vecx;
        for (int row=0;row<this->rows;row++){
            vecx.myvector.push_back(this->MyMAT[row][i]);
        }
        double n=vecx.Norm();
        My_Vec unit = UV(0,this->cols);
        My_Vec reflec_vec(this->cols);
        if(vecx.myvector[0]<0)
            reflec_vec=vecx+unit.Scalar_Mul(n);
        else
            reflec_vec=vecx-unit.Scalar_Mul(n);
        
        double normalize_ref;
        normalize_ref=reflec_vec.Norm();

    }}
        
    Matrix operator*(const Matrix& other) const{
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