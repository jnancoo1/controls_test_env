#ifndef TRANSFER_FUNCTIONS_HPP
#define TRANSFER_FUNCTIONS_HPP

#include <vector>
#include <complex>
#include <optional>
#include<Eigen/Dense>
#include "discrete_state_space.hpp"

namespace TransferFunctions {

    struct StateSpace_System {

    Eigen::MatrixXcd A;
    Eigen::MatrixXcd B;
    Eigen::MatrixXcd C;
    double D;

};

    // First Order System Parameters
    struct FirstOrderParams {
        double timeConstant;   
        double steadyStateGain;  
    };

    struct SecondOrderParams {
        double naturalFrequency; 
        double dampingRatio;     
        double steadyStateGain;  
        };

    class TransferFunction {
    private:
        std::vector<double> numerator;   
        std::vector<double> denominator; 

    public:
        // Constructor
        TransferFunction(const std::vector<double>& num, const std::vector<double>& den){
            numerator=num;
            denominator=den;

        }

        // Factory methods for common transfer functions
        static TransferFunction createFirstOrder(double K, double tau){

            std::vector<double> num = {K};
            std::vector<double> den = {tau, 1.0};
            return TransferFunction(num, den);
        };
        static TransferFunction createSecondOrder(double K, double wn, double zeta){

            std::vector<double> num = {K};
            std::vector<double> den = {1,2*zeta*wn,wn*wn};

            return TransferFunction(num, den);

        };

        const std::vector<std::complex<double>> getPoles() const{
            int n=this->denominator.size();
            std::vector<std::complex<double>> output;

            std::vector<double> normed_poles(this->denominator.size()); 

            double temp;
            if (this->denominator.size()==0||this->denominator.size()==1){
                output.push_back(0);
                return output;
            }
            else{
                //normalise
                for(int i=0;i<this->denominator.size();i++){

                    temp=this->denominator[0];

                    
                        normed_poles[i]=this->denominator[i]/temp;
                
                }
                
        Eigen::MatrixXcd Accf = Eigen::MatrixXcd::Zero(n, n);
        for(int i=0;i<n-1;i++){
            Accf(i,i+1)=1;
        }
        for(int j=0;j<n;j++){
            Accf(n-1,j)=-normed_poles[j];
        }
            Eigen::EigenSolver<Eigen::MatrixXcd> es(Accf);
            for (int i = 0; i < es.eigenvalues().size(); ++i) {
                output.push_back(es.eigenvalues()[i]);
            }
            return output;
            }
        };

        static std::vector<std::complex<double>> frequency_sweep_jw(const double& start,const double end,int n){

            std::vector<std::complex<double>> ans;
            double step=((end-start)/n);
            double temp=start;
            std::complex<double> temp2;std::complex<double> temp3(0,1);
            for(int i=0;i<n;i++){

                temp=+(i*step);//Hz
                temp2=temp*2*M_PI*temp3;
                ans.push_back(temp2);
            }
            return ans;
        }

        std::vector<double> Get_Mag_Response(TransferFunction tf,const double& start,const double end,int n){


            std::vector<std::complex<double>> freq_responses=frequency_sweep_jw(start,end,n);
            int num_size=tf.numerator.size();
            int denom_size=tf.denominator.size();
            std::complex <double> num_temp;

            std::complex <double> acc(0,0);
            std::complex <double> acc2(0,0);
            std::complex <double> acc3(0,0);
            std::vector<double> store;


        for(int iter=0;iter<n;iter++){

            for(int i=0;i<num_size;i++){

                num_temp=tf.numerator[i]*pow(freq_responses[iter],i);
                acc+=num_temp;
                }

            for(int j=0;j<denom_size;j++){

                num_temp=tf.denominator[j]*pow(freq_responses[iter],j);
                acc2+=num_temp;
                }
                acc3=acc/acc2;
                acc=0;
                acc2=0;
                store.push_back(std::abs(acc3));
            }
             return store;

            }

        std::vector<double> Get_freq_response(TransferFunction tf,const double& start,const double end,int n){


            std::vector<std::complex<double>> freq_responses=frequency_sweep_jw(start,end,n);
            int num_size=tf.numerator.size();
            int denom_size=tf.denominator.size();
            std::complex <double> num_temp;

            std::complex <double> acc(0,0);
            std::complex <double> acc2(0,0);
            std::complex <double> acc3(0,0);
            std::vector<double> store;


        for(int iter=0;iter<n;iter++){

            for(int i=0;i<num_size;i++){

                num_temp=tf.numerator[i]*pow(freq_responses[iter],i);
                acc+=num_temp;
                }

            for(int j=0;j<denom_size;j++){

                num_temp=tf.denominator[j]*pow(freq_responses[iter],j);
                acc2+=num_temp;
                }
                acc3=acc/acc2;
                acc=0;
                acc2=0;
                store.push_back(std::arg(acc3));
            }
             return store;

            }

        
        
        // Get the zeros of the transfer function
        std::vector<std::complex<double>> getZeros() const{

            int n=this->denominator.size();
            std::vector<std::complex<double>> output;

            std::vector<double> normed_zeros(this->numerator.size()); 

            double temp;
            if (this->denominator.size()==0||this->numerator.size()==1){
                output.push_back(0);
                return output;
            }
            else{
                //normalise
                for(int i=0;i<this->numerator.size();i++){

                    temp=this->numerator[0];

                    
                        normed_zeros[i]=this->numerator[i]/temp;
                
                }
                
        Eigen::MatrixXcd Accf = Eigen::MatrixXcd::Zero(n, n);
        for(int i=0;i<n-1;i++){
            Accf(i,i+1)=1;
        }
        for(int j=0;j<n;j++){
            Accf(n-1,j)=-numerator[j];
        }
            Eigen::EigenSolver<Eigen::MatrixXcd> es(Accf);
            for (int i = 0; i < es.eigenvalues().size(); ++i) {
                output.push_back(es.eigenvalues()[i]);
            }
            }
            return output;//check this over

        }
        

        StateSpace_System TF2SS(const TransferFunction tf){


            std::vector<double> normed_coeff(this->denominator.size()); 
            double temp = tf.denominator[0];
            for (int i = 0; i < tf.denominator.size(); i++) {
                normed_coeff[i] = tf.denominator[i] / temp;
            }

            
        int n=tf.denominator.size()-1;
        Eigen::MatrixXcd Accf = Eigen::MatrixXcd::Zero(n, n);
        for(int i=0;i<n-1;i++){
            Accf(i,i+1)=1;
        }
        for(int j=0;j<n;j++){
            Accf(n-1,j)=-normed_coeff[j];
        }

        Eigen::MatrixXd B(Accf.rows(), 1);
        B.setZero();
        B.row(Accf.rows()-1).setConstant(1.0);

        
        Eigen::MatrixXd C(1,Accf.rows());
        int n = tf.denominator.size();
        int m = tf.numerator.size();

        std::vector<double> num_padded(n, 0.0);
        for(int i = 0; i < m; i++) {
            num_padded[i] = tf.numerator[i] / temp; // normalize numerator by same leading denominator coeff
        }

        double bn = (m == n) ? num_padded[n - 1] : 0.0;

        for(int i = 0; i < n; i++) {
            double a_i_plus_1 = (i + 1 < n) ? normed_coeff[i + 1] : 0.0;
            C(0,i) = num_padded[i] - bn * a_i_plus_1;
        }

        
        double D = (m == n) ? bn : 0.0;


        StateSpace_System Trs;
        Trs.A=Accf;
        Trs.B=B;
        Trs.C=C; 
        Trs.D=D;
        return Trs;

        
        }


        // System identification methods
        std::optional<FirstOrderParams> identifyFirstOrder();
        std::optional<SecondOrderParams> identifySecondOrder();
    };

    // Time-domain response computations
    std::vector<double> computeStepResponse(const TransferFunction& tf, double duration, double timestep);
    std::vector<double> computeImpulseResponse(const TransferFunction& tf, double duration, double timestep);
    std::vector<double> computeRampResponse(const TransferFunction& tf, double duration, double timestep);

    // System characteristics
    double computeDCGain(const TransferFunction& tf);
    double computeBandwidth(const TransferFunction& tf);

    // First-order system characteristics
    double computeTimeConstant(const TransferFunction& tf);
    double computeSettlingTime(const TransferFunction& tf, double settlingPercentage = 0.02); // 2% by default
    double computeRiseTime(const TransferFunction& tf, double startPercentage = 0.1, double endPercentage = 0.9);
    //add max overshoot

    // Second-order system characteristics
    double computeNaturalFrequency(const TransferFunction& tf);
    double computeDampingRatio(const TransferFunction& tf);
    double computeOvershoot(const TransferFunction& tf);
    double computePeakTime(const TransferFunction& tf);
    std::vector<double> computeComplexPoleFreqDamping(const std::complex<double>& pole);

    // Performance metrics
    double computeISE(const TransferFunction& tf, double duration); // Integral Square Error
    double computeIAE(const TransferFunction& tf, double duration); // Integral Absolute Error
    double computeITAE(const TransferFunction& tf, double duration); // Integral Time Absolute Error

} 
#endif // TRANSFER_FUNCTIONS_HPP