#ifndef TRANSFER_FUNCTIONS_HPP
#define TRANSFER_FUNCTIONS_HPP

#include <vector>
#include <complex>
#include <optional>
#include<Eigen/Dense>

namespace TransferFunctions {

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

        // Evaluate the transfer function at a given frequency (s-domain)
        std::complex<double> evaluate(std::complex<double> s) const;

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

        // Get the zeros of the transfer function
        std::vector<std::complex<double>> getZeros() const;

        // Compute the frequency response at a given frequency
        std::complex<double> frequencyResponse(double omega) const;

        // System identification methods
        std::optional<FirstOrderParams> identifyFirstOrder() const;
        std::optional<SecondOrderParams> identifySecondOrder() const;
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