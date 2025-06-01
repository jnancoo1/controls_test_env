#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "discrete_state_space.hpp"
#include <Eigen/SVD>

class Analysis {
public:
    Eigen::MatrixXd compute_controllability_matrix(const Discrete_StateSpace_System& System)
    {
        int n = System.A.rows();
        int m = System.B.cols();
        Eigen::MatrixXd controllability_mat(n, n * m);

        for (int i = 0; i < n; ++i) {
            controllability_mat.block(0, i * m, n, m) = System.A.pow(i) * System.B;
        }

        return controllability_mat;
    }

    bool is_controllable(const Discrete_StateSpace_System& System) 
    {
        Eigen::MatrixXd controllability_mat = compute_controllability_matrix(System);
        Eigen::FullPivLU<Eigen::MatrixXd> lu(controllability_mat);
        return lu.rank() == System.n_states;
    }

    Eigen::MatrixXd compute_observability_matrix(const Discrete_StateSpace_System& System)
    {
        int n = System.A.rows();
        int p = System.C.rows();
        Eigen::MatrixXd observability_mat(p * n, System.A.cols());

        for (int i = 0; i < n; ++i) {
            observability_mat.block(i * p, 0, p, System.A.cols()) = System.C * System.A.pow(i);
        }

        return observability_mat;
    }

    bool is_observable(const Discrete_StateSpace_System& System) 
    {
        Eigen::MatrixXd observability_mat = compute_observability_matrix(System);
        Eigen::FullPivLU<Eigen::MatrixXd> lu(observability_mat);
        return lu.rank() == System.n_states;
    }
	bool Linear_Stability_discrete(const Discrete_StateSpace_System& System)
	{
		Eigen::VectorXcd eigenvals=System.A.eigenvalues();
		for(int i=0;i<eigenvals.size();i++){
			
			
			if(std::abs(eigenvals[i])>=1.0){
				return false;
			}
		}
	return true;}

bool Linear_Stability_cont(const Discrete_StateSpace_System& System)
	{
    Eigen::VectorXcd eigenvals = System.A.eigenvalues();
    for (int i = 0; i < eigenvals.size(); ++i) {
        if (eigenvals[i].real() >= 0.0) {
            return false;  // Unstable or marginally stable
        }
    }
    return true;
	}



    bool is_stabalizable_cont(const Discrete_StateSpace_System &System){
        Eigen::VectorXcd eigs=System.A.eigenvalues();
        int f=eigs.size();
        std::vector<bool> unstable_flag;
        int a=System.A.rows();
        int b=System.A.cols();
        for(int i=0;i<eigs.size();i++){
            if(eigs[i].real()>0){
                Eigen::MatrixXcd eye(a,b);  
                eye=Eigen::MatrixXcd::Identity(a,b);
                Eigen::MatrixXcd PBH_part_1=eigs[i]*eye-System.A;
                int n=System.B.rows();
                int m=System.B.cols();
                Eigen::MatrixXcd PBH(n,n+m);
                PBH.block(0, 0, n, n) = PBH_part_1;
                PBH.block(0, n, n, m) = System.B.cast<std::complex<double>>();
                Eigen::FullPivLU<Eigen::MatrixXcd> lu(PBH);
                if(lu.rank()<n){
                    return false;
                }


            }
        
        }
        
        return true;

    }


    bool is_detectable_cont(const Discrete_StateSpace_System &System){
        Eigen::VectorXcd eigs=System.A.eigenvalues();
        int f=eigs.size();
        std::vector<bool> unstable_flag;
        int a=System.A.rows();
        int b=System.A.cols();
        for(int i=0;i<eigs.size();i++){
            if(eigs[i].real()>0){
                Eigen::MatrixXcd eye = Eigen::MatrixXcd::Identity(a, b);
                Eigen::MatrixXcd PBH_part_1=eigs[i]*eye-System.A;
                int n=System.A.rows();
                int p=System.C.rows();
                Eigen::MatrixXcd PBH(n + p, n);
                PBH.block(0, 0, n, n) = PBH_part_1;
                PBH.block(n, 0, p, n) = System.C.cast<std::complex<double>>();
                Eigen::FullPivLU<Eigen::MatrixXcd> lu(PBH);
                if(lu.rank()<n){
                    return false;
                }


            }
        
        }
        
        return true;

    }


    bool minimality_test_cont(const Discrete_StateSpace_System &System){

        if (is_controllable(System) && is_observable(System)) {

            return true;
        }
        else{
            return false;
        }
    }

    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> controllability_decomposition(const Discrete_StateSpace_System& System){
        Eigen::MatrixXd cont_mat=compute_controllability_matrix(System);
        int n=cont_mat.rows();
        Eigen::ColPivHouseholderQR <Eigen::MatrixXd> QR(cont_mat);
        Eigen::MatrixXd Q=QR.householderQ();
        int r=QR.rank();
        Eigen::MatrixXd T(n, n);
        T << Q.leftCols(r), Q.rightCols(n - r);

        Eigen::MatrixXd Aprime=(T.inverse())*System.A*T;
        Eigen::MatrixXd Bprime=(T.inverse())*System.B;
        Eigen::MatrixXd Cprime=System.C*T;

        Eigen::MatrixXd A_cc=Aprime.topLeftCorner(r,r);
        Eigen::MatrixXd A_cu=Aprime.topRightCorner(r,n-r);
        std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> Ans={A_cc,A_cu};

        return Ans;
    }



    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> observability_decomposition(const Discrete_StateSpace_System& System){

        Eigen::MatrixXd obs_mat=compute_observability_matrix(System);
        int n=System.A.rows();
        Eigen::JacobiSVD<Eigen::MatrixXd> SVD(obs_mat);

        Eigen::MatrixXd V=SVD.matrixV();
        Eigen::MatrixXd VT=V.transpose();
        double tol = 1e-9;
        int r = (SVD.singularValues().array() > tol).count();

        Eigen::MatrixXd T(n, n);
        T << V.leftCols(r), V.rightCols(n - r);

        Eigen::MatrixXd Aprime=(T.inverse())*System.A*T;
        Eigen::MatrixXd Bprime=(T.inverse())*System.B;
        Eigen::MatrixXd Cprime=System.C*T;

        Eigen::MatrixXd A_cc=Aprime.topLeftCorner(r,r);
        Eigen::MatrixXd A_cu=Aprime.topRightCorner(r,n-r);
        std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> Ans={A_oo,A_ou};

        return Ans;


    }

		
    Eigen::MatrixXd compute_controllability_gramian(const Discrete_StateSpace_System& System); 

    Eigen::MatrixXd compute_observability_gramian(const Discrete_StateSpace_System& System);
    
    double gramian_condition_number(const Discrete_StateSpace_System& System);



	private:	
	};

#endif
