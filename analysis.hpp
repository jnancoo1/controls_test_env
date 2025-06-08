#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include <iostream>
#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include "discrete_state_space.hpp"
#include <Eigen/SVD>
#include <eigen3/unsupported/Eigen/KroneckerProduct>

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

        Eigen::MatrixXd A_oo=Aprime.topLeftCorner(r,r);
        Eigen::MatrixXd A_ou=Aprime.topRightCorner(r,n-r);
        std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> Ans={A_oo,A_ou};

        return Ans;


    }

		
    Eigen::MatrixXd compute_controllability_gramian(const Discrete_StateSpace_System& System) {
   
        Eigen::MatrixXd Q = System.B * System.B.transpose();
        int n = System.A.rows();
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
        Eigen::VectorXd vecQ = Eigen::Map<const Eigen::VectorXd>(Q.data(), Q.size());

        Eigen::MatrixXd kron1 = Eigen::kroneckerProduct(System.A, I);
        Eigen::MatrixXd kron2 = Eigen::kroneckerProduct(I, System.A);

        Eigen::VectorXd w = (kron1 + kron2).fullPivLu().solve(-vecQ);

        Eigen::MatrixXd W = Eigen::Map<Eigen::MatrixXd>(w.data(), n, n);

        return W;
}
 

    Eigen::MatrixXd compute_observability_gramian(const Discrete_StateSpace_System& System){
       
        Eigen::MatrixXd Q = (System.C) * System.C.transpose();
        int n = System.A.rows();
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
        Eigen::VectorXd vecQ = Eigen::Map<const Eigen::VectorXd>(Q.data(), Q.size());

        Eigen::MatrixXd kron1 = Eigen::kroneckerProduct(System.A, I);
        Eigen::MatrixXd kron2 = Eigen::kroneckerProduct(I, System.A);

        Eigen::VectorXd w = (kron1 + kron2).fullPivLu().solve(-vecQ);

        Eigen::MatrixXd W = Eigen::Map<Eigen::MatrixXd>(w.data(), n, n);

        return W;
    }

    std::vector<std::complex<double>> poles(const Discrete_StateSpace_System& System){

        Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver(System.A);
        Eigen::VectorXcd eigvals = eigen_solver.eigenvalues();
        std::vector<std::complex<double>> eigs(eigvals.data(), eigvals.data() + eigvals.size());

        return eigs;

    }

    std::vector<std::complex<double>> generate_z_grid(
    double r_min, double r_max,
    int r_samples,
    int theta_samples) 
    {
    std::vector<std::complex<double>> z_grid;
        z_grid.reserve(r_samples * theta_samples);

        for (int i = 0; i < r_samples; ++i) {
            double r = r_min + i * (r_max - r_min) / (r_samples - 1);
            for (int j = 0; j < theta_samples; ++j) {
                double theta = 2.0 * M_PI * j / theta_samples;
                std::complex<double> z = std::polar(r, theta);
                z_grid.push_back(z);
            }
        }
        return z_grid;
}

    std::vector<std::complex<double>> zeros(
        const Discrete_StateSpace_System& System,
        const double& r_min, const double& r_max,
        const int& r_samples,
        const int& theta_samples){

        
        int n =System.A.rows();
        int p= System.C.rows();
        int m= System.B.cols();
        Eigen::MatrixXcd R(n + p, n + m);

        std::vector<std::complex<double>> zgrid=generate_z_grid(r_min,r_max,r_samples,theta_samples);
        std::vector<std::complex<double>> zeros;
        Eigen::MatrixXcd zI_minus_A;

        Eigen::MatrixXcd eye=Eigen::MatrixXd::Identity(n, n);
        int r;
        double tol = 1e-9;
        for(int i=0;i<zgrid.size();i++){
            zI_minus_A=(zgrid[i]*eye)-(System.A.cast<std::complex<double>>);
            R.block(0, 0, n, n) = zI_minus_A;
            R.block(0, n, n, m) = -System.B.cast<std::complex<double>>();
            R.block(n, 0, p, n) = System.C.cast<std::complex<double>>();
            R.block(n, n, p, m) = System.D.cast<std::complex<double>>();

            Eigen::JacobiSVD<Eigen::MatrixXcd>SVD(R);
            r = (SVD.singularValues().array() > tol).count();
            if(r<n+p){

                zeros.push_back(zgrid[i]);
            }


        }

        return zeros;
    }

    Discrete_StateSpace_System Kalman_Decomp( const Discrete_StateSpace_System& System){


        std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> Cont,Obs;
        
        Obs=observability_decomposition(System);
        Cont=controllability_decomposition(System);
        
        Eigen::MatrixXd Obs_subspace=std::get<0>(Obs);
        Eigen::MatrixXd Cont_subspace=std::get<0>(Cont);
        Eigen::MatrixXd Diff(Obs_subspace.rows(), Obs_subspace.cols() + Cont_subspace.cols());
        Diff << Obs_subspace, -Cont_subspace;
        Eigen::MatrixXd Null_S=Diff.FullPivLU.kernel();


        int k=Obs_subspace.cols()

        Eigen::MatrixXd Null_S_S=Null_S.toprows(k);
        Eigen::MatrixXd Basis_1=Null_S_S*Obs_subspace;


    }

	private:	
	};

#endif
