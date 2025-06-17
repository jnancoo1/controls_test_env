#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include <iostream>
#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include "discrete_state_space.hpp"
#include <Eigen/SVD>
#include <eigen3/unsupported/Eigen/KroneckerProduct>

/**
 * @brief A class for analyzing properties of discrete-time state space systems
 * 
 * This class provides methods for analyzing fundamental system properties including:
 * - Controllability: ability to drive the system to any desired state
 * - Observability: ability to determine any initial state from output measurements
 * - Stability: behavior of the system's free response
 * - Stabilizability: ability to stabilize unstable modes through feedback
 *
 * For a discrete-time state space system:     * \f{eqnarray*}{
     * x[k+1] &=& Ax[k] + Bu[k] \\
     * y[k] &=& Cx[k] + Du[k]
     * \f}
 */
class Analysis {
public:
    /**
     * @brief Computes the controllability matrix for a discrete-time state space system
     * 
     * The controllability matrix is defined as:
     * \f[
     * \mathcal{C} = [B \quad AB \quad A^2B \quad \cdots \quad A^{n-1}B]
     * \f]
     * 
     * where:
     * - n is the number of states
     * - A is the state matrix
     * - B is the input matrix
     * 
     * This matrix has size n × nm, where:
     * - n is the number of states
     * - m is the number of inputs
     * 
     * @param System The discrete state space system to analyze
     * @return Eigen::MatrixXd The controllability matrix
     */
    static Eigen::MatrixXd compute_controllability_matrix(const Discrete_StateSpace_System& System)
    {
        int n = System.A.rows();
        int m = System.B.cols();
        Eigen::MatrixXd controllability_mat(n, n * m);

        for (int i = 0; i < n; ++i) {
            controllability_mat.block(0, i * m, n, m) = System.A.pow(i) * System.B;
        }

        return controllability_mat;
    }

    /**
     * @brief Checks if the system is controllable
     * 
     * A system is controllable if and only if its controllability matrix has full rank:
     * \f[
     * \text{rank}([B \quad AB \quad A^2B \quad \cdots \quad A^{n-1}B]) = n
     * \f]
     * 
     * where n is the number of states.
     * 
     * Controllability implies that there exists an input sequence that can transfer
     * the system from any initial state to any final state in finite time.
     * 
     * @param System The discrete state space system to analyze
     * @return bool True if rank(C) = n, false otherwise
     */
    bool is_controllable(const Discrete_StateSpace_System& System) 
    {
        Eigen::MatrixXd controllability_mat = compute_controllability_matrix(System);
        Eigen::FullPivLU<Eigen::MatrixXd> lu(controllability_mat);
        return lu.rank() == System.n_states;
    }

    /**
     * @brief Computes the observability matrix for a discrete-time state space system
     * 
     * The observability matrix is defined as:
     * \f[
     * \mathcal{O} = \begin{bmatrix} 
     * C \\
     * CA \\
     * CA^2 \\
     * \vdots \\
     * CA^{n-1}
     * \end{bmatrix}
     * \f]
     * 
     * This matrix has size pn × n, where:
     * - n is the number of states
     * - p is the number of outputs
     * 
     * @param System The discrete state space system to analyze
     * @return Eigen::MatrixXd The observability matrix
     */
    static Eigen::MatrixXd compute_observability_matrix(const Discrete_StateSpace_System& System)
    {
        int n = System.A.rows();
        int p = System.C.rows();
        Eigen::MatrixXd observability_mat(p * n, System.A.cols());

        for (int i = 0; i < n; ++i) {
            observability_mat.block(i * p, 0, p, System.A.cols()) = System.C * System.A.pow(i);
        }

        return observability_mat;
    }

    /**
     * @brief Checks if the system is observable
     * 
     * A system is observable if and only if its observability matrix has full rank:
     * \f[
     * \text{rank}\begin{bmatrix} 
     * C \\
     * CA \\
     * CA^2 \\
     * \vdots \\
     * CA^{n-1}
     * \end{bmatrix} = n
     * \f]
     * 
     * Observability implies that the initial state can be determined from
     * knowledge of the input and output over a finite time interval.
     * 
     * @param System The discrete state space system to analyze
     * @return bool True if rank(O) = n, false otherwise
     */
    bool is_observable(const Discrete_StateSpace_System& System) 
    {
        Eigen::MatrixXd observability_mat = compute_observability_matrix(System);
        Eigen::FullPivLU<Eigen::MatrixXd> lu(observability_mat);
        return lu.rank() == System.n_states;
    }
	
    /**
     * @brief Checks the stability of a discrete-time linear system
     * 
     * For a discrete-time system, stability requires all eigenvalues
     * to lie inside the unit circle:
     * \f[
     * |\lambda_i(A)| < 1 \quad \forall i
     * \f]
     * 
     * where \f$\lambda_i(A)\f$ are the eigenvalues of matrix A.
     * 
     * Stability types:
     * - |λ| < 1: Asymptotically stable
     * - |λ| = 1: Marginally stable
     * - |λ| > 1: Unstable
     * 
     * @param System The discrete state space system to analyze
     * @return bool True if all eigenvalues have magnitude less than 1
     */
    bool Linear_Stability_discrete(const Discrete_StateSpace_System& System)
	{
		Eigen::VectorXcd eigenvals=System.A.eigenvalues();
		for(int i=0;i<eigenvals.size();i++){
			
			
			if(std::abs(eigenvals[i])>=1.0){
				return false;
			}
		}
	return true;}

    /**
     * @brief Checks the stability of a continuous-time linear system
     * 
     * For a continuous-time system, stability requires all eigenvalues
     * to lie in the left half-plane:
     * \f[
     * \text{Re}(\lambda_i(A)) < 0 \quad \forall i
     * \f]
     * 
     * where \f$\lambda_i(A)\f$ are the eigenvalues of matrix A.
     * 
     * Stability types:
     * - \f$\text{Re}(\lambda) < 0\f$: Asymptotically stable
     * - \f$\text{Re}(\lambda) = 0\f$: Marginally stable
     * - \f$\text{Re}(\lambda) > 0\f$: Unstable
     * 
     * @param System The discrete state space system to analyze
     * @return bool True if all eigenvalues have negative real parts
     */
    static bool Linear_Stability_cont(const Discrete_StateSpace_System& System)
	{
    Eigen::VectorXcd eigenvals = System.A.eigenvalues();
    for (int i = 0; i < eigenvals.size(); ++i) {
        if (eigenvals[i].real() >= 0.0) {
            return false;  // Unstable or marginally stable
        }
    }
    return true;
	}



    /**
     * @brief Checks if the continuous-time system is stabilizable
     * 
     * A system is stabilizable if all uncontrollable modes are stable.
     * This is checked using the Popov-Belevitch-Hautus (PBH) test:
     * \f[
     * \text{rank}[sI - A \quad B] = n
     * \f]
     * 
     * for all eigenvalues s with \f$\text{Re}(s) \geq 0\f$, where:
     * - n is the number of states
     * - A is the state matrix
     * - B is the input matrix
     * 
     * If this condition is satisfied for all unstable eigenvalues,
     * then the system is stabilizable.
     * 
     * @param System The discrete state space system to analyze
     * @return bool True if the system is stabilizable
     */
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


    /**
     * @brief Checks if the continuous-time system is detectable
     * 
     * A system is detectable if all unobservable modes are stable.
     * Uses the Popov-Belevitch-Hautus (PBH) test for each unstable eigenvalue.
     * 
     * @param System The discrete state space system to analyze
     * @return bool True if the system is detectable, false otherwise
     */
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


    /**
     * @brief Performs the minimality test for continuous-time systems
     * 
     * A system is minimal if it is both controllable and observable.
     * 
     * @param System The discrete state space system to analyze
     * @return bool True if the system is minimal, false otherwise
     */
    bool minimality_test_cont(const Discrete_StateSpace_System &System){

        if (is_controllable(System) && is_observable(System)) {

            return true;
        }
        else{
            return false;
        }
    }

    /**
     * @brief Computes the controllability decomposition of a discrete-time state space system
     * 
     * This method computes the controllability matrix and performs QR decomposition
     * to find the invariant subspaces of the system.
     * 
     * @param System The discrete state space system to analyze
     * @return std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> The controllable and uncontrollable subspaces
     */
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



    /**
     * @brief Computes the observability decomposition of a discrete-time state space system
     * 
     * This method computes the observability matrix and performs SVD to find the
     * invariant subspaces of the system.
     * 
     * @param System The discrete state space system to analyze
     * @return std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> The observable and unobservable subspaces
     */
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

		
    /**
     * @brief Computes the controllability gramian for a discrete-time state space system
     * 
     * The controllability gramian is computed using the formula:
     * Wc = ∫ (Φ_A(t)B) (Φ_A(t)B)^T dt
     * where Φ_A(t) is the state transition matrix.
     * 
     * @param System The discrete state space system to analyze
     * @return Eigen::MatrixXd The controllability gramian
     */
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
 

    /**
     * @brief Computes the observability gramian for a discrete-time state space system
     * 
     * The observability gramian is computed using the formula:
     * Wo = ∫ (CΦ_A(t))^T (CΦ_A(t)) dt
     * where Φ_A(t) is the state transition matrix.
     * 
     * @param System The discrete state space system to analyze
     * @return Eigen::MatrixXd The observability gramian
     */
    Eigen::MatrixXd compute_observability_gramian(const Discrete_StateSpace_System& System){
       
        Eigen::MatrixXd Q = (System.C.transpose()) * System.C;
        int n = System.A.rows();
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
        Eigen::VectorXd vecQ = Eigen::Map<const Eigen::VectorXd>(Q.data(), Q.size());

        Eigen::MatrixXd kron1 = Eigen::kroneckerProduct(System.A.transopse(), I);
        Eigen::MatrixXd kron2 = Eigen::kroneckerProduct(I, System.A.transpose());

        Eigen::VectorXd w = (kron1 + kron2).fullPivLu().solve(-vecQ);

        Eigen::MatrixXd W = Eigen::Map<Eigen::MatrixXd>(w.data(), n, n);

        return W;
    }

    /**
     * @brief Computes the poles of the system (eigenvalues of A matrix)
     * 
     * @param System The discrete state space system to analyze
     * @return std::vector<std::complex<double>> The eigenvalues of the A matrix
     */
    std::vector<std::complex<double>> poles(const Discrete_StateSpace_System& System){

        Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver(System.A);
        Eigen::VectorXcd eigvals = eigen_solver.eigenvalues();
        std::vector<std::complex<double>> eigs(eigvals.data(), eigvals.data() + eigvals.size());

        return eigs;

    }

    /**
     * @brief Generates a grid of points in the z-plane for root locus analysis
     * 
     * @param r_min The minimum radius for the grid
     * @param r_max The maximum radius for the grid
     * @param r_samples The number of radial samples
     * @param theta_samples The number of angular samples
     * @return std::vector<std::complex<double>> The grid of points in the z-plane
     */
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

    /**
     * @brief Finds the zeros of the transfer function (poles of the closed-loop system)
     * 
     * @param System The discrete state space system to analyze
     * @param r_min The minimum radius for the grid
     * @param r_max The maximum radius for the grid
     * @param r_samples The number of radial samples
     * @param theta_samples The number of angular samples
     * @return std::vector<std::complex<double>> The zeros of the transfer function
     */
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

    /**
     * @brief Performs Kalman decomposition for state estimation
     * 
     * This method computes the observability and controllability decompositions,
     * and then combines them to form a new state space realization.
     * 
     * @param System The discrete state space system to analyze
     * @return Discrete_StateSpace_System The decomposed state space system
     */
    Discrete_StateSpace_System Kalman_Decomp( const Discrete_StateSpace_System& System){

        std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> Cont,Obs;
        
        Obs=observability_decomposition(System);
        Cont=controllability_decomposition(System);
        
        Eigen::MatrixXd Obs_subspace=std::get<0>(Obs);
        Eigen::MatrixXd Cont_subspace=std::get<0>(Cont);
        Eigen::MatrixXd Diff(Obs_subspace.rows(), Obs_subspace.cols() + Cont_subspace.cols());
        Diff << Obs_subspace, -Cont_subspace;
        Eigen::FullPivLU<Eigen::MatrixXd> lu1(Diff);
        Eigen::MatrixXd Null_S = lu1.kernel();


        int k=Obs_subspace.cols();

        Eigen::MatrixXd Null_S_S=Null_S.topRows(k);
        Eigen::MatrixXd Basis_1=Null_S_S*Obs_subspace;


        Eigen::MatrixXd Obs_ortho=(Obs_subspace.transpose()).fullPivLu().kernel();
        Eigen::MatrixXd Cont_ortho=(Cont_subspace.transpose()).fullPivLu().kernel();

        Eigen::MatrixXd Diff2(Obs_ortho.rows(), Obs_ortho.cols() + Cont_subspace.cols());
        Diff2 << Obs_ortho, -Cont_subspace;
        Eigen::MatrixXd Null_S2=Diff2.fullPivLu().kernel();

        int l=Obs_ortho.cols();

        
        Eigen::MatrixXd Null_S_S2=Null_S2.topRows(l);
        Eigen::MatrixXd Basis_2=Null_S_S2*Obs_ortho;


        Eigen::MatrixXd Diff3(Obs_subspace.rows(), Obs_subspace.cols() + Cont_ortho.cols());
        Diff3 << Obs_subspace, -Cont_ortho;
        Eigen::MatrixXd Null_S3=Diff3.fullPivLu().kernel();

                
        Eigen::MatrixXd Null_S_S3=Null_S3.topRows(k);
        Eigen::MatrixXd Basis_3=Null_S_S3*Obs_subspace;


        
        Eigen::MatrixXd Diff4(Obs_ortho.rows(), Obs_subspace.cols() + Cont_ortho.cols());
        Diff4<< Obs_ortho, -Cont_ortho;
        Eigen::MatrixXd Null_S4=Diff3.fullPivLu().kernel();

                        
        Eigen::MatrixXd Null_S_S4=Null_S4.topRows(l);
        Eigen::MatrixXd Basis_4=Null_S_S4*Obs_ortho;

        int n=Basis_1.rows();
        Eigen::MatrixXd T(n, n);
        T << Basis_1, Basis_2, Basis_3, Basis_4;

        if (T.fullPivLu().isInvertible()) {
        } 
        else {
            Eigen::MatrixXd Q = QR.householderQ();
    T = Q.leftCols(n);  // Truncate if necessary

            }

            Eigen::MatrixXd T_inv = T.inverse();


            Eigen::Matrix A_new=T_inv*System.A*T;
            Eigen::Matrix B_new=T_inv*System.B;
            Eigen::Matrix C_new=System.C*T;
            Eigen::Matrix D_new=System.D;


            Discrete_StateSpace_System Decomp;
            Decomp.A=A_new;
            Decomp.B=B_new;
            Decomp.C=C_new;
            Decomp.D=D_new;

            return Decomp;


    }

	private:	
	};

#endif
