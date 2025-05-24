

#ifndef DISCRETE_STATE_SPACE
#define DISCRETE_STATE_SPACE


#include <iostream>
#include <cmath>
#include <Eigen/Dense>



class Discrete_StateSpace_System
{
    public:

        Eigen::MatrixXd A,B,C,D;
            Discrete_StateSpace_System(const int& n_states=3,const int& n_inputs=1, const int& n_outputs=3){
                A.resize(n_states,n_states);
                B.resize(n_states,n_inputs);
                C.resize(n_outputs,n_states);
                D.resize(n_outputs,n_inputs);
                A.setZero();
                B.setZero();
                C.setZero();
                D.setZero();
        }

        Eigen::VectorXd  simulate_step (const Eigen::VectorXd& x, const Eigen::VectorXd& u){

            if(x.size()==A.rows() && u.size()==B.cols()){
                Eigen::VectorXd  X_next= A*x+B*u;

                return X_next;
            }

            throw std::invalid_argument("Invalid inputs");
        }

        Eigen::VectorXd  get_output(const Eigen::VectorXd& x, const Eigen::VectorXd& u){
        
            if(x.size()==A.rows() && u.size()==B.cols()){
                Eigen::VectorXd  output= C*x+D*u;

                return output;
            }   
                    throw std::invalid_argument("Invalid inputs");


        }
        
        Eigen::MatrixXd sim_seq_outputs(const Eigen::VectorXd& x_init,const int& no_of_steps,const Eigen::MatrixXd& u_seq,

                                        const int& n_outputs){
        Eigen::VectorXd X_next=x_init;
        Eigen::MatrixXd Out(n_outputs, no_of_steps);
        if(u_seq.cols() == no_of_steps && u_seq.rows() == B.cols()){

            for(int k=0;k<no_of_steps;k++){

                Out.col(k)=get_output(X_next,u_seq.col(k));
                X_next=simulate_step(X_next,u_seq.col(k));


            }

            return Out;
        }
        throw std::invalid_argument("Invalid inputs");


    }

        Eigen::MatrixXd sim_seq_states(const Eigen::VectorXd& x_init,const int& no_of_steps,const Eigen::MatrixXd& u_seq,
        const int& n_states){
        int i= no_of_steps;
        int checker=0;
        Eigen::VectorXd X_next=x_init;
        Eigen::MatrixXd Out(n_states, no_of_steps);
        if(u_seq.cols() == no_of_steps && u_seq.rows() == B.cols()){

            for(int k=0;k<no_of_steps;k++){

                Out.col(k)=X_next;
                X_next=simulate_step(X_next,u_seq.col(k));


            }

            return Out;
        }
        throw std::invalid_argument("Invalid inputs");

    }
    private:

};






#endif
