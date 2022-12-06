#ifndef SRC_PARALLELHYBRIDALGORITHM_H_
#define SRC_PARALLELHYBRIDALGORITHM_H_

#include <iostream>
#include "Function.h"
#include "Solution.h"
#include "Hypercube.h"
#define   IS_TESTING   false
#define   IS_TESTING1   false


using namespace std;

class ParallelHybridAlgorithmParams {
public:
	Function* f=nullptr;	// objective function

	int N=10;
    double q=10;
    double p=0.7;
    int N_max=200;
    int I_max=10;
    unsigned int h0=6;
    unsigned int hn=8;
    bool update_h0_and_hn=false;
    bool accept_sub_optimal_Pareto_solutions=false;


};


class ParallelHybridAlgorithm { 
private:
	Function *f;
	ParallelHybridAlgorithmParams *params;
  	int N;
  	double q;
  	double p;
  	int N_max;
  	int I_max;
  	vector<Solution*> *U;
  	vector<Solution*> *P;
  	unsigned int h0;
  	unsigned int hn;
  	bool update_h0_and_hn;
  	int world_size;
  	int world_rank;
  	bool print_Paret_optimal_soputions=false;
  	int total_function_evaluation_count=0;
  	double run_time=0;// run time in miliseconds

  	void optimize_master();
  	void optimize_slave();
public:

	ParallelHybridAlgorithm(ParallelHybridAlgorithmParams *params)  {
		
		this->U = new vector<Solution*>;
    	this->P = new vector<Solution*>;


    	this->f=params->f;
		this->N=params->N;
		this->q=params->q;
		this->p=params->p;
		this->N_max=params->N_max;
		this->I_max=params->I_max;
		this->h0=params->h0;
		this->hn=params->hn;
		this->update_h0_and_hn=params->update_h0_and_hn;
		this->params=params;
		//cout<<"N : "<<N<<" q : "<<q<<" p : "<<p<<" N_max : "<<N_max<<" I_max : "<<I_max<<" h0 : "<<h0<<" hn : "<<hn<<endl;

	}

	~ParallelHybridAlgorithm() {
	    for (vector<Solution*>::iterator it = U->begin(); it != U->end();it++ ) {
	        delete *it;
	    }
	    delete P;
	    delete U;
	    delete f;
	    delete params;
		
	}

	vector<Solution*> *get_U(){
		return U;
	}
  	vector<Solution*> *get_P(){
  		return P;
  	}
  	int get_total_function_evaluation_count(){
  		return total_function_evaluation_count;
  	}

  	double  get_run_time(){
  		return run_time;
  	}

	void optimize();

	void put_to_Pareto_optimal_set(Solution* s);
	void put_to_selecton_functions_Pareto_optimal_set(Solution* s, vector<Solution*> *P_sel);

};



#endif /* SRC_PARALLELHYBRIDALGORITHM_H_ */
