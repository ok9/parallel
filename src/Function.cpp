
#include "Function.h"
#include <cassert>
#include <iostream>





double* Function::evaluate(double* decision_point) {
	//cout<<"HELLO2222"<<endl;
	/*int d; //  dimension of decision space
  	double *upper_bound=nullptr; // decision  space bounds
  	double *lower_bound=nullptr;
*/
	double *decision_point_not_normalized = new double[d];
	for(int i=0;i<d;i++){
		assert(decision_point[i]<=1 && decision_point[i]>=0);
		decision_point_not_normalized[i] = lower_bound[i] + (upper_bound[i] - lower_bound[i])*decision_point[i];
	}


	/*int m; //  number of objective functions
  	double *max_values_of_objective_functions=nullptr;
  	double *min_values_of_objective_functions=nullptr; */
  	double *values_of_objective_functions = evaluate_func_values(decision_point_not_normalized);
	
	delete decision_point_not_normalized;
	return values_of_objective_functions;


	//return new double[3];
	
}





void Function::intailize(int d,  double upper_bound[], double lower_bound[]){

  	this->upper_bound=new double[d];
  	this->lower_bound=new double[d];
  	for(int i=0;i<d;i++){
  		this->upper_bound[i]=upper_bound[i];
  		this->lower_bound[i]=lower_bound[i];
  	}
  	
  	
}


