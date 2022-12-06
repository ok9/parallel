#ifndef SRC_OPTIMIZEDFUNCTION_H_
#define SRC_OPTIMIZEDFUNCTION_H_

#include <iostream>
#include <cassert>
#include "Function.h"

using namespace std;


class OptimizedFunction: public Function { 
protected:
	
	double* evaluate_func_values(double*);
  	
public:

	OptimizedFunction()  {
		d=2;//  dimension of decision space
		m=2;//  number of objective functions
		double upper_bound[]={2,3}; // decision  space bounds
	  	double lower_bound[]={-1,-4}; // decision  space bounds
		intailize(d, upper_bound, lower_bound);
	}

	~OptimizedFunction() {
		
	}
	

	

	
private:
	
  	

};









#endif /* SRC_FUNCTION_H_ */