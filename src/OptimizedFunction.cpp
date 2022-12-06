
#include "OptimizedFunction.h"
#include <cassert>
#include <iostream>

/*double upper_bound[]={2,3}; // decision  space bounds
  double lower_bound[]={-1,-4}; // decision  space bounds*/

	double* OptimizedFunction::evaluate_func_values(double* decision_point){
		//static int i=0;
		//i++;
		//cout<<"i : "<<i<<endl;
		double * func_values= new double[m];
		func_values[0]=(decision_point[0]-0.5)*(decision_point[0]-0.5)+(decision_point[1]-1)*(decision_point[1]-1);
		func_values[1]=(decision_point[0]+1)*(decision_point[0]+1)+(decision_point[1]+1)*(decision_point[1]+1);
		//func_values[0]=100*(decision_point[0]*decision_point[0]+decision_point[1]*decision_point[1]);
		//func_values[0]=decision_point[0]+decision_point[1];
		return func_values;
	}	
















