#ifndef SRC_FUNCTION_H_
#define SRC_FUNCTION_H_

#include <iostream>

using namespace std;



class Function { 
protected:
	int d; //  dimension of decision space

  	
  	int m; //  number of objective functions

  	
  	virtual double* evaluate_func_values(double*) = 0;
public:

	Function()  {
		
		
	}

	virtual ~Function() {
		//cout<<"delete ~Function"<<endl;
		delete upper_bound;
		delete lower_bound;
		delete max_values_of_objective_functions;
		delete min_values_of_objective_functions;
	}
	
	

	virtual double* evaluate(double*);
	void intailize(int d,  double upper_bound[], double lower_bound[]);

	int get_d() {
      return d;
    }
    int get_m() {
      return m;
    }


    virtual string to_string(){

			static string name ;
			name = "Function";
		

			return name;
		}

    
    double *upper_bound=nullptr; // decision  space bounds
  	double *lower_bound=nullptr;
    
    double *max_values_of_objective_functions=nullptr;
  	double *min_values_of_objective_functions=nullptr; 

};





#endif /* SRC_FUNCTION_H_ */
