#ifndef SRC_SUROGATEFUNCTION_H_
#define SRC_SUROGATEFUNCTION_H_

#include <iostream>
#include <cassert>
#include "Function.h"
#include "Solution.h"

using namespace std;


class SurogateFunction: public Function { 
protected:
	
	vector<Solution*> *_all_points;
    
    Function* _f;
    double current_value=0;
  	Solution* _current_point;

	double* evaluate_func_values(double*);
  	



private:	
	SurogateFunction()  {
		
	}

public:	

	virtual ~SurogateFunction() {
		//cout<<"delete ~SurogateFunction"<<endl;
        delete _all_points;
    }
    vector<Solution*>* get_all_points(){
        return _all_points;
    }

	double* evaluate(double* decision_point);

	SurogateFunction(Solution *s)  {
		d = s->get_f()->get_d();//  dimension of decision space
		m = 1;//  number of objective functions
		double *upper_bound=s->get_f()->upper_bound; // decision  space bounds
	  	double *lower_bound=s->get_f()->lower_bound; // decision  space bounds
		intailize(d, upper_bound, lower_bound);
		_f=s->get_f();
  		_current_point = s;
  		_current_point->set_foud_by_local_search(true);
  		_all_points = new vector<Solution*>;
	}

};




#endif /* SRC_SUROGATEFUNCTION_H_ */