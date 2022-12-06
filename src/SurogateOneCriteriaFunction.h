#ifndef SRC_SUROGATEONECRITERIAFUNCTION_H_
#define SRC_SUROGATEONECRITERIAFUNCTION_H_

#include <iostream>
#include <cassert>
#include "Function.h"
#include "Solution.h"

using namespace std;


class SurogateOneCriteriaFunction: public Function { 
protected:
	
	vector<Solution*> *_all_points;
    
    Function* _f;
    int _function_i;

  	Solution* _current_min_point;

	double* evaluate_func_values(double*);
  	



private:	
	SurogateOneCriteriaFunction()  {
		
	}

public:	

	Solution* get_current_min_point(){
		return _current_min_point;
	}

	virtual ~SurogateOneCriteriaFunction() {
		//cout<<"delete ~SurogateFunction"<<endl;
        delete _all_points;
    }
    vector<Solution*>* get_all_points(){
        return _all_points;
    }

	double* evaluate(double* decision_point);

	SurogateOneCriteriaFunction(Solution *s, int function_i)  {
		d = s->get_f()->get_d();//  dimension of decision space
		m = 1;//  number of objective functions
		double *upper_bound=s->get_f()->upper_bound; // decision  space bounds
	  	double *lower_bound=s->get_f()->lower_bound; // decision  space bounds
		intailize(d, upper_bound, lower_bound);
		_f=s->get_f();
  		_current_min_point = s;
  		//_current_min_point->set_foud_by_local_search(true);
  		_all_points = new vector<Solution*>;
  		_function_i=function_i;
	}

	virtual string to_string(){

			static string name ;
			name = "SurogateOneCriteriaFunction";
		

			return name;
		}

};




#endif /* SRC_SUROGATEFUNCTION_H_ */



