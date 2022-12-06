
#include "SurogateFunction.h"
#include <cassert>
#include <iostream>
double* SurogateFunction::evaluate(double* decision_point) {
    //cout<<"HELLO!"<<endl;
    /*int d; //  dimension of decision space
    double *upper_bound=nullptr; // decision  space bounds
    double *lower_bound=nullptr;
*/
    for(int i=0;i<d;i++){
        assert(decision_point[i]<=1 && decision_point[i]>=0);
    }


    /*int m; //  number of objective functions
    double *max_values_of_objective_functions=nullptr;
    double *min_values_of_objective_functions=nullptr; */
    double *values_of_objective_functions = evaluate_func_values(decision_point);
    
    return values_of_objective_functions;


    //return new double[3];
    
}


double* SurogateFunction::evaluate_func_values(double* decision_point){

    /*double * func_values= new double[m];
    //_current_point;
    func_values[0]=current_value;
    return func_values;*/
    Solution* p1 = new Solution(_f);
    p1->set_decision_point(decision_point);
    if(_current_point->equals(p1)){
        //cout<<"LYGU";
        delete p1;
        double *data= new double[m];
        data[0]=current_value;      
        return data;
    }else
        delete p1;        
    
    //Solution* p = new Solution(decision_point); 
    Solution* p = new Solution(_f);
    p->evaluate_function(decision_point);
    _all_points->push_back(p);
    //cout<<"---------------p->to_string()"<<p->to_string()<<endl;
    //cout<<"_all_points->size()"<<_all_points->size()<<endl;
    
    if(_current_point->is_dominated_by(p)){
        current_value--;
        _current_point=p;
        _current_point->set_foud_by_local_search(true);
    }


    double *data= new double[m];
    data[0]=current_value;      
    return data;
}

