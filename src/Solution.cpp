
#include "Solution.h"
#include <cassert>
#include <iostream>
#include <math.h> 
#include <algorithm>
#include "Hypercube.h"



void Solution::evaluate_function() {
	evaluation_of_objective_functions=f->evaluate(decision_point);
}
void Solution::evaluate_function(vector<double> loc){
	assert (decision_point==nullptr);// decision_point shoud be generated only once for solution object
	int d=f->get_d();
	decision_point=new double [d];
	for(int i=0;i<d;i++){
		decision_point[i]=loc.at(i);
	}
	evaluation_of_objective_functions=f->evaluate(decision_point);

}
void Solution::evaluate_function(double *loc){
	assert (decision_point==nullptr);// decision_point shoud be generated only once for solution object
	int d=f->get_d();
	decision_point=new double [d];
	for(int i=0;i<d;i++){
		decision_point[i]=loc[i];
	}
	evaluation_of_objective_functions=f->evaluate(decision_point);

}

void Solution::set_decision_point(vector<double> loc){
	assert (decision_point==nullptr);// decision_point shoud be generated only once for solution object
	int d=f->get_d();
	decision_point=new double [d];
	for(int i=0;i<d;i++){
		decision_point[i]=loc.at(i);
	}
}


void Solution::set_decision_point(double* loc){
	assert (decision_point==nullptr);// decision_point shoud be generated only once for solution object
	int d=f->get_d();
	decision_point=new double [d];
	for(int i=0;i<d;i++){
		decision_point[i]=loc[i]; 
	}
}


void Solution::set_decision_point_and_evaluation_of_objective_functions(double *data){
	assert (decision_point==nullptr);// decision_point shoud be generated only once for solution object
	int d=f->get_d();
	decision_point=new double [d];
	int j=0;
	for(int i=0;i<d;i++){
		decision_point[i]=data[j++]; 
	}

	assert (evaluation_of_objective_functions==nullptr);// evaluation_of_objective_functions shoud be evaluated only once for solution object	
	int m=f->get_m();
	evaluation_of_objective_functions=new double [m];
	for(int i=0;i<m;i++){
		evaluation_of_objective_functions[i]=data[j++]; 
	}

	foud_by_local_search=true;

}





void Solution::generate_decision_point(){
	assert (decision_point==nullptr);// decision_point shoud be generated only once for solution object
	int d=f->get_d();
	decision_point=new double [d];
	for(int i=0;i<d;i++){
		decision_point[i]=generate_random_number();
	}
}

void Solution::generate_decision_point_in_hypercube(Hypercube *h){
	assert (decision_point==nullptr);// decision_point shoud be generated only once for solution object
	decision_point=h->generate_decision_point_inside();
}




double Solution::generate_random_number(){

	static std::random_device rd;
    static std::default_random_engine eng(rd());
    static std::uniform_real_distribution<double> distr(0, 1);

    return distr(eng) ;

}



bool Solution::is_dominated_by(Solution *s){
	bool is_dominated=true;
	int m = f->get_m();
	for(int i=0;i < m;i++){
		if(this->evaluation_of_objective_functions[i] < s->evaluation_of_objective_functions[i]){
			is_dominated=false;
			break;
		}
	}
	
	if(is_dominated){
		bool all_equal=true;
		for(int i=0;i < m;i++){
			if(this->evaluation_of_objective_functions[i] != s->evaluation_of_objective_functions[i]){
				all_equal=false;
				break;
			}
		}
		if(all_equal){
			return false;
		}
	}
	return is_dominated;
}




bool Solution::selection_functions_is_dominated_by(Solution *s){
	bool is_dominated=true;
	
	if(-this->selection_function1_minimal_distance_to_solution_decision < -s->selection_function1_minimal_distance_to_solution_decision || 
		this->selection_function2_minimal_distance_to_Pareto_solution_function_value < s->selection_function2_minimal_distance_to_Pareto_solution_function_value){
		is_dominated=false;
	}

	if(is_dominated){
		bool all_equal=true;
		if(this->selection_function1_minimal_distance_to_solution_decision != s->selection_function1_minimal_distance_to_solution_decision || 
			this->selection_function2_minimal_distance_to_Pareto_solution_function_value != s->selection_function2_minimal_distance_to_Pareto_solution_function_value){
			all_equal=false;
		}
		
		if(all_equal){
			return false;
		}
	}
	return is_dominated;
}













void Solution::calculate_selection_functions(vector<Solution*> *U, vector<Solution*> *P){
	
	Solution *closest_solution=nullptr;
	selection_function1_minimal_distance_to_solution_decision=-1;
	for (vector<Solution*>::iterator it = U->begin(); it != U->end(); ++it) {
		double distance;
		distance= this->distance_in_decision_space(*it);
		//cout<<endl<<"this->distance_in_decision_space(*it) : "<<this->distance_in_decision_space(*it)<<endl;
		//cout<<"*it->to_string() : "<<(*it)->to_string()<<endl;
		if(selection_function1_minimal_distance_to_solution_decision<0){
			selection_function1_minimal_distance_to_solution_decision=distance;
			closest_solution=*it;
		}
		else{
			if(selection_function1_minimal_distance_to_solution_decision>distance){
				selection_function1_minimal_distance_to_solution_decision=distance;
				closest_solution=*it;
			}
		}
	}
	//cout<<"closest_solution: "<<endl;
	//cout<<closest_solution->to_string();

	selection_function2_minimal_distance_to_Pareto_solution_function_value=-1;
	for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
		double distance;
		distance= closest_solution->distance_in_objective_space(*it);
		if(selection_function2_minimal_distance_to_Pareto_solution_function_value<0){
			selection_function2_minimal_distance_to_Pareto_solution_function_value=distance;
			//cout<<"closest_solution->distance_in_objective_space(*it) : "<<closest_solution->distance_in_objective_space(*it)<<endl;
			//cout<<"(*it)->to_string() : "<<(*it)->to_string()<<endl;
		}
		else{
			if(selection_function2_minimal_distance_to_Pareto_solution_function_value>distance){
				selection_function2_minimal_distance_to_Pareto_solution_function_value=distance;
				//cout<<"closest_solution->distance_in_objective_space(*it) : "<<closest_solution->distance_in_objective_space(*it)<<endl;
				//cout<<"(*it)->to_string() : "<<(*it)->to_string()<<endl;
			}
		}
	}

	/*if(std::find(P->begin(), P->end(), closest_solution) != P->end()) {
        //cout<<"closest_solution is in P"<<endl;
        assert(selection_function2_minimal_distance_to_Pareto_solution_function_value==0);
    } else {
        //cout<<"closest_solution is not in P"<<endl;
        if(selection_function2_minimal_distance_to_Pareto_solution_function_value<0)
        	cout<<"NEGATIVE"<<selection_function2_minimal_distance_to_Pareto_solution_function_value<<endl;
        if(selection_function2_minimal_distance_to_Pareto_solution_function_value==0){
        	cout<<"ZERO"<<selection_function2_minimal_distance_to_Pareto_solution_function_value<<endl;
        	cout<<"(*std::find(P->begin(), P->end(), closest_solution))->to_string() : "<<(*std::find(P->begin(), P->end(), closest_solution))->to_string();
        	cout<<"(*std::find(P->begin(), P->end(), closest_solution)) : "<<(*std::find(P->begin(), P->end(), closest_solution));
        	cout<<"(*(P->end()))->to_string() : "<<(*(P->end()))->to_string();
        	cout<<"(*(P->end())) : "<<(*(P->end()));
        	//cout<<"(std::find(P->begin(), P->end(), closest_solution))"<<(std::find(P->begin(), P->end(), closest_solution));
        	cout<<"P->size() : "<<P->size()<<endl;
        }
        assert(selection_function2_minimal_distance_to_Pareto_solution_function_value>0);
    }*/


	assert(selection_function2_minimal_distance_to_Pareto_solution_function_value>=0);
	assert(selection_function1_minimal_distance_to_solution_decision>=0);
}

double Solution::distance_in_decision_space(Solution *s){
	double distance=0;
	int d=f->get_d();
	for(int i=0; i<d; i++){
		distance+=(this->decision_point[i]-s->decision_point[i])*(this->decision_point[i]-s->decision_point[i]);
	}
	return sqrt(distance);
}




double Solution::distance_in_objective_space(Solution *s){
	int m=f->get_m();
	double this_normalized_func_values[m];
	double s_normalized_func_values[m];
	double *max_values_of_objective_functions=this->f->max_values_of_objective_functions;
	double *min_values_of_objective_functions=this->f->min_values_of_objective_functions;
	
	bool all_not_equal=true;
	for(int i=0;i<m;i++){
		if(max_values_of_objective_functions[i]==min_values_of_objective_functions[i]){
			all_not_equal=false;
		}
	}



	if(all_not_equal)
		for(int i=0;i<m;i++){
			this_normalized_func_values[i]=(this->evaluation_of_objective_functions[i]-min_values_of_objective_functions[i])/(max_values_of_objective_functions[i]-min_values_of_objective_functions[i]);
			s_normalized_func_values[i]=(s->evaluation_of_objective_functions[i]-min_values_of_objective_functions[i])/(max_values_of_objective_functions[i]-min_values_of_objective_functions[i]);
		}
	else
		for(int i=0;i<m;i++){
			this_normalized_func_values[i]=this->evaluation_of_objective_functions[i];
			s_normalized_func_values[i]=s->evaluation_of_objective_functions[i];
		}

	double distance=0;
	for(int i=0; i<m; i++){
		distance+=(this_normalized_func_values[i]-s_normalized_func_values[i])*(this_normalized_func_values[i]-s_normalized_func_values[i]);
	}
	return sqrt(distance);
}



double Solution::distance_in_objective_space_not_normalized(Solution *s){
	int m=f->get_m();
	double this_normalized_func_values[m];
	double s_normalized_func_values[m];
	
	for(int i=0;i<m;i++){
		this_normalized_func_values[i]=this->evaluation_of_objective_functions[i];
		s_normalized_func_values[i]=s->evaluation_of_objective_functions[i];
	}

	double distance=0;
	for(int i=0; i<m; i++){
		distance+=(this_normalized_func_values[i]-s_normalized_func_values[i])*(this_normalized_func_values[i]-s_normalized_func_values[i]);
	}
	return sqrt(distance);
}




bool Solution::equals(Solution* other) {
	int d=this->f->get_d();
	for (int i = 0; i < d; i++) {
		if (decision_point[i] != other->decision_point[i]) {
			return false;
		}
	}
	return true;
}

void Solution::avg_loc(Solution& other, vector<double>& loc) {
	for (int i = 0; i < get_f()->get_d(); i++) {
		loc.push_back(get_decision_point()[i] +  4*(other.get_decision_point()[i] - get_decision_point()[i])  );
	}
}                    