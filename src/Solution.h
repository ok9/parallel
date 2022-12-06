#ifndef SRC_SOLUTION_H_
#define SRC_SOLUTION_H_

#include <iostream>
#include <random>
#include <iomanip>

#include "Function.h"


using namespace std;


class Hypercube;
class Solution { 
private:
	
  	double *decision_point=nullptr; // decision point in unit hypercube decision space
  	
  	double *evaluation_of_objective_functions=nullptr; // evaluation of objective functions

  	bool foud_by_local_search=false;

  	double selection_function1_minimal_distance_to_solution_decision=0;

  	double selection_function2_minimal_distance_to_Pareto_solution_function_value=0;



  	

  	Function *f;
public:

	double distance_in_decision_space(Solution *s);

	double distance_in_objective_space(Solution *s);
	
	double distance_in_objective_space_not_normalized(Solution *s);
	Solution(Function *f)  {
		this->f=f;
	}

	~Solution() {
		//cout<<"deleting"<<endl;
		if(	decision_point != nullptr )
			delete decision_point;
		if(	evaluation_of_objective_functions != nullptr )
			delete evaluation_of_objective_functions;		
	}
	
	void calculate_selection_functions(vector<Solution*> *U, vector<Solution*> *P);

	Function *get_f(){
		return f;
	}
	bool is_dominated_by(Solution *s);
	bool selection_functions_is_dominated_by(Solution *s);

	void evaluate_function();
	void evaluate_function(double *loc);
	void evaluate_function(vector<double> loc);
	double generate_random_number();
	void generate_decision_point();
	void set_decision_point(vector<double> loc);
	void set_decision_point(double *loc);
	void set_decision_point_and_evaluation_of_objective_functions(double *data);
	
	void generate_decision_point_in_hypercube(Hypercube *h);
	bool equals(Solution* other);

	void avg_loc(Solution& other, vector<double>& loc);

	double *get_evaluation_of_objective_functions(){
		return evaluation_of_objective_functions;
	}

	double *get_decision_point(){
		return decision_point;
	}
	bool  get_foud_by_local_search(){
		return foud_by_local_search;
	}

	void set_foud_by_local_search(bool foud_by_local_search){
		this->foud_by_local_search=foud_by_local_search;
	}

	string to_string(){

		static string name ;
		name = "decision_point = [";
		for(int i=0; decision_point != nullptr && i<f->get_d(); i++){
			name=name+" "+std::to_string(decision_point[i]);
			//cout<<"decision_point[i] : "<<std::setprecision (20)<<decision_point[i]<<endl;
		}
		name=name+"]\n";


		double *upper_bound=f->upper_bound;
        double *lower_bound=f->lower_bound;

        name =name+ "not_normalized_decision_point = [";
		for(int i=0; decision_point != nullptr && i<f->get_d(); i++){
			name=name+" "+std::to_string(lower_bound[i]+decision_point[i]*(upper_bound[i]-lower_bound[i]));
			//cout<<"decision_point[i] : "<<std::setprecision (20)<<decision_point[i]<<endl;
		}
		name=name+"]\n";


		name=name+"evaluation_of_objective_functions = [";
		for(int i=0; evaluation_of_objective_functions != nullptr &&  i<f->get_m(); i++){
			name=name+" "+std::to_string(evaluation_of_objective_functions[i]);
		}
		name=name+"]\n";
		name=name+"selection_function1_minimal_distance_to_solution_decision = "+std::to_string(selection_function1_minimal_distance_to_solution_decision)+"\n";
		name=name+"selection_function2_minimal_distance_to_Pareto_solution_function_value = "+std::to_string(selection_function2_minimal_distance_to_Pareto_solution_function_value)+"\n";
		
		/*decision_point; 
  	
  		evaluation_of_objective_functions;

  		Function *f;
*/



		return name;
	}






};



#endif /* SRC_SOLUTION_H_ */
