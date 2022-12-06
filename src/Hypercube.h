#ifndef SRC_HYPERCUBE_H_
#define SRC_HYPERCUBE_H_

#include <iostream>

#include "Solution.h"

using namespace std;



class Hypercube { 
private:

	double *boundaries_min_value=nullptr; // 
  	
  	double *boundaries_max_value=nullptr; // 
	
public:

  	Solution *s=nullptr;

  	double edge_size;

	Hypercube(Solution *s, double  edge_size)  {
		this->s=s;
		this->edge_size=edge_size;
		double *center = s->get_decision_point();
		int d=s->get_f()->get_d();

		boundaries_min_value = new double[d]; 
  		boundaries_max_value = new double[d];  

  		for(int i =0;i<d;i++){
  			//cout<<"((center[i]-edge_size*0.5)>0?(center[i]-edge_size*0.5):0)"<<((center[i]-edge_size*0.5)>0?(center[i]-edge_size*0.5):0)<<endl;
  			boundaries_min_value[i]=((center[i]-edge_size*0.5)>0?(center[i]-edge_size*0.5):0);
  			boundaries_max_value[i]=((center[i]+edge_size*0.5)<1?(center[i]+edge_size*0.5):1);
  		}

		//cout<<"Hypercube CREATE"<<endl;
	}

	~Hypercube() {
		//cout<<"deleting hypercube"<<endl;
		delete boundaries_min_value;
		delete boundaries_max_value;		
	}

	vector<Solution*> *filter_inside_points(vector<Solution*> *P);

	double *generate_decision_point_inside();

	void update_edge(double edge_size);

	string to_string(){
		static string name ;
		int d=s->get_f()->get_d();
		name = "boundaries_min_value = [";
		for(int i=0;  i<d; i++){
			name=name+" "+std::to_string(boundaries_min_value[i]);
		}
		name=name+"]\n";
		name=name+"boundaries_max_value = [";
		for(int i=0;  i<d; i++){
			name=name+" "+std::to_string(boundaries_max_value[i]);
		}
		name=name+"]\n";
		return name;
	}

	string to_string1(){

		double *upper_bound=s->get_f()->upper_bound;
		double *lower_bound=s->get_f()->lower_bound;
		//cout<<<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
		
		static string name ;
		//int d=s->get_f()->get_d();
		name = "hipercube  \n";
		
		name=name+" "+std::to_string(lower_bound[0]+(upper_bound[0]-lower_bound[0])* boundaries_min_value[0]);
		name=name+" "+std::to_string(lower_bound[1]+(upper_bound[1]-lower_bound[1])* boundaries_min_value[1]);
		
		name=name+"\n";
		name=name+" ";

		name=name+" "+std::to_string(lower_bound[0]+(upper_bound[0]-lower_bound[0])* boundaries_max_value[0]);
		name=name+" "+std::to_string(lower_bound[1]+(upper_bound[1]-lower_bound[1])* boundaries_max_value[1]);
		


		name=name+" \n";

		return name;
	}
	

};



#endif /* SRC_SOLUTION_H_ */
