
#include "Hypercube.h"
#include <cassert>
#include <iostream>

double *Hypercube::generate_decision_point_inside(){
	double *decision_point;// decision_point shoud be generated only once for solution object
	int d=s->get_f()->get_d();
	decision_point=new double [d];

	for(int i=0;i<d;i++){
		decision_point[i]=boundaries_min_value[i]+ (boundaries_max_value[i]-boundaries_min_value[i])*s->generate_random_number();
	}
	return decision_point;
}

vector<Solution*> *Hypercube::filter_inside_points(vector<Solution*> *P){

	vector<Solution*> *inside_points=new vector<Solution*>;

	for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
        int d=s->get_f()->get_d(); 
        bool all_inside=true;
  		for(int i =0;i<d;i++){
  			if( (*it)->get_decision_point()[i] < boundaries_min_value[i]  ||  (*it)->get_decision_point()[i] > boundaries_max_value[i] ) {
  				all_inside=false;
  				break;
  			}
  			
  		}
  		if(all_inside)
  			inside_points->push_back(*it);

    }

	return inside_points;
}


void Hypercube::update_edge(double edge_size){
	this->edge_size=edge_size;
	double *center = s->get_decision_point();
	int d=s->get_f()->get_d();

	for(int i =0;i<d;i++){
		//cout<<"((center[i]-edge_size*0.5)>0?(center[i]-edge_size*0.5):0)"<<((center[i]-edge_size*0.5)>0?(center[i]-edge_size*0.5):0)<<endl;
		boundaries_min_value[i]=((center[i]-edge_size*0.5)>0?(center[i]-edge_size*0.5):0);
		boundaries_max_value[i]=((center[i]+edge_size*0.5)<1?(center[i]+edge_size*0.5):1);
	}

}