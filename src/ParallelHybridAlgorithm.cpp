
#include "ParallelHybridAlgorithm.h"
#include <vector>
#include <cassert>
#include "SurogateFunction.h"
#include "SurogateOneCriteriaFunction.h"
#include "LocalSearch.h"
#include "HookeJeeves.h"
#include <mpi.h>
#include <unistd.h>
#include <time.h>


void ParallelHybridAlgorithm::put_to_Pareto_optimal_set(Solution* s){
	if(P->empty()){
		P->push_back(s);
		double *evaluation_of_objective_functions = s->get_evaluation_of_objective_functions();
		int m=s->get_f()->get_m();
		if(s->get_f()->max_values_of_objective_functions == nullptr && s->get_f()->min_values_of_objective_functions == nullptr){
			s->get_f()->max_values_of_objective_functions = new double[m];
			s->get_f()->min_values_of_objective_functions = new double[m];
			for(int i=0;i<m;i++){
				s->get_f()->max_values_of_objective_functions[i] = evaluation_of_objective_functions[i];
				s->get_f()->min_values_of_objective_functions[i] = evaluation_of_objective_functions[i];
			}
		}
	}else{
		bool dominated=false;
		for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
            if( s->is_dominated_by(*it) || s->equals(*it)){
                dominated=true;
                break;
            }
        }
        if(!dominated){
        	double *evaluation_of_objective_functions = s->get_evaluation_of_objective_functions();
        	int m=s->get_f()->get_m();
        	for(int i=0;i<m;i++){
				s->get_f()->max_values_of_objective_functions[i] = evaluation_of_objective_functions[i];
				s->get_f()->min_values_of_objective_functions[i] = evaluation_of_objective_functions[i];
			}
            for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ) {
                if((*it)->is_dominated_by(s)){
                    it=P->erase(it);   
                }else{
                	double *evaluation_of_objective_functions1 = (*it)->get_evaluation_of_objective_functions();
                	for(int i=0;i<m;i++){
						if(s->get_f()->max_values_of_objective_functions[i] < evaluation_of_objective_functions1[i])
							s->get_f()->max_values_of_objective_functions[i] = evaluation_of_objective_functions1[i];
						if(s->get_f()->min_values_of_objective_functions[i] > evaluation_of_objective_functions1[i])
							s->get_f()->min_values_of_objective_functions[i] = evaluation_of_objective_functions1[i];
					}

                    ++it;
                }
            }
            P->push_back(s);
        }
	}
}


void ParallelHybridAlgorithm::put_to_selecton_functions_Pareto_optimal_set(Solution* s, vector<Solution*> *P_sel){
	if(P_sel->empty()){
		P_sel->push_back(s);
	}else{
		bool dominated=false;
		for (vector<Solution*>::iterator it = P_sel->begin(); it != P_sel->end(); ++it) {
            if( s->selection_functions_is_dominated_by(*it) ){
                dominated=true;
                break;
            }
        }
        if(!dominated){
            for (vector<Solution*>::iterator it = P_sel->begin(); it != P_sel->end(); ) {
                if((*it)->selection_functions_is_dominated_by(s)){
                	Solution* s_to_delete;
                	s_to_delete=*it;
                    it=P_sel->erase(it); 
                    delete s_to_delete;  
                }else{
                    ++it;
                }
            }
            P_sel->push_back(s);
        }else{
        	delete s;
        }
	}
	/*for (vector<Solution*>::iterator it = P_sel->begin(); it != P_sel->end(); ++it) {
        for (vector<Solution*>::iterator it1 = P_sel->begin(); it1 != P_sel->end(); ++it1) {
            assert(!(*it)->selection_functions_is_dominated_by(*it1));
        }
    }*/

}










void ParallelHybridAlgorithm::optimize_master(){
	//cout<<"MASTER";
	vector<Solution*>* minimal_soultions = new vector<Solution*>(f->get_m());//f->get_m()
	if(IS_TESTING){
		cout<<"INITIAL RANDOM GENERATION:"<<endl;
	}
	for(int i=0;i<N;i++){
		Solution* s = new Solution(this->f);
		s->generate_decision_point();
		s->evaluate_function();
		if(IS_TESTING){
			double *upper_bound=(s)->get_f()->upper_bound;
            double *lower_bound=(s)->get_f()->lower_bound;
			cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (s)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (s)->get_decision_point()[1]<<endl;
			//cout<<s->get_evaluation_of_objective_functions()[0]<<" "<<s->get_evaluation_of_objective_functions()[1]<<endl;
		}
		U->push_back(s);
		put_to_Pareto_optimal_set(s);
	}
	int I=0;
	int N_q=N*q;


	


	while(U->size()<N_max && I<I_max ){
		I++;

		int local_evaluations=0;
		int global_evaluations=0;
		
		if(p>0){
			vector<Solution*>* old_P = new vector<Solution*>;
	        for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
	            old_P->push_back(*it);
	        }
	        if(IS_TESTING){
				cout<<"GLOBAL SEARCH (local generation): "<<endl;
				//cout<<"old_P->size() : "<<old_P->size()<<endl;
			}
			for (vector<Solution*>::iterator it = old_P->begin(); it != old_P->end(); ++it) {
				if(U->size()>=N_max)
					break;
				double edge=0.2;
				Hypercube* h=new Hypercube(*it, edge);
				vector<Solution*> *inside_hipercube_U=h->filter_inside_points(U);
				vector<Solution*> *inside_hipercube_P=h->filter_inside_points(old_P);
	    		
				while(inside_hipercube_U->size()==1){
					//cout<<"INSIDE"<<endl;
					delete inside_hipercube_U;
					delete inside_hipercube_P;
					edge+=0.2;
					h->update_edge(edge);
					inside_hipercube_U=h->filter_inside_points(U);
					inside_hipercube_P=h->filter_inside_points(old_P);
				}
				
				while(inside_hipercube_U->size()>1  && edge>=pow(2, -((int) hn))    ){
					if(IS_TESTING){
						cout<<" h->to_string1(): "<<h->to_string1()<<endl;
					}
					
					vector<Solution*> *P_sel1=new vector<Solution*>();
					for(int i=0;i<N_q;i++){
						Solution* s = new Solution(this->f);
					    s->generate_decision_point_in_hypercube(h);
					    s->calculate_selection_functions(inside_hipercube_U, inside_hipercube_P);
					    put_to_selecton_functions_Pareto_optimal_set(s, P_sel1);
					}
					
					for (vector<Solution*>::iterator it = P_sel1->begin(); it != P_sel1->end(); ++it) {
			        	(*it)->evaluate_function();
			        	if(IS_TESTING){
							double *upper_bound=(*it)->get_f()->upper_bound;
				            double *lower_bound=(*it)->get_f()->lower_bound;
							cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
							//cout<<s->get_evaluation_of_objective_functions()[0]<<" "<<s->get_evaluation_of_objective_functions()[1]<<endl;
						}
			        	U->push_back(*it);
						put_to_Pareto_optimal_set(*it);
						local_evaluations++;
			        }
					delete P_sel1;

					
					edge/=2;
					h->update_edge(edge);
					vector<Solution*> *inside_hipercube_U_new=h->filter_inside_points(inside_hipercube_U);
					vector<Solution*> *inside_hipercube_P_new=h->filter_inside_points(inside_hipercube_P);
					delete inside_hipercube_U;
					delete inside_hipercube_P;
					inside_hipercube_U=inside_hipercube_U_new;
					inside_hipercube_P=inside_hipercube_P_new;
				}
				delete h;
				delete inside_hipercube_U;
				delete inside_hipercube_P;
			}
			delete old_P;
		}
		if(IS_TESTING1)
			cout<<"After local global search U->size() : "<<U->size()<<endl;
		
		if(IS_TESTING){
			cout<<"GLOBAL SEARCH : "<<endl;
		}	
		while(global_evaluations==local_evaluations || 1-p>1.0*global_evaluations/(local_evaluations+global_evaluations)){
			if(U->size()>=N_max)
					break;
			
			vector<Solution*> *P_sel=new vector<Solution*>();
			for(int i=0;i<N_q;i++){
				Solution* s = new Solution(this->f);
			    s->generate_decision_point();
			    s->calculate_selection_functions(U, P);
			    put_to_selecton_functions_Pareto_optimal_set(s, P_sel);
			}
			
			for (vector<Solution*>::iterator it = P_sel->begin(); it != P_sel->end(); ++it) {
	        	(*it)->evaluate_function();
	        	if(IS_TESTING){
					double *upper_bound=(*it)->get_f()->upper_bound;
		            double *lower_bound=(*it)->get_f()->lower_bound;
					cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
					//cout<<s->get_evaluation_of_objective_functions()[0]<<" "<<s->get_evaluation_of_objective_functions()[1]<<endl;
				}
	        	U->push_back(*it);
				put_to_Pareto_optimal_set(*it);
				global_evaluations++;
	        }
			delete P_sel;
		}

		if(IS_TESTING1)
			cout<<"After  global search U->size() : "<<U->size()<<endl;

		if(IS_TESTING){
			cout<<"HJ LOCAL SEARCH : "<<endl;
		}




		vector<Solution*>* old_P = new vector<Solution*>;
        for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
            old_P->push_back(*it);
        }
		for (vector<Solution*>::iterator it = old_P->begin(); it != old_P->end(); ++it) {
			if(U->size()>=N_max)
					break;
	        Solution* s = *it;
	        if(s->get_foud_by_local_search())
	        	continue;
	        
	        SurogateFunction *sf = new SurogateFunction(s);
	        Solution* s1 = new Solution(sf);
	        s1->evaluate_function(s->get_decision_point());
	        
	        HookeJeevesParams *hj0 = new HookeJeevesParams();
	        hj0->_h0=h0;
	        hj0->_hn=hn;

	        if(I>1 && update_h0_and_hn){
	        	double min_dist=0;
                for (vector<Solution*>::iterator it1 = old_P->begin(); it1 != old_P->end(); ++it1) {
                    double dist=(*it)->distance_in_decision_space(*it1);
                    if(dist!=0 && min_dist==0){
                        min_dist=dist;
                    }
                    if(dist!=0 && min_dist!=0 && min_dist>dist){
                        min_dist=dist;
                    }
                }
                int new_h0 = floor(log(min_dist*2) / log(0.5))+1;
                if(new_h0<2)
                    new_h0=2;
                int new_hn;
                if(new_h0>=hn)
                    new_hn=new_h0+2;
                else
                    new_hn=hn;
                hj0->_h0=new_h0;
	        	hj0->_hn=new_hn;
	        }

	        hj0->_x0=s1;
	        hj0->_f=sf;
	        HookeJeeves *hj = new HookeJeeves(hj0);

	        hj->optimize();

	        vector<Solution*> *all_p=sf->get_all_points();
	        
	        for (vector<Solution*>::iterator it = all_p->begin(); it != all_p->end(); ++it) {
	        	if(IS_TESTING){
					double *upper_bound=(*it)->get_f()->upper_bound;
		            double *lower_bound=(*it)->get_f()->lower_bound;
					cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
					//cout<<s->get_evaluation_of_objective_functions()[0]<<" "<<s->get_evaluation_of_objective_functions()[1]<<endl;
				}
	            U->push_back(*it);
				put_to_Pareto_optimal_set(*it);
	        }

	        delete hj;
	    }

	    delete old_P;

	    if(I==1){
	    	if(IS_TESTING){
				cout<<"HJ ONE CRITERIA FIRST TIME LOCAL SEARCH : "<<endl;
			}
		    for(int i=0;i<f->get_m();i++){
		    	Solution *best=*(P->begin());
		    	for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
		            if(best->get_evaluation_of_objective_functions()[i]>(*it)->get_evaluation_of_objective_functions()[i])
		            	best=*it;
		        }
		        
		        SurogateOneCriteriaFunction *sf = new SurogateOneCriteriaFunction(best, i);
	            Solution* s1 = new Solution(sf);
	            s1->evaluate_function(best->get_decision_point());
	            
	            HookeJeevesParams *hj0 = new HookeJeevesParams();
	            hj0->_h0=h0;
	        	hj0->_hn=hn;
	            hj0->_x0=s1;
	            hj0->_f=sf;
	            HookeJeeves *hj = new HookeJeeves(hj0);
	            
	            hj->optimize();
	            minimal_soultions->at(i)=sf->get_current_min_point();
	            vector<Solution*> *all_p=sf->get_all_points();
	            
		        for (vector<Solution*>::iterator it = all_p->begin(); it != all_p->end(); ++it) {
		        	if(IS_TESTING){
						double *upper_bound=(*it)->get_f()->upper_bound;
			            double *lower_bound=(*it)->get_f()->lower_bound;
						cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
						//cout<<s->get_evaluation_of_objective_functions()[0]<<" "<<s->get_evaluation_of_objective_functions()[1]<<endl;
					}
		            U->push_back(*it);
					put_to_Pareto_optimal_set(*it);
		        }

		        delete hj;
	    	}

	    }	
	    else{
	    	if( IS_TESTING){
				cout<<"HJ ONE CRITERIA NOT FIRST TIME LOCAL SEARCH : "<<endl;
			}
	    	for(int i=0;i<f->get_m();i++){
	    		if(U->size()>=N_max)
					break;
		        
		        if(minimal_soultions->at(i)->get_evaluation_of_objective_functions()[i]  > f->min_values_of_objective_functions[i]){
		        	
		        	Solution *best=*(P->begin());
			    	for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
			            if(best->get_evaluation_of_objective_functions()[i]>(*it)->get_evaluation_of_objective_functions()[i])
			            	best=*it;
			        }
			        
			        SurogateOneCriteriaFunction *sf = new SurogateOneCriteriaFunction(best, i);
		            Solution* s1 = new Solution(sf);
		            s1->evaluate_function(best->get_decision_point());
		            
		            HookeJeevesParams *hj0 = new HookeJeevesParams();
		            hj0->_h0=h0;
	        		hj0->_hn=hn;

	        		if(I>1 && update_h0_and_hn){
	        			
			        	double min_dist=0;
		                for (vector<Solution*>::iterator it1 = P->begin(); it1 != P->end(); ++it1) {
		                    double dist=best->distance_in_decision_space(*it1);
		                    if(dist!=0 && min_dist==0){
		                        min_dist=dist;
		                    }
		                    if(dist!=0 && min_dist!=0 && min_dist>dist){
		                        min_dist=dist;
		                    }
		                }
		                
		                int new_h0 = floor(log(min_dist*2) / log(0.5))+1;
		                if(new_h0<2)
		                    new_h0=2;
		                int new_hn;
		                if(new_h0>=hn)
		                    new_hn=new_h0+2;
		                else
		                    new_hn=hn;
		                hj0->_h0=new_h0;
			        	hj0->_hn=new_hn;
			        	
			        }

		            hj0->_x0=s1;
		            hj0->_f=sf;
		            HookeJeeves *hj = new HookeJeeves(hj0);
		            
		            hj->optimize();
		            minimal_soultions->at(i)=sf->get_current_min_point();
		            vector<Solution*> *all_p=sf->get_all_points();
		            
			        for (vector<Solution*>::iterator it = all_p->begin(); it != all_p->end(); ++it) {
			        	if(IS_TESTING){
							double *upper_bound=(*it)->get_f()->upper_bound;
				            double *lower_bound=(*it)->get_f()->lower_bound;
							cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
							//cout<<s->get_evaluation_of_objective_functions()[0]<<" "<<s->get_evaluation_of_objective_functions()[1]<<endl;
						}
			            U->push_back(*it);
						put_to_Pareto_optimal_set(*it);
			        }

			        delete hj;
		        }
		        
	    	}

	    }

	    if(IS_TESTING1)
	    	cout<<"After local search U->size() : "<<U->size()<<endl;
	}
	delete minimal_soultions;
	
    //cout<<"P->size() : "<<P->size()<<endl;

    double  function_evaluation_count;
    function_evaluation_count=U->size();
    for(int i=1;i<world_size;i++){
        int size;
        MPI_Status status;
        //sleep(1);
        MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &size);
        //cout<<"status.MPI_SOURCE : "<<status.MPI_SOURCE<<" size : "<<size<<endl;
        
        double pareto_data[size];
        
        
        MPI_Recv(pareto_data, size, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        //cout<<"Master received data : "<<endl;
        //cout<<" "<<pareto_data[0]<<", ";
        //cout<<" "<<pareto_data[1]<<" "<<endl;
        function_evaluation_count+=pareto_data[1];
        double solution_data[P->at(0)->get_f()->get_m()+P->at(0)->get_f()->get_d()];
        int k=0;
        for(int j=2;j<size;j++){
        	if((j-2)!=0 && (j-2)%(P->at(0)->get_f()->get_m()+P->at(0)->get_f()->get_d())==0){
        		//cout<<"j-2 : "<<j-2<<endl;
        		Solution* s = new Solution(this->f);
        		s->set_decision_point_and_evaluation_of_objective_functions(solution_data);
        		//cout<<endl<<"recave "<<s->to_string();
        		//cout<<endl;
        		U->push_back(s);
        		if(params->accept_sub_optimal_Pareto_solutions)
        			P->push_back(s);
        		else
					put_to_Pareto_optimal_set(s);
        		k=0;
        	}
        	solution_data[k++]=pareto_data[j];
            //cout<<" "<<pareto_data[j]<<", ";
        }
        Solution* s = new Solution(this->f);
        s->set_decision_point_and_evaluation_of_objective_functions(solution_data);
        U->push_back(s);
		if(params->accept_sub_optimal_Pareto_solutions)
			P->push_back(s);
		else
			put_to_Pareto_optimal_set(s);
        //cout<<endl<<"recave "<<s->to_string();		
        //cout<<"function_evaluation_count"<<function_evaluation_count<<endl;
        //cout<<"---------------------------------"<<endl;
	    
    }
    total_function_evaluation_count=function_evaluation_count;

	return;
}











void ParallelHybridAlgorithm::optimize_slave(){
	//cout<<"SLAVE";
	vector<Solution*>* minimal_soultions = new vector<Solution*>(f->get_m());//f->get_m()
	if(IS_TESTING){
		cout<<"INITIAL RANDOM GENERATION:"<<endl;
	}
	for(int i=0;i<N;i++){
		Solution* s = new Solution(this->f);
		s->generate_decision_point();
		s->evaluate_function();
		if(IS_TESTING){
			double *upper_bound=(s)->get_f()->upper_bound;
            double *lower_bound=(s)->get_f()->lower_bound;
			cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (s)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (s)->get_decision_point()[1]<<endl;
			//cout<<s->get_evaluation_of_objective_functions()[0]<<" "<<s->get_evaluation_of_objective_functions()[1]<<endl;
		}
		U->push_back(s);
		put_to_Pareto_optimal_set(s);
	}
	int I=0;
	int N_q=N*q;
    

	while(U->size()<N_max && I<I_max ){
		I++;

		int local_evaluations=0;
		int global_evaluations=0;
		
		if(p>0){
			vector<Solution*>* old_P = new vector<Solution*>;
	        for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
	            old_P->push_back(*it);
	        }
	        if(IS_TESTING){
				cout<<"GLOBAL SEARCH (local generation): "<<endl;
				//cout<<"old_P->size() : "<<old_P->size()<<endl;
			}
			for (vector<Solution*>::iterator it = old_P->begin(); it != old_P->end(); ++it) {
				if(U->size()>=N_max)
					break;
				double edge=0.2;
				Hypercube* h=new Hypercube(*it, edge);
				vector<Solution*> *inside_hipercube_U=h->filter_inside_points(U);
				vector<Solution*> *inside_hipercube_P=h->filter_inside_points(old_P);
	    		
				while(inside_hipercube_U->size()==1){
					//cout<<"INSIDE"<<endl;
					delete inside_hipercube_U;
					delete inside_hipercube_P;
					edge+=0.2;
					h->update_edge(edge);
					inside_hipercube_U=h->filter_inside_points(U);
					inside_hipercube_P=h->filter_inside_points(old_P);
				}
				
				while(inside_hipercube_U->size()>1  && edge>=pow(2, -((int) hn))    ){
					if(IS_TESTING){
						cout<<" h->to_string1(): "<<h->to_string1()<<endl;
					}
					
					vector<Solution*> *P_sel1=new vector<Solution*>();
					for(int i=0;i<N_q;i++){
						Solution* s = new Solution(this->f);
					    s->generate_decision_point_in_hypercube(h);
					    s->calculate_selection_functions(inside_hipercube_U, inside_hipercube_P);
					    put_to_selecton_functions_Pareto_optimal_set(s, P_sel1);
					}
					
					for (vector<Solution*>::iterator it = P_sel1->begin(); it != P_sel1->end(); ++it) {
			        	(*it)->evaluate_function();
			        	if(IS_TESTING){
							double *upper_bound=(*it)->get_f()->upper_bound;
				            double *lower_bound=(*it)->get_f()->lower_bound;
							cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
							//cout<<s->get_evaluation_of_objective_functions()[0]<<" "<<s->get_evaluation_of_objective_functions()[1]<<endl;
						}
			        	U->push_back(*it);
						put_to_Pareto_optimal_set(*it);
						local_evaluations++;
			        }
					delete P_sel1;

					
					edge/=2;
					h->update_edge(edge);
					vector<Solution*> *inside_hipercube_U_new=h->filter_inside_points(inside_hipercube_U);
					vector<Solution*> *inside_hipercube_P_new=h->filter_inside_points(inside_hipercube_P);
					delete inside_hipercube_U;
					delete inside_hipercube_P;
					inside_hipercube_U=inside_hipercube_U_new;
					inside_hipercube_P=inside_hipercube_P_new;
				}
				delete h;
				delete inside_hipercube_U;
				delete inside_hipercube_P;
			}
			delete old_P;
		}
		if(IS_TESTING1)
			cout<<"After local global search U->size() : "<<U->size()<<endl;
		
		if(IS_TESTING){
			cout<<"GLOBAL SEARCH : "<<endl;
		}	
		while(global_evaluations==local_evaluations || 1-p>1.0*global_evaluations/(local_evaluations+global_evaluations)){
			if(U->size()>=N_max)
					break;
			
			vector<Solution*> *P_sel=new vector<Solution*>();
			for(int i=0;i<N_q;i++){
				Solution* s = new Solution(this->f);
			    s->generate_decision_point();
			    s->calculate_selection_functions(U, P);
			    put_to_selecton_functions_Pareto_optimal_set(s, P_sel);
			}
			
			for (vector<Solution*>::iterator it = P_sel->begin(); it != P_sel->end(); ++it) {
	        	(*it)->evaluate_function();
	        	if(IS_TESTING){
					double *upper_bound=(*it)->get_f()->upper_bound;
		            double *lower_bound=(*it)->get_f()->lower_bound;
					cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
					//cout<<s->get_evaluation_of_objective_functions()[0]<<" "<<s->get_evaluation_of_objective_functions()[1]<<endl;
				}
	        	U->push_back(*it);
				put_to_Pareto_optimal_set(*it);
				global_evaluations++;
	        }
			delete P_sel;
		}

		if(IS_TESTING1)
			cout<<"After  global search U->size() : "<<U->size()<<endl;

		if(IS_TESTING){
			cout<<"HJ LOCAL SEARCH : "<<endl;
		}




		vector<Solution*>* old_P = new vector<Solution*>;
        for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
            old_P->push_back(*it);
        }
		for (vector<Solution*>::iterator it = old_P->begin(); it != old_P->end(); ++it) {
			if(U->size()>=N_max)
					break;
	        Solution* s = *it;
	        if(s->get_foud_by_local_search())
	        	continue;
	        
	        SurogateFunction *sf = new SurogateFunction(s);
	        Solution* s1 = new Solution(sf);
	        s1->evaluate_function(s->get_decision_point());
	        
	        HookeJeevesParams *hj0 = new HookeJeevesParams();
	        hj0->_h0=h0;
	        hj0->_hn=hn;

	        if(I>1 && update_h0_and_hn){
	        	double min_dist=0;
                for (vector<Solution*>::iterator it1 = old_P->begin(); it1 != old_P->end(); ++it1) {
                    double dist=(*it)->distance_in_decision_space(*it1);
                    if(dist!=0 && min_dist==0){
                        min_dist=dist;
                    }
                    if(dist!=0 && min_dist!=0 && min_dist>dist){
                        min_dist=dist;
                    }
                }
                int new_h0 = floor(log(min_dist*2) / log(0.5))+1;
                if(new_h0<2)
                    new_h0=2;
                int new_hn;
                if(new_h0>=hn)
                    new_hn=new_h0+2;
                else
                    new_hn=hn;
                hj0->_h0=new_h0;
	        	hj0->_hn=new_hn;
	        }

	        hj0->_x0=s1;
	        hj0->_f=sf;
	        HookeJeeves *hj = new HookeJeeves(hj0);

	        hj->optimize();

	        vector<Solution*> *all_p=sf->get_all_points();
	        
	        for (vector<Solution*>::iterator it = all_p->begin(); it != all_p->end(); ++it) {
	        	if(IS_TESTING){
					double *upper_bound=(*it)->get_f()->upper_bound;
		            double *lower_bound=(*it)->get_f()->lower_bound;
					cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
					//cout<<s->get_evaluation_of_objective_functions()[0]<<" "<<s->get_evaluation_of_objective_functions()[1]<<endl;
				}
	            U->push_back(*it);
				put_to_Pareto_optimal_set(*it);
	        }

	        delete hj;
	    }

	    delete old_P;

	    if(I==1){
	    	if(IS_TESTING){
				cout<<"HJ ONE CRITERIA FIRST TIME LOCAL SEARCH : "<<endl;
			}
		    for(int i=0;i<f->get_m();i++){
		    	Solution *best=*(P->begin());
		    	for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
		            if(best->get_evaluation_of_objective_functions()[i]>(*it)->get_evaluation_of_objective_functions()[i])
		            	best=*it;
		        }
		        
		        SurogateOneCriteriaFunction *sf = new SurogateOneCriteriaFunction(best, i);
	            Solution* s1 = new Solution(sf);
	            s1->evaluate_function(best->get_decision_point());
	            
	            HookeJeevesParams *hj0 = new HookeJeevesParams();
	            hj0->_h0=h0;
	        	hj0->_hn=hn;
	            hj0->_x0=s1;
	            hj0->_f=sf;
	            HookeJeeves *hj = new HookeJeeves(hj0);
	            
	            hj->optimize();
	            minimal_soultions->at(i)=sf->get_current_min_point();
	            vector<Solution*> *all_p=sf->get_all_points();
	            
		        for (vector<Solution*>::iterator it = all_p->begin(); it != all_p->end(); ++it) {
		        	if(IS_TESTING){
						double *upper_bound=(*it)->get_f()->upper_bound;
			            double *lower_bound=(*it)->get_f()->lower_bound;
						cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
						//cout<<s->get_evaluation_of_objective_functions()[0]<<" "<<s->get_evaluation_of_objective_functions()[1]<<endl;
					}
		            U->push_back(*it);
					put_to_Pareto_optimal_set(*it);
		        }

		        delete hj;
	    	}

	    }	
	    else{
	    	if( IS_TESTING){
				cout<<"HJ ONE CRITERIA NOT FIRST TIME LOCAL SEARCH : "<<endl;
			}
	    	for(int i=0;i<f->get_m();i++){
	    		if(U->size()>=N_max)
					break;
		        
		        if(minimal_soultions->at(i)->get_evaluation_of_objective_functions()[i]  > f->min_values_of_objective_functions[i]){
		        	
		        	Solution *best=*(P->begin());
			    	for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
			            if(best->get_evaluation_of_objective_functions()[i]>(*it)->get_evaluation_of_objective_functions()[i])
			            	best=*it;
			        }
			        
			        SurogateOneCriteriaFunction *sf = new SurogateOneCriteriaFunction(best, i);
		            Solution* s1 = new Solution(sf);
		            s1->evaluate_function(best->get_decision_point());
		            
		            HookeJeevesParams *hj0 = new HookeJeevesParams();
		            hj0->_h0=h0;
	        		hj0->_hn=hn;

	        		if(I>1 && update_h0_and_hn){
	        			
			        	double min_dist=0;
		                for (vector<Solution*>::iterator it1 = P->begin(); it1 != P->end(); ++it1) {
		                    double dist=best->distance_in_decision_space(*it1);
		                    if(dist!=0 && min_dist==0){
		                        min_dist=dist;
		                    }
		                    if(dist!=0 && min_dist!=0 && min_dist>dist){
		                        min_dist=dist;
		                    }
		                }
		                
		                int new_h0 = floor(log(min_dist*2) / log(0.5))+1;
		                if(new_h0<2)
		                    new_h0=2;
		                int new_hn;
		                if(new_h0>=hn)
		                    new_hn=new_h0+2;
		                else
		                    new_hn=hn;
		                hj0->_h0=new_h0;
			        	hj0->_hn=new_hn;
			        	
			        }

		            hj0->_x0=s1;
		            hj0->_f=sf;
		            HookeJeeves *hj = new HookeJeeves(hj0);
		            
		            hj->optimize();
		            minimal_soultions->at(i)=sf->get_current_min_point();
		            vector<Solution*> *all_p=sf->get_all_points();
		            
			        for (vector<Solution*>::iterator it = all_p->begin(); it != all_p->end(); ++it) {
			        	if(IS_TESTING){
							double *upper_bound=(*it)->get_f()->upper_bound;
				            double *lower_bound=(*it)->get_f()->lower_bound;
							cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
							//cout<<s->get_evaluation_of_objective_functions()[0]<<" "<<s->get_evaluation_of_objective_functions()[1]<<endl;
						}
			            U->push_back(*it);
						put_to_Pareto_optimal_set(*it);
			        }

			        delete hj;
		        }
		        
	    	}

	    }
	    //cout<<"number : "<<number<<endl;
	    if(IS_TESTING1){
	    	cout<<"After local search U->size() : "<<U->size()<<endl;
	    }
	}
	delete minimal_soultions;
	
	int total_size=2+P->size()*(P->at(0)->get_f()->get_m()+P->at(0)->get_f()->get_d());
	
	double data_array[total_size];
	data_array[0]=world_rank;
	data_array[1]=U->size();
	int j=2;
	for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
		//cout<<endl<<"send "<<(*it)->to_string();
		for(int i=0;i<(*it)->get_f()->get_d();j++,i++){
	        data_array[j]=(*it)->get_decision_point()[i];
	    }
	    for(int i=0;i<(*it)->get_f()->get_m();j++,i++){
	        data_array[j]=(*it)->get_evaluation_of_objective_functions()[i];
	    }
	}
	
	/*if(data_array[0]==5){
        cout<<"Slave data : "<<endl;
        cout<<" "<<data_array[0]<<", ";
        cout<<" "<<data_array[1]<<" ";
        for(int j1=2;j1<j;j1++){
        	if((j1-2)%(P->at(0)->get_f()->get_m()+P->at(0)->get_f()->get_d())==0)
        		cout<<endl;
            cout<<" "<<data_array[j1]<<", ";
        }
        cout<<endl;
        cout<<"---------------------------------"<<endl;
    }*/


    MPI_Send(data_array, total_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    
	return;
}










void ParallelHybridAlgorithm::optimize() {
	
     // Start measuring time
    clock_t start = clock();

	MPI_Barrier(MPI_COMM_WORLD);
	// Get the number of processes
    
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the rank of the process
    
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    //printf("ParallelHybridAlgorithm::optimize Processor rank %d out of %d processors\n", world_rank, world_size);

    /*int delay=rand() % 60;
    cout<<"rand() % 60 :"<<delay<<endl;
    sleep(delay);*/
    
    if (world_rank == 0){ 

    	optimize_master();


    }else{
    	optimize_slave();
    }
	
	MPI_Barrier(MPI_COMM_WORLD);

	clock_t end = clock();
    run_time = double(end - start)/(CLOCKS_PER_SEC/1000);
	
	return;
}






