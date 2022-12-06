#include <cassert>
#include <iostream>
#include <fstream>

#include "HookeJeeves.h"
#include <iomanip>

Solution* HookeJeeves::eval(vector<double> &loc) {
	for (unsigned int i = 0; i < loc.size(); i++) {

		double li = 0 ;
		double ui = 1 ;
		
		if(li > loc.at(i)){
			loc.at(i)=li;
		}
		if(loc.at(i) > ui){
			loc.at(i)=ui;
		}

	}
	
	
	Solution* p = new Solution(_params->_f);
	p->evaluate_function(loc);
	
	//cout<<p->to_string();//to_string
	_x->push_back(p);//
	if (_xon->size() == 0 || _xon->back()->get_evaluation_of_objective_functions()[0] > _x->back()->get_evaluation_of_objective_functions()[0]) {
		_xon->push_back(_x->back());
	} else {
		_xon->push_back(_xon->back());
	}
	return p;
	//
}

void HookeJeeves::optimize() {

	Solution* x_base = _params->_x0;

	unsigned int i = 0;

	bool stat = true;
	while (i < _h.size() && _x->size() < _params->_N_max && stat) {
		double h = _h.at(i);
		//cout<<"h = "<<h<<endl;
		stat = search(x_base, h);
		i++;
	}
	
}

bool HookeJeeves::search(Solution* &x, double h) {
	Solution* x_base = x;
	Solution* x_center = x;
	bool stat = explore(x, x_center, h);
	//cout<<"OUT: x->to_string() : "<<x->to_string();
	while (stat && _x->size() < _params->_N_max) {
		//cout<<" x->to_string() : "<<x->to_string();
		vector<double> center_loc;
		x_base->avg_loc(*x, center_loc);
		//x_center = new Solution(center_loc);
		x_center = new Solution(_params->_f);
		x_center->set_decision_point(center_loc);
		x_base = x;
		//cout<<"WHILE"<<endl;
		//cout<<"IN WHILE1: x_center->to_string() : "<<x_center->to_string();
		stat = explore(x, x_center, h);
		//cout<<"IN WHILE2: x->to_string() : "<<x->to_string()<<" stat: "<<stat<<endl;
		delete x_center;

		if (!stat) {
			x = x_base;
			x_center = x_base;
			//cout<<"IF"<<endl;
			stat = explore(x, x_center, h);
			//cout<<"IN IF: x->to_string() : "<<x->to_string()<<" stat: "<<stat<<endl;;
		}
		if (stat) {

			if (false && fabs(x_base->get_evaluation_of_objective_functions()[0] - x->get_evaluation_of_objective_functions()[0])
					< ((HookeJeevesParams*) _params)->_tol) {
				cout<<"_params->_f->to_string()"<<_params->_f->to_string()<<endl;
				cout<<"Cant happen!!!"<<endl;
				cout<<"x_base : "<<x_base->to_string()<<endl;
				cout<<"x : "<<x->to_string()<<endl;
				exit(0);
				//cout << fabs(x_base->get_val() - x->get_val()) << " < " << _params->_tol << endl;
				return false;
			}
		}
	}
	return true;
}

bool HookeJeeves::explore(Solution* &x_base, Solution* x_center, double h) {
	// x_base - base (reference) point, an improvement over which is sought
	// x_center - center of the stencil
	// h - scale
	


	unsigned int dim = _params->_f->get_d();

	Solution* x_best = x_base;
	Solution* x_t = x_center;

	vector<double> p_loc(dim, 0);
	
	unsigned int i = 0;
	while (i < dim && _x->size() < _params->_N_max) {
		//cout<<"WORK"<<endl;
		unsigned int pos = dim - 1 - i;
		//p_loc.assign(x_t->get_loc().begin(), x_t->get_loc().end());
		
		for(unsigned int i1=0;i1<dim;i1++){
			p_loc.at(i1)=x_t->get_decision_point()[i1];
		}
		p_loc.at(pos) += h;

		Solution* p = eval(p_loc);
		assert(p);
		//if(p)
		//	cout<<"p->get_val().data[0] : "<<p->get_val().data[0]<<"  x_best->get_val().data[0] : "<< x_best->get_val().data[0]<<endl;
		if (!p || (p->get_evaluation_of_objective_functions()[0] >= x_best->get_evaluation_of_objective_functions()[0])) {
			//cout<<"WORSE OR NULL"<<endl;
			for(unsigned int i1=0;i1<dim;i1++){
				p_loc.at(i1)=x_t->get_decision_point()[i1];
			}
			p_loc.at(pos) -=  h;
			p = eval(p_loc);
		}
		assert(p);
		//if(p)
			//cout<<"p->get_evaluation_of_objective_functions()[0]: "<<std::setprecision (20)<<p->get_evaluation_of_objective_functions()[0]<<"  x_best->get_evaluation_of_objective_functions()[0] : "<< x_best->get_evaluation_of_objective_functions()[0]<<endl;		
		if (p && (p->get_evaluation_of_objective_functions()[0] < x_best->get_evaluation_of_objective_functions()[0])) {
			//cout<<"BETTER"<<endl;
			x_t = p;
			x_best = p;
		}else{
			//cout<<"WORSE OR NULL"<<endl;
		}
		
		i++;
	}

	if (!x_base->equals(x_best)) {
		x_base = x_best;
		return true;
	}
	
	return false;
}

string HookeJeeves::result() {
	

	return "OK";

}


