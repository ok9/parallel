#ifndef SRC_HOOKEJEEVES_H_
#define SRC_HOOKEJEEVES_H_

#include "Function.h"

#include "LocalSearch.h"
 

using namespace std;

class HookeJeevesParams: public LocalSearchParams {
public:
	unsigned int _h0; // largest scale exponent
	unsigned int _hn; // smallest scale exponent
	double _tol; // the difference in function values to end the search
	
	HookeJeevesParams() :
			LocalSearchParams() {
		_h0 = 6;
		_hn = 8;
		_tol = 1e-6;
	}
	virtual ~HookeJeevesParams() {
	}
};

class HookeJeeves: public LocalSearch {
public:
	vector<double> _h; // scales

	HookeJeeves(HookeJeevesParams* params) :
			LocalSearch(params) {
		assert(params->_h0<=params->_hn);
		//assert(params->_h0==6 &&  params->_hn ==8);
		for (unsigned int i = params->_h0; i <= params->_hn; i++) {
			_h.push_back(0.796875*pow(2, -((int) i)));// 0.8 dont have precise reprezentation in floating point, so number with precise representation is used: https://www.binaryconvert.com/result_double.html
														//to hold this assertion:
													   	//double x=0.796875;
    													//assert(x==4*x-x-x-x);
			//cout<<"pow(2, -((int) i)) : "<<pow(2, -((int) i))<<endl;
		}

		assert(params->_x0!=NULL);
		//vector<double> start = params->_x0->get_loc();
        //params->_lower=&(params->_f->_l);
        //params->_upper=&(params->_f->_u);
		

	}
	virtual ~HookeJeeves() {
		
	}
	Solution* eval(vector<double> &loc);
	bool explore(Solution* &x_base, Solution* xPcenter, double h);
	bool search(Solution* &x_base, double h);
	void optimize();
	string result();
};

#endif /* SRC_HOOKEJEEVES_H_ */
