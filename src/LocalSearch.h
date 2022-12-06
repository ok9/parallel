#ifndef SRC_LOCALSEARCH_H_
#define SRC_LOCALSEARCH_H_

#include <algorithm>

#include "Function.h"
#include "Solution.h"



using namespace std;

class LocalSearchParams {
public:
	Function* _f;	// objective function
	Solution* _x0;	// the starting point
	unsigned int _N_max;	// maximum number of function evaluations

	LocalSearchParams() :
			_f(0),  _x0(0), _N_max(1000000) {

	}

	virtual ~LocalSearchParams() {
		delete _f;
		delete _x0;
		//cout<<"DETELE2";
	}

	

};

class LocalSearch {
public:
	LocalSearchParams* _params;
	vector<Solution*>* _x; // unique evaluations
	vector<Solution*>* _xon; // best evaluations

	LocalSearch(LocalSearchParams* params) {
		_params = params;

		_x = new vector<Solution*>;
		_xon = new vector<Solution*>;
	}

	virtual ~LocalSearch() {
		//cout<<"DETELE1";

		for (vector<Solution*>::iterator it = _x->begin(); it != _x->end(); ++it) {
            delete *it;
        }
        delete _xon;
		delete _x;
        //cout<<"DELETE55555"<<endl;
		delete _params;
	}
};

#endif /* SRC_LOCALSEARCH_H_ */
