#include "TH1.h"
#include <iostream>
#include <vector>
#include "Axis.C"
using namespace std;

class TH1FCollection1D {
	public:
		string name;
		Axis x1, x2;
		vector<TH1F> histos;

		TH1FCollection1D(string name, Axis x1, Axis x2) {
			this->name = name;
			this->x1 = x1;
			this->x2 = x2;
			for(int k = 0; k < x2.nBins; k++) {
				histos.push_back(TH1F(Form("%s_%i", name.c_str(), k), Form("%s_%i", name.c_str(), k), x1.nBins, x1.min, x1.max));
			}
		}

};
