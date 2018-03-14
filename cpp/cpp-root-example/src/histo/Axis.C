#include <iostream>
using namespace std;

class Axis {
	public:
		int nBins;
		float min, max;

		Axis() {
		}

		Axis(int nBins, float min, float max) {
			this->nBins = nBins;
			this->min = min;
			this->max = max;
		}

		int GetBin(float value) {
			if(value < min) return -1;
			else if(value >= max) return nBins;
			else return (value - min)/((max - min)/nBins);
		}
};
