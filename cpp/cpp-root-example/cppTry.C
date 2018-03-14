#include <fstream>
#include <iostream>
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "src/histo/TH1FCollection1D.C"

void cppTry() {


	Axis QQAxis = Axis(5, 0, 5);
	Axis xAxis = Axis(20, 0, 1);
	TH1FCollection1D *hhh = new TH1FCollection1D("hhh", xAxis, QQAxis);

	for(int k = 0; k < 500; k++) {
		float randx = gRandom->Gaus(0.5, 0.1);
		float randQQ = gRandom->Uniform(5.0);
		hhh->histos[QQAxis.GetBin(randQQ)].Fill(randx);
	}

	TCanvas *can = new TCanvas("can", "can", 1000, 500);
	can->Divide(QQAxis.nBins, 1);
	for(int k = 0; k < QQAxis.nBins; k++) {
		can->cd(k+1);
		hhh->histos[k].Draw();
	}


}
