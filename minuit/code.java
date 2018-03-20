package myECcalibration;

import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;

import org.jlab.clas.physics.LorentzVector;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.utils.groups.IndexedList;
import org.jlab.io.evio.EvioDataBank;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.base.GeometryFactory;
import org.jlab.detector.geant4.v2.ECGeant4Factory;
import org.jlab.detector.geant4.v2.PCALGeant4Factory;
import org.jlab.geom.base.ConstantProvider;
import org.freehep.math.minuit.FCNBase; 
import org.freehep.math.minuit.FunctionMinimum; 
import org.freehep.math.minuit.MnMigrad; 
import org.freehep.math.minuit.MnPlot; 
import org.freehep.math.minuit.MnScan; 
import org.freehep.math.minuit.MnStrategy; 
import org.freehep.math.minuit.MnUserParameters; 
import org.freehep.math.minuit.Point; 

public class NewSingleStrip_a0a1a2 {

	static HipoDataSource reader = new HipoDataSource();
	static HipoDataSource pionreader = new HipoDataSource();
	
	static H1F hstrip = new H1F("hstrip", "hstrip", 65, 0, 65);
	static H2F ha0f = new H2F("ha0f", "ha0f", 64, 0, 64, 200, -30, 50);
	static H2F ha1f = new H2F("ha1f", "ha1f", 64, 0, 64, 200, 0, 0.05);
	static H2F ha2f = new H2F("ha2f", "ha2f", 64, 0, 64, 200, 200, 350);
	static H2F hchisqf = new H2F("hchisqf", "hchisqf", 64, 0, 64, 200, 0, 2);

	static IndexedList<H1F> histGroups = new IndexedList<H1F>(2);
	
	public static void histos() {
		for(int ihist = 0; ihist < 101; ihist++) {
			for(int istrip = 1; istrip <= 68; istrip++) {
				H1F hdeltaT = new H1F("hdeltaT", "hdeltaT", 200, -10, 10);
				hdeltaT.setTitle("Iteration " + ihist + " #DeltaT");
				hdeltaT.setTitleX("#DeltaT [ns]");
				hdeltaT.setTitleY("Counts");
				histGroups.add(hdeltaT, ihist, istrip);
			}
		}
	}

	static void processEvent(DataEvent event, int eventCounter, short setdetsector, short detECallayer, short U, short V, short W, short view, byte ECallayer, byte stripID, float slope, float intercept, float PMTx, float PMTy,
			List<Float> TDCi, List<Float> ADCi, List<Float> l2i, List<Float> l3i, List<Float> Ti, List<Float> Texpectedi, List<Float> RFTi) {
		
		int eventFinder = eventCounter - 1;
		
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter") && event.hasBank("REC::Scintillator")  && event.hasBank("REC::Event") && event.hasBank("ECAL::clusters") && event.hasBank("ECAL::adc") && event.hasBank("ECAL::tdc") && event.hasBank("RUN::rf"))
		{
			HipoDataBank part = (HipoDataBank) event.getBank("REC::Particle");
			float p = 0;
			int partcount = 0;
			int eleccount = 0;
			for(int parti = 0; parti < part.rows(); parti++)
			{
				int partid = part.getInt("pid", parti);
				if(partid == 11)
				{
					eleccount++;
					HipoDataBank detcount = (HipoDataBank) event.getBank("REC::Calorimeter");
					for(int detcounti = 0; detcounti < detcount.rows(); detcounti++)
					{
						short pindexcount = detcount.getShort("pindex", detcounti);
						byte detectorcount = detcount.getByte("detector", detcounti);
						byte sectorcount = detcount.getByte("sector", detcounti);
						if(pindexcount != parti && sectorcount == setdetsector && detectorcount == 7) partcount++;
					}
				}
			}
			HipoDataBank recev = (HipoDataBank) event.getBank("REC::Event");
			float vtime = -100;
			for(int t = 0; t < recev.rows(); t++)
			{
				float sttime = recev.getFloat("STTime", t);
				if(sttime > -100) vtime = sttime;
			}
			if(vtime > -100 && partcount == 0 && eleccount > 0)
			{
				if(partcount > 1) System.out.println("partcount: " + partcount);
				for(int i = 0; i < part.rows(); i++)
				{
					int pid = part.getInt("pid", i);
					float px = part.getFloat("px", i);
					float py = part.getFloat("py", i);
					float pz = part.getFloat("pz", i);
					p = (float) Math.sqrt(px*px+py*py+pz*pz);
					if(pid == 11)
					{
						HipoDataBank det = (HipoDataBank) event.getBank("REC::Calorimeter");
						HipoDataBank detf = (HipoDataBank) event.getBank("REC::Scintillator");
						float totEdep = 0;
						float EdepECalcluster = 0;
						float Tvertex = vtime;
						float hitx = 0;
						float hity = 0;
						float hitpath = 0;
						float hittime = 0;
						for(int j = 0; j < det.rows(); j++)
						{
							short pindex = det.getShort("pindex", j);
							byte detector = det.getByte("detector", j);
							byte sector = det.getByte("sector", j);
							byte layer = det.getByte("layer", j);
							float x = det.getFloat("x", j);
							float y = det.getFloat("y", j);
							float path = det.getFloat("path", j);
							float time = det.getFloat("time", j);
							float energy = det.getFloat("energy", j);
							if(pindex == i && detector == 7 && sector == setdetsector) totEdep = totEdep + energy;
							if(pindex == i && detector == 7 && sector == setdetsector && layer == detECallayer)
							{
								EdepECalcluster = energy;
								hitx = x;
								hity = y;
								hitpath = path;
							}
						}
						
						float xeff = (hitx+(slope*hity)-(slope*intercept))/(1+(slope*slope));
						float yeff = (slope*xeff)+intercept;
						float leff = (float) Math.sqrt(((xeff-PMTx)*(xeff-PMTx))+((yeff-PMTy)*(yeff-PMTy)));
							
						byte id = -1;
						HipoDataBank clust = (HipoDataBank) event.getBank("ECAL::clusters");
						for(int k = 0; k < clust.rows(); k++)
						{
							byte sector = clust.getByte("sector", k);
							byte layer = clust.getByte("layer", k);
							float energy = clust.getFloat("energy", k);
							if(view == U) id = clust.getByte("idU", k);
							if(view == V) id = clust.getByte("idV", k);
							if(view == W) id = clust.getByte("idW", k);
							HipoDataBank hit = (HipoDataBank) event.getBank("ECAL::hits");
							byte ECalstrip = 0;
							if(sector == setdetsector && layer == detECallayer && energy == EdepECalcluster)
							{
								for(int l = 0; l < hit.rows(); l++)
								{
									byte strip = hit.getByte("strip", l);
									byte peakid = hit.getByte("peakid", l);
									byte EChitlayer = hit.getByte("layer", l);
									if(peakid == id && EChitlayer == ECallayer && strip == stripID) ECalstrip = strip;
								}
							}
							HipoDataBank adc = (HipoDataBank) event.getBank("ECAL::adc");
							int ECalstripADC = 0;
							int ECalstripref1ADC = 0;
							int ECalstripref2ADC = 0;
							if(ECalstrip == stripID)
							{
								for(int m = 0; m < adc.rows(); m++)
								{
									byte adcsector = adc.getByte("sector", m);
									byte adclayer = adc.getByte("layer", m);
									short adcstrip = adc.getShort("component", m);
									int ADC = adc.getInt("ADC", m);
									if(adcsector == setdetsector && adclayer == ECallayer && adcstrip == stripID) ECalstripADC = ADC;
									if(adcsector == setdetsector && adclayer == ECallayer && adcstrip == stripID-1) ECalstripref1ADC = ADC;
									if(adcsector == setdetsector && adclayer == ECallayer && adcstrip == stripID+1) ECalstripref2ADC = ADC;
								}
								if(ECalstripADC > ECalstripref1ADC && ECalstripADC > ECalstripref2ADC)
								{
									HipoDataBank tdc = (HipoDataBank) event.getBank("ECAL::tdc");
									int TDCcounter = 0;
									int TDCvalue = 0;
									for(int n = 0; n < tdc.rows(); n++)
									{
										byte tdcsector = tdc.getByte("sector", n);
										byte tdclayer = tdc.getByte("layer", n);
										short tdcstrip = tdc.getShort("component", n);	
										int TDC = tdc.getInt("TDC", n);
										if(tdcsector == setdetsector && tdclayer == ECallayer && tdcstrip == stripID && TDC > 0)
										{
											TDCcounter++;
											TDCvalue = TDC;
										}
									}
									if(TDCcounter == 1 && (hitpath/30)+(leff/18.1)+Tvertex > 500)
									{
										TDCi.add((float) TDCvalue);
										ADCi.add((float) ECalstripADC);
										l2i.add(leff*leff);
										l3i.add(leff*leff*leff);
										Ti.add((float) ((hitpath/30)+(leff/18.1)+Tvertex));
										Texpectedi.add((float) ((hitpath/30)+(leff/18.1)));
										HipoDataBank runrf = (HipoDataBank) event.getBank("RUN::rf");
										for(int o = 0; o < runrf.rows(); o++)
										{
											short rfid = runrf.getShort("id", o);
											float rftime = runrf.getFloat("time", o);
											if(rfid == 1) RFTi.add(rftime);		
										}
										hstrip.fill(stripID);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	static void processPionEvent(DataEvent event, int pioneventCounter, short setdetsector, short detECallayer, short U, short V, short W, short view, byte ECallayer, byte stripID, float slope, float intercept, float PMTx, float PMTy,
			List<Float> TDCi, List<Float> ADCi, List<Float> l2i, List<Float> l3i, List<Float> Ti, List<Float> Texpectedi, List<Float> RFTi) {
		
		int eventFinder = pioneventCounter - 1;
		
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter") && event.hasBank("REC::Scintillator") && event.hasBank("REC::Event") && event.hasBank("ECAL::clusters") && event.hasBank("ECAL::adc") && event.hasBank("ECAL::tdc") && event.hasBank("RUN::rf"))
		{
			HipoDataBank part = (HipoDataBank) event.getBank("REC::Particle");
			float p = 0;
			int partcount = 0;
			int pioncount = 0;
			for(int parti = 0; parti < part.rows(); parti++)
			{
				byte partcharge = part.getByte("charge", parti);
				int partid = part.getInt("pid", parti);
				if(partcharge == -1 && partid != 11)
				{
					HipoDataBank detcount = (HipoDataBank) event.getBank("REC::Calorimeter");
					for(int detcounti = 0; detcounti < detcount.rows(); detcounti++)
					{
						short pindexcount = detcount.getShort("pindex", detcounti);
						byte detectorcount = detcount.getByte("detector", detcounti);
						byte sectorcount = detcount.getByte("sector", detcounti);
						if(pindexcount == parti && sectorcount == setdetsector) pioncount++;
						if(pindexcount != parti && sectorcount == setdetsector && detectorcount == 7) partcount++;
					}
					
				}
			}
			HipoDataBank recev = (HipoDataBank) event.getBank("REC::Event");
			if(pioncount > 0 && partcount == 0)
			{
				for(int i = 0; i < part.rows(); i++)
				{
					byte charge = part.getByte("charge", i);
					float px = part.getFloat("px", i);
					float py = part.getFloat("py", i);
					float pz = part.getFloat("pz", i);
					p = (float) Math.sqrt(px*px+py*py+pz*pz);
					float beta = (float) (p/Math.sqrt((p*p)+(0.13957*0.13957)));
					if(charge == -1 && p > 1)
					{
						HipoDataBank det = (HipoDataBank) event.getBank("REC::Calorimeter");
						HipoDataBank detf = (HipoDataBank) event.getBank("REC::Scintillator");
						float totEdep = 0;
						float EdepECalcluster = 0;
						float Tvertex = 0;
						float hitx = 0;
						float hity = 0;
						float hitpath = 0;
						float hittime = 0;
						for(int j = 0; j < det.rows(); j++)
						{
							short pindex = det.getShort("pindex", j);
							byte detector = det.getByte("detector", j);
							byte sector = det.getByte("sector", j);
							byte layer = det.getByte("layer", j);
							float x = det.getFloat("x", j);
							float y = det.getFloat("y", j);
							float path = det.getFloat("path", j);
							float time = det.getFloat("time", j);
							float energy = det.getFloat("energy", j);
							if(pindex == i && detector == 7 && sector == setdetsector) totEdep = totEdep + energy;
							if(pindex == i && detector == 7 && sector == setdetsector && layer == detECallayer)
							{
								EdepECalcluster = energy;
								hitx = x;
								hity = y;
								hitpath = path;
							}
						}
						for(int jf = 0; jf < detf.rows(); jf++)
						{
							short pindexf = detf.getShort("pindex", jf);
							byte detectorf = detf.getByte("detector", jf);
							byte sectorf = detf.getByte("sector", jf);
							byte layerf = detf.getByte("layer", jf);
							float pathf = detf.getFloat("path", jf);
							float timef = detf.getFloat("time", jf);
							if(pindexf == i && detectorf == 12 && sectorf == setdetsector && layerf == 2)
							{
								hittime = timef;
								Tvertex = timef - (pathf/(30*beta));
							}
						}
						
						float xeff = (hitx+(slope*hity)-(slope*intercept))/(1+(slope*slope));
						float yeff = (slope*xeff)+intercept;
						float leff = (float) Math.sqrt(((xeff-PMTx)*(xeff-PMTx))+((yeff-PMTy)*(yeff-PMTy)));
	
						byte id = -1;
						HipoDataBank clust = (HipoDataBank) event.getBank("ECAL::clusters");
						for(int k = 0; k < clust.rows(); k++)
						{
							byte sector = clust.getByte("sector", k);
							byte layer = clust.getByte("layer", k);
							float energy = clust.getFloat("energy", k);
							if(view == U) id = clust.getByte("idU", k);
							if(view == V) id = clust.getByte("idV", k);
							if(view == W) id = clust.getByte("idW", k);
							HipoDataBank hit = (HipoDataBank) event.getBank("ECAL::hits");
							byte ECalstrip = 0;
							if(sector == setdetsector && layer == detECallayer && energy == EdepECalcluster)
							{
								for(int l = 0; l < hit.rows(); l++)
								{
									byte strip = hit.getByte("strip", l);
									byte peakid = hit.getByte("peakid", l);
									byte EChitlayer = hit.getByte("layer", l);
									if(peakid == id && EChitlayer == ECallayer && strip == stripID) ECalstrip = strip;
								}
							}
							HipoDataBank adc = (HipoDataBank) event.getBank("ECAL::adc");
							int ECalstripADC = 0;
							int ECalstripref1ADC = 0;
							int ECalstripref2ADC = 0;
							if(ECalstrip == stripID)
							{
								for(int m = 0; m < adc.rows(); m++)
								{
									byte adcsector = adc.getByte("sector", m);
									byte adclayer = adc.getByte("layer", m);
									short adcstrip = adc.getShort("component", m);
									int ADC = adc.getInt("ADC", m);
									if(adcsector == setdetsector && adclayer == ECallayer && adcstrip == stripID) ECalstripADC = ADC;
									if(adcsector == setdetsector && adclayer == ECallayer && adcstrip == stripID-1) ECalstripref1ADC = ADC;
									if(adcsector == setdetsector && adclayer == ECallayer && adcstrip == stripID+1) ECalstripref2ADC = ADC;
								}
								if(ECalstripADC > ECalstripref1ADC && ECalstripADC > ECalstripref2ADC/* && ECalstripADC > 500*/)
								{
									HipoDataBank tdc = (HipoDataBank) event.getBank("ECAL::tdc");
									int TDCcounter = 0;
									int TDCvalue = 0;
									for(int n = 0; n < tdc.rows(); n++)
									{
										byte tdcsector = tdc.getByte("sector", n);
										byte tdclayer = tdc.getByte("layer", n);
										short tdcstrip = tdc.getShort("component", n);	
										int TDC = tdc.getInt("TDC", n);
										if(tdcsector == setdetsector && tdclayer == ECallayer && tdcstrip == stripID && TDC > 0)
										{
											TDCcounter++;
											TDCvalue = TDC;
										}
									}
									if(TDCcounter == 1 && beta > 0.99 && (hitpath/(30*beta))+(leff/18.1)+Tvertex > 500)
									{
										TDCi.add((float) TDCvalue);
										ADCi.add((float) ECalstripADC);
										l2i.add(leff*leff);
										l3i.add(leff*leff*leff);
										Ti.add((float) ((hitpath/(30*beta))+(leff/18.1)+Tvertex));
										Texpectedi.add((float) ((hitpath/(30*beta))+(leff/18.1)));
										HipoDataBank runrf = (HipoDataBank) event.getBank("RUN::rf");
										for(int o = 0; o < runrf.rows(); o++)
										{
											short rfid = runrf.getShort("id", o);
											float rftime = runrf.getFloat("time", o);
											if(rfid == 1) RFTi.add(rftime);		
										}
										hstrip.fill(stripID);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	static class FCN implements FCNBase {
		FCN(List<Float> TDCi, List<Float> ADCi, List<Float> Ti)
		{
			TDC = TDCi;
			ADC = ADCi;
			T = Ti;
		}
		public double errorDef()
		{
			return 1;
		}
		public double valueOf(double[] par)
		{
			double a0 = par[0];
			double a1 = par[1];
			double a2 = par[2];
			double chisq = 0;
			for(int n = 0; n < TDC.size(); n++)
			{
				double TDCn = TDC.get(n);
				double ADCn = ADC.get(n);
				double Tn = T.get(n);
				double delta = Tn-(-a0+(a1*TDCn)-(a2/Math.sqrt(ADCn)));
				chisq += (delta*delta/TDC.size());
			}
			return chisq;
		}
		private List<Float> TDC;
		private List<Float> ADC;
		private List<Float> T;
	}
	
	public static void main(String[] args) {
		short setdetsector = 2;
		//PCal: 1, ECInner: 4, ECOuter: 7	
		short setdetECallayer = 1;
		short U = 0;
		short V = 1;
		short W = 2;
		short setview = W;
		byte setstripID = 15;
		
		histos();
		
		short detECallayer = setdetECallayer;
		short view = setview;
		byte ECallayer = (byte) (detECallayer+view);
//		byte stripID = (byte) setstripID;
		
		ConstantProvider cp = GeometryFactory.getConstants(DetectorType.EC);
		ECGeant4Factory factory = new ECGeant4Factory(cp);
		ConstantProvider cp1 = GeometryFactory.getConstants(DetectorType.EC);
		PCALGeant4Factory factory1 = new PCALGeant4Factory(cp1);
		
		IndexedList<Double> CalConststore = new IndexedList<>(2);
		
		int stripmax = 36;
		if(ECallayer == 1) stripmax = 68;
		if(ECallayer == 2) stripmax = 62;
		if(ECallayer == 3) stripmax = 62;
		
		for(int strip = 1; strip <= stripmax; strip++)
		{
			reader.open("/Users/joshtanj/Documents/kpp_pass6_short_skim_809_810.hipo");
			pionreader.open("/Users/joshtanj/Documents/kpp_pass6_short_skim_761.hipo");
			
			byte stripID = (byte) strip;
			
			float PMTedgextot = 0;
			float PMTedgeytot = 0;
			float noPMTedgextot = 0;
			float noPMTedgeytot = 0;
			if(detECallayer != 1)
			{
				for(int a = 0; a <= detECallayer; a++)
				{
					float noPMTedgexa1 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(1).x;
					float noPMTedgexa3 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(3).x;
					float noPMTedgexa5 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(5).x;
					float noPMTedgexa7 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(7).x;
					float noPMTedgexa = (noPMTedgexa1+noPMTedgexa3+noPMTedgexa5+noPMTedgexa7)/4;
					noPMTedgextot = noPMTedgextot + noPMTedgexa;
					
					float noPMTedgeya1 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(1).y;
					float noPMTedgeya3 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(3).y;
					float noPMTedgeya5 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(5).y;
					float noPMTedgeya7 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(7).y;
					float noPMTedgeya = (noPMTedgeya1+noPMTedgeya3+noPMTedgeya5+noPMTedgeya7)/4;
					noPMTedgeytot = noPMTedgeytot + noPMTedgeya;
					
					float PMTedgexa0 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(0).x;
					float PMTedgexa2 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(2).x;
					float PMTedgexa4 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(4).x;
					float PMTedgexa6 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(6).x;
					float PMTedgexa = (PMTedgexa0+PMTedgexa2+PMTedgexa4+PMTedgexa6)/4;
					PMTedgextot = PMTedgextot + PMTedgexa;
					
					float PMTedgeya0 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(0).y;
					float PMTedgeya2 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(2).y;
					float PMTedgeya4 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(4).y;
					float PMTedgeya6 = (float) factory.getPaddle(setdetsector, view+(5*detECallayer)-19+(3*a), stripID).getVertex(6).y;
					float PMTedgeya = (PMTedgeya0+PMTedgeya2+PMTedgeya4+PMTedgeya6)/4;
					PMTedgeytot = PMTedgeytot + PMTedgeya;
				}
			}
			if(detECallayer == 1)
			{
				for(int a = 0; a <= 4; a++)
				{
					float noPMTedgexa1 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(1).x;
					float noPMTedgexa3 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(3).x;
					float noPMTedgexa5 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(5).x;
					float noPMTedgexa7 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(7).x;
					float noPMTedgexa = (noPMTedgexa1+noPMTedgexa3+noPMTedgexa5+noPMTedgexa7)/4;
					noPMTedgextot = noPMTedgextot + noPMTedgexa;
					
					float noPMTedgeya1 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(1).y;
					float noPMTedgeya3 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(3).y;
					float noPMTedgeya5 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(5).y;
					float noPMTedgeya7 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(7).y;
					float noPMTedgeya = (noPMTedgeya1+noPMTedgeya3+noPMTedgeya5+noPMTedgeya7)/4;
					noPMTedgeytot = noPMTedgeytot + noPMTedgeya;
					
					float PMTedgexa0 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(0).x;
					float PMTedgexa2 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(2).x;
					float PMTedgexa4 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(4).x;
					float PMTedgexa6 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(6).x;
					float PMTedgexa = (PMTedgexa0+PMTedgexa2+PMTedgexa4+PMTedgexa6)/4;
					PMTedgextot = PMTedgextot + PMTedgexa;
					
					float PMTedgeya0 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(0).y;
					float PMTedgeya2 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(2).y;
					float PMTedgeya4 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(4).y;
					float PMTedgeya6 = (float) factory1.getPaddle(setdetsector, view+detECallayer+(3*a), stripID).getVertex(6).y;
					float PMTedgeya = (PMTedgeya0+PMTedgeya2+PMTedgeya4+PMTedgeya6)/4;
					PMTedgeytot = PMTedgeytot + PMTedgeya;
				}
			}
			float PMTx = PMTedgextot/(detECallayer+1);
			float PMTy = PMTedgeytot/(detECallayer+1);
			float noPMTx = noPMTedgextot/(detECallayer+1);
			float noPMTy = noPMTedgeytot/(detECallayer+1);
			float slope = (PMTy-noPMTy)/(PMTx-noPMTx);
			float intercept = PMTy-(slope*PMTx);
			
			List<Float> TDCi = new ArrayList<>();
			List<Float> ADCi = new ArrayList<>();
			List<Float> l2i = new ArrayList<>();
			List<Float> l3i = new ArrayList<>();
			List<Float> Ti = new ArrayList<>();
			List<Float> Texpectedi = new ArrayList<>();
			List<Float> RFTi = new ArrayList<>();
			
			int eventCounter = 0;
			while(reader.hasEvent())
			{
				eventCounter++;
				processEvent(reader.getNextEvent(), eventCounter, setdetsector, detECallayer, U, V, W, view, ECallayer, stripID, slope, intercept, PMTx, PMTy, TDCi, ADCi, l2i, l3i, Ti, Texpectedi, RFTi);
				if(eventCounter%5000 == 0) System.out.println("e^- Event: " +eventCounter);
			}
			
			FCN theFCN = new FCN(TDCi, ADCi, Ti);
			
			MnUserParameters upar = new MnUserParameters();
			upar.add("a0", 10, 1);
			upar.add("a1", 0.02345);
			upar.add("a2", 100, 1);
			System.out.println("Initial parameters: "+upar); 
		    
		    System.out.println("start migrad");
			MnMigrad migrad = new MnMigrad(theFCN, upar); 
		    FunctionMinimum min = migrad.minimize(); 
		    if(!min.isValid()) 
		    { 
		    	//try with higher strategy 
			    	System.out.println("FM is invalid, try with strategy = 2."); 
			    	MnMigrad migrad2 = new MnMigrad(theFCN, min.userState(), new MnStrategy(2)); 
			    	min = migrad2.minimize(); 
		    }
		    System.out.println("minimum: "+min);
		    
		    IndexedList<Double> CalConst = new IndexedList<>(2);
		    
		    for(int consti = 0; consti < 3; consti++)
		    {
			    	CalConst.add(min.userParameters().value(consti), consti, 0);
		    }
		    CalConst.add(min.fval(), 3, 0);
		    
		    IndexedList<Float> TDCj = new IndexedList<>(2);
		    IndexedList<Float> ADCj = new IndexedList<>(2);
		    IndexedList<Float> Tj = new IndexedList<>(2);
		    IndexedList<Float> Texpectedj = new IndexedList<>(2);
		    IndexedList<Float> RFTj = new IndexedList<>(2);
		    List<Integer> naccept = new ArrayList<>();
		    List<Integer> nreject = new ArrayList<>();
		    
		    int size = Ti.size(); 
		    
		    for(int iteration = 0; iteration < 101; iteration++)
		    {
			    	int subindex = 0;
			    	int reject = 0;
			    	
			    	naccept.add(size);
			    	
			    	List<Float> TDCk = new ArrayList<>();
				List<Float> ADCk = new ArrayList<>();
				List<Float> Tk = new ArrayList<>();
				List<Float> Texpectedk = new ArrayList<>();
				List<Float> RFTk = new ArrayList<>();
		    	
			    for(int index = 0; index < size; index++)
			    {
				    	TDCj.add(TDCi.get(index),index, 0);
				    	ADCj.add(ADCi.get(index),index, 0);
				    	Tj.add(Ti.get(index),index, 0);
				    	Texpectedj.add(Texpectedi.get(index),index, 0);
				    	RFTj.add(RFTi.get(index),index, 0);
				    	
				    double deltaT = -CalConst.getItem(0, iteration)+(CalConst.getItem(1, iteration)*TDCj.getItem(index, iteration))-(CalConst.getItem(2, iteration)/Math.sqrt(ADCj.getItem(index, iteration)))
				    		-Tj.getItem(index, iteration);
				    
				    H1F hh = histGroups.getItem(iteration, stripID);
					hh.fill(deltaT);
			    }
			    
			    F1D f1 = new F1D("f1", "[amp]*gaus(x,[mean],[sigma])", 
						histGroups.getItem(iteration, stripID).getMean()
						-(2*histGroups.getItem(iteration, stripID).getRMS()), histGroups.getItem(iteration, stripID).getMean()-(2*histGroups.getItem(iteration, stripID).getRMS()));
				f1.setParameter(1, histGroups.getItem(iteration, stripID).getMean());
				f1.setParameter(2, histGroups.getItem(iteration, stripID).getRMS());
				DataFitter.fit(f1, histGroups.getItem(iteration, stripID), "Q");
				
				double leftcut = f1.getParameter(1)-(2*Math.abs(f1.getParameter(2)));
				double rightcut = f1.getParameter(1)+(2*Math.abs(f1.getParameter(2)));
				
				System.out.println("Iteration: " + iteration);
				System.out.println("Leftcut: " + leftcut);
				System.out.println("Rightcut: " + rightcut);
				System.out.println("Size: " + size);
			    
			    for(int index1 = 0; index1 < size; index1++)
			    {				
					double deltaT = -CalConst.getItem(0, iteration)+(CalConst.getItem(1, iteration)*TDCj.getItem(index1, iteration))-(CalConst.getItem(2, iteration)/Math.sqrt(ADCj.getItem(index1, iteration)))
							-Tj.getItem(index1, iteration);		
					
					if(deltaT > leftcut && deltaT < rightcut)
					{
						TDCj.add(TDCj.getItem(index1, iteration), subindex, iteration+1);
					    	ADCj.add(ADCj.getItem(index1, iteration), subindex, iteration+1);
					    	Tj.add(Tj.getItem(index1, iteration), subindex, iteration+1);
					    	Texpectedj.add(Texpectedj.getItem(index1, iteration), subindex, iteration+1);
					    	RFTj.add(RFTj.getItem(index1, iteration), subindex, iteration+1);
					    	subindex++;
					    	
					    	TDCk.add(TDCj.getItem(index1, iteration));
					    	ADCk.add(ADCj.getItem(index1, iteration));
					    	Tk.add(Tj.getItem(index1, iteration));
					    	Texpectedk.add(Texpectedj.getItem(index1, iteration));
					    	RFTk.add(RFTj.getItem(index1, iteration));
					}
					if(deltaT <= leftcut || deltaT >= rightcut) reject++;
			    }
			    
			    if(iteration == 0) nreject.add(reject);
			    if(iteration > 0) nreject.add(reject+nreject.get(iteration-1));
			    
			    size = subindex;
			    
			    FCN FCN = new FCN(TDCk, ADCk, Tk);
				
				MnUserParameters par = new MnUserParameters();
				par.add("a0", 10, 1);
				par.add("a1", 0.02345);
				par.add("a2", 100, 1);
				System.out.println("Initial parameters: "+par); 
			        
				MnMigrad migrad1 = new MnMigrad(FCN, par); 
			    FunctionMinimum min1 = migrad1.minimize();
			    if(!min1.isValid()) 
			    { 
				    	//try with higher strategy 
				    	System.out.println("FM is invalid, try with strategy = 2."); 
				    	MnMigrad migrad3 = new MnMigrad(FCN, min1.userState(), new MnStrategy(2)); 
				    	min1 = migrad3.minimize(); 
			    }
			    System.out.println("minimum: "+min1);
			    
			    for(int consti = 0; consti < 3; consti++)
			    {
			    		CalConst.add(min1.userParameters().value(consti), consti, iteration+1);
			    }
			    CalConst.add(min1.fval(), 3, iteration+1);
		    }
		    
		    for(int CCi = 0; CCi < 100; CCi++)
		    {
			    	if(Math.abs(CalConst.getItem(0, CCi)-CalConst.getItem(0, CCi+1)) == 0  
				    			&& Math.abs(CalConst.getItem(1, CCi)-CalConst.getItem(1, CCi+1)) == 0 
				    			&& Math.abs(CalConst.getItem(2, CCi)-CalConst.getItem(2, CCi+1)) == 0 
				    			&& Math.abs(CalConst.getItem(3, CCi)-CalConst.getItem(3, CCi+1)) == 0 && CalConst.getItem(3, CCi+1) > 0)
			    	{
			    		CalConststore.add(CalConst.getItem(0, CCi+1), 0, stripID);
			    		CalConststore.add(CalConst.getItem(1, CCi+1), 1, stripID);
			    		CalConststore.add(CalConst.getItem(2, CCi+1), 2, stripID);
			    		CalConststore.add(CalConst.getItem(3, CCi+1), 3, stripID);
			    		System.out.println("Calibration Constants for Strip: " + stripID);
			    		System.out.println("Iteration: " + CCi);
			    		break;
			    	}
		    }
		    
		    if(!CalConststore.hasItem(0, stripID) || CalConststore.getItem(3, stripID) > 10)
			{
		    		List<Float> ppTDCi = new ArrayList<>();
				List<Float> ppADCi = new ArrayList<>();
				List<Float> ppl2i = new ArrayList<>();
				List<Float> ppl3i = new ArrayList<>();
				List<Float> ppTi = new ArrayList<>();
				List<Float> ppTexpectedi = new ArrayList<>();
				List<Float> ppRFTi = new ArrayList<>();
				
				int ppeventCounter = 0;
				while(reader.hasEvent())
				{
					ppeventCounter++;
					processEvent(reader.getNextEvent(), ppeventCounter, setdetsector, detECallayer, U, V, W, view, ECallayer, stripID, slope, intercept, PMTx, PMTy, ppTDCi, ppADCi, ppl2i, ppl3i, ppTi, ppTexpectedi, ppRFTi);
					if(ppeventCounter%5000 == 0) System.out.println("e^- Event: " + ppeventCounter);
				}
				int pioneventCounter = 0;
				while(pionreader.hasEvent())
				{
					pioneventCounter++;				
					processPionEvent(pionreader.getNextEvent(), pioneventCounter, setdetsector, detECallayer, U, V, W, view, ECallayer, stripID, slope, intercept, PMTx, PMTy, ppTDCi, ppADCi, ppl2i, ppl3i, ppTi, ppTexpectedi, ppRFTi);
					if(pioneventCounter%5000 == 0) System.out.println("#pi^- Event: " + pioneventCounter);
				}

				FCN pptheFCN = new FCN(ppTDCi, ppADCi, ppTi);
				
				MnUserParameters ppupar = new MnUserParameters();
				ppupar.add("a0", 10, 1);
				ppupar.add("a1", 0.02345);
				ppupar.add("a2", 100, 1);
				System.out.println("Initial parameters: "+ppupar); 
			    
			    System.out.println("start migrad");
				MnMigrad ppmigrad = new MnMigrad(pptheFCN, ppupar); 
			    FunctionMinimum ppmin = ppmigrad.minimize(); 
			    if(!ppmin.isValid()) 
			    { 
			    	//try with higher strategy 
				    	System.out.println("FM is invalid, try with strategy = 2."); 
				    	MnMigrad migrad2 = new MnMigrad(pptheFCN, ppmin.userState(), new MnStrategy(2)); 
				    ppmin = migrad2.minimize(); 
			    }
			    System.out.println("minimum: "+ppmin);
			    
			    IndexedList<Double> ppCalConst = new IndexedList<>(2);
			    
			    for(int consti = 0; consti < 3; consti++)
			    {
				    	ppCalConst.add(ppmin.userParameters().value(consti), consti, 0);
			    }
			    ppCalConst.add(ppmin.fval(), 3, 0);
			    
			    IndexedList<Float> ppTDCj = new IndexedList<>(2);
			    IndexedList<Float> ppADCj = new IndexedList<>(2);
			    IndexedList<Float> ppTj = new IndexedList<>(2);
			    IndexedList<Float> ppTexpectedj = new IndexedList<>(2);
			    IndexedList<Float> ppRFTj = new IndexedList<>(2);
			    List<Integer> ppnaccept = new ArrayList<>();
			    List<Integer> ppnreject = new ArrayList<>();
			    
			    int ppsize = ppTi.size(); 
			    
			    for(int iteration = 0; iteration < 101; iteration++)
			    {
				    	int subindex = 0;
				    	int reject = 0;
				    	
				    	ppnaccept.add(ppsize);
				    	
				    	List<Float> TDCk = new ArrayList<>();
					List<Float> ADCk = new ArrayList<>();
					List<Float> Tk = new ArrayList<>();
					List<Float> Texpectedk = new ArrayList<>();
					List<Float> RFTk = new ArrayList<>();
			    	
				    for(int index = 0; index < ppsize; index++)
				    {
					    	ppTDCj.add(ppTDCi.get(index),index, 0);
					    	ppADCj.add(ppADCi.get(index),index, 0);
					    ppTj.add(ppTi.get(index),index, 0);
					    	ppTexpectedj.add(ppTexpectedi.get(index),index, 0);
					    	ppRFTj.add(ppRFTi.get(index),index, 0);
					    	
					    double deltaT = -ppCalConst.getItem(0, iteration)+(ppCalConst.getItem(1, iteration)*ppTDCj.getItem(index, iteration))-(ppCalConst.getItem(2, iteration)/Math.sqrt(ppADCj.getItem(index, iteration)))
					    		-ppTj.getItem(index, iteration);
					    
					    H1F hh = histGroups.getItem(iteration, stripID);
						hh.fill(deltaT);
				    }
				    
				    F1D f1 = new F1D("f1", "[amp]*gaus(x,[mean],[sigma])", 
							histGroups.getItem(iteration, stripID).getMean()
							-(2*histGroups.getItem(iteration, stripID).getRMS()), histGroups.getItem(iteration, stripID).getMean()-(2*histGroups.getItem(iteration, stripID).getRMS()));
					f1.setParameter(1, histGroups.getItem(iteration, stripID).getMean());
					f1.setParameter(2, histGroups.getItem(iteration, stripID).getRMS());
					DataFitter.fit(f1, histGroups.getItem(iteration, stripID), "Q");
					
					double leftcut = f1.getParameter(1)-(2*Math.abs(f1.getParameter(2)));
					double rightcut = f1.getParameter(1)+(2*Math.abs(f1.getParameter(2)));
					
					System.out.println("Iteration: " + iteration);
					System.out.println("Leftcut: " + leftcut);
					System.out.println("Rightcut: " + rightcut);
					System.out.println("Size: " + ppsize);
				    
				    for(int index1 = 0; index1 < ppsize; index1++)
				    {				
						double deltaT = -ppCalConst.getItem(0, iteration)+(ppCalConst.getItem(1, iteration)*ppTDCj.getItem(index1, iteration))-(ppCalConst.getItem(2, iteration)/Math.sqrt(ppADCj.getItem(index1, iteration)))
								-ppTj.getItem(index1, iteration);		
						
						if(deltaT > leftcut && deltaT < rightcut)
						{
							ppTDCj.add(ppTDCj.getItem(index1, iteration), subindex, iteration+1);
						    	ppADCj.add(ppADCj.getItem(index1, iteration), subindex, iteration+1);
						    	ppTj.add(ppTj.getItem(index1, iteration), subindex, iteration+1);
						    ppTexpectedj.add(ppTexpectedj.getItem(index1, iteration), subindex, iteration+1);
						    ppRFTj.add(ppRFTj.getItem(index1, iteration), subindex, iteration+1);
						    	subindex++;
						    	
						    	TDCk.add(ppTDCj.getItem(index1, iteration));
						    	ADCk.add(ppADCj.getItem(index1, iteration));
						    	Tk.add(ppTj.getItem(index1, iteration));
						    	Texpectedk.add(ppTexpectedj.getItem(index1, iteration));
						    	RFTk.add(ppRFTj.getItem(index1, iteration));
						}
						if(deltaT <= leftcut || deltaT >= rightcut) reject++;
				    }
				    
				    if(iteration == 0) ppnreject.add(reject);
				    if(iteration > 0) ppnreject.add(reject+ppnreject.get(iteration-1));
				    
				    ppsize = subindex;
				    
				    FCN FCN = new FCN(TDCk, ADCk, Tk);
					
					MnUserParameters par = new MnUserParameters();
					par.add("a0", 10, 1);
					par.add("a1", 0.02345);
					par.add("a2", 100, 1);
					System.out.println("Initial parameters: "+par); 
				        
					MnMigrad migrad1 = new MnMigrad(FCN, par); 
				    FunctionMinimum min1 = migrad1.minimize();
				    if(!min1.isValid()) 
				    { 
					    	//try with higher strategy 
					    	System.out.println("FM is invalid, try with strategy = 2."); 
					    	MnMigrad migrad3 = new MnMigrad(FCN, min1.userState(), new MnStrategy(2)); 
					    	min1 = migrad3.minimize(); 
				    }
				    System.out.println("minimum: "+min1);
				    
				    for(int consti = 0; consti < 3; consti++)
				    {
				    		ppCalConst.add(min1.userParameters().value(consti), consti, iteration+1);
				    }
				    ppCalConst.add(min1.fval(), 3, iteration+1);
			    }
			    
			    for(int CCi = 0; CCi < 100; CCi++)
			    {
				    	if(Math.abs(ppCalConst.getItem(0, CCi)-ppCalConst.getItem(0, CCi+1)) == 0  
					    			&& Math.abs(ppCalConst.getItem(1, CCi)-ppCalConst.getItem(1, CCi+1)) == 0 
					    			&& Math.abs(ppCalConst.getItem(2, CCi)-ppCalConst.getItem(2, CCi+1)) == 0 
					    			&& Math.abs(ppCalConst.getItem(3, CCi)-ppCalConst.getItem(3, CCi+1)) == 0 && ppCalConst.getItem(3, CCi+1) > 0)
				    	{
				    		CalConststore.add(ppCalConst.getItem(0, CCi+1), 0, stripID);
				    		CalConststore.add(ppCalConst.getItem(1, CCi+1), 1, stripID);
				    		CalConststore.add(ppCalConst.getItem(2, CCi+1), 2, stripID);
				    		CalConststore.add(ppCalConst.getItem(3, CCi+1), 3, stripID);
				    		System.out.println("Calibration Constants for Strip: " + stripID);
				    		System.out.println("Iteration: " + CCi);
				    		break;
				    	}
			    }
			}
		    pionreader.close();
		    reader.close();
		}
		
		System.out.println("Layer " + ECallayer);
		if(setview == U) System.out.println("View U");
		if(setview == V) System.out.println("View V");
		if(setview == W) System.out.println("View W");
		for(int stripi = 1; stripi <= stripmax; stripi ++)
		{
			System.out.println("Strip: " + stripi);
			System.out.println("a0: " + CalConststore.getItem(0, stripi));
			System.out.println("a1: " + CalConststore.getItem(1, stripi));
			System.out.println("a2: " + CalConststore.getItem(2, stripi));
			System.out.println("#chi^2: " + CalConststore.getItem(3, stripi));
			
			//Set strip exceptions
			if(CalConststore.hasItem(0,stripi))
			{
				ha0f.fill(stripi, CalConststore.getItem(0, stripi));
				ha1f.fill(stripi, CalConststore.getItem(1, stripi));
				ha2f.fill(stripi, CalConststore.getItem(2, stripi));				
				hchisqf.fill(stripi, CalConststore.getItem(3, stripi));				
			}
		}
		System.out.println(" ");
		System.out.println("NOTE!!!");
		for(int stripchi = 1; stripchi <= 36; stripchi ++)
		{
			if(!CalConststore.hasItem(1, stripchi))
			{
				System.out.println("Strip " + stripchi + " has no entry.");
			}
			if(CalConststore.hasItem(3, stripchi) && CalConststore.getItem(3, stripchi) > 10)
			{
				System.out.println("Strip " + stripchi + " #chi^2: " + (CalConststore.getItem(3, stripchi)));
			}
		}
			    
		JFrame frame = new JFrame("frame");
		frame.setSize(500, 500);
		EmbeddedCanvas can = new EmbeddedCanvas();
		frame.add(can);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		can.divide(1, 1);
		
		can.cd(0);
		hstrip.setTitle("Layer " + ECallayer + " View " + view + " Peaks");
		hstrip.setTitleX("Strip");
		hstrip.setTitleY("Counts");
		hstrip.setOptStat(10);
		can.getPad(0).setTitleFontSize(14);
		can.getPad(0).setAxisTitleFontSize(14);
		can.getPad(0).setAxisLabelFontSize(14);
		can.getPad(0).setStatBoxFontSize(14);
		can.draw(hstrip);
		
		JFrame frame1 = new JFrame("frame1");
		frame1.setSize(1000, 1000);
		EmbeddedCanvas can1 = new EmbeddedCanvas();
		frame1.add(can1);
		frame1.setLocationRelativeTo(null);
		frame1.setVisible(true);
		can1.divide(2, 2);
		
		can1.cd(0);
		ha0f.setTitle("a0 vs. Strip");
		ha0f.setTitleX("Strip");
		ha0f.setTitleY("a0 [ns]");
		can1.getPad(0).setTitleFontSize(14);
		can1.getPad(0).setAxisTitleFontSize(14);
		can1.getPad(0).setAxisLabelFontSize(14);
		can1.getPad(0).getAxisZ().setRange(0, 1.5);
		can1.draw(ha0f);
		can1.cd(1);
		ha1f.setTitle("a1 vs. Strip");
		ha1f.setTitleX("Strip");
		ha1f.setTitleY("a1 [ns/[TDC]]");
		can1.getPad(1).setTitleFontSize(14);
		can1.getPad(1).setAxisTitleFontSize(14);
		can1.getPad(1).setAxisLabelFontSize(14);
		can1.getPad(1).getAxisZ().setRange(0, 1.5);
		can1.draw(ha1f);
		can1.cd(2);
		ha2f.setTitle("a2 vs. Strip");
		ha2f.setTitleX("Strip");
		ha2f.setTitleY("a2 [ns[sqrt(ADC)]]");
		can1.getPad(2).setTitleFontSize(14);
		can1.getPad(2).setAxisTitleFontSize(14);
		can1.getPad(2).setAxisLabelFontSize(14);
		can1.getPad(2).getAxisZ().setRange(0, 1.5);
		can1.draw(ha2f);
		can1.cd(3);
		hchisqf.setTitle("#chi^2 vs. Strip");
		hchisqf.setTitleX("Strip");
		hchisqf.setTitleY("#chi^2");
		can1.getPad(3).setTitleFontSize(14);
		can1.getPad(3).setAxisTitleFontSize(14);
		can1.getPad(3).setAxisLabelFontSize(14);
		can1.getPad(3).getAxisZ().setRange(0, 1.5);
		can1.draw(hchisqf);
	}
}
