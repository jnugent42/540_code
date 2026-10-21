//This is a code that takes in laser data and then gets the TTS out
bool passVeto(double y[], double vetoHigh, int lenWindow){
	for (int i = 0; i<lenWindow; i++){
		if (y[i]>vetoHigh){
			std::cout << "passed the " << vetoHigh << " veto" << std::endl;
		       	return true;
		}	
	}
	return false;
	
}


void TTS(const char * PMTform, const char * face_position)  {
	//variable definition
	double pe1; //this is the number of pe of the laser file
	double volts_trigger; //this is the adc of the trigger from the file
	int lenWindow = 2000; //this is one single scope window = one trigger -> CAREFULL, it depends
	//on your sampling rate and window size (on the scope)
	int numberWindCheck = 1000 * 500;
	double x[lenWindow], y[lenWindow], trigger[lenWindow]; //this is the x and y 
	//that we want to zoom in to have an idea
	double threshold = -0.5;
	double distToLocalTrig; //for each pulse
	//int n_largeWindows = 0; // the number of large windows we have circled through
	double trigMin = 9999;
	double trigMax = -9999;
	double laserT0 = -9999.;
	double pmtT0 = -9999.;
	double laserRate = 1000; //interval in ns between two laser pulses
	std::string PMT_info = "Monitor PMT - blue splitter";

	//double TTS[numberWindCheck];
       	//double pe_max[numberWindCheck];
	double max_pe_in_window = 10;
	int offsetWind = 800; //sometimes the window starts just after the laser trigger -> that causes a mess
	//adding an offset is a simple fix
	bool laserTriggerFound = false;
	bool pmtTriggerFound = false;
	double laser_thresh = -0.5; //the drop after which we consider the laser has lased
	double offset = -(100-80.40);//from the scope
	//int count = 0;
	double lastLaserThresh = 0;
	double initialLaserThresh = 0;
   	auto hcol1 = new TH2F("QvsT","Collected charge vs hit time",100,100,-1000,100,-10,3);

	//File opening
	TFile *f = new TFile(Form("/home/a/WCTE540/PMT_data/TTSLB/%sPMT/combined_datasets/hadded_all_%s_TTSLB_%s_1ms.root", PMTform, PMTform, face_position));
	TTree *t1 = (TTree*)f->Get("PMT1data");
	t1->SetBranchAddress("pe",&pe1);
        Int_t nentries = (Int_t)t1->GetEntries();
	TTree *T = (TTree*)f->Get("Trigger");
	T->SetBranchAddress("volts",&volts_trigger);
	

	//Where we store the raw data
	TH1D *h = new TH1D("hLASnormal","pe distribution normal",100,0,10);
	//here we sdtore the TTS
	TH1D *hTTS = new TH1D("hTTS",Form("TTS: PMT thresh %.2f p.e., Laser thres %.2fV (%i triggers @ 2GHz)", threshold, laser_thresh, numberWindCheck),450, 60, 95);
	
	//That's the canvas for showing each window separately
	TCanvas *c = new TCanvas("c", "One window", 0, 0, 600, 800);
	TCanvas *cTTS = new TCanvas("cTTS", "The TTS canvas", 0, 0, 600, 800);
	
	laserT0 = -9999;
	pmtT0 = -9999;
	for (Int_t k=0; k<numberWindCheck; k++){
		//re-initialise the points at which we cross thresholds 
		laserTriggerFound = false;
		pmtTriggerFound = false;
		max_pe_in_window = 10.;
		for (Int_t i=0; i<lenWindow+1; i++) {
	                t1->GetEntry(i + k*lenWindow + offsetWind);
	                T->GetEntry(i + k*lenWindow+ offsetWind);
	                h->Fill(pe1);
			int j = i%lenWindow; 
			x[j] = j; //look at the raw output of the file
	                y[j] = -pe1;
			if (-pe1 < max_pe_in_window) max_pe_in_window = -pe1;
			trigger[j] = volts_trigger;
			//note: here we assume that we stat at 0 and then drop to -1 during the period:
			//really might not always be the case
			if (volts_trigger < laser_thresh && j>0 && trigger[j-1] > laser_thresh && !laserTriggerFound){
				laserTriggerFound = true;
			        //laserT0 = j; //-> without interpolation
				//linear interpolation part
				double a_interpolation = (trigger[j]-trigger[j-1]);
				double b_interpolation = trigger[j] - a_interpolation * j;
				laserT0 = (laser_thresh - b_interpolation) / a_interpolation;
				//if ((laserT0-lastLaserThresh) > 4 && lastLaserThresh!= -9999){
				//       	laserT0 = laserT0 - 28;
				//	std::cout << "+" << std::endl;
				//}
				//
				//else if ((laserT0-lastLaserThresh) < -4 && lastLaserThresh!= -9999){
				//       	laserT0 = laserT0 + 1;
				//	std::cout << "-" << std::endl;
				//}
				////lastLaserThresh = laserT0 * 0.5;
				//std::cout << laserT0-lastLaserThresh << std::endl;	
				//will have to change 1000 i.e the number of triggers/window
				if (k%1000 == 0) lastLaserThresh = laserT0; 
				if (k == 0) initialLaserThresh = laserT0;
				laserT0 = laserT0 - lastLaserThresh + initialLaserThresh;
			}
			if (-pe1 < threshold && j>0 && y[j-1] > threshold && !pmtTriggerFound){
				pmtTriggerFound = true;
			        //pmtT0 = j; -> without interpolation	
				//linear interpolation part
				double a_interpolation = (y[j]-y[j-1]);
				double b_interpolation = y[j] - a_interpolation * j;
				pmtT0 = (threshold - b_interpolation) / a_interpolation;
				//std::cout << pmtT0 << std::endl;
				//std::cout << a_interpolation << std::endl;
				//std::cout << b_interpolation << std::endl;
			}
			//now try to check the threshold - need to ensure we are on the way up
	//This is not done correctly and is too complicated: scrape out for now		
	//		if (y[j] > threshold && y[j] > y[j-2] && j > 0  && j< lenWindow && y[j-1] < threshold) {
	//			//do a linear interpolation
	//			double gradient = (y[j]-y[j-1])/(x[j]-x[j-1]);
	//			double offset = y[j] - x[j] * gradient;
	//			double x_thresh = (threshold - offset)/gradient;
	//			distToLocalTrig = x_thresh;
	//			std::cout << "The threshold " << threshold << "  was crossed at " << x_thresh << std::endl;
	//		       //	<< " y[i], y[i-1], y[i-2]: " << y[i] << ", " << y[i-1] << ", " << y[i-2]<< std::endl;
	//		}
	//		if (j == 0 && i>0 && distToLocalTrig > 0){
	//		//at the end of the window, we take the measure only if it passes the veto
	//			if (passVeto(y,1.1, lenWindow)) {
	//				hTTS -> Fill(distToLocalTrig*0.5 - offset);
	//				if (distToLocalTrig*0.5 -offset< trigMin) trigMin = distToLocalTrig * 0.5 - offset;	
	//				if (distToLocalTrig*0.5-offset > trigMax) trigMax = distToLocalTrig * 0.5 - offset;	
	//			}
	//		}
	//
	//		//in the raising time: we are quickly increasing no risk of jitter
	
		}//end of reading through the pe
		//some ploting
//
		//std::cout << "plotting works" << std::endl;
		//careful, the conversion in seconds here depends on your sampling rate, I use 2GHz
		if (laserTriggerFound and pmtTriggerFound) {
			double TTS = -99999.;
			std::cout << "Both trigger found for window number " << k << " : Laser threshol crossed at " << laserT0 * 0.5 << " ns and PMT threshold has been crossed at " << pmtT0 * 0.5 << "ns delay is " << (pmtT0-laserT0) * 0.5 << "ns" << std::endl;
			
			if ((pmtT0-laserT0)> 0) TTS = (pmtT0-laserT0) * 0.5;
			
			else {
				std::cout << "Careful -> looking in the previous window!!!" << std::endl;
				TTS = laserRate + (pmtT0-laserT0) * 0.5;
			}
			
			hTTS->Fill(TTS);
			//quick study of TTS vs pe
			//TTS[count]=((pmtT0-laserT0) * 0.5);
			//pe_max[count] = (max_pe_in_window);
			//count = count + 1;
			
      			hcol1->Fill(TTS,max_pe_in_window);
			//std::cout << "Trying to plot" << std::endl;
			//very useful for checnking specific windows -> just need to change where we stop
			//TGraph *g1 = new TGraph(lenWindow, x, y);
			//TGraph *g2 = new TGraph(lenWindow, x, trigger);
			//g2->SetMarkerColor(kBlue);
			//g1->SetMarkerColor(kRed);
			//c->cd();
			//g2->Draw("");
			//g1->Draw("same *");
			//c->Draw();
		}
	}
	std::cout << "Min delay (ns) " << trigMin << " mean dealy (ns) " << hTTS->GetMean(1) << " max delay (ns) " << trigMax << " max - min = " << trigMax - trigMin << " ns "<< std::endl;

	c->cd();
	//TGraph *g = new TGraph(count,TTS,pe_max);
	//g->Draw("a *");
   	hcol1->Draw("COLZ");
	hcol1->GetXaxis()->SetTitle("Delay between laser and pmt threshold crossing (ns)");
	//hcol1->SetTitle("Hit time vs. collected charge");
	hcol1->GetYaxis()->SetTitle("Max depth of PMT pulse (arbitrary.)");
   	TPad *newpad=new TPad("newpad","a transparent pad",0.04,0.22,0.87,0.90);
	newpad->SetFillStyle(4000);
	newpad->Draw();
	newpad->cd();
	TText *t = new TText(.5,.5,Form("%s", PMT_info.c_str()));
   	t->SetX(0.12107383);
   	t->SetY(0.950383);
   	t->SetTextColor(1);
   	t->SetTextFont(22);
   	t->SetTextSize(0.04);
   	t->SetTextAngle(0);
   	t->Draw();
	c->Draw();
	
	cTTS->cd();
	hTTS->Draw();
	hTTS->GetXaxis()->SetTitle("Delay between laser and pmt threshold crossing (ns)");
	hTTS->GetYaxis()->SetTitle(Form("Occurencies/0.1ns"));
	//hTTS->SetTitle(Form("PMT threshold %.2f, Laser Theshold %.2f %i triggers of 1us at2GHz \n with interpolation", threshold, laser_thresh, numberWindCheck));
	hTTS->SetTitleSize(18, "t");
	hTTS->Fit("gaus");
	gStyle->SetOptStat("n");
	gStyle->SetOptFit(1011);
	//listOfLines->Remove(tconst);
	TPad *newpad2=new TPad("newpad","a transparent pad",0.04,0.06,0.86,0.90);
        newpad2->SetFillStyle(4000);
        newpad2->Draw();
        newpad2->cd();
        TText *t2 = new TText(.5,.5,Form("%s - %s", PMTform, face_position));
        t2->SetX(0.12107383);
        t2->SetY(0.950383);
        t2->SetTextColor(1);
        t2->SetTextFont(22);
        t2->SetTextSize(0.04);
        t2->SetTextAngle(0);
        t2->Draw();
	cTTS->Draw();

	new TFile(Form("/home/a/WCTE540/PMT_data/TTSLB/%sPMT/combined_datasets/hadded_all_%s_TTS_%s_1pe_1ms_TTSHist.root", PMTform, PMTform, face_position), "RECREATE");
        //h_DR1sPD->Write();
        //h_LAS1sPD->Write();
        hTTS->Write();
        hcol1->Write();
        //h_DR1msN->Write();
        //h_DR1msC->Write();
        //h_DRtoFit->Write();
        //h_LAStoFit->Write();


}
