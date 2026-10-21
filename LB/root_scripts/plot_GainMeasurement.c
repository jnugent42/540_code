void plot_GainMeasurement(const char * PMTform, const char * face_position) {
	double pe_DR1sPD, pe_LAS1sPD, pe_LAS1msN, pe_DR1msN;
	double factor = 1.;
	double thresh = 0.002; //we need at least that number of PE/volts max to be considered a pulse
	double threshInit = 0.002; //this is the first point we consider to be part the pulse - landau
	int lenWindow = 200;
	int smallWindow = 50;
	int nPointsPrior = 5;
	double onePEpeak_min = 0.005; // what points we fit in the one pe peak
	double onePEpeak_max = 0.03;
	TCanvas *c = new TCanvas("c", "Example", 0, 0, 600, 800);
	TCanvas *c1 = new TCanvas("c1", "Example", 0, 0, 600, 800);
	
	//const char * PMTform = "Circle"; 
	//const char * face_position = "180deg154mA"; 
	
	////first file: peak detect dark rate
	//TFile *f_DR1sPD = new TFile("SquarePMTResponseDR1s1MHzPD_allPedestal.root");
   	//TTree *t_DR1sPD = (TTree*)f_DR1sPD->Get("PMT1data");
	//t_DR1sPD->SetBranchAddress("pe",&pe_DR1sPD);
	//Int_t nentries_DR1sPD = (Int_t)t_DR1sPD->GetEntries();
	//TH1D *h_DR1sPD   = new TH1D("hDR1sPD","hDR1sPD",100, -0.1, 10);
	//for (Int_t i=0; i<nentries_DR1sPD; i++) { 
	//	t_DR1sPD->GetEntry(i);
	//	h_DR1sPD->Fill(pe_DR1sPD);
	//}

	//std::cout << "Filled the first histogram" << std::endl;
	//
	//
	////second file: peak detect laser
	//TFile *f_LAS1sPD = new TFile("SquarePMTResponseLAS1s1MHzPD_allPedestal.root");
   	//TTree *t_LAS1sPD = (TTree*)f_LAS1sPD->Get("PMT1data");
	//t_LAS1sPD->SetBranchAddress("pe",&pe_LAS1sPD);
	//Int_t nentries_LAS1sPD = (Int_t)t_LAS1sPD->GetEntries();
	//TH1D *h_LAS1sPD   = new TH1D("hLAS1sPD","hLAS1sPD",100, -0.1, 10);
	//for (Int_t i=0; i<nentries_LAS1sPD; i++) { 
	//	t_LAS1sPD->GetEntry(i);
	//	h_LAS1sPD->Fill(pe_LAS1sPD);
	//}

	//std::cout << "Filled the second  histogram" << std::endl;
	
	//third file: normal 1ms dark rate
	TFile *f_LAS1msN = new TFile(Form("/home/a/WCTE540/PMT_data/TTSLB/%sPMT/combined_datasets/hadded_all_%s_TTSLB_%s_1ms.root", PMTform, PMTform, face_position));
   	TTree *t_LAS1msN = (TTree*)f_LAS1msN->Get("PMT1data");
	t_LAS1msN->SetBranchAddress("volts",&pe_LAS1msN); //can change to pe if we are sure
	Int_t nentries_LAS1msN = (Int_t)t_LAS1msN->GetEntries();
	TH1D *h_LAS1msN  = new TH1D("hLAS1msN","hLAS1msN",200, -0.001, 0.124);
	TH1D *h_LAS1msC  = new TH1D("hLAS1msC","hLAS1msC",200, -0.001, 0.124);
	
	
	TH1D *h_LAStoFit  = new TH1D("hLASpeak","hLAS50nspeak",32, 0.005, 0.025);
	TH1D *h_DRtoFit  = new TH1D("hDRpeak","hDR50nspeak",32, 0.005, 0.025);
	
	
	int n_largeWindows = 0;
	double maxPEinThisWindow = -9999;
	double x[lenWindow], y[lenWindow];
	for (Int_t i=0; i<nentries_LAS1msN; i++) { 
		t_LAS1msN->GetEntry(i);
		h_LAS1msN->Fill(pe_LAS1msN);

		//Now the landau fit
		if (i-lenWindow*n_largeWindows != 0){
        	       if (pe_LAS1msN>=maxPEinThisWindow) maxPEinThisWindow = pe_LAS1msN;
        	       //add points to the buffer of the given window size
        	       int bufferPos = i%lenWindow;
        	       x[bufferPos] = bufferPos;
        	       y[bufferPos] = pe_LAS1msN;
        	 }

        	if (i-lenWindow*n_largeWindows == 0){
			double point = -9999;
			//std::cout << "End at "<< i << std::endl;
        	        if (maxPEinThisWindow > thresh) {
        	                int count = 0;
        	                int a = 0;
        	                double x_buf[smallWindow], y_buf[smallWindow];
        	                for (int j = 0; j<=lenWindow; j++){
        	                          if(y[j+nPointsPrior] >= threshInit and a ==0){
        	                                   x_buf[0] = 0;
        	                                   y_buf[0] = y[j];
        	                                   a = 1;
        	                                   count = 1;
        	                           }
        	                           else if (a == 1 and count <= smallWindow){
        	                                   x_buf[count] = count;
        	                                   y_buf[count] = y[j];
        	                                   count += 1;
        	                           }
        	                           else if (count == smallWindow) break; 
        	                }
        	                        
			TGraph *g1 = new TGraph(smallWindow, x_buf, y_buf);
        	        TFitResultPtr r = g1->Fit("landau", "S0"); //so that it is not drawn
        	        Double_t par0   = r->Parameter(0);
        	        Double_t par1   = r->Parameter(1);
        	        Double_t par2   = r->Parameter(2);
        	        TF1 * lan = new TF1("lan", "[0] * TMath::Landau(x, [1], [2], 0)", 0, smallWindow);
			lan->SetParameters(par0, par1, par2);
        		double m = lan->GetMaximum();
			point = (max(m, maxPEinThisWindow));
			}

			else point = maxPEinThisWindow;
        	        n_largeWindows += 1;
        	        maxPEinThisWindow = -9999.;
			h_LAS1msC->Fill(point);
			if (point<onePEpeak_max && point> onePEpeak_min ) h_LAStoFit->Fill(point);
		}//end of it is the past point of the window	
	}// end of readin the third file

	//fourth file: normal 1ms Dark Rate
//	TFile *f_DR1msN = new TFile("SquarePMTResponseDR1ms4GHzN_allPedestal_gain1.37542e6.root");
//        TTree *t_DR1msN = (TTree*)f_DR1msN->Get("PMT1data");
//        t_DR1msN->SetBranchAddress("pe",&pe_DR1msN);
//        Int_t nentries_DR1msN = (Int_t)t_DR1msN->GetEntries();
//        TH1D *h_DR1msN  = new TH1D("hDR1msN","hDR1msN",100, -0.1, 10);
//        TH1D *h_DR1msC  = new TH1D("hDR1msC","hDR1msC",100, -0.1, 10);
//
//        n_largeWindows = 0;
//        maxPEinThisWindow = -9999;
//        for (Int_t i=0; i<nentries_DR1msN; i++) {
//                t_DR1msN->GetEntry(i);
//                h_DR1msN->Fill(pe_DR1msN);
//
//                //Now the landau fit
//                if (i-lenWindow*n_largeWindows != 0){
//                       if (pe_DR1msN>=maxPEinThisWindow) maxPEinThisWindow = pe_DR1msN;
//                       //add points to the buffer of the given window size
//                       int bufferPos = i%lenWindow;
//                       x[bufferPos] = bufferPos;
//                       y[bufferPos] = pe_DR1msN;
//                 }
//
//                if (i-lenWindow*n_largeWindows == 0){
//                        double point = -9999;
//                        //std::cout << "End at "<< i << std::endl;
//                        if (maxPEinThisWindow > thresh) {
//                                int count = 0;
//                                int a = 0;
//                                double x_buf[smallWindow], y_buf[smallWindow];
//                                for (int j = 0; j<=lenWindow; j++){
//                                          if(y[j+nPointsPrior] >= threshInit and a ==0){
//                                                   x_buf[0] = 0;
//                                                   y_buf[0] = y[j];
//                                                   a = 1;
//                                                   count = 1;
//                                           }
//                                           else if (a == 1 and count <= smallWindow){
//                                                   x_buf[count] = count;
//                                                   y_buf[count] = y[j];
//                                                   count += 1;
//                                           }
//                                           else if (count == smallWindow) break;
//                                }
//
//                        TGraph *g1 = new TGraph(smallWindow, x_buf, y_buf);
//                        TFitResultPtr r = g1->Fit("landau", "S0"); //so that it is not drawn
//                        Double_t par0   = r->Parameter(0);
//                        Double_t par1   = r->Parameter(1);
//                        Double_t par2   = r->Parameter(2);
//                        TF1 * lan = new TF1("lan", "[0] * TMath::Landau(x, [1], [2], 0)", 0, smallWindow);
//                        lan->SetParameters(par0, par1, par2);
//                        double m = lan->GetMaximum();
//                        point = (max(m, maxPEinThisWindow));
//                        }
//
//                        else point = maxPEinThisWindow;
//                        n_largeWindows += 1;
//                        maxPEinThisWindow = -9999.;
//                        h_DR1msC->Fill(point);
//			if (point<onePEpeak_max && point> onePEpeak_min ) h_DRtoFit->Fill(point);
//                }//end of it is the past point of the window
//        }// end of readin the third file
//
//
//


	c->cd();
	TLegend * legend = new TLegend(0.35, 0.65, 0.86, 0.82);
        //legend->AddEntry(h_DR1sPD, "Peak Detect 1s Dark Rate");
        //legend->AddEntry(h_LAS1sPD, "Peak Detect 1s Laser");
        legend->AddEntry(h_LAS1msN, "Normal 1ms Laser");
        legend->AddEntry(h_LAS1msC, "Normal Landau Corrected  1ms Laser");
        //legend->AddEntry(h_DR1msN, "Normal 1ms Dark Rate");
        //legend->AddEntry(h_DR1msC, "Normal Landau Corrected  1ms Dark Rate");
	//h_DR1sPD->Draw("HIST");	
	//h_LAS1sPD->SetLineColor(kRed);
	//h_LAS1sPD->Draw("same HIST");	
	h_LAS1msN->SetLineColor(kBlue);
	h_LAS1msN->GetXaxis()->SetTitle("volts");
	h_LAS1msN->GetYaxis()->SetTitle("Occurencies/625nV");
	h_LAS1msN->SetTitle(Form("h_LAS1msN %s %s", PMTform, face_position));
	h_LAS1msN->SetLineStyle(kDashed);
	h_LAS1msN->Draw("same HIST");
	h_LAS1msC->SetLineColor(kBlue);
	h_LAS1msC->Draw("same HIST");


	c1->cd();
	std::cout << "The Dark Rate peak " <<PMTform <<", GAIN 1.37542e6 fitted between " << onePEpeak_min << " and " << onePEpeak_max<< " is: " << std::endl; 
	//h_DRtoFit->Fit("gaus");
	std::cout << "The LASER peak " <<PMTform <<", GAIN 1.37542e6 fitted between " << onePEpeak_min << " and " << onePEpeak_max<< " is: " << std::endl; 
	h_LAStoFit->Fit("gaus");
	//h_DR1msN->SetLineColor(kRed);
	//h_DR1msN->SetLineStyle(kDashed);
	//h_DR1msN->Draw("same HIST");
	//h_DR1msC->SetLineColor(kRed);
	//h_DR1msC->Draw("same HIST");
	legend->Draw("same");	
	new TFile(Form("/home/a/WCTE540/PMT_data/TTSLB/%sPMT/combined_datasets/hadded_all_%s_TTS_%s_1pe_1ms_LandauCorrectedVolts.root", PMTform, PMTform, face_position), "RECREATE");
	//h_DR1sPD->Write();
	//h_LAS1sPD->Write();
	h_LAS1msN->Write();
	h_LAS1msC->Write();
	//h_DR1msN->Write();
	//h_DR1msC->Write();
	//h_DRtoFit->Write();
	h_LAStoFit->Write();	
}
/*
	
	TFile *f2 = new TFile("SquarePMTResponseLAS1msN_1sNEGto0_tryOneFifth.root");
        TTree *t2 = (TTree*)f2->Get("PMT1data");
        t2->SetBranchAddress("pe",&pe2);
        Int_t nentries2 = (Int_t)t2->GetEntries();
        //TH1D *h2   = new TH1D("h - 1000ms","pe distribution PD pos 1000 * 1ms",100,0,10);
        TH1D *h3  = new TH1D("h3 - laser corrected","pe distribution PD pos 1000 * 1ms - MAX over 200samples window",100,0,10);
        TH1D *h4  = new TH1D("h4 - laser uncorrected","laser uncorrected",100,0,10);
	double maxPEinThisWindow = -9999.;
	int n_largeWindows = 0;
	double thresh = 0.005;
	int n_pe = 0;
	int lenWindow = 2000;
	int smallWindow = 50;
	std::cout << nentries2;
	double x[lenWindow], y[lenWindow];
        for (Int_t i=0; i<nentries2; i++) {
                t2->GetEntry(i);
		//std::cout << n_largeWindows << std::endl;
                h4->Fill(pe2);
		//if (pe2 <= 1.4 and pe2 >= 0.5) n_highRate_pePeak +=1;
		//if (pe2 <= 0.05) n_highRate_pedestal +=1;
		//correct for the different rates
		//std::cout << i << " " << nentries2 << std::endl;
		if (i-lenWindow*n_largeWindows != 0){
			if (pe2>=maxPEinThisWindow) maxPEinThisWindow = pe2;
			//add points to the buffer of the given window size
			int bufferPos = i%lenWindow;
			//std::cout << bufferPos << std::endl;
				//(int) - n_largeWindows * 4000;
			x[bufferPos] = bufferPos;
			y[bufferPos] = pe2;
		}
		
		//std::cout << maxPEinThisWindow << std::endl;
		
		if (i-lenWindow*n_largeWindows == 0){

			if (maxPEinThisWindow > thresh) {
				std::cout << "Indeed: " << maxPEinThisWindow << std::endl;
				int count = 0;
				int a = 0;
				double x_buf[smallWindow], y_buf[smallWindow];
				for (int j = 0; j<=4000; j++){
					if(y[j+4] >= thresh and a ==0){
						std::cout << j << std::endl; 
						x_buf[0] = 0;
						y_buf[0] = y[j];
						a = 1;
						count = 1;
					}
					else if (a == 1 and count <= smallWindow){
						x_buf[count] = count;
						y_buf[count] = y[j];
						count += 1;
					}
					else if (count == smallWindow) break; //assume there is only one hit in this window
				}
				std::cout << "YUP" << x_buf[25] << " " << y[25]<< " " << x[25] << " "<< y_buf[25] << " " << y_buf[26]<< std::endl;
				//extract a cleaner maximum
				TGraph *g1 = new TGraph(smallWindow, x_buf, y_buf);
        			TFitResultPtr r = g1->Fit("landau", "S0"); //so that it is not drawn
				Double_t par0   = r->Parameter(0);
        			Double_t par1   = r->Parameter(1);
        			Double_t par2   = r->Parameter(2);
				TF1 * lan = new TF1("lan", "[0] * TMath::Landau(x, [1], [2], 0)", 0, smallWindow);
				lan->SetParameters(par0, par1, par2);
        			double m = lan->GetMaximum();
				h3->Fill(max(m, maxPEinThisWindow));
				std::cout << "MAx here is " << m << std::endl;
			        g1->Draw();
				lan->Draw("same");
				n_pe += 1;	
			}

			else h3->Fill(maxPEinThisWindow);
			n_largeWindows += 1;
			maxPEinThisWindow = -9999.;
		}
        }

	TLegend * legend = new TLegend(0.35, 0.65, 0.86, 0.82);
	legend->AddEntry(h, "Peak Detect 1s Dark Rate");
	legend->AddEntry(h3, "Fitted corrected LASER 1ms - max (fit max, datapoints)");
	c->cd(1);
	
	TFile *f3 = new TFile("comparisonFile_full.root");
        TH1D *h2 = (TH1D*)f3->Get("h3 - h corrected");
	h2->SetLineColor(kRed);
	h3->SetLineColor(kGreen);
	h4->SetLineColor(kBlack);
	h4->Scale(5);
	h3->Scale(5);
	//h3->Scale(20);
	legend->AddEntry(h2, "Fitted corrected Dark Rate 1ms - fit max");
	legend->AddEntry(h4, "Not corrected LASER 1ms");
	//h2->Scale(factor/h2->GetEntries());
	//h->Scale(factor/h->GetEntries());
	h->GetXaxis()->SetTitle("P.E.");
	h2->Draw("HIST");
	h3->Draw("same HIST"); 
	h4->Draw("same HIST"); 
	TH1D * h_ref = (TH1D *)h->Clone(); 
	TH1D * h2_ref = (TH1D *)h2->Clone(); 
	legend->Draw("same");
	h->Draw("same HIST");
	std::cout << "Full number of pe with the DR is " << n_pe << std::endl; 

// 	Now, check the ratio between the two 
//	c->cd(2);
//	TLegend * legend2 = new TLegend(0.35, 0.65, 0.86, 0.82);
//	//legend2->AddEntry(h_ref, "Peak detect 1 x 1s / peak detect 1,000 x 1ms");
//	//h2->Scale(factor * h->GetEntries());
//	//g1->Draw("A*")
//	TH1D * h_test = (TH1D *) h->Clone();
//	TH1D * hlas_test = (TH1D *) h3->Clone();
//	legend2->AddEntry(h_test, "(PD 1s DR - N fitted DR) / N fitted DR");
//	legend2->AddEntry(hlas_test, "(N fitted Laser - N fitted DR) / N fitted DR");
//	h_test->Add(h2, -1);
//	h_test->Divide(h2); 
//	h_test->Draw(" ");
//	hlas_test->Add(h2, -1);
//	hlas_test->Divide(h2); 
//	h_test->GetYaxis()->SetTitle("test - refDR / ref DR");
//	hlas_test->Draw("same");
//	legend2->Draw("same");
//
//	new TFile("comparisonFile_laser.root", "RECREATE");
//	h3->Write();
//	h->Write();
//	h2->Write();
//	//g1->Write();
}

*/
