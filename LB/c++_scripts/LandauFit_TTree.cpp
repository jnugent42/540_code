#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TText.h"
#include "TStyle.h"


// Landau function with 2 different constants (left and right of peak)
double Landau2C(double *x, double *par){
	double N = par[0];
	double landau = TMath::Landau(x[0], par[1], par[2], kFALSE);
	double left = par[3];
	double right = par[4];

	if(x[0] < par[1]){
		return N * landau + left;
	} else{
		return N * landau + right;
	}
}


//this is a code to landau fit 1 channel of PMT data
void LandauFit(TFile * f_LAS1msN, const char* data_name, TH1D * h_LAS1msN, TH1D * h_LAS1msC, TH1D * h_LAS1msQ, TH1D * h_chi2, const char* path_to_output_file){
	double pe_DR1sPD, pe_LAS1sPD, pe_LAS1msN, pe_DR1msN;
	double thresh = 0.009; //we need at least that number of PE/volts max to be considered a pulse
	double threshInit = thresh; //this is the first point we consider to be part the pulse - landau
	int lenWindow = 2000; //for every smaller window of xx points we look for a pulse
	int smallWindow = 100; //when we have found a pulse(above threh.) we will fit it over xx points 
	int nPointsPrior = 40; //We offset the fitting window by xx points
	int Fitted = 0; //Total fitted waveforms
	int NotFitted = 0; // Potentially bad fit
	int DarkPulse = 0; // Potentially dark noise
	int reps = 100; // Total 1ms Windows

	double point = -9999;

	std::cout << "For " << data_name << " : " << std::endl;

	//check how many points we have in tottal
	TTree *t_LAS1msN = (TTree*)f_LAS1msN->Get(data_name);
	t_LAS1msN->SetBranchAddress("volts",&pe_LAS1msN); //can change to pe if we are sure
	Int_t nentries_LAS1msN = t_LAS1msN->GetEntries(); 
	std::cout << "Total entries = " << nentries_LAS1msN << std::endl;
	std::cout << "Total Laser Trigger = " << nentries_LAS1msN/lenWindow << std::endl; // 2Gsample/s / LaserFreq

	//Int_t nentries_LAS1msN = 2000000 * reps ; (Int_t)t_LAS1msN->GetEntries();//total = (samples in 1 window * number of windows) 
	//std::cout << "Total entries = " << nentries_LAS1msN << std::endl;
	//std::cout << "Total Laser Trigger = " << reps*100 << std::endl;

	double normal;
	double corrected;
	TTree tree(Form("%s", data_name),Form("%s_corrected",data_name));
	TTree treeNormal(Form("%s", data_name),Form("%s_normal",data_name));
	treeNormal.Branch("normal",&normal);
	tree.Branch("corrected",&corrected);

	//Prepare for looping	
	int n_largeWindows = 0;
	double maxPEinThisWindow = -9999; 
	int maxPos = 0;
	double x[lenWindow], y[lenWindow]; //create an array x,y, for the given sub-window
	for(int k = 0; k < lenWindow; k++){
		x[k] = 0.0;
		y[k] = 0.0;
	}
	for (Int_t i=0; i<nentries_LAS1msN; i++) {
		//this reads in the ith entry of the root file and updates pe_LAS1msN with the value
		//of the voltage correcponsing to this point along the waveform
		t_LAS1msN->GetEntry(i);

		//store the uncorrected voltages for future reference
		h_LAS1msN->Fill(pe_LAS1msN);
		normal = pe_LAS1msN;

		//Now the landau fit
		if (i-lenWindow*n_largeWindows != 0){

			//store the maximum value and its position in this window
			if (pe_LAS1msN >= maxPEinThisWindow){ 
				maxPEinThisWindow = pe_LAS1msN;
				maxPos = (i%lenWindow)*0.5; // time after trigger
			}
			//add points to the buffer of the given window size
			int bufferPos = i%lenWindow;
			x[bufferPos] = bufferPos;
			y[bufferPos] = pe_LAS1msN;
		}
		else if (i-lenWindow*n_largeWindows == 0){
			//this is the end of the window
			point = 0; //-9999;

			//if the max is above threshold then it is a pulse
			if (maxPEinThisWindow > thresh) {
				int count = 0;
				int a = 0; // just to prevent re-doing the same code

				if(maxPos < 30 || maxPos>100){
					DarkPulse += 1; //outside of TTS range for PMT
				}

				double x_buf[smallWindow], y_buf[smallWindow];
				for (int g = 0; g<smallWindow; g++){
					x_buf[g] = 0.0;
					y_buf[g] = 0.0;
				}

				for (int j = 1; j<lenWindow; j++){
					//find j such that nPointsPrior after that we are above thres
					//if(y[j+nPointsPrior] >= threshInit && a == 0 && (j+nPointsPrior) < lenWindow){
					if(y[j+nPointsPrior] == maxPEinThisWindow && a == 0 && (j+nPointsPrior) < lenWindow){
						//x_buf[0] = 0;
						x_buf[0] = j*0.5; // to get time after trigger
						y_buf[0] = y[j];
						a = 1;
						count = 1;
					}

					//after we have found the begining of the window, add points
					//to it until it is finished 
					else if (a == 1 && count < smallWindow){
						//x_buf[count] = count;
						x_buf[count] = j*0.5; // to get time after trigger
						y_buf[count] = y[j];
						count += 1;
					}
					
					else if (count == smallWindow) break;
				}
				
				if(a == 1){
					//now that we have the zoomed in array of points around the peak, fit it with 
					//a Landau distribution
					TGraph *g1 = new TGraph(smallWindow, x_buf, y_buf);
					g1->SetTitle(";Time after trigger (ns); Volts");
					
					/*
					TF1 * lan = new TF1("lan", "[0] * TMath::Landau(x, [1], [2], 0) + [3]", 0, smallWindow);
					lan->SetParameters(0.01, maxPos, 0.5, 0);
					lan->SetParLimits(0, 0, 1);
					lan->SetParLimits(1, maxPos-15, maxPos+15);
					lan->SetParLimits(2, 0, 1);
					lan->SetParLimits(3, 0, 0.005);
					*/

					TF1 *lan = new TF1("lan", Landau2C, x_buf[0], x_buf[99], 5);
					lan->SetParameters(0.01, maxPos, 0.5, 0, 0);
					lan->SetParLimits(0, 0, 1);
					lan->SetParLimits(1, maxPos-15, maxPos+15);
					lan->SetParLimits(2, 0, 1);
					//lan->SetParLimits(3, -0.001, 0.005);
					//lan->SetParLimits(4, -0.001, 0.005);
					
					g1->Fit(lan, "QS0"); //so that it is not drawn
					double m = lan->GetMaximum();
					double area = lan->Integral(maxPos-5, maxPos+5, 1e-9); // assume all pulses have similar width
					h_LAS1msQ->Fill(area);

					//then we take either the data point maximum over the xx points or the maximum
					//of the landau fit, this limits the issue that come with finite sampling frequnecy 
					point = (std::max(m, maxPEinThisWindow));
					Fitted += 1;

					//Write the Waveforms
					//g1->Write(Form("%s_Graph_Window_%d",data_name,n_largeWindows-1));
					//lan->Write(Form("%s_LanFit_Window_%d",data_name,n_largeWindows-1));

					//See goodness of Fit
					double chi2_val = lan-> GetChisquare();
					h_chi2->Fill(chi2_val);	
					if(chi2_val > 0.00050 || m >= 1.2*maxPEinThisWindow || m <= 0.8*maxPEinThisWindow){
						NotFitted +=1;
						//g1->Write(Form("%s_Graph_Window_%d",data_name,n_largeWindows-1));
						//lan->Write(Form("%s_LanFit_Window_%d",data_name,n_largeWindows-1));
					}
				}
			}
			//if we haven't found any pulse anyways then we save the max of the window
			else{
			       	point = maxPEinThisWindow;
			}

			//prepare for the next window
			n_largeWindows += 1;
			corrected = point;
			maxPEinThisWindow = -9999.;

			//Store the point in the histogram
			h_LAS1msC->Fill(point);

		}//end of it is the past point of the window
	}// end of readin the third file
	std::cout << "Total fitted waveforms = " << Fitted << std::endl;
	std::cout << "Poorly fitted waveforms = " << NotFitted << std::endl;
	std::cout << "Potential Dark Noise = " << DarkPulse << std::endl;

}//end of the LandauFit function

int main(int argc, char **argv) {
	char * PMTform = NULL;
	char * face_position = NULL;
	char * path_to_input_file = NULL;
	const char * path_to_output_file = "../data/output_LandauFit.root";
	char option;
	int triggerChannel = 0;
	//the following piece of code reads in some user input
	while ((option = getopt(argc, argv, "f:p:r:t:b:o:e:i:h")) != -1 )
	{
		switch(option)
		{
			case 'f':
				PMTform = optarg;
				break;
			case 'p':
				face_position = optarg;
				break;
			case 'i':
				path_to_input_file = optarg;
				break;
			case 'o':
				path_to_output_file = optarg;
				break;
			case 't':
				triggerChannel = std::stoi(optarg);
				break;
			default:
				return 0;
		}
	}


	TCanvas *c = new TCanvas("c", "Example", 0, 0, 600, 800);
	TCanvas *c1 = new TCanvas("c1", "Example", 0, 0, 600, 800);

	//read in the data taken with the normal aquisition mode	
	TFile *f_LAS1msN = new TFile(Form("%s", path_to_input_file),"READ");

	//Create the histograms that we are going to fill with the data and with the landau corrected
	TH1D *h_LAS1N  = new TH1D("PMT1_Normal", "PMT1_Normal", 100, 0.00, 0.08);
	TH1D *h_LAS1C  = new TH1D("PMT1_Corrected","PMT1_Corrected", 100, 0.00, 0.08);
	TH1D *h_LAS1Q  = new TH1D("PMT1_Charge","PMT1_Charge", 200, 0, 1);
	TH1D *h_LAS1Chi2  = new TH1D("PMT1_Chi2","PMT1_Chi2", 100, -0.00005, 0.0006);

	TH1D *h_LAS2N  = new TH1D("PMT2_Normal", "PMT2_Normal", 120, 0.0, 0.08);
	TH1D *h_LAS2C  = new TH1D("PMT2_Corrected","PMT2_Corrected", 100, 0.00, 0.08);
	TH1D *h_LAS2Q  = new TH1D("PMT2_Charge","PMT2_Charge", 200, 0, 1);
	TH1D *h_LAS2Chi2  = new TH1D("PMT2_Chi2","PMT2_Chi2", 100, -0.00005, 0.0006);

	TH1D *h_LAS3N  = new TH1D("PMT3_Normal", "PMT3_Normal", 500, 0, 0.5);
	TH1D *h_LAS3C  = new TH1D("PMT3_Corrected","PMT3_Corrected", 500, 0, 0.5);
	TH1D *h_LAS3Q  = new TH1D("PMT3_Charge","PMT3_Charge", 400, 0, 2);
	TH1D *h_LAS3Chi2  = new TH1D("PMT3_Chi2","PMT3_Chi2", 100, -0.00005, 0.0006);

	TH1D *h_LAS4N  = new TH1D("PMT4_Normal", "PMT4_Normal", 100, 0.00, 0.08);
	TH1D *h_LAS4C  = new TH1D("PMT4_Corrected","PMT4_Corrected", 100, 0.00, 0.08);
	TH1D *h_LAS4Q  = new TH1D("PMT4_Charge","PMT4_Charge", 200, 0, 2);
	TH1D *h_LAS4Chi2  = new TH1D("PMT4_Chi2","PMT4_Chi2", 100, -0.00005, 0.0006);

	TFile* o = new TFile(Form("%s",path_to_output_file), "RECREATE");
	//do the conversion

	LandauFit(f_LAS1msN, "PMT1data", h_LAS1N, h_LAS1C, h_LAS1Q, h_LAS1Chi2, path_to_output_file);
	LandauFit(f_LAS1msN, "PMT2data", h_LAS2N, h_LAS2C, h_LAS2Q, h_LAS2Chi2, path_to_output_file);
	LandauFit(f_LAS1msN, "PMT3data", h_LAS3N, h_LAS3C, h_LAS3Q, h_LAS3Chi2, path_to_output_file);
	//LandauFit(f_LAS1msN, "PMT4data", h_LAS4N, h_LAS4C, h_LAS4Q, h_LAS4Chi2, path_to_output_file);

	TF1 * gaus = new TF1("gaus", "gaus", 0, 0.1);

	//write the two histograms to the output file
	h_LAS1N->Write();
	h_LAS1C->Write();
	h_LAS1Q->Write();
	h_LAS1Chi2->Write();

	h_LAS2N->Write();
	h_LAS2C->Write();
	h_LAS2Q->Write();
	h_LAS2Chi2->Write();

	h_LAS3Q->Fit("gaus");
	h_LAS3Q->Draw();
	gaus->Draw("same");
	h_LAS3N->Write();
	h_LAS3C->Write();
	h_LAS3Q->Write();
	h_LAS3Chi2->Write();

	h_LAS4N->Write();
	h_LAS4C->Write();
	h_LAS4Q->Write();
	h_LAS4Chi2->Write();

	o->Close();
	f_LAS1msN->Close();

}
