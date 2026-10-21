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



//this is a code to landau fit 1 channel of PMT data
void LandauFit(TFile * f_LAS1msN, char* data_name, TH1D * h_LAS1msN, TH1D * h_LAS1msC){
	double pe_DR1sPD, pe_LAS1sPD, pe_LAS1msN, pe_DR1msN;
        double factor = 1.;
        double thresh = 0.002; //we need at least that number of PE/volts max to be considered a pulse
        double threshInit = 0.002; //this is the first point we consider to be part the pulse - landau
        int lenWindow = 200; //for every smaller window of 200 points we look for a pulse
        int smallWindow = 50; //when we have found a pulse(above threh.) we will fit it over 50 points 
        int nPointsPrior = 5; //We offset the fitting window by 5 points
        double onePEpeak_min = 0.005; // what points we fit in the one pe peak
        double onePEpeak_max = 0.03;

	TTree *t_LAS1msN = (TTree*)f_LAS1msN->Get(data_name);
        t_LAS1msN->SetBranchAddress("volts",&pe_LAS1msN); //can change to pe if we are sure
        Int_t nentries_LAS1msN = (Int_t)t_LAS1msN->GetEntries(); 
	//check how many points we have in tottal


	int n_largeWindows = 0;
        double maxPEinThisWindow = -9999; 
        double x[lenWindow], y[lenWindow]; //create an array x,y, for the given sub-window
        for (Int_t i=0; i<nentries_LAS1msN; i++) {
		//this reads in the ith entry of the root file and updates pe_LAS1msN with the value
		//of the voltage correcponsing to this point along the waveform
                t_LAS1msN->GetEntry(i);
		//store the uncorrected voltages for future reference
                h_LAS1msN->Fill(pe_LAS1msN);

                //Now the landau fit
                if (i-lenWindow*n_largeWindows != 0){
		       //store the maximum value in this window
                       if (pe_LAS1msN>=maxPEinThisWindow) maxPEinThisWindow = pe_LAS1msN;
                       //add points to the buffer of the given window size
                       int bufferPos = i%lenWindow;
                       x[bufferPos] = bufferPos;
                       y[bufferPos] = pe_LAS1msN;
                 }
		
		if (i-lenWindow*n_largeWindows == 0){
			//this is the end of the window
                        double point = -9999;
			//if the max is above threshold then it is a pulse
                        if (maxPEinThisWindow > thresh) {
                                int count = 0;
                                int a = 0; // just to prevent re-doing the same code
                                double x_buf[smallWindow], y_buf[smallWindow];

                                for (int j = 0; j<=lenWindow; j++){
					  //find j such that nPointsPrior after that we are above thres
                                          if(y[j+nPointsPrior] >= threshInit and a ==0){
                                                   x_buf[0] = 0;
                                                   y_buf[0] = y[j];
                                                   a = 1;
                                                   count = 1;
                                           }
					  //after we have found the begining of the window, add points
					  //to it until it is finished 
                                           else if (a == 1 and count <= smallWindow){
                                                   x_buf[count] = count;
                                                   y_buf[count] = y[j];
                                                   count += 1;
                                           }
                                           else if (count == smallWindow) break;
                                }

			//now that we have the zoomed in array of points around the peak, fit it with 
			//a Landau distribution
                        TGraph *g1 = new TGraph(smallWindow, x_buf, y_buf);
                        TFitResultPtr r = g1->Fit("landau", "S0"); //so that it is not drawn
                        Double_t par0   = r->Parameter(0);
                        Double_t par1   = r->Parameter(1);
                        Double_t par2   = r->Parameter(2);
                        TF1 * lan = new TF1("lan", "[0] * TMath::Landau(x, [1], [2], 0)", 0, smallWindow);
                        lan->SetParameters(par0, par1, par2);
                        double m = lan->GetMaximum();
			//then we take either the data point maximum over the 200 points or the maximum
			//of the landau fit, this limits the issue that come with finite sampling 
			//frequnecy 
                        point = (std::max(m, maxPEinThisWindow));
                        }
			//if we haven't found any pulse anyways then we save the max of the window
                        else point = maxPEinThisWindow;
                        //prepare for the name window
			n_largeWindows += 1;
                        maxPEinThisWindow = -9999.;
			//Store the point in the histogram
                        h_LAS1msC->Fill(point);
                }//end of it is the past point of the window
        }// end of readin the third file

}//end of the LandauFit function

int main(int argc, char **argv) {
	char * PMTform = NULL;
        char * face_position = NULL;
        char * path_to_input_file = NULL;
        char * path_to_output_file = "../Data/output_LandauFit.root";
        char option;
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
			default:
                                return 0;
                }
        }


	TCanvas *c = new TCanvas("c", "Example", 0, 0, 600, 800);
	TCanvas *c1 = new TCanvas("c1", "Example", 0, 0, 600, 800);

	//read in the data taken with the normal aquisition mode	
	TFile *f_LAS1msN = new TFile(Form("%s", path_to_input_file),"READ");
	//read the information of the first PMT branch, will need to loop over the 
	//other PMTs if more than one is used at any one time
   	//TTree *t_LAS1msN = (TTree*)f_LAS1msN->Get("PMT1data");
	
	//Create the histograms that we are going to fill with the data and with the landau corrected
	TH1D *h_LAS1N  = new TH1D("PMT1_Normal", "PMT1_Normal",120, -0.001, 0.2);
	TH1D *h_LAS1C  = new TH1D("PMT1_Corrected","PMT1_Corrected",120, -0.001, 0.2); 

	//do the conversion
	LandauFit(f_LAS1msN, "PMT1data", h_LAS1N, h_LAS1C);

        TFile* o = new TFile(Form("%s",path_to_output_file), "RECREATE");
	//write the two histograms to the output file
	h_LAS1N->Write();
	h_LAS1C->Write();
	o->Close();
	f_LAS1msN->Close();
		
}
