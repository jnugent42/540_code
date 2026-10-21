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
void LandauFit(TFile * f_LAS1msN, char* data_name, char* path_to_output_file){
	double pe_LAS1msN;
	int reps = 100; // Total 1ms Windows

	TTree *t_LAS1msN = (TTree*)f_LAS1msN->Get(data_name);
	t_LAS1msN->SetBranchAddress("volts",&pe_LAS1msN); //can change to pe if we are sure
	//Int_t nentries_LAS1msN = t_LAS1msN->GetEntries(); 
	Int_t nentries_LAS1msN = 2000000 * reps ; (Int_t)t_LAS1msN->GetEntries();//total = (samples in 1 window * number of windows) 
	//std::cout << "Total entries = " << nentries_LAS1msN << std::endl;
	//std::cout << "Total Laser Trigger = " << reps*100 << std::endl;

	double normal;
	double corrected;
	TTree tree(Form("%s", data_name),Form("%s_corrected",data_name));
	TTree treeNormal(Form("%s", data_name),Form("%s_normal",data_name));
	treeNormal.Branch("normal",&normal);
	tree.Branch("corrected",&corrected);

	int start_window = 0;
	std::cout << "Enter starting window to look at:" << std::endl;
	std::cin >> start_window;
	//std::cout<< "";

	int end_window = 200;
	std::cout << "Enter Final window to look at:" << std::endl;
	std::cin >> end_window;
	//std::cout<< "";

	int n = 2000;
	//int start_window = 50;
	//int end_window = 90;
	int offset = 0;
	double vol[n];
	double pos[n];

	for(int i=start_window; i<=end_window; i++){
		for(int j=i*n; j<=(i+1)*n; j++){
			t_LAS1msN->GetEntry(j+offset);
			vol[j%n] = pe_LAS1msN * 1000;
			pos[j%n] = j%n * 0.5;
		}
		TGraph *g_wave = new TGraph(n, pos, vol);
		g_wave->SetTitle("; Time after trigger (ns); mV");
		g_wave->Write(Form("%s_Graph_Window_%d",data_name,i));

	}
	
}

int main(int argc, char **argv) {
	char * PMTform = NULL;
	char * face_position = NULL;
	char * path_to_input_file = NULL;
	char * path_to_output_file = "../data/output_LandauFit.root";
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
	TFile* o = new TFile(Form("%s",path_to_output_file), "RECREATE");
	
	//do the conversion

	LandauFit(f_LAS1msN, "PMT1data", path_to_output_file);
	LandauFit(f_LAS1msN, "PMT2data", path_to_output_file);
	LandauFit(f_LAS1msN, "PMT3data", path_to_output_file);
	//LandauFit(f_LAS1msN, "PMT4data", path_to_output_file);

	o->Close();
	f_LAS1msN->Close();

}
