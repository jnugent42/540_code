{
	TCanvas *c = new TCanvas("c", "Example", 0, 0, 600, 800);
	TCanvas *c1 = new TCanvas("c1", "Example", 0, 0, 600, 800);
	
	TFile *f_input = new TFile("/home/a/WCTE540/PMT_data/Aug31/DR_FT_LightsON/aug31_DR_FT_lightsON.root","READ");

	//TH1D *pmt_1 = (TH1D*)f_input->Get("PMT1_Corrected");
	TH1D *pmt_1 = (TH1D*)f_input->Get("PMT1data");

	//TTree *t_LAS1msN = (TTree*)f_LAS1msN->Get("PMT1data");
	
	TH1D *h_LAS1C  = new TH1D("PMT1_Corrected","PMT1_Corrected",120, -0.001, 0.059); //0.124, 200 bins


	
        TFile* o = new TFile("", "RECREATE");
	h_LAS1C->Write();
	o->Close();

	VoltageFit(h_LAS1C);
		
}



void VoltageFit(TH1D * h_LAS1msC) {
  double total = 0;
  // gets the center value of each of the bins, then multiplies
  //by number of counts to get the total voltage 
  for(int x = 0; x < h_LAS1msC->GetNbinsX(); x++) {
    for(int y = 0; y < h_LAS1msC->GetNbinsX(); y++) {
      double newTotal =  h_LAS1msC->GetXaxis()->GetBinCenter(x);// * h_LAS1msC->GetBinContent(x, y);
      total = total + newTotal;
      cout << total << endl;
      cout <<  h_LAS1msC->GetBinContent(x) << endl;
    }
  }
  // outside loop
  cout << "final total: " << total << endl;
  
}
