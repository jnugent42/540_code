{
  TFile *f1 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/SquarePMT/combined_datasets/hadded_all_Square_TTS_120mA_1pe_1ms_LandauCorrectedVolts.root");
  f1->ls();
  TH1D *h1f; f1->GetObject("hLAS1msC", h1f);
  h1f->SetLineColor(1);
  
  TFile *f2 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/SquarePMT/combined_datasets/hadded_all_Square_TTS_Test130mA_1pe_1ms_LandauCorrectedVolts.root");
  f2->ls();
  TH1D *h1s; f2->GetObject("hLAS1msC", h1s);
  
  TFile *f3 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/SquarePMT/combined_datasets/hadded_all_Square_TTS_Test114mA_1pe_1ms_LandauCorrectedVolts.root");
  f3->ls();
  TH1D *h1t; f3->GetObject("hLAS1msC", h1t);
  
  TFile *f4 = TFile::Open("/home/a/WCTE540/PMT_data/TTSLB/SquarePMT/combined_datasets/hadded_all_Square_TTS_180degLOSCAB_1pe_1ms_LandauCorrectedVolts.root");
  f4->ls();
  TH1D *h14; f4->GetObject("hLAS1msC", h14);
  
  TFile *f5 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/SquarePMT/combined_datasets/hadded_all_Square_TTS_centre_1pe_1ms_LandauCorrectedVolts.root");
  f5->ls();
  TH1D *h15; f5->GetObject("hLAS1msC", h15);
  
  TFile *f6 = TFile::Open("/home/a/WCTE540/PMT_data/TTSLB/SquarePMT/combined_datasets/hadded_all_Square_TTS_180deg150mA_1pe_1ms_LandauCorrectedVolts.root");
  f6->ls();
  TH1D *h16; f6->GetObject("hLAS1msC", h16);
  
  TLegend * legend = new TLegend(0.35, 0.65, 0.86, 0.82);
  legend->AddEntry(h1f, "Square PMT - Fibre + Att. 134mA");
  legend->AddEntry(h1s, "Square PMT - Fibre + Att. 130mA");
  legend->AddEntry(h1t, "Square PMT - Fibre + Att. 114mA");
  //legend->AddEntry(h16, "Square PMT - LB 150mA");
  //legend->AddEntry(h14, "Square PMT - LB 154mA");
  legend->AddEntry(h15, "Square PMT - Fibre + Att. 124mA (initial)");

  h1s->SetTitle("Landau - corrected spectra");
  h1f->GetXaxis()->SetTitle("Volts");
  h1f->GetYaxis()->SetTitle("Occurencies");
  h1s->SetLineColor(4);
  h1t->SetLineColor(3);
  h14->SetLineColor(3);
  h15->SetLineColor(6);
  h16->SetLineColor(14);
  h1f->Draw("hist");
  h1s->Draw("hist same");
  h1t->Draw("hist same");
  //h14->Draw("hist same");
  h15->Draw("hist same");
  //h16->Draw("hist same");
  legend->Draw();
}
