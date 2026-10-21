{
  TFile *f1 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/SquarePMT/combined_datasets/hadded_all_Square_TTS_testWithMonitor3_1pe_1ms_LandauCorrectedVolts_new.root");
  f1->ls();
  TH1D *h1f; f1->GetObject("hLAS1msC", h1f);
  h1f->SetLineColor(kMagenta);
  
  TFile *f2 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/LosangePMT/combined_datasets/hadded_all_Losange_TTS_testWithMonitor4_1pe_1ms_LandauCorrectedVolts_new.root");
  f2->ls();
  TH1D *h1s; f2->GetObject("hLAS1msC", h1s);
  
  TFile *f3 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/CirclePMT/combined_datasets/hadded_all_Circle_TTS_testWithMonitor4_1pe_1ms_LandauCorrectedVolts_new.root");
  f3->ls();
  TH1D *h1t; f3->GetObject("hLAS1msC", h1t);
  
  TFile *f4 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/TrianglePMT/combined_datasets/hadded_all_Triangle_TTS_testWithMonitor4_1pe_1ms_LandauCorrectedVolts_new.root");
  f4->ls();
  TH1D *h14; f4->GetObject("hLAS1msC", h14);
  
  TLegend * legend = new TLegend(0.35, 0.65, 0.86, 0.82);
  legend->AddEntry(h1f, "Square PMT");
  legend->AddEntry(h1s, "Losange PMT");
  legend->AddEntry(h1t, "Circle PMT");
  legend->AddEntry(h14, "Triangle PMT");

  

  h1f->SetTitle("Landau - corrected spectra");
  h1f->GetXaxis()->SetTitle("Volts");
  h1f->GetYaxis()->SetTitle("Occurencies");
  
  h1s->SetLineColor(9);
  h1t->SetLineColor(1);
  h14->SetLineColor(3);

  h1f->Draw("hist");
  h1s->Draw("hist same");
  h1t->Draw("hist same");
  h14->Draw("hist same");
  legend->Draw();
}
