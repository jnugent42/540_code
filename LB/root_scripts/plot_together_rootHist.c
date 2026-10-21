{
  //TFile *f1 = TFile::Open("/home/a/WCTE540/PMT_data/TTSLB/LosangePMT/combined_datasets/hadded_all_Losange_TTS_180deg-try3_1pe_1ms_LandauCorrectedVolts.root");
  //TFile *f1 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/LosangePMT/combined_datasets/hadded_all_Losange_TTS_Try120mA_1pe_1ms_LandauCorrectedVolts.root");
  TFile *f1 = TFile::Open("/home/a/WCTE540/PMT_data/TTSLB/LosangePMT/combined_datasets/hadded_all_Losange_TTS_360deg154mA_1pe_1ms_LandauCorrectedVolts.root");
  f1->ls();
  TH1D *h1f; f1->GetObject("hLAS1msC", h1f);
  h1f->SetLineColor(kMagenta);
  
  TFile *f2 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/LosangePMT/hadded_all_Losange_TTS_testWithMonitor6_1pe_1ms_LandauCorrectedVolts_new.root");
  f2->ls();
  TH1D *h1s; f2->GetObject("hLAS1msC", h1s);
  
  TFile *f3 = TFile::Open("/home/a/WCTE540/PMT_data/TTSLB/LosangePMT/combined_datasets/hadded_all_Losange_TTS_270deg154mA_1pe_1ms_LandauCorrectedVolts.root");
  f3->ls();
  TH1D *h1t; f3->GetObject("hLAS1msC", h1t);
  
  TFile *f4 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/LosangePMT/combined_datasets/hadded_all_Losange_TTS_left_1pe_1ms_LandauCorrectedVolts_new.root");
  //TFile *f4 = TFile::Open("/home/a/WCTE540/PMT_data/TTSLB/LosangePMT/combined_datasets/hadded_all_Losange_TTS_180deg-try6_1pe_1ms_LandauCorrectedVolts.root");
  f4->ls();
  TH1D *h14; f4->GetObject("hLAS1msC", h14);
  
  TFile *f5 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/LosangePMT/combined_datasets/hadded_all_Losange_TTS_centre_1pe_1ms_LandauCorrectedVolts_new.root");
  f5->ls();
  TH1D *h15; f5->GetObject("hLAS1msC", h15);
  
  TFile *f6 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/LosangePMT/combined_datasets/hadded_all_Losange_TTS_right_1pe_1ms_LandauCorrectedVolts_new.root");
 //TFile *f6 = TFile::Open("/home/a/WCTE540/PMT_data/TTSLB/LosangePMT/combined_datasets/hadded_all_Losange_TTS_180deg150mA_1pe_1ms_LandauCorrectedVolts.root");
  f6->ls();
  TH1D *h16; f6->GetObject("hLAS1msC", h16);

  TFile *f7 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/LosangePMT/combined_datasets/hadded_all_Losange_TTS_top_1pe_1ms_LandauCorrectedVolts_new.root");
  //TFile *f7 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/LosangePMT/combined_datasets/hadded_all_Losange_TTS_Try120mA_1pe_1ms_LandauCorrectedVolts.root");
  f7->ls();
  TH1D *h17; f7->GetObject("hLAS1msC", h17);
  
  TFile *f8 = TFile::Open("/home/a/WCTE540/PMT_data/TTSLB/LosangePMT/combined_datasets/hadded_all_Losange_TTS_180deg154mA_1pe_1ms_LandauCorrectedVolts.root");
  f8->ls();
  TH1D *h18; f8->GetObject("hLAS1msC", h18);
  

  TFile *f9 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/LosangePMT/combined_datasets/hadded_all_Losange_TTS_bottom_1pe_1ms_LandauCorrectedVolts_new.root");
  //TFile *f9 = TFile::Open("/home/a/WCTE540/PMT_data/TTS/LosangePMT/combined_datasets/hadded_all_Losange_TTS_Try123mA_1pe_1ms_LandauCorrectedVolts.root");
  f9->ls();
  TH1D *h19; f9->GetObject("hLAS1msC", h19);
  
  TLegend * legend = new TLegend(0.35, 0.65, 0.86, 0.82);
  legend->AddEntry(h1f, "Losange PMT - LB 154mA 360deg");
  legend->AddEntry(h1s, "Losange PMT - LB 154mA 90deg");
  legend->AddEntry(h1t, "Losange PMT - LB 154mA 270deg");
  legend->AddEntry(h18, "Losange PMT - LB 154mA 180deg");
  //legend->AddEntry(h14, "Losange PMT - LB 154mA");
  //legend->AddEntry(h15, "Losange PMT - Fibre + Att. 124mA (initial)");
  //legend->AddEntry(h1f, "Losange PMT - Fibre + Att. 120mA");
  //legend->AddEntry(h18, "Losange PMT - Fibre + Att. 128mA");
  //legend->AddEntry(h19, "Losange PMT - Fibre + Att. 123mA");

  
  legend->AddEntry(h14, "Losange PMT - Fibre + Att. 124mA (left)");
  legend->AddEntry(h15, "Losange PMT - Fibre + Att. 124mA (centre)");
  legend->AddEntry(h16, "Losange PMT - Fibre + Att. 124mA (right)");
  legend->AddEntry(h17, "Losange PMT - Fibre + Att. 124mA (top)");
  legend->AddEntry(h19, "Losange PMT - Fibre + Att. 124mA (bottom)");

  h1f->SetTitle("Landau - corrected spectra");
  h1f->GetXaxis()->SetTitle("Volts");
  h1f->GetYaxis()->SetTitle("Occurencies");
  
  h1s->SetLineColor(9);
  h1t->SetLineColor(1);
  h14->SetLineColor(3);
  //h17->SetLineStyle(2);
  h15->SetLineColor(4);
  //h15->SetLineStyle(2);
  h16->SetLineColor(14);
  h17->SetLineColor(1);
  h18->SetLineColor(7);
  //h18->SetLineStyle(2);
  h19->SetLineColor(2);
  //h19->SetLineStyle(2);
  
  //h1f->Scale(1.0/h1f->GetEntries());
  //h1s->Scale(1.0/h1s->GetEntries());
  //h1t->Scale(1.0/h1t->GetEntries());
  //h14->Scale(1.0/h14->GetEntries());
  //h15->Scale(1.0/h15->GetEntries());
  //h16->Scale(1.0/h16->GetEntries());
  //h17->Scale(1.0/h17->GetEntries());
  //h18->Scale(1.0/h18->GetEntries());
  //h19->Scale(1.0/h19->GetEntries());

  h1f->Draw("hist");
  h1s->Draw("hist same");
  h1t->Draw("hist same");
  h15->Draw("hist same");
  h14->Draw("hist same");
  h16->Draw("hist same");
  h17->Draw("hist same");
  h18->Draw("hist same");
  h19->Draw("hist same");
  legend->Draw();
}
