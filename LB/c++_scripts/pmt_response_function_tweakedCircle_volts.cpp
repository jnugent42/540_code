/// PMT Response Function
/// Based on:
/// E.H. Bellamy, et. al., NIM A 339 (1994) 468.
/// and "Fitting Single Photo-electron peak," Qing He, Princeton, 2010

#include "pmt_response_function.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TText.h"

const double pi = acos(-1.0);



/// factorial using TMath::Factorial
double factorial( unsigned n ){
  return TMath::Factorial( n );
}


/// Poisson term
double poisson_term( double mu, unsigned n ){
  return std::pow( mu, n ) * std::exp( -mu ) / factorial(n);
}

/// Gaussian part of signal for n'th photoelectron
/// This gaussian is defined only for x>0 and normalized accordingly
/// such that
///  integra( 0, infty, gaussian ) = 1
double gaussian( double x, double Qn, double sigman ){
  //double norm = 1.0 / ( sigman * std::sqrt( 2*pi ) );
  double norm = sqrt(2/pi) / ( sigman * (1 + std::erf( Qn / ( sqrt(2)*sigman ) ) ) );

  return norm * std::exp( -std::pow( (x-Qn)/(sqrt(2)*sigman), 2 ) );
}

double sign( double x ){
  if (x > 0.) return 1.;
  if (x < 0.) return -1.;
  return 0.;
}

/// PMT response Sreal(x) from
/// E.H. Bellamy, et. al., NIM A 339 (1994) 468.
///
/// *** Exclude Pedestal term. ***
/// x[0] is charge
/// p[0] is N, overall normalization
/// p[1] is Q1, charge of single photoelectron
/// p[2] is sigma1, width of single photoelectron gaussian
/// p[3] is mu, the poisson mean number of photoelectrons
/// p[4] is w, the fraction of the background signal that is exponential
/// p[5] is alpha, the exponential constant
double pmtresponse( double * xx, double * p ){
  double x     = xx[0];
  double N     = p[0];
  double Q1    = p[1];
  double s1    = p[2];
  double mu    = p[3];
  double w     = p[4];
  double alpha = p[5];

  /// how many terms to include?
  /// go to 3 sigma above mean
  double approx_sigma = std::sqrt( mu );
  unsigned nmax = std::max(2, int(approx_sigma*3+mu) );

  double signal = 0;
  for (unsigned n=1; n<=nmax; ++n){
    double Qn = n * Q1;
    double sigman = std::sqrt(n) * s1;
    double poisson = poisson_term( mu, n );
    double exp_term = alpha / 2 * std::exp( -alpha * ( x - Qn - alpha*sigman*sigman/2 ) );
    double gaus_term = gaussian( x, Qn, sigman );
    double erf1 = std::erf( fabs(n*Q1 + sigman*sigman*alpha) / ( sigman*sqrt(2.) ) );
    double erfarg = x - Qn - sigman*sigman*alpha;
    double erf2 = sign(erfarg)*std::erf( fabs(erfarg) / ( sigman*sqrt(2.) ) );
    double IGn_term = exp_term*( erf1 + erf2 );
    signal += poisson*( (1-w)*gaus_term + w*IGn_term );
  }
  return N*signal;
}

///  Background response
/// Just the background portionof the function
double pmtbackgroundresponse( double * xx, double *p ){
  double x     = xx[0];
  double N     = p[0];
  double Q1    = p[1];
  double s1    = p[2];
  double mu    = p[3];
  double w     = p[4];
  double alpha = p[5];

  /// how many terms to include?
  /// go to 3 sigma above mean
  double approx_sigma = std::sqrt( mu );
  //unsigned nmax = std::max(2, int(approx_sigma*3+mu) );
  unsigned nmax = 2;

  double signal = 0;
  for (unsigned n=1; n<=nmax; n++){
    double Qn = n * Q1;
    double sigman = std::sqrt(n) * s1;
    double poisson = poisson_term( mu, n );
    double exp_term = alpha / 2 * std::exp( -alpha * ( x - Qn - alpha*sigman*sigman/2 ) );
    double gaus_term = gaussian( x, Qn, sigman );
    double erf1 = std::erf( fabs(n*Q1 + sigman*sigman*alpha) / ( sigman*sqrt(2.) ) );
    double erfarg = x - Qn - sigman*sigman*alpha;
    double erf2 = sign(erfarg)*std::erf( fabs(erfarg) / ( sigman*sqrt(2.) ) );
    double IGn_term = exp_term*( erf1 + erf2 );
    signal += poisson* w * IGn_term ;
  }
  return N*signal;
}

/// N'th photoelectron ideal response
/// make seventh parameter N (the photoelectron number)
double pmtG_n( double * xx, double *p ){
  double x     = xx[0];
  double N     = p[0];
  double Q1    = p[1];
  double s1    = p[2];
  double mu    = p[3];
  double w     = p[4];
  double alpha = p[5];
  double n     = p[6];

  double signal = 0;

  double Qn = n * Q1;
  double sigman = std::sqrt(n) * s1;
  double poisson = poisson_term( mu, n );
  double exp_term = alpha / 2 * std::exp( -alpha * ( x - Qn - alpha*sigman*sigman/2 ) );
  double gaus_term = gaussian( x, Qn, sigman );
  signal += poisson*( (1-w)*gaus_term );

  return N*signal;
}

/// Return a vector of functions that sums to be approxametely the Bellamy function
/// First term is the background, then 1 pe, 2 pe, ... etc
std::vector< TF1* > get_pmt_response_components( double * p ){
  double pars[ 7 ];
  for ( unsigned i = 0; i < 6; ++i){
    pars[i] = p[i];
  }
  pars[6]=0.;

  std::vector< TF1* > fcns;
  TF1* bg = new TF1("Background", pmtbackgroundresponse, 0., 5000., 6 );
  bg->SetParameters( pars );
  //bg->SetLineColor( 52 );
  fcns.push_back( bg );

  double mu=p[3];
  double approx_sigma = std::sqrt( mu );
  unsigned nmax = std::max(5, int(approx_sigma*3+mu) );

  for (unsigned n=1; n<=nmax; ++n ){
    std::ostringstream fname;
    fname <<"fnpe_"<<n;
    TF1* ftmp = new TF1( fname.str().c_str(), pmtG_n, 0., 5000., 7 );
    pars[6] = n;
    ftmp->SetNpx(1000);
    ftmp->SetLineWidth(3);
    ftmp->SetParameters( pars );
    ftmp->SetLineColor( 52 + 5*n );
    fcns.push_back( ftmp );
  }
  return fcns;
}



/////****************************************************************
///// Repeat all of above functions but including pedestal
/////****************************************************************
///
/// PMT response Sreal(x) from
/// E.H. Bellamy, et. al., NIM A 339 (1994) 468.
///
/// *** INclude Pedestal term. ***
/// x[0] is charge
/// p[0] is N, overall normalization
/// p[1] is Q0, charge of pedestal
/// p[2] is sigma0, width of pedestal
/// p[3] is Q1, charge of single photoelectron
/// p[4] is sigma1, width of single photoelectron gaussian
/// p[5] is mu, the poisson mean number of photoelectrons
/// p[6] is w, the fraction of the background signal that is exponential
/// p[7] is alpha, the exponential constant
double pmtresponseped( double * xx, double * p ){
  double x     = xx[0];
  double N     = p[0];
  double Q0    = p[1];
  double s0    = p[2];
  double Q1    = p[3];
  double s1    = p[4];
  double mu    = p[5];
  double w     = p[6];
  double alpha = p[7];
  double A2    = p[8];
  double A3    = p[9];

  /// how many terms to include?
  /// go to 3 sigma above mean
  double approx_sigma = std::sqrt( mu );
  //unsigned nmax = std::max(5, int(approx_sigma*3+mu) ); //number of pe? //change 1 for 3
  unsigned nmax = 2;

  double signal = 0;
  //A2 = 1.;
  //A3 = 1.;

  // Use approx equation from eq. 10 of Bellamy
  double exp0term = w * alpha * std::exp( -alpha * ( x - Q0 ) );
  if ( x < Q0 ) exp0term = 0.;
  double gaus0term = (1-w) * gaussian( x, Q0, s0 );
  signal = poisson_term( mu, 0 ) * ( gaus0term + exp0term );


  for (unsigned n=1; n<=nmax; n++){
    double gaus_term = gaussian( x, Q0 + n*Q1 + w/alpha, std::sqrt(n) * s1 );
    if (n==2) signal += A2 * poisson_term(mu,n) * gaus_term ;
    //if (n==3) signal += TMath::Abs(A3)* poisson_term(mu,n) * gaus_term ;

    else signal += 1 * poisson_term(mu,n) * gaus_term ;
    //signal += poisson_term(mu,n) * gaus_term ;
  }
  return N*signal;
}

///  Background response *** with pedesal ***
/// Just the background portion of the function, excludes pedestal term
double pmtbackgroundresponseped( double * xx, double *p ){
  double x     = xx[0];
  double N     = p[0];
  double Q0    = p[1];
  double s0    = p[2];
  double Q1    = p[3];
  double s1    = p[4];
  double mu    = p[5];
  double w     = p[6];
  double alpha = p[7];

  double signal = 0;

  // Use approx equation from eq. 10 of Bellamy
  double exp0term = w * alpha * std::exp( -alpha * ( x - Q0 ) );
  if ( x < Q0 ) exp0term = 0.;
  double gaus0term = (1-w) * gaussian( x, Q0, s0 );
  signal = poisson_term( mu, 0 ) * ( gaus0term + exp0term );

  return N*signal;
}


/// N'th photoelectron ideal response
/// make ninth parameter n (the photoelectron number)
/// 0 photoelectron term is 1pe
double pmtG_n_ped( double * xx, double *p ){
  double x     = xx[0];
  double N     = p[0];
  double Q0    = p[1];
  double s0    = p[2];
  double Q1    = p[3];
  double s1    = p[4];
  double mu    = p[5];
  double w     = p[6];
  double alpha = p[7];
  unsigned n   = unsigned(p[10]);
  double A2    = p[8];
  double A3    = p[9];  

  double signal = 0;

  //A2 = 1; 
  //A3 = 1;
  double gaus_term = gaussian( x, Q0 + n*Q1 + w/alpha, std::sqrt(n) * s1 );
  //signal += poisson_term(mu,n) * gaus_term ;
  if (n==1) signal += poisson_term(mu,n) * gaus_term ;
  //else if (n==3) {
	//  signal += TMath::Abs(A3) * poisson_term(mu,n) * gaus_term ;
	  //std::cout << "Position of the 3 pe peak :" << Q0 + n*Q1 + w/alpha << std::endl;
  //}
  else if (n==2){
	  signal += A2 * poisson_term(mu,n) * gaus_term ;
	  //std::cout << "Position of the " << n << " pe peak " << Q0 + n*Q1 + w/alpha << " n: " << n << std::endl;
  }

  return N*signal;
}

/// Return a vector of functions that sums to be approxametely the Bellamy function
/// First term is the background, then 0 pe, 1 pe, 2 pe, ... etc
std::vector< TF1* > get_pmt_response_components_ped( double * p ){
  double pars[ 11 ];
  double peToVolts = 0.01;
  for ( unsigned i = 0; i < 10; ++i){
    pars[i] = p[i];
  }
  pars[10]=0.;

  std::vector< TF1* > fcns;
  TF1* bg = new TF1("Background", pmtbackgroundresponseped, -0.1 * peToVolts, 6 * peToVolts, 10 );
  bg->SetParameters( pars );
  //bg->SetLineColor( 52 );
  std::cout << "Integral BG " << " " << bg->Integral(-0.1 * peToVolts, 6 * peToVolts) * pars[0] << std::endl;
  fcns.push_back( bg );

  double mu=p[5];
  double approx_sigma = std::sqrt( mu );
  //unsigned nmax = std::max(5, int(approx_sigma*3+mu) );
  unsigned nmax = 2;
  for (unsigned n=1; n<=nmax; n++ ){
    std::ostringstream fname;
    fname <<"fnpe_"<<n;
    TF1* ftmp = new TF1( fname.str().c_str(), pmtG_n_ped, -0.1 * peToVolts, 6 * peToVolts, 11 );
    std::cout << " n " << n << " mu " << mu << " approx sigma " << approx_sigma << std::endl; 
    pars[10] = n;
    ftmp->SetNpx(1000);
    ftmp->SetLineWidth(3);
    ftmp->SetParameters( pars );
    ftmp->SetLineColor( 52 + 5*n );
    std::cout << "Integral " << n << " " << ftmp->Integral(-0.1 * peToVolts, 6 * peToVolts) << std::endl;
    fcns.push_back( ftmp );
  }
  return fcns;
}


PMTResponsePed * PMTResponsePed::fInstance = nullptr;
TH1D* PMTResponsePed::fPDF = nullptr;
double PMTResponsePed::fWid = 0;
TFile* PMTResponsePed::fin = nullptr;

void PMTResponsePed::set_pedestal( const std::string& fname, const std::string& histname ){
  if ( fInstance == nullptr ) fInstance = new PMTResponsePed( 0.0 );
  fWid = 0.0;
  if ( fin == nullptr ) fin = new TFile( fname.c_str(), "read" );
  if ( !fin ) {
    std::cout<<"Could not find background root file: "<<fname<<std::endl;
    exit(0);
  }
  fPDF = (TH1D*) fin->Get( histname.c_str() );
  if ( !fPDF ){
    std::cout<<"Could not find histogram: "<<histname<<" in root file: "<<fname<<std::endl;
    exit(0);
  }
}


/*
Model 1:
x[0] = charge
p[0] = normalization (count)
p[1] = q1 single pe mean charge
p[2] = s1 single pe charge width
p[3] = mu mean number of pe
p[4] = w weight of exponential bg (0-1)
p[5] = alpha exponential constant bg
Change background to be part of pedestal signal!
 */
double model1( double * x, double *p ){
  double N = p[0];
  double q1 = p[1];
  double s1 = p[2];
  double mu = p[3];
  double w = p[4];
  double alpha = p[5];

  double retval = 0.0;

  retval = w * alpha * exp( -alpha*x[0] );

  retval += (1-w)*poisson_term( mu, 0 ) * PMTResponsePed::get_prob_density( x[0] );

  for ( unsigned npe =1; npe<5; ++npe ){
    retval += (1-w)*poisson_term( mu, npe ) * gaussian( x[0], npe*q1, sqrt( npe )*s1 )
      ;
  }

  return retval * N;
}


double model1bg( double * x, double *p ){
  double N = p[0];
  double q1 = p[1];
  double s1 = p[2];
  double mu = p[3];
  double w = p[4];
  double alpha = p[5];

  double retval = 0.0;

  //double binwid = PMTResponsePed_BinWid::GetBinWid();
  //double expnorm = alpha / (1.0 - exp( -alpha*binwid ) );
  //if ( x[0] < binwid ) retval = (1-w)*poisson_term( mu, 0 )/binwid;


  retval = w * alpha * exp( -alpha*x[0] );
  retval += (1-w)* poisson_term( mu, 0 ) * PMTResponsePed::get_prob_density( x[0] );


  return retval * N;
}


double model1npe( double * x, double *p ){
  double N = p[0];
  double q1 = p[1];
  double s1 = p[2];
  double mu = p[3];
  double w = p[4];
  double alpha = p[5];

  unsigned npe = p[6];

  double retval = 0.0;

  retval += (1-w) * poisson_term( mu, npe ) * gaussian( x[0], npe*q1, sqrt( npe ) * s1 );

  return retval * N;
}



/// Return a vector of functions that sums to be approxametely model1
/// First term is 0 pe + bg, 1 pe, 2 pe, ... etc
std::vector< TF1* > get_model1_components( double * p ){
  double pars[ 7 ];
  for ( unsigned i = 0; i < 7; ++i){
    pars[i] = p[i];
  }
  pars[6]=0.;

  std::vector< TF1* > fcns;
  TF1* bg = new TF1("Background", model1bg, 0., 0.065, 6 );
  bg->SetParameters( pars );
  //bg->SetLineColor( 52 );
  fcns.push_back( bg );

  for (unsigned n=1; n<6; ++n ){
    std::ostringstream fname;
    fname <<"fnpe_"<<n;
    TF1* ftmp = new TF1( fname.str().c_str(), model1npe, 0., 0.065, 7 );
    pars[6] = n;
    ftmp->SetLineWidth(3);
    ftmp->SetParameters( pars );
    ftmp->SetLineColor( 52 + 5*n );
    fcns.push_back( ftmp );
  }
  for ( TF1 * fcn : fcns ) fcn->SetNpx(1000);
  return fcns;
}




/// Test case code to test the fitting ***with pedestal***
int main(void){
//void pmt_response_function_tweakedCircle_volts(){
    
    const char * PMTform = "None";
    const char * face_position = "None";
    const char * TTS =  "None";
    
  //const char * PMTform = "Losange";
  //const char * face_position = "bottom";
  // try it out
  //TFile* fin = new TFile( "/vols/t2k/users/yia20/WCTE540/Analysis/Dark_Rate_540/05-07-04-2022/PMT_Square/peak_detect_2-0kV/root_files/PMT_Square_2-0kV_ADC_output.root", "read" );

  
  //TFile* fin = new TFile( Form("/home/a/WCTE540/PMT_data/%s/%sPMT/combined_datasets/hadded_all_%s_TTS_%s_1pe_1ms_LandauCorrectedVolts.root", TTS, PMTform, PMTform, face_position), "read" );

  TFile* fin = new TFile("../data/landau_run0023_triggerChannel2_laserPower155_Ch1Square120_Ch2Trigger-9999_Ch3Losange30_Ch4Triangle0_phi180.root");

  //Gain_TrianglePMT_1ms_0.200limit_pe1.17831e6.root
  //TTree *tree = (TTree*)fin->Get("testTree");
  TH1D * h = (TH1D*)fin->Get("PMT3_Corrected");
  std::cout << "Got the file" << std::endl;
  double gain = 1.068e6;
  double peToVolts = 0.01; //we now have the thing in number of volts -> better to calculate the gain directly 
  //gain * 50 * 1.6e-19 / (0.8e-9);
  std::cout << "The pe to volt conversion is " << peToVolts << " and the original gain is " << gain << std::endl;
 // int nevent = tree->GetEntries();

 // double volts;

 // TBranch *branch = tree->GetBranch("volts");
 // branch->SetAddress(&volts);

 // TH1D *h = new TH1D("h", "Volts", 220, 0., 0.069 );

 // for(int i = 0 ; i < nevent ; i++)
 //   tree->GetEntry(i);
 //   h->Fill(volts);

 // }
  //TH1D* h = (TH1D*) fin->Get("hADC_RO");//get_sum_histogram( fin, "hsum_4525" );

  //TH1D* h=(TH1D*)fin->Get("hqall");
  h->SetLineWidth(3);
  //h->Rebin(2);

  //PMTResponsePed::set_binwid( h->GetBinWidth(1) );
  //PMTResponsePed_BinWid::SetBinWid( h->GetBinWidth(1) );

  TCanvas *c=new TCanvas();
  c->cd();
  h->Draw("E1");
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.5);

  double Nfix = h->Integral( 0, h->GetNbinsX()+1, "width" );

  TF1* ff = new TF1( "pmt_response", pmtresponseped, -0.1 * peToVolts, 6 * peToVolts, 10 );

  double N0 = h->Integral(0,5);
  double Nrest = h->Integral(5,100); // beyond bin 15 dominated by exp background
  double height_onepe = h->Integral(5,1000); //4th bin is just above DR
  double mufix = Nrest / N0;
  mufix = log( mufix + 1 ); // correction factor to go from estimated mu, to true
  //double mufix = 5;
  std::cout << "Nfix = "<<Nfix<<" Nrest="<<Nrest<<" / N0 = "<<N0<<" = "<<mufix<<std::endl;
  std::cout << "mufix = " << mufix << std::endl;

  ff->SetNpx(1000);
  //Best values are square, circle, losange, triangle
  double best_Q0[4]= {5.27777e-4, 4.94150e-4, 5.04368e-4, 4.95821e-4};
  double best_s0[4] = {2.72169e-4, 2.7621e-4, 2.72831e-4, 2.77949e-4};
  double best_Q1[4]={1.57682e-2, 9.9694e-3, 1.10866e-2, 1.22630e-2};
  double best_s1[4] = {4.13079e-3, 3.52959e-3, 3.55860e-3, 3.34752e-3};
  double best_w[4] = {3.76437e-3, 3.42151e-5, -1.27453e-5, 2.75938e-3};
  double best_a[4] = {1.07855e2, -3.12277e-2, 1.16969e-2, 1.34141e2};

  int PMTnumber = 1; //-1 free fit, 0 Square fixed, 1 Circle, 2 Losange. 3 Tirangle

  //dec ent-ish guess for square
  ff->SetParNames("N", "Q_{0}", "#sigma_{0}", "Q_{1}", "#sigma_{1}", "#mu", "w", "#alpha", "A2", "A3" );
  ff->SetParameters( Nfix, 0.0005, 0.00021, 0.013, 0.001, mufix, 0.00005 * peToVolts, 0.045 * peToVolts, 1, 0  );
  ff->FixParameter( 0, Nfix );
  ff->FixParameter( 1, ff->GetParameter(1) );
  ff->FixParameter( 2, ff->GetParameter(2));
//   ff->FixParameter( 3, ff->GetParameter(3) );
//   ff->FixParameter( 8, ff->GetParameter(8) ) ;
//   ff->FixParameter( 9, ff->GetParameter(9) );

//   if (PMTnumber!=-1){
// // 	ff->FixParameter( 1, best_Q0[PMTnumber]);
// // 	ff->FixParameter( 2, best_s0[PMTnumber]);
// // 	ff->FixParameter( 3, best_Q1[PMTnumber]);
// // 	ff->FixParameter( 4, best_s1[PMTnumber]);
// // 	ff->FixParameter( 6, best_w[PMTnumber]);
// // 	ff->FixParameter( 7, best_a[PMTnumber]);
// // 	ff->FixParameter( 8, 0);
//   }

  h->Fit(ff, "", "", 0.01 * peToVolts, 6 * peToVolts ); //background at the end?

  for ( unsigned ipar=0 ; ipar < 10; ++ipar ) ff->ReleaseParameter(ipar);
  //ff->FixParameter( 0, ff->GetParameter(0) );
  //ff->FixParameter( 3, ff->GetParameter(3) ) ;
  //ff->FixParameter( 2, ff->GetParameter(2) ) ;
  //ff->FixParameter( 1, ff->GetParameter(1) ) ;
//   ff->FixParameter( 8, ff->GetParameter(8) ) ;
/*
  if (PMTnumber!=-1){
	ff->FixParameter( 1, best_Q0[PMTnumber]);
	ff->FixParameter( 2, best_s0[PMTnumber]);
	ff->FixParameter( 3, best_Q1[PMTnumber]);
	ff->FixParameter( 4, best_s1[PMTnumber]);
	ff->FixParameter( 6, best_w[PMTnumber]);
	ff->FixParameter( 7, best_a[PMTnumber]);
	ff->FixParameter( 8, 0);
  }*/

//   ff->FixParameter( 9, ff->GetParameter(9) );
  //ff->FixParameter( 5, ff->GetParameter(5) );
  //ff->FixParameter( 2, ff->GetParameter(2) );
  //ff->FixParameter( 6, ff->GetParameter(6) );
  //ff->FixParameter( 5, 0.068 );
  //ff->FixParameter(1, 0);

  ff->SetLineWidth(3);
  ff->SetLineColor(kRed+2);


  TFitResultPtr r = h->Fit(ff, "S", "",-0.2 * peToVolts, 6 * peToVolts ); //whole fit? here in  pe
  //h->Draw();
  //ff->Draw("same");
  //ff->FitResult();

  //ff->SetParameter(1,400);
  //ff->SetParameter(2,80);

  //h->Draw("same");
  //ff->SetParameter(0,2.85828e+05);
  //ff->SetParameter(1,-2.27603e-02);
  //ff->SetParameter(5,0.3);
  //ff->SetParameter(7,10.02);

  //std::cout << *ff->GetParameters() << std::endl; 
//These are all of the n-pe responses
  std::vector< TF1* > fcmp = get_pmt_response_components_ped( ff->GetParameters() );
  //std::vector< TF1* > fcmp_err = ff->GetErrors();

  int count = 0;
  TCanvas *c2=new TCanvas();
  TF1 * one_pe = new TF1("1pe", "gaus", -0.001, 0.012);
  c->cd();
  for ( TF1* f : fcmp ){
    f->Draw("same");
    if (count == 1){
	  one_pe->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
//	  height_onepe = f->Integral(-0.2 * peToVolts, 6 * peToVolts);
	  std::cout << f->GetParameter(3) << " " << f->GetParameter(4) << " " << f->GetParameter(5) << std::endl;
    }
    count += 1;
  }
  
  one_pe->GetXaxis()->SetTitle("volts");
  c2->cd();
  one_pe->Draw("AL");
  std::cout << "height one pe " << height_onepe << std::endl; 
  c2->Draw();
  c->cd();

  double Q0 = ff->GetParameter(1);
  double Q0_err = r->ParError(1);
  double Q1 = ff->GetParameter(3);
  double Q1_err = r->ParError(3);
  double chi2 = r->Chi2();
  double Q_err = TMath::Sqrt(Q0_err*Q0_err + Q1_err*Q1_err);
  double R = 50.;
  double Tr = 0.8e-9;
  double e = 1.6e-19;
  double gain_guess = (((Q1-Q0) * Tr)/R)/e ;
  double gain_error = ((Q_err * Tr)/R)/e;

  std::cout << PMTform << " " << face_position  << " Q0 = " << Q0 << " +/- " << Q0_err << " volts Q1 = " << Q1 << " +/- " << Q1_err << " volts gain = " << gain_guess << " +/- " << gain_error << std::endl;
  
  TPad *newpad=new TPad("newpad","a transparent pad",0.04,0.22,0.87,0.90);
  newpad->SetFillStyle(4000);
  newpad->Draw();
  newpad->cd();
  TText *t = new TText(.5,.5,Form("%s %s - gain %.3e +/- %.3e \n chi2 = %.f", PMTform, face_position, gain_guess, gain_error, chi2 ));
  t->SetX(0.12107383);
  t->SetY(0.950383);
  TText *t2 = new TText(.5,.5,Form("Counts above pedestal %.4e", height_onepe ));
  t2->SetX(0.12107383);
  t2->SetY(0.900383);
  t->SetTextColor(1);
  t->SetTextFont(22);
  t->SetTextSize(0.05);
  t->SetTextAngle(0);
  t->Draw();
  t2->SetTextColor(1);
  t2->SetTextFont(22);
  t2->SetTextSize(0.05);
  t2->SetTextAngle(0);
  t2->Draw();

  h->GetXaxis()->SetRangeUser(-0.1 * peToVolts, 6 * peToVolts);
  c->SetLogy();
  c->SaveAs(Form("./%s_%s_Working.pdf", PMTform, face_position));


  return 0;
}
