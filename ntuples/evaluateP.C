#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "Math/MinimizerOptions.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TEfficiency.h"

using namespace std;

TString year_;
TString process_;
TString dir_ = "./";
TString suffix_;
double min_, max_;
int nBins_=50;
pair<double, double> CountEventsWithFit(TH1 *hist);

void evaluateP(TString file = "./ntuBsMC2017.root",  TString cutEvt_ = "", TString cut_ = "", TString suffix = "")
{
    gErrorIgnoreLevel = kWarning;
    TFile *f = new TFile(file);
    TTree *t = (TTree*)f->Get("PDsecondTree");
    if(file.Contains("16")) year_ = "2016";
    if(file.Contains("17")) year_ = "2017";
    if(file.Contains("Bs")){
        process_ = "BsJPsiPhi";
        min_ = 5.25;
        max_ = 5.50;
    }
    if(file.Contains("Bu")){
        process_ = "BuJPsiK";
        min_ = 5.10;
        max_ = 5.50;
    }

    if(file.Contains("MC")) process_ = process_ + "MC";
    if(file.Contains("Data")) process_ = process_ + "Data";

    TH1F *ssB       = new TH1F( "ssB", "ssB", nBins_, min_, max_ );
    TH1F *ssB_TT    = new TH1F( "ssB_TT", "ssB_TT", nBins_, min_, max_ );
    TH1F *ssB_RT    = new TH1F( "ssB_RT", "ssB_RT", nBins_, min_, max_ );
    TH1F *ssB_WT    = new TH1F( "ssB_WT", "ssB_WT", nBins_, min_, max_ );

    TString cut = "1";
    //cut += "&&!isnan(muoDxy)&&!isnan(muoJetDFprob)&&!isinf(muoJetEnergyRatio)&&!isinf(muoConeEnergyRatio)";
    TString cutEvt = "hltJpsiMu&&ssbIsTight";

    if(cutEvt_ != "") cutEvt = cutEvt + "&&" + cutEvt_;
    if(cut_ != "")    cut = cut + "&&" + cut_;

    TString base =  "evtWeight*((" + cutEvt;
    TString cutTT = base + ")&&(" + cut + ")&&osMuon)";
    TString cutRT = base + ")&&(" + cut + ")&&osMuon&&osMuonTag==1)";
    TString cutWT = base + ")&&(" + cut + ")&&osMuon&&osMuonTag==0)";

    t->Project("ssB", "ssbMass", base + "))" );
    t->Project("ssB_RT", "ssbMass", cutRT );
    t->Project("ssB_WT", "ssbMass", cutWT );
    t->Project("ssB_TT", "ssbMass", cutTT );

    pair<double, double> nTot;
    pair<double, double> nTT;
    pair<double, double> nRT;
    pair<double, double> nWT;
    pair<double, double> eff;
    pair<double, double> w;
    pair<double, double> power;

    nTot.first = ssB->Integral();
    nRT.first = ssB_RT->Integral();
    nWT.first = ssB_WT->Integral();
    nTT.first = ssB_TT->Integral();

    nTot.second = sqrt(nTot.first);
    nRT.second = sqrt(nRT.first);
    nWT.second = sqrt(nWT.first);
    nTT.second = sqrt(nTT.first);

    if(file.Contains("Data")){
        nTot = CountEventsWithFit(ssB);
        nRT = CountEventsWithFit(ssB_RT);
        nWT = CountEventsWithFit(ssB_WT);        
        nTT = CountEventsWithFit(ssB_TT);        
    }

    double tagged = nRT.first+nWT.first;

    eff.first = (nTT.first)/(nTot.first);
    w.first = nWT.first/(nTT.first);
    power.first = eff.first*pow((1-2*w.first), 2);

    if(file.Contains("Data")){
        double TT2 = pow(nTT.first,2);
        double WT2 = pow(nWT.first,2);
        double Tot2 = pow(nTot.first,2);
        double sTT2 = pow(nTT.second,2);
        double sWT2 = pow(nWT.second,2);
        double sTot2 = pow(nTot.second,2);

        eff.second  = sqrt( sTT2/Tot2 + sTot2*TT2/pow(Tot2,2) );
        w.second    = sqrt( sWT2/TT2 + sTT2*WT2/pow(TT2,2) );
        power.second = sqrt( sTT2*pow((TT2-4*WT2),2)/(Tot2*pow(TT2,2)) 
            + 16*sWT2*pow(nTT.first-2*nWT.first,2)/(Tot2*TT2)
            + sTot2*pow(nTT.first-2*nWT.first,4)/(pow(Tot2,2)*TT2) );
    }else{
        double cl = 0.6827;
        double alpha = 0.5;
        eff.second = (TEfficiency::Bayesian(nTot.first,tagged,cl,alpha,alpha,1)-TEfficiency::Bayesian(nTot.first,tagged,cl,alpha,alpha,0))/2;
        w.second = (TEfficiency::Bayesian(tagged,nWT.first,cl,alpha,alpha,1)-TEfficiency::Bayesian(tagged,nWT.first,cl,alpha,alpha,0))/2;
        power.second = sqrt(16*pow(eff.first,2)*pow(1-2*w.first,2)*pow(w.second,2) + pow(1-2*w.first,4)*pow(eff.second,2));
    }

    cout<<endl;

    cout<<"Bs = "<<nTot.first<<" +- "<<(int)nTot.second<<endl;
    cout<<std::setprecision(4);
    cout<<"Eff = "<<100*eff.first<<" +- "<<100*eff.second<<" %"<<" ["<<100*eff.second/eff.first<<"%]"<<endl;
    cout<<"Mistag = "<<100*w.first<<" +- "<<100*w.second<<" %"<<" ["<<100*w.second/w.first<<"%]"<<endl;
    cout<<"Power = "<<100*power.first<<" +- "<<100*power.second<<" %"<<" ["<<100*power.second/power.first<<"%]"<<endl;

    cout<<endl;

    return;

}

pair<double, double> CountEventsWithFit(TH1 *hist){

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );

    double mean = 5.3663;
    double sigma = 0.015;
    if(process_.Contains("BsJPsiPhi")) mean = 5.3663;
    if(process_.Contains("BuJPsiK"))   mean = 5.2793;
    if(process_.Contains("BdJPsiKx"))  mean = 5.2796;

    TString sgnDef = "[1]*TMath::Gaus(x, [0], [4], true)";
    sgnDef +=       "+[2]*TMath::Gaus(x, [0], [5], true)";
    sgnDef +=       "+[3]*TMath::Gaus(x, [0], [6], true)";

    TString bkgDef = "[7]+[8]*TMath::Erfc([9]*(x-[10]))";
    if(hist->GetEntries()<=250) bkgDef = "[7]";

    TString funcDef = sgnDef + "+" + bkgDef;

    TF1 *func = new TF1("func", funcDef, min_, max_);

    //SIGNAL
    double limit = hist->GetEntries()*hist->GetBinWidth(1);

    func->SetParameter(0, mean);
    func->SetParameter(1, limit/3);
    func->SetParameter(2, limit/3);
    func->SetParameter(3, limit/3);
    func->SetParameter(4, sigma);
    func->SetParameter(5, sigma);
    func->SetParameter(6, sigma);
    func->SetParLimits(0, mean-sigma, mean+sigma);
    func->SetParLimits(1, 0, limit);
    func->SetParLimits(2, 0, limit);
    func->SetParLimits(3, 0, limit);
    func->SetParLimits(4, sigma/2, sigma*2);
    func->SetParLimits(5, sigma/2, sigma*2);
    func->SetParLimits(6, sigma/2, sigma*2);

    //BKG
    double bkgHeight = 0.;
    int nBinsBkgEst = 7;
    for(int i=0; i<nBinsBkgEst; i++)
        bkgHeight += hist->GetBinContent(hist->GetNbinsX()-i);
    bkgHeight /= nBinsBkgEst;

    func->SetParameter(7, bkgHeight);
    if(hist->GetEntries()>250){
        func->SetParameter(8, hist->GetBinContent(1)/2);
        func->SetParameter(9, 10);
        func->SetParameter(10, 5);
        func->SetParLimits(7, 0, bkgHeight*2);
        func->SetParLimits(8, 0, hist->GetBinContent(1));
        //func->SetParLimits(9, 10, 1e3);
        //func->SetParLimits(10, 5.0, mean);
    }

    TFitResultPtr r = hist->Fit("func","MRLS");
    TF1 *fit = hist->GetFunction("func");
/*
    auto c1 = new TCanvas();
    func->SetNpx(2000);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(.75);
    hist->Draw("PE");
    hist->SetMinimum(0.1);

    TF1 *f1 = new TF1("f1","[0]*TMath::Gaus(x, [1], [2], true)", min_, max_);
    TF1 *f2 = new TF1("f2","[0]*TMath::Gaus(x, [1], [2], true)", min_, max_);
    TF1 *f3 = new TF1("f3","[0]*TMath::Gaus(x, [1], [2], true)", min_, max_);
    TF1 *f4 = new TF1("f4","[0]+[1]*TMath::Erfc([2]*(x-[3]))", min_, max_);
    if(hist->GetEntries()<=250) f4 = new TF1("f4","[0]", min_, max_);

    f1->SetParameters(fit->GetParameter(1),fit->GetParameter(0),fit->GetParameter(4));
    f2->SetParameters(fit->GetParameter(2),fit->GetParameter(0),fit->GetParameter(5));
    f3->SetParameters(fit->GetParameter(3),fit->GetParameter(0),fit->GetParameter(6));
    if(hist->GetEntries()<=250) f4->SetParameter(0,fit->GetParameter(7));
    else f4->SetParameters(fit->GetParameter(7),fit->GetParameter(8),fit->GetParameter(9),fit->GetParameter(10));

    f1->SetLineColor(kBlue);
    f2->SetLineColor(kViolet);
    f3->SetLineColor(kAzure);
    f4->SetLineColor(kOrange);
    f1->SetLineStyle(2);
    f2->SetLineStyle(2);
    f3->SetLineStyle(2);
    f4->SetLineStyle(2);

    f1->Draw("same");
    f2->Draw("same");
    f3->Draw("same");
    f4->Draw("same");
    TString name = hist->GetName(); 

    c1->Print(dir_ + name + "_" + year_ + suffix_ + ".pdf");
    c1->Print(dir_ + name + "_" + year_ + suffix_ + ".png");
*/
    double nEvt = fit->GetParameter(1);
    nEvt += fit->GetParameter(2);
    nEvt += fit->GetParameter(3);
    nEvt/=hist->GetBinWidth(1);

    TMatrixDSym cov = r->GetCovarianceMatrix();
    double errN = sqrt(cov(1,1)+cov(2,2)+cov(3,3)+2*(cov(1,2)+cov(1,3)+cov(2,3)))/hist->GetBinWidth(1);

    return make_pair(nEvt, errN);
}
