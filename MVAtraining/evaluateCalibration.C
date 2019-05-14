#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/PyMethodBase.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "Math/MinimizerOptions.h"
#include "TFitResult.h"

using namespace std;

TString process_;
TString dirPath_;
double min_, max_;
int nBins_ = 50;

pair<double, double> CountEventsWithFit(TH1 *hist, TString name);

void evaluateCalibration(TString file_ = "./ntuples/ntuBsMC2017.root"
    , TString method_ = "DNNOsMuonHLTJpsiMu_31025"
    , bool useTightSelection_ = false
    , int nEvents_ = -1
    )
{
    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptFit();

    cout<<"----- Parameters -----"<<endl;
    cout<<"file_ = "<<file_<<endl;
    cout<<"method_ = "<<method_<<endl;
    cout<<"useTightSelection_ = "<<useTightSelection_<<endl;
    cout<<"nEvents_ = "<<nEvents_<<endl;

    cout<<endl<<"----- BEGIN CODE"<<endl;
    TString path = "/lustre/cmswork/abragagn/BPH/BTag/OSMuon/src/PDAnalysis/Ntu/bin/";
    TString methodpath = "MVAtraining/dataset/weights/TMVAClassification_";

    auto *f = new TFile(path + file_);
    auto *t = (TTree*)f->Get("PDsecondTree");

    bool isData = false;

    if(file_.Contains("Bs")){
        process_ = "BsJPsiPhi";
        min_ = 5.25;
        max_ = 5.50;
        dirPath_ = "./Bs";
    }
    if(file_.Contains("Bu")){
        process_ = "BuJPsiK";
        min_ = 5.10;
        max_ = 5.50;
        dirPath_ = "./Bu";
    }
    if(file_.Contains("MC")){
        process_ = process_ + "MC";
        dirPath_ += "MC";
    }
    if(file_.Contains("Data")){
        process_ = process_ + "Data";
        dirPath_ += "Data";
        isData = true;
    }

    // CALIBRATION
    auto *f2 = new TFile("OSMuonTaggerCalibration" + process_ + ".root");
    f2->cd();
    TF1 *fCal = (TF1*)f2->Get("osMuonCal");
    f2->Close();
    delete f2;
    f->cd();

    cout<<"----- FILE OPEN"<<endl;

    double totalP = 0.;
    double totalPerr = 0;

    //MVA VARIABLES
    float muoPt;
    float muoEta;
    float muoDxy;
    float muoExy;
    float muoDz;
    float muoEz;
    float muoSoftMvaValue;
    float muoDrB;
    float muoPFIso;
    float muoConeCleanPt;
    float muoConeCleanPtRel;
    float muoConeCleanDr;
    float muoConeCleanEnergyRatio;
    float muoConeCleanQ;
    //TAGGING VARIABLES
    int osMuon, osMuonTag, osMuonCharge, ssbLund;
    //EVENT VARIABLES
    float ssbMass;
    int evtWeight, hltJpsiMu, ssbIsTight;

    //BOOKING
    t->SetBranchAddress("muoPt", &muoPt);
    t->SetBranchAddress("muoEta", &muoEta);
    t->SetBranchAddress("muoDxy", &muoDxy);
    t->SetBranchAddress("muoExy", &muoExy);
    t->SetBranchAddress("muoDz", &muoDz);
    t->SetBranchAddress("muoEz", &muoEz);
    t->SetBranchAddress("muoSoftMvaValue", &muoSoftMvaValue);
    t->SetBranchAddress("muoDrB", &muoDrB);
    t->SetBranchAddress("muoPFIso", &muoPFIso);
    t->SetBranchAddress("muoConeCleanPt", &muoConeCleanPt);
    t->SetBranchAddress("muoConeCleanPtRel", &muoConeCleanPtRel);
    t->SetBranchAddress("muoConeCleanDr", &muoConeCleanDr);
    t->SetBranchAddress("muoConeCleanEnergyRatio", &muoConeCleanEnergyRatio);
    t->SetBranchAddress("muoConeCleanQ", &muoConeCleanQ);
    t->SetBranchAddress("osMuon", &osMuon);
    t->SetBranchAddress("osMuonTag", &osMuonTag);
    t->SetBranchAddress("evtWeight", &evtWeight);
    t->SetBranchAddress("muoCharge", &osMuonCharge);
    t->SetBranchAddress("ssbLund", &ssbLund);
    t->SetBranchAddress("ssbMass", &ssbMass);
    t->SetBranchAddress("hltJpsiMu", &hltJpsiMu);
    t->SetBranchAddress("ssbIsTight", &ssbIsTight);

    //COMPUTE MVA
    TMVA::PyMethodBase::PyInitialize();
    TMVA::Reader reader("!Color:Silent");
    double mvaValue = -1.;

    reader.AddVariable("muoPt", &muoPt);
    reader.AddVariable("muoEta", &muoEta);
    reader.AddVariable("muoDxy", &muoDxy);
    reader.AddVariable("muoExy", &muoExy);
    reader.AddVariable("muoDz", &muoDz);
    reader.AddVariable("muoEz", &muoEz);
    reader.AddVariable("muoSoftMvaValue", &muoSoftMvaValue);
    reader.AddVariable("muoDrB", &muoDrB);
    reader.AddVariable("muoPFIso", &muoPFIso);
    reader.AddVariable("muoConeCleanPt", &muoConeCleanPt);
    reader.AddVariable("muoConeCleanPtRel", &muoConeCleanPtRel);
    reader.AddVariable("muoConeCleanDr", &muoConeCleanDr);
    reader.AddVariable("muoConeCleanEnergyRatio", &muoConeCleanEnergyRatio);
    reader.AddVariable("muoConeCleanQ", &muoConeCleanQ);
    reader.BookMVA( method_, path + methodpath + method_ + ".weights.xml" );

    auto *hMass    = new TH1F( "hMass","hMass", nBins_, min_, max_ );
    auto *hMassRT  = new TH1F( "hMassRT","hMassRT", nBins_, min_, max_ );
    auto *hMassWT  = new TH1F( "hMassWT","hMassWT", nBins_, min_, max_ );
    auto *hMassTagged = new TH1F( "hMassTagged","hMassTagged", nBins_, min_, max_ );

    cout<<"----- BOOKING COMPLETED"<<endl;

    //EVENT LOOP
    cout<<"----- BEGIN LOOP"<<endl;

    if(nEvents_ == -1) nEvents_ = t->GetEntries();

    for(int i=0; i<nEvents_; ++i){
        if(i%100000==0) cout<<"----- at event "<<i<<endl;

        t->GetEntry(i);

        //EVENT SELECTION
        if(!hltJpsiMu) continue;
        if(useTightSelection_ && !ssbIsTight) continue;
        hMass->Fill(ssbMass, evtWeight);

        //MUON SELECTION
        if(!osMuon) continue;
        if(muoSoftMvaValue<=0.35) continue;

        //NAN INF VARIABLES BUG FIX
        bool nanFlag = isnan(muoDxy); // || isnan(muoJetDFprob);
        bool infFlag = isinf(muoConeCleanEnergyRatio);
        if(nanFlag || infFlag) continue;

        //COMPLEX VARIABLE COMPUTATION
        // None at the moment

        //TAGGING
        mvaValue = reader.EvaluateMVA(method_);
        double evtW = fCal->Eval(1-mvaValue); // calibration
        double evtWerr = sqrt(pow(fCal->GetParError(0),2)+pow((fCal->GetParError(1))*(1-mvaValue),2));
        totalP += pow(1.-2.*evtW, 2)*evtWeight;
        double errp = sqrt( pow(8*evtW-4,2)*pow(evtWerr,2) );
        totalPerr += pow(evtWeight*errp,2);

        int evtTag = -1*osMuonCharge;
        bool isTagRight = TMath::Sign(1, ssbLund) == evtTag;
        if(isTagRight!=osMuonTag) cout<<"!!!! isTagRight != osMuonTag"<<endl;

        hMassTagged->Fill(ssbMass, evtWeight);

        if(isTagRight){
            hMassRT->Fill(ssbMass, evtWeight);
        }else{
            hMassWT->Fill(ssbMass, evtWeight);
        }

    }

    cout<<"----- EVENTS LOOP ENDED"<<endl;

    // PERFORMANCE OUTPUT
    int nRT = hMassRT->Integral(); //Integral() takes in consideration event weights
    int nWT = hMassWT->Integral();
    int nTagged = hMassTagged->Integral();
    int nTot = hMass->Integral();

    if(isData){ //for data fit mass
        nRT = CountEventsWithFit(hMassRT, "histTotMass_RT").first;
        nWT = CountEventsWithFit(hMassWT, "histTotMass_WT").first;
        nTagged = CountEventsWithFit(hMassTagged, "histTotMass_Tagged").first;
        nTot = CountEventsWithFit(hMass, "histTotMass").first;
    }

    double effBase = (double)(nRT+nWT)/nTot;
    double wBase = (double)nWT/(nRT+nWT);
    double pBase = effBase*pow(1.-2.*wBase,2);

    cout<<"nTot "<<nTot<<endl;
    cout<<"nRT "<<nRT<<endl;
    cout<<"nWT "<<nWT<<endl;

    cout<<endl;
    cout<<"Base efficiency = "<<100*effBase<<"%"<<endl;
    cout<<"Base mistag = "<<100*wBase<<"%"<<endl;
    cout<<"Base power = "<<100*pBase<<"%"<<endl;

    double park = totalP;
    totalP /= (double)nTot;
    totalPerr /= pow((double)nTot,2);
    totalPerr += nTot*pow(park,2)/pow((double)nTot,4);
    totalPerr = sqrt(totalPerr);

    cout<<endl;
    cout<<"Per-event-mistag power = "<<100.*totalP<<" +- "<<100.*totalPerr<<" % (+"<<100*(totalP - pBase)/pBase<<"%)"<<endl;
    cout<<endl;

    f->Close();
    delete f;
    return;

}

// ----------------------------------------------
pair<double, double> CountEventsWithFit(TH1 *hist, TString name = "hist"){

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
    }

    func->SetNpx(2000);

    auto c5 = new TCanvas();
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(.75);
    TFitResultPtr r = hist->Fit("func","MRLS");
    hist->Draw("PE");
    hist->SetMinimum(0);

    TF1 *fit = hist->GetFunction("func");
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
    c5->Print(dirPath_ + "/" + name + ".pdf");

    double nEvt = fit->GetParameter(1);
    nEvt += fit->GetParameter(2);
    nEvt += fit->GetParameter(3);
    nEvt/=hist->GetBinWidth(1);

    TMatrixDSym cov = r->GetCovarianceMatrix();
    double errN = sqrt(cov(1,1)+cov(2,2)+cov(3,3)+2*(cov(1,2)+cov(1,3)+cov(2,3)))/hist->GetBinWidth(1);

    return make_pair(nEvt, errN);
}