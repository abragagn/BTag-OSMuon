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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "Math/MinimizerOptions.h"
#include "TKDE.h"

using namespace std;

TString process_;
float min_, max_;
int nBins_ = 50;

float CountEventsWithFit(TH1 *hist, TString name);

void fitMVA(TString file_ = "ntuBsMC2017.root"
            , TString method_ = "DNNOsMuonHLTJpsiMu_test241"
            , TString mode_ = "CREATE"
            , bool addMva_ = false      //currently deprecated
            , bool readMva_ = false     //currently deprecated
            , bool useTightSelection = false
            , int nEvents_ = -1)
{
    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptFit();
    TString path = "/lustre/cmswork/abragagn/BPH/BTag/OSMuon/src/PDAnalysis/Ntu/bin/ntuples/";

    cout<<"----- Parameters -----"<<endl;
    cout<<"file_ "<<file_<<endl;
    cout<<"method_ "<<method_<<endl;
    cout<<"mode_ "<<mode_<<endl;
    cout<<"addMva_ "<<addMva_<<endl;
    cout<<"readMva_ "<<readMva_<<endl;
    cout<<"nEvents_ "<<nEvents_<<endl;

    if(mode_ != "CREATE" && mode_ != "USE"){
        cout<<"WRONG MODE_"<<endl;
        return;
    }

    if(addMva_ && readMva_){
        cout<<"CAN'T READ AND WRITE MVA AT THE SAME TIME"<<endl;
        return;
    }

    cout<<endl<<"----- BEGIN CODE"<<endl;

    auto *f = new TFile(path + file_);
    auto *t = (TTree*)f->Get("PDsecondTree");
    cout<<"----- FILE READ"<<endl;

//----------DECLARE STUFF USED IN USE MODE----------
    TString HLT;
    if(method_.Contains("JpsiMu")) HLT = "HltJpsiMu";
    if(method_.Contains("JpsiTrkTrk")) HLT = "HltJpsiTrkTrk";
    if(method_.Contains("JpsiTrk") && !method_.Contains("JpsiTrkTrk")) HLT = "HltJpsiTrk";

    if(file_.Contains("Bs")){
        process_ = "BsJPsiPhi";
        min_ = 5.25;
        max_ = 5.50;
    }
    if(file_.Contains("Bu")){
        process_ = "BuJPsiK";
        min_ = 5.10;
        max_ = 5.50;
    }

    if(file_.Contains("MC")) process_ = process_ + "MC";
    if(file_.Contains("Data")) process_ = process_ + "Data";

    cout<<"HLT "<<HLT<<endl;

    int nCat;
    float *catEdgeL;
    float *catEdgeR;
    float *catMistag;
    
    TF1 *perEvtW;
    TString perEvtWFormula;
    float avgW;

    float totPCut = 0.;
    float totPCat = 0.;
    float totPFit = 0.;
    float totPKde = 0.;

    TGraph *g_pdfW;
    TGraph *g_pdfW_extended;

    float evtWcat = -1.;
    float evtWfit = -1.;
    float evtWkde = -1.;
    float evtWkde_ext = -1.;

    int nBinCheck = 20;
    auto *wCalc = new float[nBinCheck];
    TH1F **hMassRT[4];
    TH1F **hMassWT[4];
    for(int i=0; i<4; ++i){
        hMassRT[i] = new TH1F*[nBinCheck];
        hMassWT[i] = new TH1F*[nBinCheck];
    }
    
    float pass = 1./nBinCheck;  //1. = 1.-0. (mva max - mva min)

//----------READ INPUT FILE AND BOOK METHODS----------
    if(mode_ == "USE")
    {
        cout<<"----- USE METHOD MODE"<<endl;

        std::ifstream ifs("OSMuonTagger" + HLT + "Categories.txt", std::ifstream::in);
        if(!ifs.is_open()) return;

        ifs >> method_;
        //CAT
        ifs >> nCat;
        catEdgeL  = new float[nCat];
        catEdgeR  = new float[nCat];
        catMistag = new float[nCat];
        for(int i=0; i<nCat; ++i){
            ifs >> catEdgeL[i];
            ifs >> catEdgeR[i];
            ifs >> catMistag[i];
        }
        cout<<"Input categories [L-edge R-edge mistag]"<<endl;
        for(int i=0; i<nCat; ++i)
            cout<<catEdgeL[i]<<" "<<catEdgeR[i]<<" "<<catMistag[i]<<endl;

        //FIT
        ifs >> perEvtWFormula;
        int nPar;
        ifs >> avgW;
        ifs >> nPar;
        cout<<endl<<"Input function"<<endl;
        cout<<perEvtWFormula<<endl;
        cout<<"avgW = "<<avgW<<endl;
        perEvtW = new TF1("perEvtW", perEvtWFormula, 0., 1.);
        for(int i=0;i<nPar;++i){
            float par;
            ifs >> par;
            perEvtW->SetParameter(i, par);
            cout<<"p"<<i<<" = "<<par<<endl;
        }

        ifs.close();

        auto *f2 = new TFile("OSMuonTagger" + HLT + "KDE.root");
        f2->cd();
        g_pdfW = (TGraph*)f2 ->Get("pdfW");
        g_pdfW_extended = (TGraph*)f2 ->Get("pdfW_extended");
        f2->Close();
        delete f2;
        f->cd();

        for(int j=0;j<nBinCheck;++j){
            wCalc[j] = (float)j*pass + pass/2;
            for(int i=0; i<4; ++i){
                hMassRT[i][j] = new TH1F(TString::Format("mbRT%i%i", i, j),"", nBins_, min_, max_);
                hMassWT[i][j] = new TH1F(TString::Format("mbWT%i%i", i, j),"", nBins_, min_, max_);
            }
        } 

        cout<<"----- METHOD BOOKED"<<endl;      
    }

//----------DECLARE VARIABLES----------
    cout<<"----- VARIABLE DECLARATION"<<endl;

    //MVA VARIABLES
    float muoPt;
    float absmuoEta;
    float muoDxy;
    float absmuoDz;
    float muoSoftMvaValue;
    float muoDrB;
    float muoPFIso;
    float muoJetConePt;
    float muoJetConePtRel;
    float muoJetConeDr;
    float muoJetConeEnergyRatio;
    float muoJetDFprob;
    float muoJetConeSize;
    float muoJetConeQ;
    float muoCharge;
    //OTHER VARIABLES
    float muoEta;
    float muoDz;
    float muoJetPt, muoConePt;
    float muoJetPtRel, muoConePtRel;
    float muoJetDr, muoConeDr;
    float muoJetEnergyRatio, muoConeEnergyRatio;
    int   muoJetSize, muoConeSize;
    float muoJetQ, muoConeQ;

    //TAGGING VARIABLES
    int osMuon, osMuonTag, osMuonCharge, ssbLund;
    //EVENT VARIABLES
    float ssbMass;
    int evtWeight, hltJpsiMu, hltJpsiTrkTrk, hltJpsiTrk, ssbIsTight;

    t->SetBranchAddress("muoPt", &muoPt);
    t->SetBranchAddress("muoEta", &muoEta);
    t->SetBranchAddress("muoDxy", &muoDxy);
    t->SetBranchAddress("muoDz", &muoDz);
    t->SetBranchAddress("muoSoftMvaValue", &muoSoftMvaValue);
    t->SetBranchAddress("muoDrB", &muoDrB);
    t->SetBranchAddress("muoPFIso", &muoPFIso);
    t->SetBranchAddress("muoJetDFprob", &muoJetDFprob);
    t->SetBranchAddress("muoJetPt", &muoJetPt);
    t->SetBranchAddress("muoJetPtRel", &muoJetPtRel);
    t->SetBranchAddress("muoJetDr", &muoJetDr);
    t->SetBranchAddress("muoJetEnergyRatio", &muoJetEnergyRatio);
    t->SetBranchAddress("muoJetSize", &muoJetSize);
    t->SetBranchAddress("muoJetQ", &muoJetQ);
    t->SetBranchAddress("muoConePt", &muoConePt);
    t->SetBranchAddress("muoConePtRel", &muoConePtRel);
    t->SetBranchAddress("muoConeDr", &muoConeDr);
    t->SetBranchAddress("muoConeEnergyRatio", &muoConeEnergyRatio);
    t->SetBranchAddress("muoConeSize", &muoConeSize);
    t->SetBranchAddress("muoConeQ", &muoConeQ);
    t->SetBranchAddress("osMuon", &osMuon);
    t->SetBranchAddress("osMuonTag", &osMuonTag);
    t->SetBranchAddress("evtWeight", &evtWeight);
    t->SetBranchAddress("muoCharge", &osMuonCharge);
    t->SetBranchAddress("ssbLund", &ssbLund);
    t->SetBranchAddress("ssbMass", &ssbMass);
    t->SetBranchAddress("hltJpsiMu", &hltJpsiMu);
    t->SetBranchAddress("hltJpsiTrkTrk", &hltJpsiTrkTrk);
    t->SetBranchAddress("hltJpsiTrk", &hltJpsiTrk);
    t->SetBranchAddress("ssbIsTight", &ssbIsTight);

//----------COMPUTE MVA----------
    cout<<"----- BEGIN MVA SETUP"<<endl;

    TMVA::Reader reader("!Color:Silent");
    float mvaValue = -1.;

    TMVA::PyMethodBase::PyInitialize();

    reader.AddVariable("muoPt", &muoPt);
    reader.AddVariable("abs_muoEta := fabs(muoEta)", &absmuoEta);
    reader.AddVariable("muoDxy", &muoDxy);
    reader.AddVariable("abs_muoDz := fabs(muoDz)", &absmuoDz);
    reader.AddVariable("muoSoftMvaValue", &muoSoftMvaValue);
    reader.AddVariable("muoDrB", &muoDrB);
    reader.AddVariable("muoPFIso", &muoPFIso);
    reader.AddVariable("muoJetConePt := muoJetPt != -1 ? muoJetPt : muoConePt", &muoJetConePt);
    reader.AddVariable("muoJetConePtRel := muoJetPt != -1 ? muoJetPtRel : muoConePtRel", &muoJetConePtRel);
    reader.AddVariable("muoJetConeDr := muoJetPt != -1 ? muoJetDr : muoConeDr", &muoJetConeDr);
    reader.AddVariable("muoJetConeEnergyRatio := muoJetPt != -1 ? muoJetEnergyRatio : muoConeEnergyRatio", &muoJetConeEnergyRatio);
    reader.AddVariable("muoJetDFprob", &muoJetDFprob);
    reader.AddVariable("muoJetConeSize := muoJetPt != -1 ? muoJetSize : muoConeSize", &muoJetConeSize);
    reader.AddVariable("muoJetConeQ := muoJetPt != -1 ? muoJetQ : muoConeQ", &muoJetConeQ);
    reader.BookMVA( method_, path + "MVAtraining/dataset/weights/TMVAClassification_" + method_ + ".weights.xml" );

    int nBinsMva = 1000;
    auto *mva    = new TH1F( "mva",    "mva",    nBinsMva, 0.0, 1.0 );
    auto *mva_RT = new TH1F( "mva_RT", "mva_RT", nBinsMva, 0.0, 1.0 );
    auto *mva_WT = new TH1F( "mva_WT", "mva_WT", nBinsMva, 0.0, 1.0 );
    auto *evtMistag = new TH1F( "evtMistag", "evtMistag", 50, 0.0, 1.0 );
    mva_WT->SetLineColor(kRed);

    vector<double> vKDERT;
    vector<double> vKDEWT;

    //EVENT LOOP
    cout<<"----- BEGIN LOOP"<<endl;

    if(nEvents_ == -1) nEvents_ = t->GetEntries();
    int nEventsRead = 0;
    for(int i=0; i<nEvents_; ++i){
        if(i%100000==0) cout<<"----- at event "<<i<<endl;
        t->GetEntry(i);

        //PRESELECTION EVENT
        if(HLT == "HltJpsiMu" && !hltJpsiMu) continue;
        if(HLT == "HltJpsiTrkTrk" && !hltJpsiTrkTrk) continue;
        if(HLT == "HltJpsiTrk" && !hltJpsiTrk) continue;
        if(useTightSelection && !ssbIsTight) continue;

        nEventsRead++;

        //PRESELECTION MUON
        if(!osMuon) continue;
        if((fabs(muoEta)<1.2 && muoSoftMvaValue<=0.891)) continue;
        if((fabs(muoEta)>=1.2 && muoSoftMvaValue<=0.8925)) continue;

        //CHECK NAN INF VARIABLES, THIS IS A BUG FIX
        bool nanFlag = isnan(muoDxy) || isnan(muoJetDFprob);
        bool infFlag = isinf(muoJetEnergyRatio) || isinf(muoConeEnergyRatio);
        if(nanFlag || infFlag) continue;

        absmuoEta = fabs(muoEta);
        absmuoDz = fabs(muoDz);
        muoJetConePt = muoJetPt != -1 ? muoJetPt : muoConePt;
        muoJetConePtRel = muoJetPt != -1 ? muoJetPtRel : muoConePtRel;
        muoJetConeDr = muoJetPt != -1 ? muoJetDr : muoConeDr;
        muoJetConeEnergyRatio = muoJetPt != -1 ? muoJetEnergyRatio : muoConeEnergyRatio;
        muoJetConeSize = muoJetPt != -1 ? (float)muoJetSize : (float)muoConeSize;
        muoJetConeQ = muoJetPt != -1 ? muoJetQ : muoConeQ;

        mvaValue = reader.EvaluateMVA(method_);

        int evtTag = -1*osMuonCharge;

        if(mode_=="USE"){
            if( mvaValue < catEdgeL[0] ){ //EVENTS WITH MVA<MVA_MIN USED IN CREATING THE METHOD
                cout<<"-----Undercat "<<mvaValue<<" at "<<i<<endl;
                evtWcat = catMistag[0];
                evtWfit = perEvtW->Eval(catEdgeL[0]);
            }else if( mvaValue >= catEdgeR[nCat-1] ){ //EVENTS WITH MVA>=MVA_MAX USED IN CREATING THE METHOD
                cout<<"-----Overcat "<<mvaValue<<" at "<<i<<endl;
                evtWcat = catMistag[nCat-1];
                evtWfit = perEvtW->Eval(catEdgeR[nCat-1]);
            }else{
                for(int j=0; j<nCat; ++j){ //EVENTS WITH MVA IN [MVA_MIN, MVA_MAX] USED IN CREATING THE METHOD
                    if(( mvaValue >= catEdgeL[j]) && (mvaValue < catEdgeR[j]) )
                    { 
                        evtWcat = catMistag[j]; //CAT
                        break; 
                    }
                }
                evtWfit = perEvtW->Eval(mvaValue);  //FIT
            }
            evtWkde = g_pdfW->Eval(mvaValue);   //KDE
            evtWkde_ext = g_pdfW_extended->Eval(mvaValue);  //KDE EXTENDED

            totPCat += pow(1.-2.*evtWcat, 2)*evtWeight;
            totPKde += pow(1.-2.*evtWkde ,2)*evtWeight;

            evtMistag->Fill(evtWkde, evtWeight);

            for(int j=0;j<nBinCheck;++j){
                if( (evtWcat>=(float)j*pass) && (evtWcat<((float)j*pass+pass)) ){
                    if(TMath::Sign(1, ssbLund) == evtTag) hMassRT[0][j]->Fill(ssbMass, evtWeight);
                    if(TMath::Sign(1, ssbLund) != evtTag) hMassWT[0][j]->Fill(ssbMass, evtWeight);
                    break;
                }
            }
            for(int j=0;j<nBinCheck;++j){
                if( (evtWfit>=(float)j*pass) && (evtWfit<((float)j*pass+pass)) ){
                    if(TMath::Sign(1, ssbLund) == evtTag) hMassRT[1][j]->Fill(ssbMass, evtWeight);
                    if(TMath::Sign(1, ssbLund) != evtTag) hMassWT[1][j]->Fill(ssbMass, evtWeight);
                    break;
                }
            }
            for(int j=0;j<nBinCheck;++j){             
                if( (evtWkde>=(float)j*pass) && (evtWkde<((float)j*pass+pass)) ){
                    if(TMath::Sign(1, ssbLund) == evtTag) hMassRT[2][j]->Fill(ssbMass, evtWeight);
                    if(TMath::Sign(1, ssbLund) != evtTag) hMassWT[2][j]->Fill(ssbMass, evtWeight);
                    break;
                }
            }                
            for(int j=0;j<nBinCheck;++j){
                if( (evtWkde_ext>=(float)j*pass) && (evtWkde_ext<((float)j*pass+pass)) ){
                    if(TMath::Sign(1, ssbLund) == evtTag) hMassRT[3][j]->Fill(ssbMass, evtWeight);
                    if(TMath::Sign(1, ssbLund) != evtTag) hMassWT[3][j]->Fill(ssbMass, evtWeight);
                }
            }

        }

        mva->Fill(mvaValue, evtWeight);
        if(osMuonTag == 1){
            mva_RT->Fill(mvaValue, evtWeight);
            vKDERT.push_back(mvaValue);
            if(evtWeight==2) vKDERT.push_back(mvaValue); //avoid using weights in kde
        }
        if(osMuonTag == 0){
            mva_WT->Fill(mvaValue, evtWeight);
            vKDEWT.push_back(mvaValue);
            if(evtWeight==2) vKDEWT.push_back(mvaValue);
        }
    }

    int nRT = mva_RT->Integral(); //Integral() take in consideration weight, GetEntries() not
    int nWT = mva_WT->Integral();

    cout<<"----- MVA HISTOGRAMS FILLED"<<endl;
    cout<<"nRT "<<nRT<<endl;
    cout<<"nWT "<<nWT<<endl;

//----------USE----------
    if(mode_ == "USE"){

        for(int i=0; i<4; ++i){
            for(int j=0;j<nBinCheck;++j){
                cout<<i<<" "<<j<<" --- "<<hMassRT[i][j]->Integral()<<" "<<hMassWT[i][j]->Integral()<<endl;
            }
            cout<<endl;
        }

        for(int i=0; i<4; ++i){
            cout<<"TYPE "<<i<<endl;

            vector<float> vX;
            vector<float> vY;
            vector<float> vEY;

            int minEntries = 5;

            int rebinThr = 1000;
            if(file_.Contains("Data")) rebinThr = 5000;

            for(int j=0;j<nBinCheck;++j){
                float nMassRT = hMassRT[i][j]->Integral();
                float nMassWT = hMassWT[i][j]->Integral();
                if(nMassRT < minEntries && nMassWT < minEntries ) continue;

                if(nMassRT<rebinThr) hMassRT[i][j]->Rebin();
                if(nMassWT<rebinThr) hMassWT[i][j]->Rebin();

                if(nMassRT >= minEntries) nMassRT = CountEventsWithFit(hMassRT[i][j], TString::Format("hMassRT%i%i", i, j));
                if(nMassWT >= minEntries) nMassWT = CountEventsWithFit(hMassWT[i][j], TString::Format("hMassWT%i%i", i, j));

                vX.push_back( wCalc[j] );
                vY.push_back( nMassWT/(nMassWT+nMassRT) );
                if(nMassWT == 0) vEY.push_back(1/nMassRT);
                else if(nMassRT == 0) vEY.push_back(1/nMassWT);
                else vEY.push_back(sqrt((nMassWT * nMassRT)/pow((nMassWT + nMassRT),3) ));

                cout<<"BIN "<<j<<", wCalc "<<wCalc[j]<<", nRT "<<nMassRT<<", nWT "<<nMassWT<<", wMeas "<<nMassWT/(nMassWT+nMassRT);
                cout<<" +- "<<sqrt((nMassWT * nMassRT)/pow((nMassRT),3) )<<endl;
            }

            auto *grW = new TGraphErrors(vX.size(),&vX[0],&vY[0],0,&vEY[0]);
            auto *fCheck = new TF1("fCheck","[0]+[1]*(x-[2])",0.,1.);
            fCheck->FixParameter(2, avgW);
            fCheck->SetParameter(0, avgW);
            fCheck->SetParameter(1, 1);

            grW->Fit("fCheck");
            auto *myfunc = grW->GetFunction("fCheck");

            vector<float> wResY;
            vector<float> wResEY;

            for (unsigned int j=0;j<vX.size();++j) { 
                wResY.push_back((vY[j] - myfunc->Eval(vX[j]))/vEY[j]);
                wResEY.push_back(1.);
            } 

            auto *grWres = new TGraphErrors(vX.size(),&vX[0],&wResY[0],0,&wResEY[0]);

            auto *c30 = new TCanvas();
            c30->Divide(1,2);
            c30->cd(1);
            grW->SetMarkerStyle(20);
            grW->SetMarkerSize(.5);
            grW->SetMaximum(1.);
            grW->SetMinimum(0);
            grW->Draw("APE");
            c30->cd(2);
            grWres->SetMarkerStyle(20);
            grWres->SetMarkerSize(.5);
            grWres->Draw("APE");
            auto *y0_ = new TF1("","0.",0.,1.);
            auto *y1_ = new TF1("","1.",0.,1.);
            auto *y1__ = new TF1("","-1.",0.,1.);
            auto *y2_ = new TF1("","2.",0.,1.);
            auto *y2__ = new TF1("","-2.",0.,1.);
            auto *y3_ = new TF1("","3.",0.,1.);
            auto *y3__ = new TF1("","-3.",0.,1.);
            y0_->SetLineColor(kBlack);
            y1_->SetLineColor(kGreen);
            y1__->SetLineColor(kGreen);
            y2_->SetLineColor(kOrange);
            y2__->SetLineColor(kOrange);
            y3_->SetLineColor(kRed);
            y3__->SetLineColor(kRed);
            y1_->SetLineWidth(1);
            y1__->SetLineWidth(1);
            y2_->SetLineWidth(1);
            y2__->SetLineWidth(1);
            y3_->SetLineWidth(1);
            y3__->SetLineWidth(1);
            y1_->SetLineStyle(2);
            y1__->SetLineStyle(2);
            y2_->SetLineStyle(2);
            y2__->SetLineStyle(2);
            y3_->SetLineStyle(2);
            y3__->SetLineStyle(2);
            y0_->Draw("SAME");
            y1_->Draw("SAME");
            y1__->Draw("SAME");
            y2_->Draw("SAME");
            y2__->Draw("SAME");
            y3_->Draw("SAME");
            y3__->Draw("SAME");

            c30->Print("validation" + process_ + HLT + TString::Format("_TYPE%i.pdf", i));

            float p0 = myfunc->GetParameter(0);
            float p1 = myfunc->GetParameter(0);

            cout<<myfunc->GetParameter(0)-myfunc->GetParameter(1)*avgW<<" +- "<<sqrt(pow(myfunc->GetParError(0),2) + pow(avgW,2)*pow(myfunc->GetParError(1),2) )<<endl<<endl;

        }

        auto *c32 = new TCanvas();
        evtMistag->DrawCopy("hist");
        c32->Print("perEvtMistagHist" + process_ + HLT + ".pdf");
    }

//----------CREATE----------
    if(mode_ == "CREATE")
    {
        //----------KDE----------
        sort(vKDERT.begin(), vKDERT.end()); //not really necessary
        sort(vKDEWT.begin(), vKDEWT.end());
        float xMin = min(vKDERT[0], vKDEWT[0]);
        float xMax = max(vKDERT[vKDERT.size()-1], vKDEWT[vKDEWT.size()-1]);
        float rho = .5;
        cout<<"rho "<<rho<<endl;
        cout<<"xMin "<<xMin<<endl;
        cout<<"xMax "<<xMax<<endl;

        auto *kdeRT = new TKDE(vKDERT.size(), &vKDERT[0], 0., 1., 
            "KernelType:Gaussian;Iteration:Adaptive;Mirror:NoMirror;Binning:RelaxedBinning", rho);
        auto *kdeWT = new TKDE(vKDEWT.size(), &vKDEWT[0], 0., 1., 
            "KernelType:Gaussian;Iteration:Adaptive;Mirror:NoMirror;Binning:RelaxedBinning", rho);

        //plotting
        auto *pdfRT = new TF1("pdfRT",kdeRT,0.,1.,0);
        auto *pdpdfWT = new TF1("pdpdfWT",kdeWT,0.,1.,0);

        pdfRT->SetLineColor(kBlue);
        pdfRT->SetNpx(2000);
        pdpdfWT->SetLineColor(kRed);
        pdpdfWT->SetNpx(2000);

        auto *c10 = new TCanvas();
        c10->Divide(2,2);
        c10->cd(1);
        mva_RT->DrawCopy("hist");
        c10->cd(2);
        mva_WT->DrawCopy("hist");
        c10->cd(3);
        pdfRT->Draw();
        c10->cd(4);
        pdpdfWT->Draw();

        //Lambda functions, NEED all the objects called inside
        auto pdfW = new TF1("pdfW",
            [&](double *x, double *p)
            { 
                double yRT = vKDERT.size()*kdeRT->GetValue(x[0]); 
                double yWT = vKDEWT.size()*kdeWT->GetValue(x[0]); 
                return yWT/(yWT+yRT);
            }, 
            0.,1.,0);

        //this function force pdfW to be monotonically decreasing/constant near the boundaries
        auto pdfW_extended = new TF1("pdfW_extended",
            [&](double *x, double *p)
            {
                if(x[0]<pdfW->GetMaximumX())
                    return pdfW->GetMaximum();
                if(x[0]>=pdfW->GetMinimumX())
                    return pdfW->GetMinimum();
                return pdfW->Eval(x[0]);
            }, 
            0.,1.,0);

        pdfW->SetNpx(100); //number of points to build the TGraph 
        pdfW_extended->SetNpx(100);
        auto *out_pdfW = new TGraph(pdfW);
        auto *out_pdfW_extended = new TGraph(pdfW_extended);

        auto *c16 = new TCanvas();    
        pdfW_extended->SetMarkerStyle(20);
        pdfW_extended->SetLineColor(kBlue);
        pdfW_extended->SetMarkerSize(.5);
        pdfW_extended->SetNpx(100);
        pdfW_extended->DrawClone("PL");
        pdfW->SetMarkerStyle(20);
        pdfW->SetMarkerSize(1.);
        pdfW->SetNpx(100);
        pdfW->DrawClone("P SAME");

        //----------CATEGORIES----------
        nCat = 25;
        int nTotTagged = nRT + nWT;
        int catSize = nTotTagged / nCat;
        cout<<endl<<"nTotTagged "<<nTotTagged<<endl<<"catSize "<<catSize<<endl<<endl;

        int   cTot    = 0;
        float cCenter = 0;
        int   lastCat = 0;
        
        float   *catEdge    = new float[nCat];
        float   *catCenter  = new float[nCat];
        int     *catRT      = new int[nCat];
        int     *catWT      = new int[nCat];
        float   *catEff     = new float[nCat];
        float   *catW       = new float[nCat];
        float   *catP       = new float[nCat];

        float   *vexl = new float[nCat];
        float   *vexh = new float[nCat];
        float   *vey = new float[nCat];

        for(int i=0; i<nCat; ++i){
            catRT[i] = 0;
            catWT[i] = 0;
            catCenter[i] = 0;
        }

        int cCat = 0;
        for(int i = 1; i<=nBinsMva; ++i){
            cCenter += mva->GetBinCenter(i)*mva->GetBinContent(i);
            catRT[cCat] += mva_RT->GetBinContent(i);
            catWT[cCat] += mva_WT->GetBinContent(i);
            cTot        += mva->GetBinContent(i);
            //cout<<"bin "<<i<<", "<<mva_RT->GetBinContent(i) + mva_WT->GetBinContent(i);
            //cout<<" ("<<cTot<<")"<<" ["<<hTot<<"]"<<endl;
            if(cTot >= catSize){
                catEdge[cCat] = mva->GetXaxis()->GetBinLowEdge(i+1);
                catCenter[cCat] = (float)cCenter/(float)cTot;
                cout<<" -----> cat "<<cCat<<", "<<cTot<<" [-"<<catEdge[cCat]<<"]"<<endl;
                cTot = 0;
                cCenter = 0;
                cCat++;
            }
            if(i==nBinsMva){
                catCenter[cCat] = (float)cCenter/(float)cTot;
                catEdge[cCat] = 1.;
                lastCat = cCat;
                cout<<" -----> cat "<<cCat<<", "<<cTot<<" ["<<catEdge[cCat]<<"]"<<endl;
            }
        }

        cout<<endl;

        nCat = lastCat + 1;
        float totEff = 0;
        float totW = 0;
        float avgD = 0;
        float totP = 0;

        for(int i=0; i<nCat; ++i)
        {
            catEff[i] = (float)(catWT[i] + catRT[i]) / (float)nEventsRead;
            totEff += catEff[i];

            catW[i] = (float)catWT[i] / (float)(catWT[i] + catRT[i]);
            avgD += catEff[i]*(1-2*catW[i]);

            catP[i] = catEff[i]*pow(1-2*catW[i],2);
            totP += catP[i];

            vexh[i] = catEdge[i] - catCenter[i];
            if(i!=0) vexl[i] = catCenter[i] - catEdge[i-1];
            else     vexl[i] = catCenter[i] - 0.;

            vey[i] = sqrt(((float)catWT[i] * (float)catRT[i])/pow(((float)catWT[i] + (float)catRT[i]),3) );

            cout<<"CAT "<<i+1;
            if(i!=0) cout<<" ["<<catEdge[i-1]<<" - "<<catEdge[i]<<" ] ";
            else     cout<<" ["<<0.<<" - "<<catEdge[i]<<" ] ";
            cout<<" - eff = "<<100*catEff[i]<<"%, w = "<<100*catW[i];
            cout<<"%, P = "<<100*catP[i]<<"%"<<endl;
            //cout<<" - "<<catCenter[i]<<" - ";
            //cout<<catRT[i] + catWT[i]<<endl;
        }

        avgD /= totEff;
        avgW = (1-avgD)/2;

        //FUNCTION FIT
        TGraphAsymmErrors *grAsy = new TGraphAsymmErrors(nCat,catCenter,catW,vexl,vexh,vey,vey);
        auto *gr = new TGraph(nCat,catCenter,catW);

        perEvtWFormula = "[0]+[1]*(1-TMath::Erf([2]+[3]*x))";
        auto *fitErf = new TF1("fitErf", perEvtWFormula, 0., 1.);
        fitErf->SetParameter(0, catW[nCat-1]);
        fitErf->SetParameter(1, (catW[0]-catW[nCat-1])/2);
        fitErf->SetParLimits(0, 0., 0.3);
        fitErf->SetParLimits(1, 0.5/2, 1./2);

        auto *fitTanH = new TF1("fitTanH", "[0]+[1]*(1-TMath::TanH([2]+[3]*x))", 0., 1.);
        auto *fitATan = new TF1("fitATan", "[0]+[1]*(1-TMath::ATan([2]+[3]*x))", 0., 1.);

        gr->Fit("fitErf","MELR");

        //----------OUTPUT STREAM----------
        ofstream ofs;
        ofs.open ("OSMuonTagger" + HLT + "Categories.txt");
        //CAT
        ofs<<method_<<endl;
        ofs<<nCat<<endl;
        for(int i=0; i<nCat; ++i)
        {
            if(i!=0) ofs<<catEdge[i-1];
            else     ofs<<0.;
            ofs<<" "<<catEdge[i]<<" ";
            ofs<<catW[i]<<endl;
        }
        //FUNCTION
        ofs<<perEvtWFormula<<endl;
        ofs<<avgW<<endl;
        ofs<<fitErf->GetNumberFreeParameters()<<endl;
        for(int j=0;j<fitErf->GetNumberFreeParameters();++j)
            ofs<<fitErf->GetParameter(j)<<endl;
        ofs.close();

        auto *fo = new TFile("OSMuonTagger" + HLT + "KDE.root", "RECREATE");
        fo->cd();
        out_pdfW->Write();
        out_pdfW_extended->Write();
        fo->Close();
        f->cd();

        //----------PRINT----------
        cout<<endl<<"Cat Eff = "<<100*totEff<<"%"<<endl;
        cout<<"Cat W = "<<100*avgW<<"%"<<endl;
        cout<<"Cat P = "<<100*totP<<"%"<<endl;

        auto *c2 = new TCanvas();
        if(mva_WT->GetMaximum()>mva_RT->GetMaximum()){
            mva_WT->DrawCopy("hist");
            mva_RT->DrawCopy("hist same");            
        }else{
            mva_RT->DrawCopy("hist");
            mva_WT->DrawCopy("hist same");                   
        }


        auto *c3 = new TCanvas();
        pdfW_extended->DrawClone("");
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(.5);
        gr->SetTitle("per-event mistag");
        gr->GetXaxis()->SetTitle("MVA output");
        gr->SetMinimum(0.);
        gr->Draw("PE SAME");
        grAsy->Draw("E SAME");
        pdfW->SetLineColor(kGreen);
        pdfW->DrawClone("P SAME");

        c3->Update();

        c10->Print("kdeDistributions" + process_ + HLT + ".pdf");
        c2->Print("mvaDistributions" + process_ + HLT + ".pdf");
        c3->Print("perEventW" + process_ + HLT + ".pdf");
    }

    cout<<endl<<"NoCat Eff = "<<100.*(float)(nRT+nWT)/nEventsRead<<"%"<<endl;
    cout<<"NoCat W = "<<100.*(float)nWT/(nRT+nWT)<<"%"<<endl;
    cout<<"NoCat P = "<<100.*((float)(nRT+nWT)/nEventsRead)*pow(1.-2.*((float)nWT/(nRT+nWT)),2)<<"%"<<endl;

    cout<<endl<<"Cat P = "<<100.*totPCat/(float)nEventsRead<<"%"<<endl;
    cout<<"KDE P = "<<100.*totPKde/(float)nEventsRead<<"%"<<endl;

    f->Close();
    delete f;
    return;

}

float CountEventsWithFit(TH1 *hist, TString name = "hist"){

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );

    float mean = 5.3663;
    float sigma = 0.015;
    if(process_.Contains("BsJPsiPhi")) mean = 5.3663;
    if(process_.Contains("BuJPsiK"))   mean = 5.2793;
    if(process_.Contains("BdJPsiKx"))  mean = 5.2796;

    TString sgnDef = "[1]*TMath::Gaus(x, [0], [4], true)";
    sgnDef +=       "+[2]*TMath::Gaus(x, [0], [5], true)";
    sgnDef +=       "+[3]*TMath::Gaus(x, [0], [6], true)";

    TString bkgDef;

    if(process_.Contains("MC")) bkgDef = "[7]+[8]*x";
    if(process_.Contains("Data")) bkgDef = "[7]+[8]*TMath::Erfc([9]*(x-[10]))";

    TString funcDef = sgnDef + "+" + bkgDef;

    TF1 *func = new TF1("func", funcDef, min_, max_);

    float limit = hist->GetEntries()*hist->GetBinWidth(1);

    float bkgHeight = 0;
    for(int i=0; i<5; i++){
        bkgHeight += hist->GetBinContent(i+i);
        bkgHeight += hist->GetBinContent(hist->GetNbinsX()-i);
    }

    bkgHeight /= 10;

    func->SetParameter(0, mean);
    func->SetParameter(1, limit/3);
    func->SetParameter(2, limit/3);
    func->SetParameter(3, limit/3);
    func->SetParameter(4, sigma);
    func->SetParameter(5, sigma);
    func->SetParameter(6, sigma);
    func->SetParameter(7, bkgHeight);

    func->SetParLimits(0, mean-sigma, mean+sigma);
    func->SetParLimits(1, 0, hist->GetEntries());
    func->SetParLimits(2, 0, hist->GetEntries());
    func->SetParLimits(3, 0, hist->GetEntries());
    func->SetParLimits(4, sigma/2, sigma*2);
    func->SetParLimits(5, sigma/2, sigma*2);
    func->SetParLimits(6, sigma/2, sigma*2);
    func->SetParLimits(7, bkgHeight/2, bkgHeight*2);

    if(process_.Contains("MC")){
        func->SetParameter(8, 0);
    }

    if(process_.Contains("Data")){
        func->SetParameter(7, bkgHeight);
        func->SetParameter(8, 1);
        func->SetParameter(9, 20);
        func->SetParameter(10, 5.10);
        func->SetParLimits(8, 0, bkgHeight);
        func->SetParLimits(9, 10, 1e3);
        func->SetParLimits(10, 5.0, mean);
    }

    func->SetNpx(2000);

    auto c5 = new TCanvas();
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(.75);
    hist->Fit("func","MRLQ");
    hist->Draw("PE");
    hist->SetMinimum(0.1);

    TF1 *fit = hist->GetFunction("func");
    TF1 *f1 = new TF1("f1","[0]*TMath::Gaus(x, [1], [2], true)", min_, max_);
    TF1 *f2 = new TF1("f2","[0]*TMath::Gaus(x, [1], [2], true)", min_, max_);
    TF1 *f3 = new TF1("f3","[0]*TMath::Gaus(x, [1], [2], true)", min_, max_);
    TF1 *f4;
    TF1 *f5;

    f1->SetParameters(fit->GetParameter(1),fit->GetParameter(0),fit->GetParameter(4));
    f2->SetParameters(fit->GetParameter(2),fit->GetParameter(0),fit->GetParameter(5));
    f3->SetParameters(fit->GetParameter(3),fit->GetParameter(0),fit->GetParameter(6));

    if(process_.Contains("MC")){
        f4 = new TF1("f4","[0]+[1]*x", min_, max_);
        f4->SetParameter(fit->GetParameter(7), fit->GetParameter(8));
    }
    if(process_.Contains("Data")){
        f4 = new TF1("f4","[0]", min_, max_);
        f5 = new TF1("f5","[0]*TMath::Erfc([1]*(x-[2]))", min_, max_);
        f4->SetParameter(0, fit->GetParameter(7));
        f5->SetParameters(fit->GetParameter(8),fit->GetParameter(9),fit->GetParameter(10));
    }

    f1->SetLineColor(kBlue);
    f2->SetLineColor(kViolet);
    f3->SetLineColor(kAzure);
    f4->SetLineColor(kOrange);
    if(process_.Contains("Data")) f5->SetLineColor(kGreen);
    f1->SetLineStyle(2);
    f2->SetLineStyle(2);
    f3->SetLineStyle(2);
    f4->SetLineStyle(2);
    if(process_.Contains("Data")) f5->SetLineStyle(2);

    f1->Draw("same");
    f2->Draw("same");
    f3->Draw("same");
    f4->Draw("same");
    if(process_.Contains("Data")) f5->Draw("same");
    c5->Print(name + ".pdf");

    float nEvt = fit->GetParameter(1);
    nEvt += fit->GetParameter(2);
    nEvt += fit->GetParameter(3);
    nEvt/=hist->GetBinWidth(1);

    return nEvt;

}