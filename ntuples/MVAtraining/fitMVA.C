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
TString dirPath_;
float min_, max_;
int nBins_ = 50;

pair<float, float> CountEventsWithFit(TH1 *hist, TString name);

void fitMVA(TString file_ = "ntuBsMC2017.root"
            , TString method_ = "DNNOsMuonHLTJpsiMu_test241"
            , TString mode_ = "CREATE"
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
    cout<<"nEvents_ "<<nEvents_<<endl;

    if(mode_ != "CREATE" && mode_ != "USE"){
        cout<<"WRONG MODE_"<<endl;
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
    }

    cout<<"HLT "<<HLT<<endl;

    int nCat;
    float *catEdgeL;
    float *catEdgeR;
    float *catMistag;
    float *catMistagErr;
    
    TF1 *perEvtW;
    TString perEvtWFormula;
    float avgW;

    float totPCut = 0.;
    float totPCat = 0.;
    float totPFit = 0.;
    float totPKde = 0.;

    TGraph *g_pdfW;
    TGraph *g_pdfW_extended;

    float evtWcat[2] = {-1, -1};
    float evtWfit[2] = {-1, -1};
    float evtWkde[2] = {-1, -1};
    float evtWkde_ext[2] = {-1, -1};

    int nBinCheck = 30;
    TH1F **hMassRT[4];
    TH1F **hMassWT[4];

    float **wCalc = new float*[4];
    for(int i=0; i<4; ++i){
        hMassRT[i] = new TH1F*[nBinCheck];
        hMassWT[i] = new TH1F*[nBinCheck];
        wCalc[i] = new float[nBinCheck];
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
        catEdgeL     = new float[nCat];
        catEdgeR     = new float[nCat];
        catMistag    = new float[nCat];
        catMistagErr = new float[nCat];
        for(int i=0; i<nCat; ++i){
            ifs >> catEdgeL[i];
            ifs >> catEdgeR[i];
            ifs >> catMistag[i];
            ifs >> catMistagErr[i];
        }
        cout<<"Input categories [L-edge R-edge mistag error]"<<endl;
        for(int i=0; i<nCat; ++i)
            cout<<catEdgeL[i]<<" "<<catEdgeR[i]<<" "<<catMistag[i]<<" "<<catMistagErr[i]<<endl;

        ifs.close();

        //FIT
        auto *f2 = new TFile("OSMuonTagger" + HLT + "Functions.root");
        f2->cd();
        perEvtW = (TF1*)f2->Get("fitErf")
        g_pdfW = (TGraph*)f2 ->Get("pdfW");
        g_pdfW_extended = (TGraph*)f2 ->Get("pdfW_extended");
        f2->Close();
        delete f2;
        f->cd();

        for(int j=0;j<nBinCheck;++j){
            for(int i=0; i<4; ++i){
                hMassRT[i][j] = new TH1F(TString::Format("mbRT%i%i", i, j),"", nBins_, min_, max_);
                hMassWT[i][j] = new TH1F(TString::Format("mbWT%i%i", i, j),"", nBins_, min_, max_);
                wCalc[i][j] = 0.;
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
    auto *htotMass    = new TH1F( "htotMass","htotMass", nBins_, min_, max_ );
    auto *htotMass_RT = new TH1F( "htotMass_RT","htotMass_RT", nBins_, min_, max_ );
    auto *htotMass_WT = new TH1F( "htotMass_WT","htotMass_WT", nBins_, min_, max_ );
    auto *htotMass_Tag = new TH1F( "htotMass_Tag","htotMass_Tag", nBins_, min_, max_ );

    auto *evtMistagSignal = new TH1F( "evtMistagSignal", "evtMistagSignal", 1000, 0.0, 1.0 );
    auto *evtMistagSignal_bis = new TH1F( "evtMistagSignal_bis", "evtMistagSignal_bis", 1000, 0.0, 1.0 );
    auto *evtMistagSignal_Cat = new TH1F( "evtMistagSignal_Cat", "evtMistagSignal_Cat", 1000, 0.0, 1.0 );
    auto *evtMistagSignal_Fit = new TH1F( "evtMistagSignal_Fit", "evtMistagSignal_Fit", 1000, 0.0, 1.0 );
    auto *evtMistagSide = new TH1F( "evtMistagSide", "evtMistagSide", 1000, 0.0, 1.0 );
    float nSignal=0;
    float nSide=0;    

    mva_WT->SetLineColor(kRed);

    vector<double> vKDERT;
    vector<double> vKDEWT;

    //EVENT LOOP
    cout<<"----- BEGIN LOOP"<<endl;

    if(nEvents_ == -1) nEvents_ = t->GetEntries();
    int nEvtsRead = 0;
    for(int i=0; i<nEvents_; ++i){
        if(i%100000==0) cout<<"----- at event "<<i<<endl;
        t->GetEntry(i);

        //PRESELECTION EVENT
        if(HLT == "HltJpsiMu" && !hltJpsiMu) continue;
        if(HLT == "HltJpsiTrkTrk" && !hltJpsiTrkTrk) continue;
        if(HLT == "HltJpsiTrk" && !hltJpsiTrk) continue;
        if(useTightSelection && !ssbIsTight) continue;

        nEvtsRead = nEvtsRead + (int)evtWeight;
        htotMass->Fill(ssbMass, evtWeight);

        //PRESELECTION MUON
        if(!osMuon) continue;
        if((fabs(muoEta)<1.2 && muoSoftMvaValue<=0.891)) continue;
        if((fabs(muoEta)>=1.2 && muoSoftMvaValue<=0.8925)) continue;

        //CHECK NAN INF VARIABLES, THIS IS A BUG FIX
        bool nanFlag = isnan(muoDxy) || isnan(muoJetDFprob);
        bool infFlag = isinf(muoJetEnergyRatio) || isinf(muoConeEnergyRatio);
        if(nanFlag || infFlag) continue;

        //SWITCH VARIABLE COMPUTATION
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
                evtWcat[0] = catMistag[0];
                evtWcat[1] = catMistagErr[0];
                evtWfit[0] = perEvtW->Eval(catEdgeL[0]);
            }else if( mvaValue >= catEdgeR[nCat-1] ){ //EVENTS WITH MVA>=MVA_MAX USED IN CREATING THE METHOD
                cout<<"-----Overcat "<<mvaValue<<" at "<<i<<endl;
                evtWcat[0] = catMistag[nCat-1];
                evtWcat[1] = catMistagErr[nCat-1];
                evtWfit[0] = perEvtW->Eval(catEdgeR[nCat-1]);
            }else{
                for(int j=0; j<nCat; ++j){ //EVENTS WITH MVA IN [MVA_MIN, MVA_MAX] USED IN CREATING THE METHOD
                    if(( mvaValue >= catEdgeL[j]) && (mvaValue < catEdgeR[j]) )
                    { 
                        evtWcat[0] = catMistag[j]; //CAT
                        evtWcat[1] = catMistagErr[j];
                        break; 
                    }
                }
                evtWfit[0] = perEvtW->Eval(mvaValue);  //FIT
            }
            evtWkde[0] = g_pdfW->Eval(mvaValue, 0, "S");   //KDE
            evtWkde_ext[0] = g_pdfW_extended->Eval(mvaValue, 0, "S");  //KDE EXTENDED

            totPCat += pow(1.-2.*evtWcat[0], 2)*evtWeight;
            totPKde += pow(1.-2.*evtWkde[0] ,2)*evtWeight;

            if( (ssbMass>5.18) && (ssbMass<5.37) ){
                evtMistagSignal->Fill(evtWkde_ext[0], evtWeight);
                evtMistagSignal_bis->Fill(evtWkde[0], evtWeight);
                evtMistagSignal_Cat->Fill(evtWcat[0], evtWeight);
                evtMistagSignal_Fit->Fill(evtWfit[0], evtWeight);
            }
            else evtMistagSide->Fill(evtWkde_ext[0], evtWeight);


            for(int j=0;j<nBinCheck;++j){
                if( (evtWcat[0]>=(float)j*pass) && (evtWcat[0]<((float)j*pass+pass)) ){
                    if(TMath::Sign(1, ssbLund) == evtTag) hMassRT[0][j]->Fill(ssbMass, evtWeight);
                    if(TMath::Sign(1, ssbLund) != evtTag) hMassWT[0][j]->Fill(ssbMass, evtWeight);
                    wCalc[0][j] += evtWcat[0]*evtWeight;
                    break;
                }
            }
            for(int j=0;j<nBinCheck;++j){
                if( (evtWfit[0]>=(float)j*pass) && (evtWfit[0]<((float)j*pass+pass)) ){
                    if(TMath::Sign(1, ssbLund) == evtTag) hMassRT[1][j]->Fill(ssbMass, evtWeight);
                    if(TMath::Sign(1, ssbLund) != evtTag) hMassWT[1][j]->Fill(ssbMass, evtWeight);
                    wCalc[1][j] += evtWfit[0]*evtWeight;
                    break;
                }
            }
            for(int j=0;j<nBinCheck;++j){             
                if( (evtWkde[0]>=(float)j*pass) && (evtWkde[0]<((float)j*pass+pass)) ){
                    if(TMath::Sign(1, ssbLund) == evtTag) hMassRT[2][j]->Fill(ssbMass, evtWeight);
                    if(TMath::Sign(1, ssbLund) != evtTag) hMassWT[2][j]->Fill(ssbMass, evtWeight);
                    wCalc[2][j] += evtWkde[0]*evtWeight;
                    break;
                }
            }
            for(int j=0;j<nBinCheck;++j){
                if( (evtWkde_ext[0]>=(float)j*pass) && (evtWkde_ext[0]<((float)j*pass+pass)) ){
                    if(TMath::Sign(1, ssbLund) == evtTag) hMassRT[3][j]->Fill(ssbMass, evtWeight);
                    if(TMath::Sign(1, ssbLund) != evtTag) hMassWT[3][j]->Fill(ssbMass, evtWeight);
                    wCalc[3][j] += evtWkde_ext[0]*evtWeight;
                }
            }

        }

        mva->Fill(mvaValue, evtWeight);
        htotMass_Tag->Fill(ssbMass, evtWeight);
        if(osMuonTag == 1){
            mva_RT->Fill(mvaValue, evtWeight);
            vKDERT.push_back(mvaValue);
            if(evtWeight==2) vKDERT.push_back(mvaValue); //avoid using weights in kde
            htotMass_RT->Fill(ssbMass, evtWeight);
        }
        if(osMuonTag == 0){
            mva_WT->Fill(mvaValue, evtWeight);
            vKDEWT.push_back(mvaValue);
            if(evtWeight==2) vKDEWT.push_back(mvaValue);
            htotMass_WT->Fill(ssbMass, evtWeight);
        }
    }

    int nRT = mva_RT->Integral(); //Integral() take in consideration weight, GetEntries() not
    int nWT = mva_WT->Integral();

    if(file_.Contains("Data")){ //for data fit mass
        nRT = CountEventsWithFit(htotMass_RT, "histTotMass_RT").first;
        nWT = CountEventsWithFit(htotMass_WT, "histTotMass_WT").first;
        nEvtsRead = CountEventsWithFit(htotMass, "histTotMass").first;
    }

    cout<<"----- MVA HISTOGRAMS FILLED"<<endl;
    cout<<"nRT "<<nRT<<endl;
    cout<<"nWT "<<nWT<<endl;

    cout<<endl<<"Base Eff = "<<100.*(float)(nRT+nWT)/nEvtsRead<<"%"<<endl;
    cout<<"Base W = "<<100.*(float)nWT/(nRT+nWT)<<"%"<<endl;
    cout<<"Base P = "<<100.*((float)(nRT+nWT)/nEvtsRead)*pow(1.-2.*((float)nWT/(nRT+nWT)),2)<<"%"<<endl;

    if(file_.Contains("MC")){
        cout<<endl<<"Cat P = "<<100.*totPCat/(float)nEvtsRead<<"%"<<endl;
        cout<<"KDE P = "<<100.*totPKde/(float)nEvtsRead<<"%"<<endl;
    }else if(file_.Contains("Data")){
        totPCat = 0.;
        totPKde = 0.;
    }

//----------VALIDATION OF THE METHOD----------
    if(mode_ == "USE"){
        for(int i=0; i<4; ++i){
            for(int j=0;j<nBinCheck;++j){
                wCalc[i][j] /= (hMassRT[i][j]->Integral() + hMassWT[i][j]->Integral());
            }
        }

        for(int i=0; i<4; ++i){
            cout<<"TYPE "<<i<<endl;

            vector<float> vX;
            vector<float> vY;
            vector<float> vEX;
            vector<float> vEY;

            int minEntries = 0;
            int rebinThr = 5000;

            for(int j=0;j<nBinCheck;++j){
                pair<float, float> nMassRT; // .first = nEvt; .second = sigma(nEvt)
                pair<float, float> nMassWT;

                nMassRT.first = hMassRT[i][j]->Integral();
                nMassWT.first = hMassWT[i][j]->Integral();
                nMassRT.second = sqrt(nMassRT.first);
                nMassWT.second = sqrt(nMassWT.first);

                if(nMassRT.first <= minEntries && nMassWT.first <= minEntries ) continue;

                if(process_.Contains("Data")){
                    if(nMassRT.first<rebinThr) hMassRT[i][j]->Rebin();
                    if(nMassWT.first<rebinThr) hMassWT[i][j]->Rebin();
                    if(nMassRT.first >= minEntries) nMassRT = CountEventsWithFit(hMassRT[i][j], TString::Format("hMassRT%i%i", i, j));
                    if(nMassWT.first >= minEntries) nMassWT = CountEventsWithFit(hMassWT[i][j], TString::Format("hMassWT%i%i", i, j));
                }

                vX.push_back( wCalc[i][j] );
                vY.push_back( nMassWT.first/(nMassWT.first + nMassRT.first) ); // MEASURED MISTAG

                float errMistag, errRT, errWT;

                if(nMassWT.first == 0){
                    errWT = 1;
                    errRT = nMassRT.second;
                }
                else if(nMassRT.first == 0){
                    errRT = 1;
                    errWT = nMassWT.second;
                }else{
                    errRT = nMassRT.second;
                    errWT = nMassWT.second;
                }

                errMistag =  sqrt(pow(nMassWT.first,2)*pow(errRT,2) + pow(nMassRT.first,2)*pow(errWT,2)  )/pow(nMassRT.first + nMassWT.first, 2);
                vEY.push_back(errMistag);

                if(file_.Contains("Data") && i==0){ //category tagging performance evaluation
                    totPCat += (nMassRT.first + nMassWT.first)*pow(1-2*wCalc[i][j],2);
                }

                if(file_.Contains("Data") && i==3){ //KDE tagging performance evaluation
                    totPKde += (nMassRT.first + nMassWT.first)*pow(1-2*wCalc[i][j],2);
                }

                cout<<"BIN "<<j<<", wCalc "<<wCalc[i][j]<<", nRT "<<nMassRT.first<<"+-"<<errRT<<", nWT "<<nMassWT.first<<"+-"<<errWT;
                cout<<", wMeas "<<nMassWT.first/(nMassWT.first + nMassRT.first);
                cout<<" +- "<<errMistag<<endl;

            }

            auto *grW = new TGraphErrors(vX.size(),&vX[0],&vY[0],0,&vEY[0]);
            //auto *fCheck = new TF1("fCheck","[0]+[1]*(x-[2])",0.,1.); //alternative form from previous analysis
            //fCheck->FixParameter(2, avgW);
            auto *fCheck = new TF1("fCheck","[0]+[1]*x",0.,1.);
            fCheck->SetParameters(avgW, 1);

            grW->Fit("fCheck");
            auto *myfunc = grW->GetFunction("fCheck");
            auto *fitDraw = new TF1("fitDraw","[0]+[1]*x",0.,.5);
            fitDraw->SetParameters(myfunc->GetParameter(0), myfunc->GetParameter(1));

            vector<float> wResY;
            vector<float> wResEY;

            for (unsigned int j=0;j<vX.size();++j) { 
                wResY.push_back((vY[j] - myfunc->Eval(vX[j]))/vEY[j]);
                wResEY.push_back(1.);
            } 

            auto *grWres = new TGraphErrors(vX.size(),&vX[0],&wResY[0],0,&wResEY[0]);

            auto *c30 = new TCanvas("c30","c30",1000,1600);
            c30->Divide(1,2);
            c30->cd(1);
            grW->SetMarkerStyle(20);
            grW->SetMarkerSize(1);
            grW->SetMaximum(1.);
            grW->SetMinimum(0);
            grW->GetXaxis()->SetLimits(0.0,0.9);
            grW->SetTitle("");
            grW->GetXaxis()->SetTitle("w calc.");
            grW->GetYaxis()->SetTitle("w meas.");

            grW->Draw("APZ");
            fitDraw->Draw("same");
            c30->cd(2);
            grWres->SetMarkerStyle(20);
            grWres->SetMarkerSize(1);
            grWres->GetXaxis()->SetLimits(0.0,0.9);
            grWres->Draw("APZ");
            grWres->SetTitle("");
            grWres->GetYaxis()->SetTitle("# std dev");
            auto *y0_ = new TF1("","0.",0.,1.);
            y0_->SetLineColor(kBlack);
            y0_->Draw("SAME");

            c30->Print(dirPath_ + "/" + "validation" + process_ + HLT + TString::Format("_TYPE%i.pdf", i));

            float p0 = myfunc->GetParameter(0);
            float p1 = myfunc->GetParameter(0);
        }

        auto *c32 = new TCanvas();
        c32->Divide(2,2);
        c32->cd(1);
        g_pdfW_extended->Draw("APC");
        c32->cd(2);
        evtMistagSignal_Fit->DrawCopy("hist");
        c32->cd(3);
        evtMistagSignal_bis->DrawCopy("hist");
        c32->cd(4);
        evtMistagSignal->DrawCopy("hist");
        c32->Print(dirPath_ + "/" + "perEvtMistagHist" + process_ + HLT + ".pdf");
    }

    if(file_.Contains("Data")){
        cout<<endl<<"Cat P = "<<100.*totPCat/(float)nEvtsRead<<"%"<<endl;
        cout<<"KDE P = "<<100.*totPKde/(float)nEvtsRead<<"%"<<endl;
    }
//----------CREATE----------
    if(mode_ == "CREATE")
    {
        //----------KDE----------
        sort(vKDERT.begin(), vKDERT.end()); //not really necessary
        sort(vKDEWT.begin(), vKDEWT.end());
        float xMin = min(vKDERT[0], vKDEWT[0]);
        float xMax = max(vKDERT[vKDERT.size()-1], vKDEWT[vKDEWT.size()-1]);
        float rho = 1.;
        cout<<"rho "<<rho<<endl;
        cout<<"xMin "<<xMin<<endl;
        cout<<"xMax "<<xMax<<endl;

        auto *kdeRT = new TKDE(vKDERT.size(), &vKDERT[0], 0., 1., 
            "KernelType:Gaussian;Iteration:Adaptive;Mirror:NoMirror;Binning:RelaxedBinning", rho);
        auto *kdeWT = new TKDE(vKDEWT.size(), &vKDEWT[0], 0., 1., 
            "KernelType:Gaussian;Iteration:Adaptive;Mirror:NoMirror;Binning:RelaxedBinning", rho);

        cout<<"----- KDE COMPUTED"<<endl;

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

        cout<<"----- KDE LAMBDA FUNCTIONS COMPUTED"<<endl;

        int nPoints = 300;//number of points to build the TGraph 

        pdfW->SetNpx(nPoints); 
        pdfW_extended->SetNpx(nPoints);
        auto *out_pdfW = new TGraph(pdfW);
        auto *out_pdfW_extended = new TGraph(pdfW_extended);

        auto *c16 = new TCanvas();    
        pdfW_extended->SetMarkerStyle(20);
        pdfW_extended->SetLineColor(kBlue);
        pdfW_extended->SetMarkerSize(.5);
        pdfW_extended->DrawClone("PL");
        pdfW->SetMarkerStyle(20);
        pdfW->SetMarkerSize(1.);
        pdfW->DrawClone("P SAME");

        cout<<"----- KDE TGRAPHS COMPUTED"<<endl;

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
        float   *catWerr    = new float[nCat];
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
            catEff[i] = (float)(catWT[i] + catRT[i]) / (float)nEvtsRead;
            totEff += catEff[i];

            catW[i] = (float)catWT[i] / (float)(catWT[i] + catRT[i]);
            catWerr[i] = sqrt(pow((float)catWT[i],2)*pow(sqrt((float)catRT[i]),2) + pow((float)catRT[i],2)*pow(sqrt((float)catWT[i]),2)  )/pow((float)catRT[i] + (float)catWT[i], 2);
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
            ofs<<catW[i]<<" "<<catWerr[i]<<endl;
        }

        ofs.close();

        //FUNCTIONS
        auto *fo = new TFile("OSMuonTagger" + HLT + "Functions.root", "RECREATE");
        fo->cd();
        out_pdfW->Write();
        out_pdfW_extended->Write();
        fitErf->Write();
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
        gr->GetYXaxis()->SetTitle("per-event mistag");
        gr->GetXaxis()->SetTitle("MVA output");
        gr->SetMinimum(0.);
        gr->Draw("APE SAME");
        grAsy->Draw("E SAME");
        pdfW->SetLineColor(kGreen);
        pdfW->DrawClone("P SAME");

        c3->Update();

        c10->Print("kdeDistributions" + process_ + HLT + ".pdf");
        c2->Print("mvaDistributions" + process_ + HLT + ".pdf");
        c3->Print("perEventW" + process_ + HLT + ".pdf");
    }

    f->Close();
    delete f;
    return;

}

pair<float, float> CountEventsWithFit(TH1 *hist, TString name = "hist"){

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );

    float mean = 5.3663;
    float sigma = 0.015;
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
    float limit = hist->GetEntries()*hist->GetBinWidth(1);

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
    float bkgHeight = 0.;
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

    float nEvt = fit->GetParameter(1);
    nEvt += fit->GetParameter(2);
    nEvt += fit->GetParameter(3);
    nEvt/=hist->GetBinWidth(1);

    TMatrixDSym cov = r->GetCovarianceMatrix();
    float errN = sqrt(cov(1,1)+cov(2,2)+cov(3,3)+2*(cov(1,2)+cov(1,3)+cov(2,3)))/hist->GetBinWidth(1);

    return make_pair(nEvt, errN);
}