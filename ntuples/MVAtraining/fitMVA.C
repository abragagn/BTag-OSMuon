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

int fitMVA(TString file_ = "../ntuBsMC2017.root"
            , TString method_ = "DNNOsMuon2017test231"
            , TString mode_ = "CREATE"
            , bool addMva_ = false
            , bool readMva_ = false
            , int nEvents_ = -1)
{
    cout<<"file_ "<<file_<<endl;
    cout<<"method_ "<<method_<<endl;
    cout<<"mode_ "<<mode_<<endl;
    cout<<"addMva_ "<<addMva_<<endl;
    cout<<"readMva_ "<<readMva_<<endl;
    cout<<"nEvents_ "<<nEvents_<<endl;


    gErrorIgnoreLevel = kWarning;
    if(mode_ != "CREATE" && mode_ != "ADD"){
        cout<<"WRONG MODE_"<<endl;
        return 0;
    }

    if(addMva_ && readMva_){
        cout<<"CAN'T READ AND WRITE MVA AT THE SAME TIME"<<endl;
        return 0;
    }

    TString fileopt = "READ";
    if(mode_ == "ADD") fileopt = "UPDATE";
    if(mode_ == "CREATE" && addMva_) fileopt = "UPDATE";

    cout<<"fileopt "<<fileopt<<endl;
    cout<<endl<<"----- BEGIN CODE"<<endl;
    auto *f = new TFile(file_, fileopt);
    auto *t = (TTree*)f->Get("PDsecondTree");
    cout<<"----- FILE READ"<<endl;

//----------DELETING BRANCHES----------
    cout<<"----- DELETING BRANCHES TO WRITE"<<endl;
    if(addMva_ && (t->GetBranch(method_))) {
        t->GetListOfBranches()->Remove(t->GetBranch(method_));
        t->Write();
    }
    if(mode_ == "USE"){
        if(t->GetBranch("evtMistagCat"))    t->GetListOfBranches()->Remove("evtMistagCat");
        if(t->GetBranch("evtMistagFit"))    t->GetListOfBranches()->Remove("evtMistagFit");
        if(t->GetBranch("evtMistagKde"))    t->GetListOfBranches()->Remove("evtMistagKde");
        if(t->GetBranch("evtMistagKdeExt")) t->GetListOfBranches()->Remove("evtMistagKdeExt");
        t->Write();
    }


//----------DECLARE STUFF----------
    float *catEdgeL;
    float *catEdgeR;
    float *catMistag;
    int nCat;
    TString bestFit;
    TF1 *perEvtW;
    float avgW;
    float totPCut = 0.;
    float totPCat = 0.;
    float totPFunc = 0.;
    TBranch *bMva;
    TGraph *g_pdfW;
    TGraph *g_pdfW_extended;
    TBranch *bCat;
    TBranch *bFit;
    TBranch *bKde;
    TBranch *bKde_ext;
    float evtWcat = -1.;
    float evtWfit = -1.;
    float evtWkde = -1.;
    float evtWkde_ext = -1.;

//----------READ INPUT FILE----------
    if(mode_=="ADD")
    {
        cout<<"----- ADD MODE_"<<endl;

        std::ifstream ifs("OSMuonTaggerCategories.txt",std::ifstream::in);
        if(!ifs.is_open()) return 0;

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
        cout<<"Input categories "<<"L R W"<<endl;
        for(int i=0; i<nCat; ++i)
            cout<<catEdgeL[i]<<" "<<catEdgeR[i]<<" "<<catMistag[i]<<endl;

        //FUNC
        ifs >> bestFit;
        int nPar;
        ifs >> avgW;
        ifs >> nPar;
        cout<<endl<<"Input function"<<endl;
        cout<<bestFit<<endl;
        cout<<"avgW = "<<avgW<<endl;
        perEvtW = new TF1("perEvtW", bestFit, 0., 1.);
        for(int i=0;i<nPar;++i){
            float par;
            ifs >> par;
            perEvtW->SetParameter(i, par);
            cout<<"p"<<i<<" = "<<par<<endl;
        }

        ifs.close();

        auto *f2 = new TFile("OSMuonTaggerKDE.root");
        f2->cd();
        g_pdfW = (TGraph*)f2 ->Get("pdfW");
        g_pdfW_extended = (TGraph*)f2 ->Get("pdfW_extended");
        f2->Close();
        delete f2;
        f->cd();

        if(!(t->GetBranch("evtMistagCat"))) bCat = t->Branch("evtMistagCat",&evtWcat,"evtMistagCat/F");
        if(!(t->GetBranch("evtMistagFit"))) bFit = t->Branch("evtMistagFit",&evtWfit,"evtMistagFit/F");
        if(!(t->GetBranch("evtMistagKde"))) bKde = t->Branch("evtMistagKde",&evtWkde,"evtMistagKde/F");
        if(!(t->GetBranch("evtMistagKdeExt"))) bKde_ext = t->Branch("evtMistagKdeExt",&evtWkde_ext,"evtMistagKdeExt/F");
    }

//----------DECLARE VARIABLES----------
    cout<<"----- VARIABLE DECLARATION"<<endl;

    //MVA VARS
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
    //LOW LEVEL VARS
    float muoEta;
    float muoDz;
    float muoJetPt;
    float muoJetPtRel;
    float muoJetDr;
    float muoJetEnergyRatio;
    int muoJetSize;
    float muoJetQ;
    float muoConePt;
    float muoConePtRel;
    float muoConeDr;
    float muoConeEnergyRatio;
    int muoConeSize;
    float muoConeQ;
    //TAG VARS
    int osMuon, osMuonTag, osMuonCharge, ssbLund;
    float evtWeight, ssbMass;

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

//----------COMPUTE MVA----------
    cout<<"----- BEGIN MVA SETUP"<<endl;

    TMVA::Reader reader("!Color:Silent");
    float mvaValue = -1.;
    if(addMva_ && !(t->GetBranch(method_)) bMva = t->Branch(method_, &mvaValue, method_ + "/F");

    if(readMva_){
        t->SetBranchAddress(method_, &mvaValue);
    }else{
        TMVA::PyMethodBase::PyInitialize();

        reader.AddVariable("muoPt", &muoPt);
        reader.AddVariable("abs_muoEta := abs(muoEta)", &absmuoEta);
        reader.AddVariable("muoDxy", &muoDxy);
        reader.AddVariable("abs_muoDz := abs(muoDz)", &absmuoDz);
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
        reader.BookMVA( method_, "dataset/weights/TMVAClassification_" + method_ + ".weights.xml" );
    }

    int nBinsMva = 1000;
    auto *mva    = new TH1F( "mva",    "mva",    nBinsMva, 0.0, 1.0 );
    auto *mva_RT = new TH1F( "mva_RT", "mva_RT", nBinsMva, 0.0, 1.0 );
    auto *mva_WT = new TH1F( "mva_WT", "mva_WT", nBinsMva, 0.0, 1.0 );

    vector<double> vKDERT;
    vector<double> vKDEWT;

    //EVENT LOOP
    cout<<"----- BEGIN LOOP"<<endl;
    if(nEvents_ == -1) nEvents_ = t->GetEntries();
    for(int i=0; i<nEvents_; ++i){
        if(i%100000==0) cout<<"----- at event "<<i<<endl;
        t->GetEntry(i);
        if(!osMuon) continue;
        bool nanFlag = isnan(muoDxy) || isnan(muoJetDFprob) 
                    || isinf(muoJetEnergyRatio) || isinf(muoConeEnergyRatio);
        if(nanFlag) continue;

        if(!readMva_){
            absmuoEta = abs(muoEta);
            absmuoDz = abs(muoDz);
            muoJetConePt = muoJetPt != -1 ? muoJetPt : muoConePt;
            muoJetConePtRel = muoJetPt != -1 ? muoJetPtRel : muoConePtRel;
            muoJetConeDr = muoJetPt != -1 ? muoJetDr : muoConeDr;
            muoJetConeEnergyRatio = muoJetPt != -1 ? muoJetEnergyRatio : muoConeEnergyRatio;
            muoJetConeSize = muoJetPt != -1 ? (float)muoJetSize : (float)muoConeSize;
            muoJetConeQ = muoJetPt != -1 ? muoJetQ : muoConeQ;

            mvaValue = reader.EvaluateMVA(method_);
        }
        if(addMva_) bMva->Fill();

        int evtTagCat = -1*osMuonCharge;
        int evtTagFunc = -1*osMuonCharge;

        if(mode_=="ADD"){

            if( mvaValue < catEdgeL[0] ){
                cout<<"-----Undercat "<<mvaValue<<" at "<<i<<endl;
                evtWcat = catMistag[0];
                evtWfit = perEvtW->Eval(catEdgeL[0]);
            }else if( mvaValue >= catEdgeR[nCat-1] ){
                cout<<"-----Overcat "<<mvaValue<<" at "<<i<<endl;
                evtWcat = catMistag[nCat-1];
                evtWfit = perEvtW->Eval(catEdgeR[nCat-1]);

            }else{
                for(int j=0; j<nCat; ++j){
                    if(( mvaValue >= catEdgeL[j]) 
                        && (mvaValue < catEdgeR[j]) )
                    { 
                        evtWcat = catMistag[j]; break; 
                    }
                }

                evtWfit = perEvtW->Eval(mvaValue);
            }
            evtWkde = g_pdfW->Eval(mvaValue);
            evtWkde_ext = g_pdfW_extended->Eval(mvaValue);

            bCat->Fill();
            bFit->Fill();
            bKde->Fill();
            bKde_ext->Fill();

            totPCat += 1./(float)nEvents_*pow(1.-2.*evtWcat, 2)*evtWeight;
            totPFunc += 1./(float)nEvents_*pow(1.-2.*evtWkde ,2)*evtWeight;
        }

        mva->Fill(mvaValue, evtWeight);
        if(osMuonTag == 1){
            mva_RT->Fill(mvaValue, evtWeight);
            vKDERT.push_back(mvaValue);
            if(evtWeight==2) vKDERT.push_back(mvaValue);
        }
        if(osMuonTag == 0){
            mva_WT->Fill(mvaValue, evtWeight);
            vKDEWT.push_back(mvaValue);
            if(evtWeight==2) vKDEWT.push_back(mvaValue);
        }
    }

    int nRT = mva_RT->Integral();
    int nWT = mva_WT->Integral();

    cout<<"----- MVA HISTOGRAMS FILLED"<<endl;
    cout<<"nRT "<<nRT<<endl;
    cout<<"nWT "<<nWT<<endl;

//----------CREATE----------
    if(mode_ == "CREATE")
    {
        //----------KDE----------
        sort(vKDERT.begin(), vKDERT.end());
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

        auto *pdfRT = new TF1("pdfRT",kdeRT,0.,1.,0);
        auto *pdpdfWT = new TF1("pdpdfWT",kdeWT,0.,1.,0);

        pdfRT->SetLineColor(kBlue);
        pdfRT->SetNpx(2000);
        pdpdfWT->SetLineColor(kRed);
        pdpdfWT->SetNpx(2000);

        auto *c10 = new TCanvas();
        c10->Divide(2,2);
        c10->cd(1);
        mva_RT->Draw("hist");
        c10->cd(2);
        mva_WT->Draw("hist");
        c10->cd(3);
        pdfRT->Draw();
        c10->cd(4);
        pdpdfWT->Draw();

        auto pdfW = new TF1("pdfW",
            [&](double *x, double *p)
            { 
                double yRT = vKDERT.size()*kdeRT->GetValue(x[0]); 
                double yWT = vKDEWT.size()*kdeWT->GetValue(x[0]); 
                return yWT/(yWT+yRT);
            }, 
            0.,1.,0);

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

        pdfW_extended->SetNpx(100);
        pdfW->SetNpx(100);
        auto *out_pdfW = new TGraph(pdfW);
        auto *out_pdfW_extended = new TGraph(pdfW_extended);

        auto *c16 = new TCanvas();    
        pdfW_extended->SetMarkerStyle(20);
        pdfW_extended->SetLineColor(kBlue);
        pdfW_extended->SetMarkerSize(.5);
        pdfW_extended->SetNpx(50);
        pdfW_extended->DrawClone("PL");
        pdfW->SetMarkerStyle(20);
        pdfW->SetMarkerSize(1.);
        pdfW->SetNpx(50);
        pdfW->DrawClone("P SAME");

        //----------CATEGORIES----------
        nCat = 25;
        int nTotTagged = nRT + nWT;
        int catSize = nTotTagged / nCat;
        cout<<endl<<"nTotTagged "<<nTotTagged<<endl<<"catSize "<<catSize<<endl<<endl;

        int cTot = 0;
        float cCenter = 0;
        int lastCat = 0;
        
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
            catEff[i] = (float)(catWT[i] + catRT[i]) / (float)nEvents_;
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

        bestFit = "[0]+[1]*(1-TMath::Erf([2]+[3]*x))";
        auto *fitErf = new TF1("fitErf", bestFit, 0., 1.);
        fitErf->SetParameter(0, catW[nCat-1]);
        fitErf->SetParameter(1, (catW[0]-catW[nCat-1])/2);
        fitErf->SetParLimits(0, 0., 0.3);
        fitErf->SetParLimits(1, 0.5/2, 1./2);

        auto *fitTanH = new TF1("fitTanH", "[0]+[1]*(1-TMath::TanH([2]+[3]*x))", 0., 1.);
        auto *fitATan = new TF1("fitATan", "[0]+[1]*(1-TMath::ATan([2]+[3]*x))", 0., 1.);

        gr->Fit("fitErf","MELR");

        //----------OUTPUT STREAM----------
        ofstream ofs;
        ofs.open ("OSMuonTaggerCategories.txt");
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
        ofs<<bestFit<<endl;
        ofs<<avgW<<endl;
        ofs<<fitErf->GetNumberFreeParameters()<<endl;
        for(int j=0;j<fitErf->GetNumberFreeParameters();++j)
            ofs<<fitErf->GetParameter(j)<<endl;
        ofs.close();

        auto *fo = new TFile("OSMuonTaggerKDE.root", "RECREATE");
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
        mva_WT->Draw("hist");
        mva_RT->Draw("hist same");
        mva_WT->SetLineColor(kRed);

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

        c10->Print("kdeDistributions.pdf");
        c2->Print("mvaDistributions.pdf");
        c3->Print("perEventW.pdf");
    }

    cout<<endl<<"NoCat Eff = "<<100.*(float)(nRT+nWT)/nEvents_<<"%"<<endl;
    cout<<"NoCat W = "<<100.*(float)nWT/(nRT+nWT)<<"%"<<endl;
    cout<<"NoCat P = "<<100.*((float)(nRT+nWT)/nEvents_)*pow(1.-2.*((float)nWT/(nRT+nWT)),2)<<"%"<<endl;

    cout<<endl<<"Cat P = "<<100.*totPCat<<"%"<<endl;
    cout<<"Func P = "<<100.*totPFunc<<"%"<<endl;


    if(fileopt == "UPDATE")  t->Write();
    f->Close();
    delete f;
    return 0;

}
