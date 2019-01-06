#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <sstream>

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
#include "TTree.h"
#include "TFile.h"
#include "TKDE.h"
#include "TROOT.h"

using namespace std;

TString year_;
TString process_;
TString dir_ = "./";
TString suffix_;
float min_ = 5.1;
float max_ = 5.5;
int nBins_ = 50;

void setGvars(TString filename)
{
    if(filename.Contains("16")) year_ = "2016";
    if(filename.Contains("17")) year_ = "2017";
    if(filename.Contains("Bs")) process_ = "BsJPsiPhi";
    TString dirPath( filename(0, filename.Index("ntu")) ); 
    dir_ = dirPath;
}

int fitMVA(TString file = "../BsMC/ntuBsMC2017.root"
            , TString method = "DNNOsMuon2017test231"
            , TString mode = "CREATE")
{
    cout<<"----- BEGIN CODE"<<endl;

    gErrorIgnoreLevel = kWarning;
    if(mode != "CREATE" && mode != "USE" && mode != "CHECK"){
        cout<<"WRONG MODE"<<endl;
        return 0;
    }

    TFile *f = new TFile(file);
    TTree *t = (TTree*)f ->Get("PDsecondTree");
    cout<<"----- FILE OPEN"<<endl;
    setGvars(file);

//----------DEFINE STUFF----------
    float *catEdgeL;
    float *catEdgeR;
    float *catMistag;
    int nCat;
    TString bestFit;
    TF1 *perEvtW;
    float totPCut = 0.;
    float totPCat = 0.;
    float totPFunc = 0.;
    float avgW;
    float *wCalc;
    float *wRT;
    float *wWT;
    int nBinCheck = 20;
    float step;
    TH1 *hCheck;

//----------READ INPUT FILE----------
    if(mode=="USE" || mode=="CHECK")
    {
        std::ifstream ifs("OSMuonTaggerDefinition.txt",std::ifstream::in);
        if(!ifs.is_open()) return 0;

        ifs >> method;

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

        if(mode == "CHECK"){

            wCalc = new float[nBinCheck];
            wRT   = new float[nBinCheck];
            wWT   = new float[nBinCheck];

            step = 0.5/nBinCheck;

            for(int j=0;j<nBinCheck;++j){
                wCalc[j] = (float)j*step + step/2;
                wRT[j] = 0;
                wWT[j] = 0;
            }
        }
    }


    cout<<"----- BEGIN MVA SETUP"<<endl;
    
//----------COMPUTE MVA----------
    TMVA::Reader reader("!Color:Silent");
    TMVA::PyMethodBase::PyInitialize();

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
    reader.BookMVA( method, "dataset/weights/TMVAClassification_" + method + ".weights.xml" );

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

    //OTHER VARS
    int osMuon, osMuonTag, osMuonCharge, ssbLund;
    float evtWeight;

    t->SetBranchAddress("osMuon", &osMuon);
    t->SetBranchAddress("osMuonTag", &osMuonTag);
    t->SetBranchAddress("evtWeight", &evtWeight);
    t->SetBranchAddress("muoCharge", &osMuonCharge);
    t->SetBranchAddress("ssbLund", &ssbLund);

    int nBinsMva = 1000;
    TH1F *mva    = new TH1F( "mva",    "mva",    nBinsMva, 0.0, 1.0 );
    TH1F *mva_RT = new TH1F( "mva_RT", "mva_RT", nBinsMva, 0.0, 1.0 );
    TH1F *mva_WT = new TH1F( "mva_WT", "mva_WT", nBinsMva, 0.0, 1.0 );

    vector<double> vKDERT;
    vector<double> vKDEWT;

    //EVENT LOOP
    cout<<"----- BEGIN LOOP"<<endl;
    int nEvents = 1000000;//t->GetEntries();
    for(int i=0; i<nEvents; ++i){
        //if(i!=15285) continue;
        if(i%100000==0) cout<<"----- at event "<<i<<endl;
        t->GetEntry(i);
        if(!osMuon) continue;
        bool nanFlag = isnan(muoDxy) || isnan(muoJetDFprob) 
                    || isinf(muoJetEnergyRatio) || isinf(muoConeEnergyRatio);
        if(nanFlag) continue;

        absmuoEta = abs(muoEta);
        absmuoDz = abs(muoDz);
        muoJetConePt = muoJetPt != -1 ? muoJetPt : muoConePt;
        muoJetConePtRel = muoJetPt != -1 ? muoJetPtRel : muoConePtRel;
        muoJetConeDr = muoJetPt != -1 ? muoJetDr : muoConeDr;
        muoJetConeEnergyRatio = muoJetPt != -1 ? muoJetEnergyRatio : muoConeEnergyRatio;
        muoJetConeSize = muoJetPt != -1 ? (float)muoJetSize : (float)muoConeSize;
        muoJetConeQ = muoJetPt != -1 ? muoJetQ : muoConeQ;

        float mvaValue = reader.EvaluateMVA(method);

        float evtWCat = -1.;
        float evtWFunc = -1.;
        int evtTagCat = -1*osMuonCharge;
        int evtTagFunc = -1*osMuonCharge;

        if(mode=="USE" || mode=="CHECK"){

            if( mvaValue < catEdgeL[0] ){
                cout<<"Undercat "<<mvaValue<<" at "<<i<<endl;
                evtWCat = catMistag[0];
                evtWFunc = perEvtW->Eval(catEdgeL[0]);
            }else if( mvaValue >= catEdgeR[nCat-1] ){
                cout<<"-----Overcat "<<mvaValue<<" at "<<i<<endl;
                evtWCat = catMistag[nCat-1];
                evtWFunc = perEvtW->Eval(catEdgeR[nCat-1]);

            }else{
                for(int j=0; j<nCat; ++j){
                    if(( mvaValue >= catEdgeL[j]) 
                        && (mvaValue < catEdgeR[j]) )
                    { 
                        evtWCat = catMistag[j]; break; 
                    }
                }

                evtWFunc = perEvtW->Eval(mvaValue);
            }
            if(evtWCat>0.5){ evtWCat = 1 - evtWCat; evtTagCat *= -1;}
            totPCat += 1./(float)nEvents*pow(1.-2.*evtWCat, 2)*evtWeight;

            if(evtWFunc>0.5){ evtTagFunc *= -1; evtWFunc = 1 - evtWFunc; }
            totPFunc += 1./(float)nEvents*pow(1.-2.*evtWFunc ,2)*evtWeight;

            if(mode == "CHECK"){
                for(int j=0;j<nBinCheck;++j){
                    if( (evtWFunc>=(float)j*step) && (evtWFunc<((float)j*step+step)) ){
                        if(TMath::Sign(1, ssbLund) == evtTagFunc) wRT[j] += evtWeight;
                        if(TMath::Sign(1, ssbLund) != evtTagFunc) wWT[j] += evtWeight;
                        break;
                    }
                }
            }

        }

        mva->Fill(mvaValue, evtWeight);
        if(osMuonTag == 1){
            mva_RT->Fill(mvaValue, evtWeight);
            vKDERT.push_back(mvaValue);
            if(evtWeight==2.) vKDERT.push_back(mvaValue);
        }
        if(osMuonTag == 0){
            mva_WT->Fill(mvaValue, evtWeight);
            vKDEWT.push_back(mvaValue);
            if(evtWeight==2.) vKDEWT.push_back(mvaValue);
        }
    }

    if(mode == "CHECK"){
        vector<float> wX;
        vector<float> wY;
        vector<float> weY;

        for(int j=0;j<nBinCheck;++j){
            if(wRT[j] == 0 && wWT[j] == 0 ) continue;

            wX.push_back(wCalc[j]);
            wY.push_back(wWT[j]/(wWT[j]+wRT[j]));
            if(wWT[j] == 0) weY.push_back(1/wRT[j]);
            else if(wRT[j] == 0) weY.push_back(1/wWT[j]);
            else weY.push_back(sqrt((wWT[j] * wRT[j])/pow((wWT[j] + wRT[j]),3) ));

            cout<<"BIN "<<j<<", wCalc "<<wCalc[j]<<", wRT "<<wRT[j]<<", wWT "<<wWT[j]<<", wMeas "<<wWT[j]/(wWT[j]+wRT[j]);
            cout<<" +- "<<sqrt((wWT[j] * wRT[j])/pow((wWT[j] + wRT[j]),3) )<<endl;
        }

        TGraph *grW = new TGraphErrors(wX.size(),&wX[0],&wY[0],0,&weY[0]);
        TF1 *fCheck = new TF1("fCheck","[0]+[1]*(x-[2])",0.,.5);
        fCheck->FixParameter(2, avgW);
        fCheck->SetParameter(0, avgW);
        fCheck->SetParameter(1, 1);

        grW->Fit("fCheck");
        TF1 *myfunc = grW->GetFunction("fCheck");

        vector<float> wResY;
        vector<float> wResEY;

        for (unsigned int j=0;j<wX.size();++j) { 
            wResY.push_back((wY[j] - myfunc->Eval(wX[j]))/weY[j]);
            wResEY.push_back(1.);
        } 

        TGraph *grWres = new TGraphErrors(wX.size(),&wX[0],&wResY[0],0,&wResEY[0]);

        TCanvas *c30 = new TCanvas();
        c30->Divide(1,2);
        c30->cd(1);
        grW->SetMarkerStyle(20);
        grW->SetMarkerSize(.5);
        grW->Draw("APE");
        c30->cd(2);
        grWres->SetMarkerStyle(20);
        grWres->SetMarkerSize(.5);
        grWres->Draw("APE");

        //c30->Print("check.jpg");

        float p0 = myfunc->GetParameter(0);
        float p1 = myfunc->GetParameter(0);

        cout<<myfunc->GetParameter(0)-myfunc->GetParameter(1)<<" +- "<<sqrt(pow(myfunc->GetParError(0),2) + pow(avgW,2)*pow(myfunc->GetParError(1),2) )<<endl;

    }


    int nRT = mva_RT->Integral();
    int nWT = mva_WT->Integral();
/*
    float tailFraction = 0.001;
    sort(vKDERT.begin(), vKDERT.end());
    sort(vKDEWT.begin(), vKDEWT.end());
    for(int i=0;i<int(nRT*tailFraction/2);++i) vKDERT.pop_back(); 
    for(int i=0;i<int(nWT*tailFraction/2);++i) vKDEWT.pop_back();
    vKDERT.erase(vKDERT.begin(), vKDERT.begin()+(int)(nRT*tailFraction/2));
    vKDEWT.erase(vKDEWT.begin(), vKDEWT.begin()+(int)(nWT*tailFraction/2));
*/
    sort(vKDERT.begin(), vKDERT.end());
    sort(vKDEWT.begin(), vKDEWT.end());

    float xMin = min(vKDERT[0], vKDEWT[0]);
    float xMax = max(vKDERT[vKDERT.size()-1], vKDEWT[vKDEWT.size()-1]);

    cout<<"----- MVA HISTOGRAMS FILLED"<<endl;

//----------CREATE----------
    if(mode == "CREATE")
    {
        //----------KDE----------
        float rho = .5;

        //mva_RT->GetXaxis()->SetRangeUser(xMin, xMax);
        //mva_WT->GetXaxis()->SetRangeUser(xMin, xMax);

        cout<<"nRT "<<nRT<<endl;
        cout<<"nWT "<<nWT<<endl;
        cout<<"xMin "<<xMin<<endl;
        cout<<"xMax "<<xMax<<endl;
        cout<<"rho "<<rho<<endl;

        TKDE *kdeRT = new TKDE(vKDERT.size(), &vKDERT[0], 0., 1., "", rho);
        TKDE *kdeWT = new TKDE(vKDEWT.size(), &vKDEWT[0], 0., 1., "", rho);

        TF1 *pdfRT = new TF1("pdfRT",kdeRT,0.,1.,0);
        TF1 *pdfWT = new TF1("pdfWT",kdeWT,0.,1.,0);

        pdfRT->SetLineColor(kBlue);
        pdfRT->SetNpx(2000);
        pdfWT->SetLineColor(kRed);
        pdfWT->SetNpx(2000);

        TCanvas *c10 = new TCanvas();
        c10->Divide(2,2);
        c10->cd(1);
        mva_RT->Draw("hist");
        c10->cd(2);
        mva_WT->Draw("hist");
        c10->cd(3);
        pdfRT->Draw();
        c10->cd(4);
        pdfWT->Draw();

        auto fW = new TF1("fW",
            [&](double *x, double *p)
            { 
                double yRT = vKDERT.size()*kdeRT->GetValue(x[0]); 
                double yWT = vKDEWT.size()*kdeWT->GetValue(x[0]); 
                return yWT/(yWT+yRT);
            }, 
            0.,1.,0);

        auto fW_extended = new TF1("fW_extended",
            [&](double *x, double *p)
            {
                if(x[0]<fW->GetMaximumX())
                    return fW->GetMaximum();
                if(x[0]>=fW->GetMinimumX())
                    return fW->GetMinimum();
                return fW->Eval(x[0]);
            }, 
            0.,1.,0);

        TCanvas *c16 = new TCanvas();    
        fW_extended->SetMarkerStyle(20);
        fW_extended->SetLineColor(kBlue);
        fW_extended->SetMarkerSize(.5);
        fW_extended->SetNpx(40);
        fW_extended->DrawClone("PL");
        fW->SetMarkerStyle(20);
        fW->SetMarkerSize(1.);
        fW->SetNpx(50);
        fW->DrawClone("P SAME");

        //----------CATEGORIES----------
    /*
        TString base =  "evtWeight*((";
        TString cutRT = base + ")&&osMuon&&osMuonTag==1)";
        TString cutWT = base + ")&&osMuon&&osMuonTag==0)";

        t->Project("mva", "osMuonTagMvaValue", cutRT + "||" + cutWT);
        t->Project("mva_RT", "osMuonTagMvaValue", cutRT);
        t->Project("mva_WT", "osMuonTagMvaValue", cutWT);
    */

    /*
        TH1F *mB       = new TH1F( "mB", "mB", nBins_, min_, max_ );
        TH1F *mB_RT    = new TH1F( "mB_RT", "mB_RT", nBins_, min_, max_ );
        TH1F *mB_WT    = new TH1F( "mB_WT", "mB_WT", nBins_, min_, max_ );
        t->Project("mB", "mBMass", base + "))");
        t->Project("mB_RT", "mBMass", cutRT );
        t->Project("mB_WT", "mBMass", cutWT );
    */

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
            catEff[i] = (float)(catWT[i] + catRT[i]) / (float)nEvents;
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
        TGraph *gr = new TGraph(nCat,catCenter,catW);

        TF1 *fitErf = new TF1("fitErf","[0]+[1]*(1-TMath::Erf([2]+[3]*x))",0.,1.);
        fitErf->SetParameter(0, catW[nCat-1]);
        fitErf->SetParameter(1, (catW[0]-catW[nCat-1])/2);
        fitErf->SetParLimits(0, 0., 0.3);
        fitErf->SetParLimits(1, 0.5/2, 1./2);

        TF1 *fitErfLim = new TF1("fitErfLim","[0]+[1]*(1-TMath::Erf([2]+[3]*x))",0.,1.);
        fitErfLim->SetParameter(0, catW[nCat-1]);
        fitErfLim->SetParameter(1, (catW[0]-catW[nCat-1])/2);
        fitErfLim->SetParLimits(0, 0., catW[nCat-1]);
        fitErfLim->SetParLimits(1, 0., (catW[0]-catW[nCat-1])/2);

        TF1 *fitTanH = new TF1("fitTanH","[0]+[1]*(1-TMath::TanH([2]+[3]*x))",0.,1.);
        TF1 *fitATan = new TF1("fitATan","[0]+[1]*(1-TMath::ATan([2]+[3]*x))",0.,1.);

        gr->Fit("fitErf","MELR");
        //gr->Fit("fitErfLim","MELR+");

        bestFit = "[0]+[1]*(1-TMath::Erf([2]+[3]*x))";

//----------OUTPUT STREAM----------
        ofstream ofs;
        ofs.open ("OSMuonTaggerDefinition.txt");
        //CAT
        ofs<<method<<endl;
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

/*        TFile *fo = new TFile("OSMuonTaggerDefinition.root");
        fo-cd();
        kdeRT->Write();
        fo->Close();
        f->cd()*/
//----------PRINT----------


        cout<<endl<<"Cat Eff = "<<100*totEff<<"%"<<endl;
        cout<<"Cat W = "<<100*avgW<<"%"<<endl;
        cout<<"Cat P = "<<100*totP<<"%"<<endl;

        TCanvas *c2 = new TCanvas();
        mva_WT->Draw("hist");
        mva_RT->Draw("hist same");
        mva_WT->SetLineColor(kRed);

        TCanvas *c3 = new TCanvas();
        fW_extended->DrawClone("");
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(.5);
        gr->SetTitle("per-event mistag");
        gr->GetXaxis()->SetTitle("MVA output");
        gr->SetMinimum(0.);
        gr->Draw("PE SAME");
        grAsy->Draw("E SAME");
        fW->SetLineColor(kGreen);
        fW->DrawClone("P SAME");

        c3->Update();

        //c10->Print("kde.jpg");
        //c2->Print("mva.jpg");
        //c3->Print("perEventW.jpg");
    }

    cout<<endl<<"NoCat Eff = "<<100.*(float)(nRT+nWT)/nEvents<<"%"<<endl;
    cout<<"NoCat W = "<<100.*(float)nWT/(nRT+nWT)<<"%"<<endl;
    cout<<"NoCat P = "<<100.*((float)(nRT+nWT)/nEvents)*pow(1.-2.*((float)nWT/(nRT+nWT)),2)<<"%"<<endl;

    cout<<endl<<"Cat P = "<<100.*totPCat<<"%"<<endl;
    cout<<"Func P = "<<100.*totPFunc<<"%"<<endl;

    return 0;

}

int main( int argc, char** argv )
{
    cout<<"----- BEGIN MAIN"<<endl;
    return fitMVA(argv[0],argv[1],argv[2]); 
}