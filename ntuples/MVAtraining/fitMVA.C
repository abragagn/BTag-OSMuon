TString year_;
TString process_;
TString dir_ = "./";
TString suffix_;
float min_ = 5.1;
float max_ = 5.5;
int nBins_=50;

void setGvars(TString filename)
{
    if(filename.Contains("16")) year_ = "2016";
    if(filename.Contains("17")) year_ = "2017";
    if(filename.Contains("Bs")) process_ = "BsJPsiPhi";
    TString dirPath( filename(0, filename.Index("ntu")) ); 
    dir_ = dirPath;
}

void fitMVA(TString file = "../BsMC/ntuBsMC2017.root", TString method = "BDTOsMuon2017test999")
{
    gErrorIgnoreLevel = kWarning;
    TFile *f = new TFile(file);
    TTree *t = (TTree*)f ->Get("PDsecondTree");

    setGvars(file);

    TMVA::Reader reader("!Color:Silent");
    TMVA::PyMethodBase::PyInitialize();

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

    int osMuon, osMuonTag;
    float evtWeight;

    t->SetBranchAddress("osMuon", &osMuon);
    t->SetBranchAddress("osMuonTag", &osMuonTag);
    t->SetBranchAddress("evtWeight", &evtWeight);

    int nBinsMva = 500;
    TH1F *mva    = new TH1F( "mva",    "mva",    nBinsMva, 0.0, 1.0 );
    TH1F *mva_RT = new TH1F( "mva_RT", "mva_RT", nBinsMva, 0.0, 1.0 );
    TH1F *mva_WT = new TH1F( "mva_WT", "mva_WT", nBinsMva, 0.0, 1.0 );

    int nEvents = t->GetEntries();

    std::vector<double> vKDERT;
    std::vector<double> vKDERT_w;
    std::vector<double> vKDEWT;
    std::vector<double> vKDEWT_w;

    for(int i=0; i<nEvents; ++i){
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

        mva->Fill(mvaValue, evtWeight);
        if(osMuonTag == 1){
            mva_RT->Fill(mvaValue, evtWeight);
            vKDERT.push_back(mvaValue);
            vKDERT_w.push_back(evtWeight);
        }
        if(osMuonTag == 0){
            mva_WT->Fill(mvaValue, evtWeight);
            vKDEWT.push_back(mvaValue);
            vKDEWT_w.push_back(evtWeight);
        }
    }

    cout<<"----- MVA HISTOGRAMS FILLED"<<endl;

    //----------KDE----------
    int nRT = vKDERT.size();
    int nWT = vKDEWT.size();

    auto itMaxRT = max_element(std::begin(vKDERT), std::end(vKDERT));
    auto itMinRT = min_element(std::begin(vKDERT), std::end(vKDERT));
    auto itMaxWT = max_element(std::begin(vKDEWT), std::end(vKDEWT));
    auto itMinWT = min_element(std::begin(vKDEWT), std::end(vKDEWT));

    float xMin = *itMinWT < *itMinRT ? *itMinWT : *itMinRT;
    float xMax = *itMaxRT > *itMaxWT ? *itMaxRT : *itMaxWT;

    float rhoRT = 0.1; // = pow(nRT,-1./5);
    float rhoWT = 0.1; // = pow(nWT,-1./5);

    mva_RT->GetXaxis()->SetRangeUser(xMin, xMax);
    mva_WT->GetXaxis()->SetRangeUser(xMin, xMax);

    cout<<"nRT "<<nRT<<endl;
    cout<<"nWT "<<nWT<<endl;
    cout<<"xMin "<<xMin<<endl;
    cout<<"xMax "<<xMax<<endl;    
    cout<<"rhoRT "<<rhoRT<<endl;
    cout<<"rhoWT "<<rhoWT<<endl;

    TKDE::EMirror mirror = TKDE::kMirrorBoth;

    TKDE *kdeRT = new TKDE(nRT, &vKDERT[0],&vKDERT_w[0],xMin,xMax,
        "KernelType:Epanechnikov;Iteration:Adaptive;Mirror:noMirror;Binning:RelaxedBinning", rhoRT);
    TKDE *kdeWT = new TKDE(nWT, &vKDEWT[0],&vKDEWT_w[0],xMin,xMax,
        "KernelType:Epanechnikov;Iteration:Adaptive;Mirror:noMirror;Binning:RelaxedBinning", rhoWT);

    TF1 *pdfRT = new TF1("pdfRT",kdeRT,xMin,xMax,0);
    TF1 *pdfWT = new TF1("pdfWT",kdeWT,xMin,xMax,0);

    pdfRT->SetLineColor(kBlue);
    pdfRT->SetNpx(2000);
    pdfWT->SetLineColor(kRed);
    pdfWT->SetNpx(2000);

/*    auto fRT = new TF1("fRT",
        [&](double *x, double *p) { 
            return mva_RT->Integral()*(*kdeRT)(x[0]); }, 
        xMin,xMax,0);

    auto fWT = new TF1("fRT",
        [&](double *x, double *p) { 
            return mva_WT->Integral()*(*kdeWT)(x[0]); }, 
        xMin,xMax,0);
*/

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
        [&](double *x, double *p) { 
            double yRT = nRT*(*kdeRT)(x[0]); 
            double yWT = nWT*(*kdeWT)(x[0]); return yWT/(yWT+yRT); }, 
        xMin,xMax,0);

    TCanvas *c12 = new TCanvas();
    fW->SetMarkerStyle(20);
    fW->SetMarkerSize(1.);
    fW->SetNpx(25);
    fW->DrawClone("P");

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


    int nCat = 25;
    int nTot = nRT + nWT;
    int catSize = nTot / nCat;
    cout<<endl<<nTot<<endl<<catSize<<endl<<endl;

    int cTot = 0;
    int hTot = 0;
    int cCenter = 0;
    
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
        cTot = cTot + mva_RT->GetBinContent(i) + mva_WT->GetBinContent(i);
        hTot = hTot  + mva->GetBinContent(i);
        cout<<"bin "<<i<<", "<<mva_RT->GetBinContent(i) + mva_WT->GetBinContent(i);
        cout<<" ("<<cTot<<")"<<" ["<<hTot<<"]"<<endl;
        if(cTot >= catSize){
            catEdge[cCat] = mva->GetXaxis()->GetBinLowEdge(i+1);
            catCenter[cCat] = (float)cCenter/(float)cTot;
            cout<<" -----> cat "<<cCat<<", "<<cTot<<" ["<<catEdge[cCat]<<"]"<<endl;
            cTot = 0;
            cCenter = 0;
            cCat++;
        }
        if(i==nBinsMva){
            catCenter[cCat] = (float)cCenter/(float)cTot;
            catEdge[cCat] = xMax;
        }
    }

    cout<<endl;

    for(int i=0; i<nCat; ++i)
    {
        catW[i] = (float)catWT[i] / (float)(catWT[i] + catRT[i]);
        vexh[i] = catEdge[i] - catCenter[i];
        if(i!=0) vexl[i] = catCenter[i] - catEdge[i-1];
        else     vexl[i] = catCenter[i] - xMin;

        vey[i] = sqrt(((float)catWT[i] * (float)catRT[i])/pow(((float)catWT[i] + (float)catRT[i]),3) );


        cout<<"CAT "<<i+1;
        if(i!=0) cout<<" ["<<catEdge[i-1]<<" - "<<catEdge[i]<<" ] ";
        else     cout<<" ["<<xMin<<" - "<<catEdge[i]<<" ] ";
        cout<<" - "<<catCenter[i]<<" - ";
        cout<<catRT[i] + catWT[i]<<endl;
    }


    TGraph* gr = new TGraphAsymmErrors(nCat,catCenter,catW,vexl,vexh,vey,vey);
    TCanvas *c1 = new TCanvas();
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(.5);
    gr->SetTitle("per-event mistag");
    gr->GetXaxis()->SetTitle("MVA output");
    gr->Draw("APE");

    TCanvas *c2 = new TCanvas();
    mva_WT->Draw("hist");
    mva_RT->Draw("hist same");
    mva_WT->SetLineColor(kRed);

    TCanvas *c3 = new TCanvas();
    gr->Draw("APE");
    fW->SetLineColor(kGreen);
    fW->DrawClone("P SAME");
}
