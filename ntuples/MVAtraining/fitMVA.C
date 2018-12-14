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

void fitMVA(TString file = "../BsMC/ntuBsMC2017.root", TString method = "DNNOsMuon2017test231")
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

    int nBinsMva = 1000;
    TH1F *mva   = new TH1F( "mva", "mva", nBinsMva, 0.0, 1.0 );
    TH1F *mva_RT   = new TH1F( "mva_RT", "mva_RT", nBinsMva, 0.0, 1.0 );
    TH1F *mva_WT   = new TH1F( "mva_WT", "mva_WT", nBinsMva, 0.0, 1.0 );

    int nEvents = t->GetEntries();

    for(int i=0; i<nEvents; ++i){
        if(i%100000==0) cout<<"----- at event "<<i<<endl;
        t->GetEntry(i);
        if(!osMuon) continue;

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
        if(osMuonTag == 1) mva_RT->Fill(mvaValue, evtWeight);
        if(osMuonTag == 0) mva_WT->Fill(mvaValue, evtWeight);
    }

    cout<<"----- MVA HISTOGRAMS FILLED"<<endl;

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

    int nPoints = 20;
    int catSize = mva->GetEntries() / nPoints;

    cout<<catSize<<endl<<endl;

    int cTot = 0;
    int cCenter = 0;
    
    float   *catEdge    = new float[nPoints+1];
    float   *catCenter  = new float[nPoints];
    int     *catRT      = new int[nPoints];
    int     *catWT      = new int[nPoints];
    float   *catEff     = new float[nPoints];
    float   *catW       = new float[nPoints];
    float   *catP       = new float[nPoints];

    catEdge[nPoints] = 1;
    catEdge[0] = 0;
    for(int i=0; i<nPoints; ++i){
        catRT[i] = 0;
        catWT[i] = 0;
        catCenter[i] = 0;
    }

    int cCat = 0;
    for(int i = 1; i<=nBinsMva; ++i){
        cTot += mva->GetBinContent(i);
        cCenter += mva->GetBinCenter(i)*mva->GetBinContent(i);
        catRT[cCat] += mva_RT->GetBinContent(i);
        catWT[cCat] += mva_WT->GetBinContent(i);
        if(cTot >= catSize){
            catEdge[cCat+1] = mva->GetXaxis()->GetBinLowEdge(i+1);
            catCenter[cCat] = (float)cCenter/(float)cTot;
            cTot = 0;
            cCenter = 0;
            cCat++;
        }
    }

    for(int i=0; i<nPoints; ++i)
    {
        catW[i] = (float)catWT[i] / (float)(catWT[i] + catRT[i]);
        cout<<"CAT "<<i+1<<" ["<<catEdge[i]<<" - "<<catEdge[i+1]<<" ] ";
        cout<<" - "<<catCenter[i]<<" - ";
        cout<<catRT[i]<<" "<<catWT[i]<<endl;
    }


    TGraph* gr = new TGraph(nPoints,catCenter,catW);
    TCanvas *c1 = new TCanvas();
    gr->Draw("AP*");

}
