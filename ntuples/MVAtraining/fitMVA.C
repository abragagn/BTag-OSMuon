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

void fitMVA(TString file = "../BsMC/ntuBsMC2017.root", TString method = "")
{
    gErrorIgnoreLevel = kWarning;
    TFile *f = new TFile(file);
    TTree *t = (TTree*)f ->Get("PDsecondTree");

    setGvars(file);

    TMVA::Reader reader("!Color:Silent");
    TMVA::PyMethodBase::PyInitialize();

    float muoPt_;
    float absmuoEta_;
    float muoDxy_;
    float absmuoDz_;
    float muoSoftMvaValue_;
    float muoDrB_;
    float muoPFIso_;
    float muoJetConePt_;
    float muoJetConePtRel_;
    float muoJetConeDr_;
    float muoJetConeEnergyRatio_;
    float muoJetDFprob;
    float muoJetDFprob_;
    float muoJetConeSize_;
    float muoJetConeQ_;
    float muoCharge_;

    reader.AddVariable( "muoPt", &muoPt_);
    reader.AddVariable( "abs_muoEta", &absmuoEta_);
    reader.AddVariable( "muoDxy", &muoDxy_);
    reader.AddVariable( "abs_muoDz", &absmuoDz_);
    reader.AddVariable( "muoSoftMvaValue", &muoSoftMvaValue_);
    reader.AddVariable( "muoDrB", &muoDrB_);
    reader.AddVariable( "muoPFIso", &muoPFIso_);
    reader.AddVariable( "muoJetConePt", &muoJetConePt_);
    reader.AddVariable( "muoJetConePtRel", &muoJetConePtRel_);
    reader.AddVariable( "muoJetConeDr", &muoJetConeDr_);
    reader.AddVariable( "muoJetConeEnergyRatio", &muoJetConeEnergyRatio_);
    reader.AddVariable( "muoJetDFprob", &muoJetDFprob_);
    reader.AddVariable( "muoJetConeSize", &muoJetConeSize_);
    reader.AddVariable( "muoJetConeConeQ", &muoJetConeQ_);
    reader.BookMVA( method, "dataset/weights/" + method + ".weights.xml" );

    t->SetBranchAddress("muoPt", &muoPt_);
    t->SetBranchAddress("abs(muoEta)", &absmuoEta_);
    t->SetBranchAddress("muoDxy", &muoDxy_);
    t->SetBranchAddress("abs(muoDz)", &absmuoDz_);
    t->SetBranchAddress("muoSoftMvaValue", &muoSoftMvaValue_);
    t->SetBranchAddress("muoPFIso", &muoPFIso_);
    t->SetBranchAddress("muoJetPt != -1 ? muoJetPt : muoConePt", &muoJetConePt_);
    t->SetBranchAddress("muoJetPt != -1 ? muoJetPtRel : muoConePtRel", &muoJetConePtRel_);
    t->SetBranchAddress("muoJetPt != -1 ? muoJetDr : muoConeDr", &muoJetConeDr_);
    t->SetBranchAddress("muoJetPt != -1 ? muoJetEnergyRatio : muoConeEnergyRatio", &muoJetConeEnergyRatio_);
    t->SetBranchAddress("muoJetDFprob", &muoJetDFprob_);
    t->SetBranchAddress("muoJetPt != -1 ? muoJetSize : muoConeSize", &muoJetConeSize_);
    t->SetBranchAddress("muoJetPt != -1 ? muoJetQ : muoConeQ", &muoJetConeQ_);

    int osMuon_, osMuonTag_, evtWeight_;

    t->SetBranchAddress("osMuon", &osMuon_);
    t->SetBranchAddress("osMuonTag", &osMuonTag_);
    t->SetBranchAddress("evtWeight", &evtWeight_);

    int nBinsMva = 1000;
    TH1F *mva   = new TH1F( "mva", "mva", nBinsMva, 0.0, 1.0 );
    TH1F *mva_RT   = new TH1F( "mva_RT", "mva_RT", nBinsMva, 0.0, 1.0 );
    TH1F *mva_WT   = new TH1F( "mva_WT", "mva_WT", nBinsMva, 0.0, 1.0 );

    int nEvents = t->GetEntries();

    for(int i =0; i<nEvents; ++i){
        t->GetEntry(i);
        if(!osMuon_) continue;
        float mvaValue = reader.EvaluateMVA(method);
        mva->Fill(mvaValue, evtWeight_);
        if(osMuonTag_ == 1) mva_RT->Fill(mvaValue, evtWeight_);
        if(osMuonTag_ == 0) mva_WT->Fill(mvaValue, evtWeight_);
    }

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
