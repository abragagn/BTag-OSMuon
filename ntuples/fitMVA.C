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

void fitMVA(TString file = "./BsMC/ntuBsMC2016.root", TString cutEvt = "")
{
    gErrorIgnoreLevel = kWarning;
    TFile *f = new TFile(file);
    TTree *t = (TTree*)f ->Get("PDsecondTree");

    setGvars(file);

    TH1F *mB       = new TH1F( "mB", "mB", nBins_, min_, max_ );
    TH1F *mB_RT    = new TH1F( "mB_RT", "mB_RT", nBins_, min_, max_ );
    TH1F *mB_WT    = new TH1F( "mB_WT", "mB_WT", nBins_, min_, max_ );

    int nBinsMva = 1000;

    TH1F *mva   = new TH1F( "mva", "mva", nBinsMva, 0.0, 1.0 );
    TH1F *mva_RT   = new TH1F( "mva_RT", "mva_RT", nBinsMva, 0.0, 1.0 );
    TH1F *mva_WT   = new TH1F( "mva_WT", "mva_WT", nBinsMva, 0.0, 1.0 );

    if(cutEvt=="")  cutEvt = "1";
    TString base =  "evtWeight*((" + cutEvt;
    TString cutRT = base + ")&&osMuon&&osMuonTag==1)";
    TString cutWT = base + ")&&osMuon&&osMuonTag==0)";

    //t->Project("mB", "mBMass", base + "))");
    //t->Project("mB_RT", "mBMass", cutRT );
    //t->Project("mB_WT", "mBMass", cutWT );

    t->Project("mva", "osMuonTagMvaValue", cutRT + "||" + cutWT);
    t->Project("mva_RT", "osMuonTagMvaValue", cutRT);
    t->Project("mva_WT", "osMuonTagMvaValue", cutWT);

    int nPoints = 20;
    int catSize = mva->GetEntries() / nPoints;

    cout<<catSize<<endl<<endl;

    int cTot = 0;
    
    float   *catEdges   = new float[nPoints+1];
    int     *catRT      = new int[nPoints];
    int     *catWT      = new int[nPoints];
    float   *catEff     = new float[nPoints];
    float   *catW       = new float[nPoints];
    float   *catP       = new float[nPoints];

    catEdges[nPoints] = 1;
    catEdges[0] = 0;
    for(int i=0; i<nPoints; ++i){
        catRT[i] = 0;
        catWT[i] = 0;
    }

    int cCat = 0;
    for(int i = 1; i<=nBinsMva; ++i){
        cTot += mva->GetBinContent(i);
        catRT[cCat] += mva_RT->GetBinContent(i);
        catWT[cCat] += mva_WT->GetBinContent(i);
        if(cTot >= catSize){
            catEdges[cCat+1] = mva->GetXaxis()->GetBinLowEdge(i+1);
            cTot = 0;
            cCat++;
        }
    }

    for(int i=0; i<nPoints; ++i)
    {
        catW[i] = (float)catWT[i] / (float)(catWT[i] + catRT[i]);
        cout<<"CAT "<<i+1<<" ["<<catEdges[i]<<" - "<<catEdges[i+1]<<" ] ";
        cout<<catRT[i]<<" "<<catWT[i]<<endl;
    }



    TGraph* gr = new TGraph(nPoints,catEdgesW0,catW);
    TCanvas *c1 = new TCanvas();
    gr->Draw("AP*");

}
