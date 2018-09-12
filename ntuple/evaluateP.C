TString year;
float CountEventsWithFit(TH1 *hist, TString process);

void evaluateP(TString file = "./2016/ntuBs2016.root", TCut cut = "")
{
    TFile *f = new TFile(file);
    TTree *t = (TTree*)f ->Get("PDsecondTree");
    if(file.Contains("16")) year = "2016";
    if(file.Contains("17")) year = "2017";

    float min = 5.25;
    float max = 5.50;
    int nBins = 250;

    TH1F *ssB       = new TH1F( "ssB", "ssB", nBins, min, max );
    TH1F *ssB_RT    = new TH1F( "ssB_RT", "ssB_RT", nBins, min, max );
    TH1F *ssB_WT    = new TH1F( "ssB_WT", "ssB_WT", nBins, min, max );

    TCut base = cut + "osMuon==1";
    TCut cutRT = base + "osMuonTag==1";
    TCut cutWT = base + "osMuonTag==0";

    t->Project("ssB", "ssbMass", "");
    t->Project("ssB_RT", "ssbMass", cutRT );
    t->Project("ssB_WT", "ssbMass", cutWT );

    float nTot = CountEventsWithFit(ssB, "BsJPsiPhi");
    float nRT = CountEventsWithFit(ssB_RT, "BsJPsiPhi");
    float nWT = CountEventsWithFit(ssB_WT, "BsJPsiPhi");

    float eff = (nRT+nWT)/nTot;
    float w = nWT/(nRT+nWT);

    cout<<"Eff = "<<100*eff<<"%"<<endl;
    cout<<"Mistag = "<<100*w<<"%"<<endl;
    cout<<"Power = "<< 100*eff*pow((1-2*w), 2)<<"%"<<endl;

    return;

}

float CountEventsWithFit(TH1 *hist, TString process){

    float mean=5.3663;
    if(process=="BsJPsiPhi")   mean=5.3663;
    if(process=="BuJPsiK")     mean=5.2793;
    if(process=="BdJPsiKx")    mean=5.2796;

    float sigma = 0.015;

    TString sgnDef = "[1]*exp(-0.5*((x-[0])/[4])**2)";
    sgnDef +=       "+[2]*exp(-0.5*((x-[0])/[5])**2)";
    sgnDef +=       "+[3]*exp(-0.5*((x-[0])/[6])**2)";
    TString bkgDef = "[7]*exp([8]*x)";
    TString funcDef = sgnDef + "+" + bkgDef;

    TF1 *func = new TF1("func", funcDef, 5.25, 5.5);

    func->FixParameter(0, mean);

    func->SetParameter(1, 1);
    func->SetParameter(2, 1);
    func->SetParameter(3, 1);

    func->SetParameter(4, sigma);
    func->SetParameter(5, sigma);
    func->SetParameter(6, sigma);

    func->SetParameter(7, 1);
    func->SetParameter(8, -1);

    func->SetParLimits(0, mean-sigma, mean+sigma);

    func->SetParLimits(1, 0, 1e9);
    func->SetParLimits(2, 0, 1e9);
    func->SetParLimits(3, 0, 1e9);

    func->SetParLimits(4, sigma/3, sigma*3);
    func->SetParLimits(5, sigma/3, sigma*3);
    func->SetParLimits(6, sigma/3, sigma*3);

    func->SetParLimits(7, 0, 1e9);
    func->SetParLimits(8, -1e9, 0);

    auto c1 = new TCanvas();
    hist->Draw("");
    hist->Fit("func","MRLQ");

    TString name = hist->GetName(); 

    c1->Print(name + "_" + year + ".pdf");

    TF1 *fit = hist->GetFunction("func");

    double a[3] = {
        fit->GetParameter(1),
        fit->GetParameter(2),
        fit->GetParameter(3)
    };

    double s[3] = {
        fit->GetParameter(4),
        fit->GetParameter(5),
        fit->GetParameter(6)
    };

    float nEvt = TMath::Sqrt(2*TMath::Pi())*(a[0]*s[0]+a[1]*s[1]+a[2]*s[2]);
    nEvt/=hist->GetBinWidth(0);

    return nEvt;

}