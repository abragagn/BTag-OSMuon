TString year;
TString process;
TString dir = "./";
TString suffix_;
float min_, max_;
int nBins_;
float CountEventsWithFit(TH1 *hist, TString process);

void evaluateP(TString file = "./BsMC/ntuBsMC2016.root", TCut cut = "", TString suffix = "")
{
    TFile *f = new TFile(file);
    TTree *t = (TTree*)f ->Get("PDsecondTree");
    if(file.Contains("16")) year = "2016";
    if(file.Contains("17")) year = "2017";
    if(file.Contains("Bs")){
        process = "BsJPsiPhi";
        min_ = 5.00;
        max_ = 5.50;
        nBins_ = 250;
    }
    if(file.Contains("Bp")){
        process = "BuJPsiK";
        min_ = 5.00;
        max_ = 5.50;
        nBins_ = 250;
    }

    TString dirPath( file(0, file.Index("ntu")) ); 
    dir = dirPath;
    suffix_ = suffix;

    TH1F *ssB       = new TH1F( "ssB", "ssB", nBins_, min_, max_ );
    TH1F *ssB_RT    = new TH1F( "ssB_RT", "ssB_RT", nBins_, min_, max_ );
    TH1F *ssB_WT    = new TH1F( "ssB_WT", "ssB_WT", nBins_, min_, max_ );

    TCut base = cut + "osMuon==1";
    TCut cutRT = base + "osMuonTag==1";
    TCut cutWT = base + "osMuonTag==0";

    t->Project("ssB", "ssbMass", "");
    t->Project("ssB_RT", "ssbMass", cutRT );
    t->Project("ssB_WT", "ssbMass", cutWT );

    float nTot = CountEventsWithFit(ssB, process);
    float nRT = CountEventsWithFit(ssB_RT, process);
    float nWT = CountEventsWithFit(ssB_WT, process);

    float eff = (nRT+nWT)/nTot;
    float w = nWT/(nRT+nWT);

    cout<<"Eff = "<<100*eff<<"%"<<endl;
    cout<<"Mistag = "<<100*w<<"%"<<endl;
    cout<<"Power = "<< 100*eff*pow((1-2*w), 2)<<"%"<<endl;

    return;

}

float CountEventsWithFit(TH1 *hist, TString process){

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );

    float mean=5.3663;
    if(process=="BsJPsiPhi")   mean=5.3663;
    if(process=="BuJPsiK")     mean=5.2793;
    if(process=="BdJPsiKx")    mean=5.2796;

    float sigma = 0.015;

    TString sgnDef = "[1]*TMath::Gaus(x, [0], [4], true)";
    sgnDef +=       "+[2]*TMath::Gaus(x, [0], [5], true)";
    sgnDef +=       "+[3]*TMath::Gaus(x, [0], [6], true)";
    TString bkgDef = "[7]+[8]*TMath::Erfc([9]*(x-[10]))";
    TString funcDef = sgnDef + "+" + bkgDef;

    TF1 *func = new TF1("func", funcDef, min_, max_);

    func->SetParameter(0, mean);

    func->SetParameter(1, 1);
    func->SetParameter(2, 1);
    func->SetParameter(3, 1);

    func->SetParameter(4, sigma);
    func->SetParameter(5, sigma);
    func->SetParameter(6, sigma);

    func->SetParameter(7, hist->GetBinContent(nBins_-1));
    func->SetParameter(8, 1);
    func->SetParameter(9, 20);
    func->SetParameter(10, 5.10);

    func->SetParLimits(0, mean-sigma, mean+sigma);

    func->SetParLimits(1, 0, hist->GetEntries());
    func->SetParLimits(2, 0, hist->GetEntries());
    func->SetParLimits(3, 0, hist->GetEntries());

    func->SetParLimits(4, 0, sigma*3);
    func->SetParLimits(5, 0, sigma*3);
    func->SetParLimits(6, 0, sigma*3);

    func->SetParLimits(7, 0, hist->GetBinContent(nBins_-1)*1.5);
    func->SetParLimits(8, 0, hist->GetBinContent(nBins_-1));
    func->SetParLimits(9, 10, 1e3);
    func->SetParLimits(10, 5.0, mean);

    auto c1 = new TCanvas();
    hist->Draw("");
    hist->SetMinimum(0.1);
    hist->Fit("func","MRL");

    TF1 *fit = hist->GetFunction("func");

    TF1 *f1 = new TF1("f1","[0]*TMath::Gaus(x, [1], [2], true)", min_, max_);
    TF1 *f2 = new TF1("f2","[0]*TMath::Gaus(x, [1], [2], true)", min_, max_);
    TF1 *f3 = new TF1("f3","[0]*TMath::Gaus(x, [1], [2], true)", min_, max_);
    TF1 *f4 = new TF1("f4","[0]", min_, max_);
    TF1 *f5 = new TF1("f5","[0]*TMath::Erfc([1]*(x-[2]))", min_, max_);

    f1->SetParameters(fit->GetParameter(1),fit->GetParameter(0),fit->GetParameter(4));
    f2->SetParameters(fit->GetParameter(2),fit->GetParameter(0),fit->GetParameter(5));
    f3->SetParameters(fit->GetParameter(3),fit->GetParameter(0),fit->GetParameter(6));
    f4->SetParameter(0, fit->GetParameter(7));
    f5->SetParameters(fit->GetParameter(8),fit->GetParameter(9),fit->GetParameter(10));

    f1->SetLineColor(kBlue);
    f2->SetLineColor(kViolet);
    f3->SetLineColor(kAzure);
    f4->SetLineColor(kOrange);
    f5->SetLineColor(kGreen);
    f1->SetLineStyle(2);
    f2->SetLineStyle(2);
    f3->SetLineStyle(2);
    f4->SetLineStyle(2);
    f5->SetLineStyle(2);

    f1->Draw("same");
    f2->Draw("same");
    f3->Draw("same");
    f4->Draw("same");
    f5->Draw("same");

    TString name = hist->GetName(); 

    c1->Print(dir + name + "_" + year + suffix_ + ".pdf");

    float nEvt = fit->GetParameter(1);
    nEvt += fit->GetParameter(2);
    nEvt += fit->GetParameter(3);

    nEvt/=hist->GetBinWidth(0);

    return nEvt;

}