void customFOM( TString var = "muoPt", TString cutDir = "down", TString addCut = "", bool verbose = false )
{

    if((cutDir != "down") && (cutDir != "up")) return;

    gStyle->SetOptTitle(1); gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);   gStyle->SetStatBorderSize(0);
    gStyle->SetStatX(.49);  gStyle->SetStatY(.89);

    TFile *f = new TFile("ntuBs2016.root");
    TTree *t = (TTree*)f ->Get("PDsecondTree");

    int nBins=100;
    float min=t->GetMinimum(var);
    float max=t->GetMaximum(var);

    TH1F *Tot = new TH1F( "Tot", "Tot", 100, 5.0, 6.0 );
    TH1F *Sgn = new TH1F( "Sgn", "Sgn", nBins, min, max );
    TH1F *Bkg = new TH1F( "Bkg", "Bkg", nBins, min, max );

    TCut generalCuts = "osMuon==1&&((fabs(muoEta)<1.2&&muoSoftMvaValue>0.26)||((fabs(muoEta)>=1.2&&muoSoftMvaValue>0.50)))";
    generalCuts += addCut;
    TCut sgnDef = generalCuts + "osMuonTag==1";
    TCut bkgDef = generalCuts + "osMuonTag==0";


    t->Project("Tot", "ssbMass", "");
    t->Project("Sgn", var, sgnDef);
    t->Project("Bkg", var, bkgDef);

    int nTot = Tot->GetEntries();
    float nSgn=0, nBkg=0;

    float *x = new float[nBins];
    float *yE = new float[nBins];
    float *yW = new float[nBins];
    float *yP = new float[nBins];

    for(int i=0; i<nBins; ++i){
        x[i]=0;
        yE[i]=0;
        yW[i]=0;
        yP[i]=0;
    }

    float bestCut = 0;
    float bestP = 0, bestE = 0, bestW = 0.5;

    if(verbose) cout<<endl<<"nTot = "<<nTot<<endl;
    if(verbose) cout<<"x[bin]  nSgn  nBkg  yE[bin]  yW[bin]  yP[bin]"<<endl<<endl;

    if(cutDir == "down"){
        for(int bin=nBins; bin>0; --bin){

            float varCut = Sgn->GetBinLowEdge(bin);

            nSgn += Sgn->GetBinContent(bin);
            nBkg += Bkg->GetBinContent(bin);
            if(nSgn + nBkg == 0) continue;
            float eff = (nSgn + nBkg)/ nTot;
            float w = nBkg / (nSgn + nBkg);
            float p = eff*pow(1-2*w,2);
            x[bin-1] = varCut;
            yE[bin-1] = 100*eff ;
            yW[bin-1] = 100*w;
            yP[bin-1] = 100*p;

            if(verbose) cout<<x[bin-1]<<"    "<<nSgn<<"    "<<nBkg<<"    ";
            if(verbose) cout<<yE[bin-1]<<"    "<<yW[bin-1]<<"    "<<yP[bin-1]<<endl;

            if(p > bestP){
                bestCut = varCut;
                bestP = p;
                bestE = eff;
                bestW = w;
            }

        }       
    }else{
        for(int bin=0; bin<=nBins; ++bin){

            float varCut = Sgn->GetBinLowEdge(bin+1);

            nSgn += Sgn->GetBinContent(bin);
            nBkg += Bkg->GetBinContent(bin);
            if(nSgn + nBkg == 0) continue;
            float eff = (nSgn + nBkg)/ nTot;
            float w = nBkg / (nSgn + nBkg);
            float p = eff*pow(1-2*w,2);
            x[bin] = varCut;
            yE[bin] = 100*eff ;
            yW[bin] = 100*w;
            yP[bin] = 100*p;

            if(verbose) cout<<x[bin]<<"    "<<nSgn<<"    "<<nBkg<<"    ";
            if(verbose) cout<<yE[bin]<<"    "<<yW[bin]<<"    "<<yP[bin]<<endl;

            if(p > bestP){
                bestCut = varCut;
                bestP = p;
                bestE = eff;
                bestW = w;
            }

        }           
    }


    TGraph* grP = new TGraph(nBins, x, yP);
    TGraph* grE = new TGraph(nBins, x, yE);
    TGraph* grW = new TGraph(nBins, x, yW);

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 1000);
    c1->Divide(1,3);
    TString opt = "AL";
    c1->cd(1); grP->SetTitle("Power"); grP->Draw(opt); 
    c1->cd(2); grE->SetTitle("Efficiency"); grE->Draw(opt);
    c1->cd(3); grW->SetTitle("Mistag"); grW->Draw(opt);

    c1->Print("FOM_" + var + ".png");
    c1->DrawClone();

    cout<<"best cut = "<<bestCut<<endl;
    cout<<"(P="<<100*bestP<<"%, eff="<<100*bestE<<"%, w="<<100*bestW<<"%)"<<endl;

    delete Tot;
    delete Sgn;
    delete Bkg;
    delete grP;
    delete grE;
    delete grW;
    delete c1;
    delete t;
    delete f;

    return;

}
