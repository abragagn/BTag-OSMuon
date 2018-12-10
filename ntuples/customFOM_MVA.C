float calculateFOM( TString eta = "B", TString addCut = "", bool verbose = false )
{

	gStyle->SetOptTitle(1); gStyle->SetOptStat(0);
	gStyle->SetOptFit(0); 	gStyle->SetStatBorderSize(0);
	gStyle->SetStatX(.49); 	gStyle->SetStatY(.89);

	TFile *f = new TFile("./BsMC/ntuBsMC2017.root");
	TTree *t = (TTree*)f ->Get("PDsecondTree");

	int nBins=50;
	float min=0.86;
	float max=0.92;

	TH1F *Tot = new TH1F( "Tot", "Tot", 100, 5.0, 6.0 );
	TH1F *Sgn = new TH1F( "Sgn", "Sgn", nBins, min, max );
	TH1F *Bkg = new TH1F( "Bkg", "Bkg", nBins, min, max );

	TCut generalCuts;
	if(eta == "B") generalCuts = "abs(muoEta)<1.2&&osMuon==1";
	if(eta == "E") generalCuts = "abs(muoEta)>=1.2&&osMuon==1";
	generalCuts += addCut;
	TCut sgnDef = generalCuts + "osMuonTag==1";
	TCut bkgDef = generalCuts + "osMuonTag==0";

	TString var = "muoSoftMvaValue";

	t->Project("Tot", "ssbMass", "");
	t->Project("Sgn", var, sgnDef);
	t->Project("Bkg", var, bkgDef);

	int nTot = Tot->GetEntries();
	float nSgn=0, nBkg=0;

	float *x = new float[nBins];
	float *yE = new float[nBins];
	float *yW = new float[nBins];
	float *yP = new float[nBins];
	float bestCut = 0;
	float bestP = 0, bestE = 0, bestW = 0.5;

	if(verbose) cout<<endl<<"nTot = "<<nTot<<endl;
	if(verbose) cout<<"x[bin]  nSgn  nBkg  yE[bin]  yW[bin]  yP[bin]"<<endl<<endl;

	for(int bin=nBins; bin>0; --bin){

		float varCut = Sgn->GetBinLowEdge(bin);

		nSgn += Sgn->GetBinContent(bin);
		nBkg += Bkg->GetBinContent(bin);
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

	TGraph* grP = new TGraph(nBins, x, yP);
	TGraph* grE = new TGraph(nBins, x, yE);
	TGraph* grW = new TGraph(nBins, x, yW);

	TCanvas *c1 = new TCanvas("c1", "c1", 600, 1000);
	c1->Divide(1,3);
	TString opt = "AC";
	c1->cd(1); grP->SetTitle("Power"); grP->Draw(opt); 
	c1->cd(2); grE->SetTitle("Efficiency"); grE->Draw(opt);
	c1->cd(3); grW->SetTitle("Mistag"); grW->Draw(opt);

	c1->Print("FOM_" + eta + ".png");

	cout<<"best cut = "<<bestCut<<endl;
	cout<<"(P="<<100*bestP<<"%, eff="<<100*bestE<<"%, w="<<100*bestW<<"%)"<<endl;

	return bestP;

}

void customFOM_MVA(TString addCut = "", bool verbose = false)
{
	cout<<"    BARREL"<<endl;
	float pB = calculateFOM("B", addCut, verbose);
	cout<<endl<<"    ENDCAP"<<endl;
	float pE = calculateFOM("E", addCut, verbose);

	cout<<endl<<"Total P = "<<100*(pB+pE)<<endl;

	return;
}

