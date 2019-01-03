TString year_;
TString process_;
TString dir_ = "./";
TString suffix_;
float min_ = 5.1;
float max_ = 5.5;
int nBins_=50;

float *catEdgeL;
float *catEdgeR;
float *catMistag;
int nCat;
TString bestFit;
TF1 *perEvtW;

void setGvars(TString filename)
{
    if(filename.Contains("16")) year_ = "2016";
    if(filename.Contains("17")) year_ = "2017";
    if(filename.Contains("Bs")) process_ = "BsJPsiPhi";
    TString dirPath( filename(0, filename.Index("ntu")) ); 
    dir_ = dirPath;
}

void fitMVA(TString file = "../BsMC/ntuBsMC2017.root"
            , TString method = "BDTOsMuon2017test999"
            , TString mode = "CREATE")
{
    gErrorIgnoreLevel = kWarning;
    if(mode != "CREATE" && mode != "CAT" && mode != "FUNCTION"){
        cout<<"WRONG MODE"<<endl;
        return false;
    }

    TFile *f = new TFile(file);
    TTree *t = (TTree*)f ->Get("PDsecondTree");

    setGvars(file);



//----------READ FILE----------
    if(mode == "CAT")
    {
        std::ifstream ifs("OSMuonTaggerDefinition.txt",std::ifstream::in);
        if(!ifs.is_open()) return;

        ifs >> nCat;
        catEdgeL = new float[nCat];
        catEdgeR = new float[nCat];
        catMistag = new float[nCat];

        for(int i=0; i<nCat; ++i){
            ifs >> catEdgeL[i];
            ifs >> catEdgeR[i];
            ifs >> catMistag[i];
        }

        ifs.close();
    }


    if(mode == "FUNCTION")
    {
        std::ifstream ifs("OSMuonTaggerDefinition.txt",std::ifstream::in);
        if(!ifs.is_open()) return;

        //CAT
        ifs >> nCat;
        catEdgeL = new float[nCat];
        catEdgeR = new float[nCat];
        catMistag = new float[nCat];

        for(int i=0; i<nCat; ++i){
            ifs >> catEdgeL[i];
            ifs >> catEdgeR[i];
            ifs >> catMistag[i];
        }
        //FUNC
        ifs >> bestFit;
        int nPar;
        ifs >> nPar;
        perEvtW = new TF1("perEvtW", bestFit, catEdgeL[0], catEdgeR[nCat]);
        for(int i=0;i<nPar;++i){
            float par;
            ifs >> par;
            perEvtW->SetParameter(i, par);
        }

        ifs.close();
    }

    
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
    TH1F *mva    = new TH1F( "mva",    "mva",    nBinsMva, 0.0, 1.0 );
    TH1F *mva_RT = new TH1F( "mva_RT", "mva_RT", nBinsMva, 0.0, 1.0 );
    TH1F *mva_WT = new TH1F( "mva_WT", "mva_WT", nBinsMva, 0.0, 1.0 );

    std::vector<double> vKDERT;
    std::vector<double> vKDERT_w;
    std::vector<double> vKDEWT;
    std::vector<double> vKDEWT_w;

    //EVENT LOOP
    int nEvents = t->GetEntries();
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

        int evtCat = -1;
        float evtW = -1;
        if(mode=="CAT"){
            if( mvaValue < catEdgeL[0] ){
                evtCat = -1;
                evtW = catMistag[0];
            }else if( mvaValue >= catEdgeR[nCat] ){
                evtCat = -1;
                evtW = catMistag[nCat];
            }else{
                for(int j=0; j<nCat; ++j){
                    if(( mvaValue >= catEdgeL[j]) 
                        && (mvaValue < catEdgeR[j]) )
                    {
                        evtCat = j;
                        evtW = catMistag[j];
                        break;
                    }
                }
            }
        }

        if(mode=="FUNCTION") evtW = perEvtW->Eval(mvaValue);

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

    int nRT = vKDERT.size();
    int nWT = vKDEWT.size();

    auto itMaxRT = max_element(std::begin(vKDERT), std::end(vKDERT));
    auto itMinRT = min_element(std::begin(vKDERT), std::end(vKDERT));
    auto itMaxWT = max_element(std::begin(vKDEWT), std::end(vKDEWT));
    auto itMinWT = min_element(std::begin(vKDEWT), std::end(vKDEWT));

    float xMin = *itMinWT < *itMinRT ? *itMinWT : *itMinRT;
    float xMax = *itMaxRT > *itMaxWT ? *itMaxRT : *itMaxWT;

    cout<<"----- MVA HISTOGRAMS FILLED"<<endl;

//----------CREATE CATEGORIES----------
    if(mode == "CREATE")
    {

        //----------KDE----------
        float rhoRT = pow(nRT,-1./5);
        float rhoWT = pow(nWT,-1./5);

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


        nCat = 25;
        int nTotTagged = nRT + nWT;
        int catSize = nTotTagged / nCat;
        cout<<endl<<nTotTagged<<endl<<catSize<<endl<<endl;

        int cTot = 0;
        int hTot = 0;
        int cCenter = 0;
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
            cTot = cTot + mva_RT->GetBinContent(i) + mva_WT->GetBinContent(i);
            hTot = hTot  + mva->GetBinContent(i);
            //cout<<"bin "<<i<<", "<<mva_RT->GetBinContent(i) + mva_WT->GetBinContent(i);
            //cout<<" ("<<cTot<<")"<<" ["<<hTot<<"]"<<endl;
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
            else     vexl[i] = catCenter[i] - xMin;

            vey[i] = sqrt(((float)catWT[i] * (float)catRT[i])/pow(((float)catWT[i] + (float)catRT[i]),3) );

            cout<<"CAT "<<i+1;
            if(i!=0) cout<<" ["<<catEdge[i-1]<<" - "<<catEdge[i]<<" ] ";
            else     cout<<" ["<<xMin<<" - "<<catEdge[i]<<" ] ";
            cout<<" - eff = "<<100*catEff[i]<<"%, w = "<<100*catW[i];
            cout<<"%, P = "<<100*catP[i]<<"%"<<endl;
            //cout<<" - "<<catCenter[i]<<" - ";
            //cout<<catRT[i] + catWT[i]<<endl;
        }

        //FUNCTION FIT
        TGraph *gr = new TGraphAsymmErrors(nCat,catCenter,catW,vexl,vexh,vey,vey);

        TF1 *fitErf = new TF1("fitErf","[0]+[1]*(1-TMath::Erf([2]+[3]*x))",xMin,xMax);
        fitErf->SetParameter(0, catW[nCat-1]);
        fitErf->SetParameter(1, (catW[0]-catW[nCat-1])/2);

        TF1 *fitErfLim = new TF1("fitErfLim","[0]+[1]*(1-TMath::Erf([2]+[3]*x))",xMin,xMax);
        fitErfLim->SetParameter(0, catW[nCat-1]);
        fitErfLim->SetParameter(1, (catW[0]-catW[nCat-1])/2);
        fitErfLim->SetParLimits(0, 0., catW[nCat-1]);
        fitErfLim->SetParLimits(1, 0., (catW[0]-catW[nCat-1])/2);

        TF1 *fitTanH = new TF1("fitTanH","[0]+[1]*(1-TMath::TanH([2]+[3]*x))",xMin,xMax);
        bestFit = "[0]+[1]*(1-TMath::ATan([2]+[3]*x))";
        TF1 *fitATan = new TF1("fitATan","[0]+[1]*(1-TMath::ATan([2]+[3]*x))",xMin,xMax);

        gr->Fit("fitErf","MELR");
        //gr->Fit("fitErfLim","MELR+");
        gr->Fit("fitTanH","MELR+");
        gr->Fit("fitATan","MELR+");

        cout<<"Chi2"<<endl;
        cout<<"fitErf "<<fitErf->GetChisquare()<<endl;
        cout<<"fitTanH "<<fitTanH->GetChisquare()<<endl;
        cout<<"fitATan "<<fitATan->GetChisquare()<<endl;

        //OUTPUT STREAM
        ofstream ofs;
        ofs.open ("OSMuonTaggerDefinition.txt");
        //CAT
        ofs<<nCat<<endl;
        for(int i=0; i<nCat; ++i)
        {
            if(i!=0) ofs<<catEdge[i-1];
            else     ofs<<xMin;
            ofs<<" "<<catEdge[i]<<" ";
            ofs<<catW[i]<<endl;

        }
        //FUNCTION
        ofs<<bestFit<<endl;
        ofs<<fitATan->GetNumberFreeParameters()<<endl;
        for(int j=0;j<fitATan->GetNumberFreeParameters();++j)
        {
            ofs<<fitATan->GetParameter(j)<<endl;
        }
        ofs.close();

        avgD /= totEff;
        float avgW = (1-avgD)/2;

        cout<<endl<<"Cat Eff = "<<100*totEff<<"%"<<endl;
        cout<<"Cat W = "<<100*avgW<<"%"<<endl;
        cout<<"Cat P = "<<100*totP<<"%"<<endl;

        TCanvas *c2 = new TCanvas();
        mva_WT->Draw("hist");
        mva_RT->Draw("hist same");
        mva_WT->SetLineColor(kRed);

        TCanvas *c3 = new TCanvas();
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(.5);
        gr->SetTitle("per-event mistag");
        gr->GetXaxis()->SetTitle("MVA output");
        gr->Draw("APE");
        fW->SetNpx(nCat);
        fW->SetLineColor(kGreen);
        fW->DrawClone("P SAME");
        c3->BuildLegend();
    }

    cout<<endl<<"NoCat Eff = "<<100.*(float)(nRT+nWT)/nEvents<<"%"<<endl;
    cout<<"NoCat W = "<<100.*(float)nWT/(nRT+nWT)<<"%"<<endl;
    cout<<"NoCat P = "<<100.*((float)(nRT+nWT)/nEvents)*pow(1.-2.*((float)nWT/(nRT+nWT)),2)<<"%"<<endl;


}
