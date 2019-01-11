#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>

#include "PDAnalyzer.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TFile.h"

// additional features
#include "PDSecondNtupleWriter.h"
#include "PDMuonVar.cc"
#include "PDSoftMuonMvaEstimator.cc"
#include "AlbertoUtil.cc"
#include "OSMuonMvaTag.cc"

using namespace std;

/*
pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2016Lists/BsToJpsiPhi_BMuonFilter_2016_DCAP.list hist.root -v outputFile ntu.root -v histoMode RECREATE -v use_gen t -n 10000
*/
PDAnalyzer::PDAnalyzer() {

    std::cout << "new PDAnalyzer" << std::endl;

    // default values can be set in the analyzer class contructor

    setUserParameter( "verbose", "f" );

    setUserParameter( "process", "BsJPsiPhi" );
    setUserParameter( "useHLT", "false" );

    setUserParameter( "outputFile", "ntu.root" );

    setUserParameter( "muonIdWpBarrel", "0.00" ); 
    setUserParameter( "muonIdWpEndcap", "0.00" ); 

    setUserParameter( "muonMvaMethod",      "BDTMuonID2017woIPwIso" ); 
    setUserParameter( "osMuonTagMvaMethod", "BDTOsMuon2016Jet" ); 

    setUserParameter( "ptCut", "40.0" ); //needed for paolo's code for unknow reasons

}


PDAnalyzer::~PDAnalyzer() {
}



void PDAnalyzer::beginJob() {

    PDAnalyzerUtil::beginJob();

    // user parameters are retrieved as strings by using their names;
    // numeric parameters ( int, float or whatever ) can be directly set
    // by passing the corresponding variable,
    // e.g. getUserParameter( "name", x )

    getUserParameter( "verbose", verbose );

    getUserParameter( "process", process );
    getUserParameter( "useHLT", useHLT );

    getUserParameter( "outputFile", outputFile );

    getUserParameter( "muonIdWpBarrel", muonIdWpBarrel ); 
    getUserParameter( "muonIdWpEndcap", muonIdWpEndcap ); 

    getUserParameter( "muonMvaMethod", muonMvaMethod );
    getUserParameter( "osMuonTagMvaMethod", osMuonTagMvaMethod );

    getUserParameter( "ptCut", ptCut ); //needed for paolo's code for unknow reasons

//  additional features
    tWriter = new PDSecondNtupleWriter; // second ntuple
    tWriter->open( getUserParameter("outputFile"), "RECREATE" ); // second ntuple

    setOsMuonCuts(muonIdWpBarrel, muonIdWpEndcap, 1. );

    inizializeMuonMvaReader( muonMvaMethod );
    //inizializeOSMuonMvaTagReader( osMuonTagMvaMethod );

    if(process=="BsJPsiPhi") SetBsMassRange(5.20, 5.50);
    if(process=="BuJPsiK") SetBuMassRange(5.1, 5.50);

    return;

}


void PDAnalyzer::book() {

    // putting "autoSavedObject" in front of the histo creation 
    // it's automatically marked for saving on file; the option 
    // is uneffective when not using the full utility

    float min = 5.0;
    float max = 5.5;
    float nbin = 250;


    autoSavedObject =
    hmass_ssB       = new TH1D( "hmass_ssB", "hmass_ssB", nbin, min, max );

    autoSavedObject =
    hmass_ssB_os    = new TH1D( "hmass_ssB_os", "hmass_ssB_os", nbin, min, max );

    autoSavedObject =
    hmass_ssB_osWT  = new TH1D( "hmass_ssB_osWT", "hmass_ssB_osWT", nbin, min, max );

    autoSavedObject =
    hmass_ssB_osRT  = new TH1D( "hmass_ssB_osRT", "hmass_ssB_osRT", nbin, min, max );

    autoSavedObject =
    hmass_ssB_osCC  = new TH1D( "hmass_ssB_osCC", "hmass_ssB_osCC", nbin, min, max );
    
    autoSavedObject =
    hmass_ssB_osRC  = new TH1D( "hmass_ssB_osRC", "hmass_ssB_osRC", nbin, min, max );
    
    autoSavedObject =
    hmass_ssB_osWC  = new TH1D( "hmass_ssB_osWC", "hmass_ssB_osWC", nbin, min, max );

    autoSavedObject =
    hTest   = new TH1D( "hTest", "hTest", 5, 0, 5 );

    autoSavedObject =
    hTest2   = new TH1D( "hTest2", "hTest2", 500, -0.3, 0.3 );

    return;

}


void PDAnalyzer::reset() {

    autoReset();
    return;
}


bool PDAnalyzer::analyze( int entry, int event_file, int event_tot ) {

    if ( verbose ) {
        cout << " +++++++++++++++++++++++++++ " << endl;
        cout << "entry: "
             << entry << " " << event_file << " " << event_tot << endl;
        cout << "run: " << runNumber << " , "
             << "evt: " << eventNumber << endl;
    }
    else {
        if ( (!(event_tot%10) && event_tot<100 ) || 
        (!(event_tot %100) && event_tot<1000 ) || 
        (!(event_tot %1000)&& event_tot<10000 ) || 
        (!(event_tot %10000) && event_tot<100000 ) || 
        (!(event_tot %100000) && event_tot<1000000 ) || 
        (!(event_tot %1000000) && event_tot<10000000 ) )
            cout << " == at event " << event_file << " " << event_tot << endl;
    }

// additional features
    computeMuonVar();
    tWriter->Reset();
    convSpheCart(jetPt, jetEta, jetPhi, jetPx, jetPy, jetPz);
    convSpheCart(muoPt, muoEta, muoPhi, muoPx, muoPy, muoPz);
    convSpheCart(trkPt, trkEta, trkPhi, trkPx, trkPy, trkPz);
    convSpheCart(pfcPt, pfcEta, pfcPhi, pfcPx, pfcPy, pfcPz);

    if( !((process=="BsJPsiPhi")||(process=="BuJPsiK")) ) {
        cout<<"!$!#$@$% PROCESS NAME WRONG"<<endl;
        return false;
    }

//------------------------------------------------HLT---------------------------------------

    bool jpsimu = false;
    bool jpsitktk = false;
    bool jpsitk = false;

    if(hlt(PDEnumString::HLT_Dimuon0_Jpsi3p5_Muon2_v)||hlt(PDEnumString::HLT_Dimuon0_Jpsi_Muon_v)) jpsimu = true;
    if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrkTrk_Displaced_v)) jpsitktk =  true;
    if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrk_Displaced_v)) jpsitk = true;

    if( jpsimu ) SetJpsiMuCut();
    if( !jpsimu ) SetJpsiTrktrkCut();


    if(useHLT && process=="BsJPsiPhi" && !(jpsimu || jpsitktk)) return false;
    if(useHLT && process=="BuJPsiK" && !(jpsimu || jpsitk)) return false;


//------------------------------------------------SEARCH FOR SS---------------------------------------

    int iSsB = GetCandidate(process);
    if(iSsB<0) return false;

    bool isTight = false;
    int iSsBtight = GetTightCandidate(process);
    if(iSsBtight>=0){
        isTight = true;
        iSsB = iSsBtight;
    }

    int iJPsi = (subVtxFromSV(iSsB)).at(0);
    vector <int> tkJpsi = tracksFromSV(iJPsi);
    vector <int> tkSsB = tracksFromSV(iSsB);

    TLorentzVector tB = GetTLorentzVecFromJpsiX(iSsB);

    //generation information

    vector <int> ListLongLivedB;
    vector <int> ListB;
    int genBindex = -1;
    int ssBLund = 0;
    int tagMix = -1;
    float evtWeight = 1;

    if(use_gen){
        for( unsigned int i=0 ; i<genId->size() ; ++i ){
            if( IsB(i) ) ListB.push_back(i);
            if(TagMixStatus( i ) == 2) continue;
            unsigned int Code = abs(genId->at(i));
            if( Code == 511 || Code == 521 || Code == 531 || Code == 541 || Code == 5122 ) ListLongLivedB.push_back(i);
        }

        genBindex = GetClosestGenLongLivedB( tB.Eta(), tB.Phi(), tB.Pt(), &ListLongLivedB);
        if(genBindex<0) return false;

        ssBLund = genId->at(genBindex);
        if((process=="BsJPsiPhi") && (abs(ssBLund)!=531)) return false;
        if((process=="BuJPsiK") && (abs(ssBLund)!=521)) return false;

        tagMix = TagMixStatus( genBindex );
        if(tagMix == 2) return false;
        if(tagMix == 1) ssBLund*=-1;

        for(auto it:ListLongLivedB){
            if(it == genBindex) continue;
            if(abs(genId->at(it)) == abs(ssBLund)) evtWeight = 2;
        }


    }else{
        if(process=="BsJPsiPhi") ssBLund = ((double)rand() / (RAND_MAX)) < 0.5 ? +531 : -531; //this should not be used
        if(process=="BuJPsiK"){
            for( unsigned int i=0; i<tkSsB.size(); ++i ){
            if( tkSsB[i] == tkJpsi[0] || tkSsB[i] == tkJpsi[1] ) continue;
            ssBLund = trkCharge->at(tkSsB[i]) > 0 ? +521 : -521;
            }
        }
    }


    int iSsPV = GetBestPV(iSsB, tB);
    if(iSsPV < 0) return false;

    setSsForTag(iSsB, iSsPV);

    //FILLING SS
    (tWriter->ssbPt) = tB.Pt();
    (tWriter->ssbEta) = tB.Eta();
    (tWriter->ssbPhi) = tB.Phi();
    (tWriter->ssbMass) = svtMass->at(iSsB);
    (tWriter->ssbIsTight) = isTight;

    (tWriter->ssbLxy) = GetCt2D(tB, iSsB) / (MassBs/tB.Pt());
    (tWriter->ssbCt2D) = GetCt2D(tB, iSsB);
    (tWriter->ssbCt2DErr) = GetCt2DErr(tB, iSsB, iSsPV);
    (tWriter->ssbCt2DSigmaUnit) = GetCt2D(tB, iSsB, iSsPV)/GetCt2DErr(tB, iSsB, iSsPV);
    (tWriter->ssbCt3D) = GetCt3D(tB, iSsB, iSsPV);
    (tWriter->ssbCt3DErr) = GetCt3DErr(tB, iSsB, iSsPV);
    (tWriter->ssbCt3DSigmaUnit) = GetCt3D(tB, iSsB, iSsPV)/GetCt3DErr(tB, iSsB, iSsPV);

    (tWriter->ssbSVT) = iSsB;
    (tWriter->ssbPVT) = iSsPV;
    
    (tWriter->ssbLund) = ssBLund;

    (tWriter->hltJpsiMu) = jpsimu;
    (tWriter->hltJpsiTrkTrk) = jpsitktk;
    (tWriter->hltJpsiTrk) = jpsitk;
    
    (tWriter->evtWeight) = evtWeight;
    (tWriter->evtNb) = ListLongLivedB.size();

    hmass_ssB->Fill(svtMass->at(iSsB), evtWeight);
    
//-----------------------------------------OPPOSITE SIDE-----------------------------------------

    int bestMuIndex = getOsMuon();
    int tagDecision = getOsMuonTag();

    if( tagDecision == 0 ){
        (tWriter->osMuon) = 0 ;
        (tWriter->osMuonTag) = -1 ;
        (tWriter->osMuonChargeInfo) = -1;
        (tWriter->evtNumber)= event_tot ;
        tWriter->fill();
        return true;
    }

    (tWriter->osMuon) = 1;
    //(tWriter->osMuonTagMvaValue) = getOsMuonTagMvaValue();
    //(tWriter->osMuonTagMistag) = 

    hmass_ssB_os->Fill(svtMass->at(iSsB), evtWeight);

    if( TMath::Sign(1, ssBLund) == tagDecision ){ 
        hmass_ssB_osRT->Fill(svtMass->at(iSsB), evtWeight);
        (tWriter->osMuonTag) = 1 ;
    }

    if( TMath::Sign(1, ssBLund) != tagDecision ){
        hmass_ssB_osWT->Fill(svtMass->at(iSsB), evtWeight);
        (tWriter->osMuonTag) = 0 ;
    }

    //COMPLEX TAGGING VARIABLES
    //INDICES
    int iMuon = bestMuIndex;
    int itkmu = muonTrack( iMuon, PDEnumString::muInner );
    //GEN INFO
    int genMuIndex = -1;
    int muoLund=0, muoAncestor=-1; 
    
    genMuIndex = GetClosestGen( muoEta->at(iMuon), muoPhi->at(iMuon), muoPt->at(iMuon) );
    if( genMuIndex >= 0 ) {
        muoLund = genId->at(genMuIndex);
        muoAncestor = GetAncestor( genMuIndex, &ListB ); 
    }

    //Cone
    float kappa = 1;
    float drCharge = 0.4;

    //JET variables
    int iJet = trkJet->at(itkmu);
    if(iJet<0 && trkPFC->at(itkmu)>=0) iJet=pfcJet->at(trkPFC->at(itkmu));  
    TVector3 vMu(muoPx->at(iMuon), muoPy->at(iMuon), muoPz->at(iMuon));

    float muoJetPtRel = -1;
    float muoJetDr = -1;
    float muoJetEnergyRatio = -1;
    float muoJetCSV = -1;
    float muoJetDFprob = -1;
    float muoJetSize = -1;
    float muoJetQ = -1;
    float muoJetPt = -1;

    if(iJet>=0){
        vector <int> jet_pfcs = pfCandFromJet( iJet );
        TVector3 vJet(jetPx->at(iJet), jetPy->at(iJet), jetPz->at(iJet));
        muoJetPt = jetPt->at(iJet);
        muoJetDr = deltaR(jetEta->at(iJet), jetPhi->at(iJet), muoEta->at( iMuon ), muoPhi->at(iMuon));
        muoJetEnergyRatio = muoE->at(iMuon) / jetE->at(iJet);
        vJet -= vMu;
        muoJetPtRel = muoPt->at( iMuon ) * (vMu.Unit() * vJet.Unit());
        muoJetSize = jet_pfcs.size();
        muoJetQ = GetJetCharge(iJet, kappa);
        muoJetQ *= trkCharge->at(itkmu); 
        muoJetCSV = jetCSV->at(iJet);
        muoJetDFprob = GetJetProbb(iJet);
    }

    //SVT variables

    int osSvt = GetBestSvtFromTrack(itkmu);

    float muoSvtPtRel = -1;
    float muoSvtDr = -1;
    float muoSvtEnergyRatio = -1;
    float muoSvtCSV = -1;
    float muoSvtDFprob = -1;
    float muoSvtSize = -1;
    float muoSvtQ = -1;
    float muoSvtPt = -1;

    if(osSvt>=0){
        vector <int> tkSvt = tracksFromSV(osSvt);
        TVector3 vSvt;
        for(auto it:tkSvt){
            TVector3 v;
            v.SetXYZ(trkPx->at(it), trkPy->at(it), trkPz->at(it));
            vSvt += v;
        }
        muoSvtPt = vSvt.Pt();
        muoSvtDr = deltaR(vSvt.Eta(), vSvt.Phi(), muoEta->at( iMuon ), muoPhi->at(iMuon));
        muoSvtEnergyRatio = -1;
        muoSvtCSV = -1;
        muoSvtDFprob = -1;
        vSvt -= vMu;
        muoSvtPtRel = muoPt->at( iMuon ) * (vMu.Unit() * vSvt.Unit());
        muoSvtSize = svtNTracks->at(osSvt);
        muoSvtQ = GetSvtCharge(osSvt, kappa);
        muoSvtQ *= trkCharge->at(itkmu);
    }

    //CONE variables

    float muoConePtRel = -1;
    float muoConeDr = -1;
    float muoConeEnergyRatio = -1;
    float muoConeCSV = -1;
    float muoConeDFprob = -1;
    float muoConeSize = 0;
    float muoConeQ = -1;
    float muoConePt = -1;

    TLorentzVector tCone, tMu;
    tCone.SetPtEtaPhiM(0.,0.,0.,0.);
    tMu.SetPtEtaPhiM(muoPt->at( iMuon ), muoEta->at( iMuon ), muoPhi->at( iMuon ), MassMu);
    float qCone=0, ptCone=0;

    for(int ipf=0; ipf<nPF; ++ipf){
        float pfpfc = pfcPt->at(ipf);
        float etapfc = pfcEta->at(ipf);
        if( deltaR(etapfc, pfcPhi->at( ipf ), muoEta->at( iMuon ), muoPhi->at( iMuon )) > drCharge) continue;
        if(std::find(tkSsB.begin(), tkSsB.end(), pfcTrk->at(ipf)) != tkSsB.end()) continue;
        if(pfpfc < 0.2) continue;
        if(fabs(etapfc) > 2.5) continue;        
        TLorentzVector a;
        a.SetPxPyPzE(pfcPx->at(ipf), pfcPy->at(ipf), pfcPz->at(ipf), pfcE->at(ipf));
        tCone += a;
        ++muoConeSize;
        qCone += pfcCharge->at(ipf) * pow(pfpfc, kappa);
        ptCone += pow(pfpfc, kappa);
    }

    if(ptCone != 0) qCone /= ptCone;
    else qCone = 1;
    qCone *= trkCharge->at(itkmu);

    muoConePt = tCone.Pt();
    muoConeDr = deltaR(tCone.Eta(), tCone.Phi(), muoEta->at( iMuon ), muoPhi->at(iMuon));
    muoConeEnergyRatio = muoE->at(iMuon) / tCone.E();
    muoConeCSV = -1;
    muoConeDFprob = -1;
    tCone -= tMu;
    muoConePtRel = muoPt->at( iMuon ) * (tMu.Vect().Unit() * tCone.Vect().Unit());
    muoConeQ = qCone;

    bool debugJet = false;

    if(debugJet && muoAncestor>=0 && iJet>=0){

        vector <int> jet_pfcs = pfCandFromJet( iJet );
        cout<<endl;
        printDaughterTree(muoAncestor, "");
        cout<<endl;

        cout<<"dXY = "<<GetSignedDxy(iMuon, iSsPV)<<endl;
        cout<<"osB: "<<genPt->at( muoAncestor )<<" "<<genEta->at( muoAncestor )<<" "<<genPhi->at( muoAncestor )<<endl;
        cout<<"jet: "<<jetPt->at( iJet )<<" "<<jetEta->at( iJet )<<" "<<jetPhi->at( iJet )<<endl;
        cout<<"muo: "<<muoPt->at( iMuon )<<" "<<muoEta->at( iMuon )<<" "<<muoPhi->at( iMuon )<<endl;

        
        for(int it:jet_pfcs){
            int g = GetClosestGen( pfcEta->at(it), pfcPhi->at(it), pfcPt->at(it) );
            if(g>=0) cout<<genId->at(g)<<" "; else cout<<"noGen ";
        }
        cout<<endl;

        hTest->Fill(GetSignedDxy(iMuon, iSsPV));
        hTest2->Fill( fabs(dXY( itkmu, pvtX->at(iSsPV), pvtY->at(iSsPV) )) * dSign(itkmu, jetPx->at( iJet ), jetPy->at( iJet ), pvtX->at(iSsPV), pvtY->at(iSsPV)) );
    }

    //------------------------------------------------FILLING------------------------------------------------
    (tWriter->muoPt) = muoPt->at( iMuon );
    (tWriter->muoEta) = muoEta->at( iMuon );
    (tWriter->muoPhi) = muoPhi->at(iMuon);
    (tWriter->muoCharge) = trkCharge->at(itkmu);

    (tWriter->muoDxy) = GetSignedDxy(iMuon, iSsPV);
    (tWriter->muoDz) = dZ(itkmu, iSsPV);
    (tWriter->muoExy) = trkExy->at(itkmu);
    (tWriter->muoEz) = trkEz->at(itkmu);

    (tWriter->muoSoftMvaValue) = computeMuonMva(iMuon);

    (tWriter->muoLund) = muoLund;
    (tWriter->muoAncestor) = muoAncestor; 

    (tWriter->muoDrB) = deltaR(tB.Eta(), tB.Phi(), muoEta->at(iMuon), muoPhi->at(iMuon));
    (tWriter->muoPFIso) = GetMuoPFiso(iMuon);

    (tWriter->muoJetPt) = muoJetPt;
    (tWriter->muoJetPtRel) = muoJetPtRel;
    (tWriter->muoJetDr) = muoJetDr;
    (tWriter->muoJetEnergyRatio) = muoJetEnergyRatio;
    (tWriter->muoJetQ) = muoJetQ;
    (tWriter->muoJetCSV) = muoJetCSV;
    (tWriter->muoJetDFprob) = muoJetDFprob;
    (tWriter->muoJetSize) = muoJetSize;

    (tWriter->muoHowMany) = getNosMuons();

    (tWriter->muoSvtPt)  = muoSvtPt;
    (tWriter->muoSvtPtRel) = muoSvtPtRel;
    (tWriter->muoSvtDr) = muoSvtDr;
    (tWriter->muoSvtEnergyRatio) = muoSvtEnergyRatio;
    (tWriter->muoSvtQ) = muoSvtQ;
    (tWriter->muoSvtCSV) = muoSvtCSV;
    (tWriter->muoSvtDFprob) = muoSvtDFprob;
    (tWriter->muoSvtSize) = muoSvtSize;

    (tWriter->muoConePt) = muoConePt;
    (tWriter->muoConePtRel) = muoConePtRel;
    (tWriter->muoConeDr) = muoConeDr;
    (tWriter->muoConeEnergyRatio) = muoConeEnergyRatio;
    (tWriter->muoConeQ) = muoConeQ;
    (tWriter->muoConeCSV) = muoConeCSV;
    (tWriter->muoConeDFprob) = muoConeDFprob;
    (tWriter->muoConeSize) = muoConeSize;

    //------------------------------------------------TAG------------------------------------------------

    //CHARGE CORRELATION 
    if( muoAncestor >=0 ){
        if( TMath::Sign(1, ssBLund) == -1*trkCharge->at(itkmu) ){
            hmass_ssB_osCC->Fill(svtMass->at(iSsB), evtWeight);
            (tWriter->osMuonChargeInfo) = 1 ;
        }else{
            hmass_ssB_osWC->Fill(svtMass->at(iSsB), evtWeight);
            (tWriter->osMuonChargeInfo) = 0 ;
        }
    }else{
        hmass_ssB_osRC->Fill(svtMass->at(iSsB), evtWeight);
        (tWriter->osMuonChargeInfo) = 2 ;
    }

    (tWriter->evtNumber)=( event_tot );
    tWriter->fill();

    return true;

}

void PDAnalyzer::endJob() {

// additional features
    tWriter->close();   // second ntuple

    cout<<endl;

    cout<<"-----TAG RESULTS-----"<<endl;

    float eff = CountEventsWithFit(hmass_ssB_os, process) / CountEventsWithFit(hmass_ssB, process);
    float w  = CountEventsWithFit(hmass_ssB_osWT, process) / CountEventsWithFit(hmass_ssB_os, process);
    float power = eff*pow(1-2*w, 2);

    float tot = CountEventsWithFit(hmass_ssB_osCC, process ) + CountEventsWithFit(hmass_ssB_osWC, process) + CountEventsWithFit(hmass_ssB_osRC, process);

    cout<<"CC   WC  RC"<<endl;
    cout<< CountEventsWithFit(hmass_ssB_osCC, process)/tot<<" "<< CountEventsWithFit(hmass_ssB_osWC, process)/tot<<" "<< CountEventsWithFit(hmass_ssB_osRC, process)/tot<<endl<<endl;

    cout<<"#B    eff%    w%    P%"<<endl;
    cout<< CountEventsWithFit(hmass_ssB, process)<<" "<<eff*100<<" "<<w*100<<" "<<power*100<<endl;

    return;
}


void PDAnalyzer::save() {
#   if UTIL_USE == FULL
    // explicit saving not necessary for "autoSavedObjects"
    autoSave();
#elif UTIL_USE == BARE
    // explicit save histos when not using the full utility

#endif

    return;
}

// ======MY FUNCTIONS===============================================================================
int PDAnalyzer::GetBestSvtFromTrack(int trkIndex )
{
    vector <int> svtList = sVtsWithTrack( trkIndex );
    int index = -1;
    float bestChi2 = 1e9;

    for(auto it : svtList){
        if( svtChi2->at(it)/svtNDOF->at(it) > bestChi2 ) continue;
        index = it;
        bestChi2 = svtChi2->at(it)/svtNDOF->at(it);
    }

    return index;
}
float PDAnalyzer::GetSvtCharge(int iSvt, float kappa)
{

    float QSvt = 0;
    float ptSvt = 0;

    vector <int> list = tracksFromSV(iSvt);

    for(int it:list){

       float pt = trkPt->at(it);

       if(pt<0.2) continue;
       if(fabs(trkEta->at(it))>2.5) continue;

       QSvt += trkCharge->at(it) * pow(pt, kappa);
       ptSvt += pow(pt, kappa);

    }

    QSvt /= ptSvt;

    return QSvt; 

}
