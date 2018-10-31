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

using namespace std;

/*
pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2016Lists/BsToJpsiPhi_BMuonFilter_2016_DCAP.list hist.root -v outputFile ntu.root -v histoMode RECREATE -v use_gen t -n 10000
*/
PDAnalyzer::PDAnalyzer() {

    std::cout << "new PDAnalyzer" << std::endl;

    // default values can be set in the analyzer class contructor

    setUserParameter( "verbose", "f" );

    setUserParameter( "process", "BsJPsiPhi" );
    setUserParameter( "useTightSel", "f" );

    setUserParameter( "minPtMuon", "2." );
    setUserParameter( "maxEtaMuon", "2.4" );

    setUserParameter( "outputFile", "ntu.root" );

    setUserParameter( "muonIdWpBarrel", "0.00" ); 
    setUserParameter( "muonIdWpEndcap", "0.00" ); 

    setUserParameter( "muoDzCut", "1." ); 
    setUserParameter( "muoPFIsoCut", "5" ); 

    setUserParameter( "mvaMethod", "DNNGlobal2016woIPwIso" ); 

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

    getUserParameter( "minPtMuon", minPtMuon );
    getUserParameter( "maxEtaMuon", maxEtaMuon );

    getUserParameter( "outputFile", outputFile );

    getUserParameter( "muonIdWpBarrel", muonIdWpBarrel ); 
    getUserParameter( "muonIdWpEndcap", muonIdWpEndcap ); 

    getUserParameter( "muoDzCut", muoDzCut ); 
    getUserParameter( "muoPFIsoCut", muoPFIsoCut ); 

    getUserParameter( "mvaMethod", mvaMethod );

    getUserParameter( "ptCut", ptCut ); //needed for paolo's code for unknow reasons

/// additional features
    tWriter = new PDSecondNtupleWriter; // second ntuple
    tWriter->open( getUserParameter("outputFile"), "RECREATE" ); // second ntuple

    setupReader( mvaMethod );

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

    if(useHLT){
        if(process=="BsJPsiPhi"){
            if(hlt(PDEnumString::HLT_Dimuon0_Jpsi3p5_Muon2_v)||hlt(PDEnumString::HLT_Dimuon0_Jpsi_Muon_v)) jpsimu = true;
            if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrkTrk_Displaced_v)) jpsitktk =  true;
            if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrk_Displaced_v)) jpsitk = true;
            if( !(jpsimu || jpsitktk || jpsitk) ) return false;
        }else if(process=="BuJPsiK"){
            if(hlt(PDEnumString::HLT_Dimuon0_Jpsi3p5_Muon2_v)||hlt(PDEnumString::HLT_Dimuon0_Jpsi_Muon_v)) jpsimu = true;
            if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrk_Displaced_v)) jpsitk = true;
            if( !(jpsimu || jpsitk) ) return false;
        }
    }


//------------------------------------------------SEARCH FOR SS---------------------------------------

    int iSsB = GetCandidate(process);
    if(iSsB<0) return false;

    bool isTight = false;
    int iSsBtight = GetTightCandidate(process, -0.1, -0.1);
    if(iSsBtight >= 0){
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

        if(ListLongLivedB.size()!=2) return false; //only evts with two b hadrons            

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

    hmass_ssB->Fill(svtMass->at(iSsB), evtWeight);
    
//-----------------------------------------OPPOSITE SIDE-----------------------------------------

//-----------------------------------------SELECTION-----------------------------------------

    int bestMuIndex=-1;
    float bestMuPt = minPtMuon;
    int nMuonsSel = 0;

    for(int iMuon = 0; iMuon < nMuons; ++iMuon ){

        int itkmu = muonTrack( iMuon, PDEnumString::muInner );
        if(itkmu<0) continue;

        if(std::find(tkSsB.begin(), tkSsB.end(), itkmu) != tkSsB.end()) continue;

        if(muoPt->at( iMuon )<minPtMuon) continue;
        if(abs(muoEta->at( iMuon ))>maxEtaMuon) continue;
       
        if(!isMvaMuon(iMuon, muonIdWpBarrel, muonIdWpEndcap)) continue;

        if(abs(dZ(itkmu, iSsPV)) > muoDzCut) continue;
        if(GetMuoPFiso(iMuon) > muoPFIsoCut)  continue;

        ++nMuonsSel;

        if(muoPt->at( iMuon ) > bestMuPt){
            bestMuPt = muoPt->at( iMuon );
            bestMuIndex = iMuon;
        }
    }


//-----------------------------------------TAGGING VARIABLES-----------------------------------------

    if( bestMuIndex < 0 ){

        (tWriter->osMuon) = 0 ;
        (tWriter->osMuonTag) = -1 ;
        (tWriter->osMuonChargeInfo) = -1;
        (tWriter->evtNumber)= event_tot ;
        tWriter->fill();

        return true;
    }


    (tWriter->osMuon) = 1;
    hmass_ssB_os->Fill(svtMass->at(iSsB), evtWeight);

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


    //COMPLEX TAGGING VARIABLES

    //Muon Cone Charge
    float kappa = 1;
    float drCharge = 0.4;
    float qCone=0, ptCone=0;

    for(int ipf = 0; ipf<nPF; ++ipf){

        float pfpfc = pfcPt->at(ipf);
        float etapfc = pfcEta->at(ipf);

        if( deltaR(etapfc, pfcPhi->at(ipf), muoEta->at(iMuon), muoPhi->at(iMuon)) > drCharge ) continue;

        if(std::find(tkSsB.begin(), tkSsB.end(), pfcTrk->at(ipf)) != tkSsB.end()) continue;

        if(pfpfc < 0.2) continue;
        if(abs(etapfc) > 2.5) continue;

        qCone += pfcCharge->at(ipf) * pow(pfpfc, kappa);
        ptCone += pow(pfpfc, kappa);

    }

    qCone /= ptCone;
    qCone *= trkCharge->at(itkmu); //removing dependency from the muon charge sign
    
    //Jet variables
    int iJet = trkJet->at(itkmu);
    if(iJet<0 && trkPFC->at(itkmu)>=0) iJet=pfcJet->at(trkPFC->at(itkmu));  

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

        TVector3 vjet(jetPx->at(iJet), jetPy->at(iJet), jetPz->at(iJet));
        TVector3 vmu(muoPx->at(iMuon), muoPy->at(iMuon), muoPz->at(iMuon));

        muoJetPt = jetPt->at(iJet);

        muoJetDr = deltaR(jetEta->at(iJet), jetPhi->at(iJet), muoEta->at( iMuon ), muoPhi->at(iMuon));
        
        muoJetEnergyRatio = muoE->at(iMuon) / jetE->at(iJet);

        vjet -= vmu;
        muoJetPtRel = muoPt->at( iMuon ) * (vmu.Unit() * vjet.Unit());

        muoJetSize = jet_pfcs.size();

        muoJetQ = GetJetCharge(iJet, kappa);
        muoJetQ *= trkCharge->at(itkmu); 

        muoJetCSV = jetCSV->at(iJet);

        muoJetDFprob = GetJetProbb(iJet);

    } 

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
        hTest2->Fill( abs(dXY( itkmu, pvtX->at(iSsPV), pvtY->at(iSsPV) )) * dSign(itkmu, jetPx->at( iJet ), jetPy->at( iJet ), pvtX->at(iSsPV), pvtY->at(iSsPV)) );
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

    (tWriter->muoSoftMvaValue) = computeMva(iMuon);

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
    (tWriter->muoQCone) = qCone;

    (tWriter->muoHowMany) = nMuonsSel;


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

    //TAG
    if( TMath::Sign(1, ssBLund) == -1*trkCharge->at(itkmu) ){ 
        hmass_ssB_osRT->Fill(svtMass->at(iSsB), evtWeight);
        (tWriter->osMuonTag) = 1 ;
    }

    if( TMath::Sign(1, ssBLund) != -1*trkCharge->at(itkmu) ){
        hmass_ssB_osWT->Fill(svtMass->at(iSsB), evtWeight);
        (tWriter->osMuonTag) = 0 ;
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
