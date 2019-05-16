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
pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2017Lists/BsToJpsiPhi_2017_DCAP.list hist.root -v outputFile ntu.root -v histoMode RECREATE -v use_gen t -n 10000
*/
PDAnalyzer::PDAnalyzer() {

    std::cout << "new PDAnalyzer" << std::endl;

    // default values can be set in the analyzer class contructor

    setUserParameter( "verbose", "f" );

    setUserParameter( "process", "BsJPsiPhi" );
    setUserParameter( "useHLT", "true" );
    setUserParameter( "writeVars", "true" );

    setUserParameter( "outputFile", "ntu.root" );

    setUserParameter( "muonIdWp", "0.35" ); 

    setUserParameter( "muonMvaMethod",      "DNNMuonIDFull2017woIPwIso" ); 
    setUserParameter( "osMuonTagMvaMethod", "DNNOsMuonHLTJpsiMu_test241" ); 

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
    getUserParameter( "writeVars", writeVars );

    getUserParameter( "outputFile", outputFile );

    getUserParameter( "muonIdWp", muonIdWp ); 

    getUserParameter( "muonMvaMethod", muonMvaMethod );
    getUserParameter( "osMuonTagMvaMethod", osMuonTagMvaMethod );

    getUserParameter( "ptCut", ptCut ); //needed for paolo's code for unknow reasons

//  additional features
    tWriter = new PDSecondNtupleWriter; // second ntuple
    tWriter->open( getUserParameter("outputFile"), "RECREATE" ); // second ntuple

    setOsMuonMvaCut(muonIdWp);
    setOsMuonDzCut( 1.0 );

    inizializeMuonMvaReader();
    inizializeOSMuonMvaTagReader( osMuonTagMvaMethod );
    bool osInit = inizializeOSMuonMvaMistagMethods();
    if(!osInit) cout<<"METHOD NOT INIZIALIZATED. ABORT"<<endl;

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
    hTest   = new TH1D( "hTest", "hTest", 100, 0, 1 );

    autoSavedObject =
    hTest2   = new TH1D( "hTest2", "hTest2", 100, 0, 1 );

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
    inizializeTagVariables();
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

    if( !jpsimu ) return false;
    SetJpsiMuCut();

//------------------------------------------------SEARCH FOR SS---------------------------------------

    int ssbSVT = GetCandidate(process);
    if(ssbSVT<0) return false;

    bool isTight = false;
    int ssbSVTtight = GetTightCandidate(process);
    if(ssbSVTtight>=0){
        isTight = true;
        ssbSVT = ssbSVTtight;
    }

    int iJPsi = (subVtxFromSV(ssbSVT)).at(0);
    vector <int> tkJpsi = tracksFromSV(iJPsi);
    vector <int> tkSsB = tracksFromSV(ssbSVT);

    TLorentzVector tB = GetTLorentzVecFromJpsiX(ssbSVT);

    //generation information

    vector <int> ListLongLivedB;
    vector <int> ListB;
    int genBindex = -1;
    int ssBLund = 0;
    int tagMix = -1;
    float evtWeight = 1;

    if(use_gen){
        for( uint i=0 ; i<genId->size() ; ++i ){
            if(TagMixStatus( i ) == 2) continue;
            if( IsB(i) ) ListB.push_back(i);
            uint Code = abs(genId->at(i));
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
            for( auto it:tkSsB ){
                if( it == tkJpsi[0] || it == tkJpsi[1] ) continue;
                ssBLund = trkCharge->at(it) > 0 ? +521 : -521;
            }
        }
    }


    int ssbPVT = GetBestPV(ssbSVT, tB);
    if(ssbPVT < 0) return false;

    setSsForTag(ssbSVT, ssbPVT);

    //FILLING SS
    (tWriter->ssbPt) = tB.Pt();
    (tWriter->ssbEta) = tB.Eta();
    (tWriter->ssbPhi) = tB.Phi();
    (tWriter->ssbMass) = svtMass->at(ssbSVT);
    (tWriter->ssbIsTight) = isTight;

    (tWriter->ssbLxy) = GetCt2D(tB, ssbSVT) / (MassBs/tB.Pt());
    (tWriter->ssbCt2D) = GetCt2D(tB, ssbSVT);
    (tWriter->ssbCt2DErr) = GetCt2DErr(tB, ssbSVT, ssbPVT);
    (tWriter->ssbCt2DSigmaUnit) = GetCt2D(tB, ssbSVT, ssbPVT)/GetCt2DErr(tB, ssbSVT, ssbPVT);
    (tWriter->ssbCt3D) = GetCt3D(tB, ssbSVT, ssbPVT);
    (tWriter->ssbCt3DErr) = GetCt3DErr(tB, ssbSVT, ssbPVT);
    (tWriter->ssbCt3DSigmaUnit) = GetCt3D(tB, ssbSVT, ssbPVT)/GetCt3DErr(tB, ssbSVT, ssbPVT);

    (tWriter->ssbSVT) = ssbSVT;
    (tWriter->ssbPVT) = ssbPVT;
    
    (tWriter->ssbLund) = ssBLund;

    (tWriter->hltJpsiMu) = jpsimu;
    (tWriter->hltJpsiTrkTrk) = jpsitktk;
    (tWriter->hltJpsiTrk) = jpsitk;
    
    (tWriter->evtWeight) = evtWeight;
    (tWriter->evtNb) = ListLongLivedB.size();

    hmass_ssB->Fill(svtMass->at(ssbSVT), evtWeight);
    
//-----------------------------------------OPPOSITE SIDE-----------------------------------------

    int bestMuIndex = getOsMuon();
    int tagDecision = getOsMuonTag();

    if( tagDecision == 0 ){
        (tWriter->osMuon) = 0;
        (tWriter->osMuonTag) = -1;
        (tWriter->osMuonChargeInfo) = -1;
        (tWriter->evtNumber) = event_tot;
        tWriter->fill();
        return true;
    }

    (tWriter->osMuon) = 1;
    float osMuonTagMvaValue = getOsMuonTagMvaValue();
    pair<float,float> osMuonTagMistag = getOsMuonTagMistagProb(2);
    (tWriter->osMuonTagMvaValue) = osMuonTagMvaValue;
    (tWriter->osMuonTagMistag) = osMuonTagMistag.first;

    hTest->Fill(osMuonTagMistag.first);

    hmass_ssB_os->Fill(svtMass->at(ssbSVT), evtWeight);

    if( TMath::Sign(1, ssBLund) == tagDecision ){ 
        hmass_ssB_osRT->Fill(svtMass->at(ssbSVT), evtWeight);
        (tWriter->osMuonTag) = 1 ;
    }

    if( TMath::Sign(1, ssBLund) != tagDecision ){
        hmass_ssB_osWT->Fill(svtMass->at(ssbSVT), evtWeight);
        (tWriter->osMuonTag) = 0 ;
    }

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
    if(writeVars){
        float kappa = 1;
        float drCone = 0.4;
        bool osMuonJet = false;

        //JET variables
        int iJet = trkJet->at(itkmu);
        if(iJet<0 && trkPFC->at(itkmu)>=0) iJet=pfcJet->at(trkPFC->at(itkmu));  
        TVector3 vMu(muoPx->at(iMuon), muoPy->at(iMuon), muoPz->at(iMuon));

        float muoJetPtRel = -1;
        float muoJetDr = -1;
        float muoJetEnergyRatio = -1;
        int   muoJetSize = -1;
        float muoJetQ = -1;
        float muoJetPt = -1;

        float muoJetCSV = -1;
        float muoJetDFprob = -1;
        int   muoJetNDau = -1;
        float muoJetNHF = -1;
        float muoJetNEF = -1;
        float muoJetCHF = -1;
        float muoJetCEF = -1;
        int   muoJetNCH = -1;

        if(iJet>=0){
            osMuonJet = true;
            vector <int> jet_pfcs = pfCandFromJet( iJet );
            vector <int> jet_tks = tracksFromJet( iJet );

            vector <int> pfToRemove;
            vector <int> pfToKeep;
            for(auto it:jet_pfcs){
                int tk = pfcTrk->at(it);
                if(tk<0) continue;
                if(fabs(dZ(tk, ssbPVT))>=1.0) pfToRemove.push_back(it);
            }

            TLorentzVector tJet(0., 0., 0., 0.);
            for(auto it:jet_pfcs){
                if ( std::find(pfToRemove.begin(), pfToRemove.end(), it) != pfToRemove.end() ) continue;
                pfToKeep.push_back(it);
                TLorentzVector tTPf;
                tTPf.SetPxPyPzE(pfcPx->at(it),pfcPy->at(it),pfcPz->at(it),pfcE->at(it));
                tJet += tTPf;
            }

            TVector3 vJet(tJet.Px(), tJet.Py(), tJet.Pz());
            muoJetPt = tJet.Pt();
            muoJetDr = deltaR(tJet.Eta(), tJet.Phi(), muoEta->at( iMuon ), muoPhi->at(iMuon));
            if(tJet.E() != 0) muoJetEnergyRatio = muoE->at(iMuon) / tJet.E();
            else muoJetEnergyRatio = 1;
            vJet -= vMu;
            muoJetPtRel = muoPt->at( iMuon ) * (vMu.Unit() * vJet.Unit());
            muoJetSize = pfToKeep.size();
            muoJetQ = GetListPfcCharge(&pfToKeep, kappa) * trkCharge->at(itkmu);

            muoJetCSV = jetCSV->at(iJet);
            muoJetDFprob = GetJetProbb(iJet);
            muoJetNDau = jetNDau->at(iJet);
            muoJetNHF = jetNHF->at(iJet);
            muoJetNEF = jetNEF->at(iJet);
            muoJetCHF = jetCHF->at(iJet);
            muoJetCEF = jetCEF->at(iJet);
            muoJetNCH = jetNCH->at(iJet);
        }

        //CONE variables
        float muoConePtRel = -1;
        float muoConeDr = -1;
        float muoConeEnergyRatio = -1;
        int   muoConeSize = 0;
        float muoConeQ = -1;
        float muoConePt = -1;
        float muoConeNF = 0;
        float muoConeCF = 0;
        int   muoConeNCH = 0;

        TLorentzVector tCone(0.,0.,0.,0.), tMu;
        tMu.SetPtEtaPhiM(muoPt->at(iMuon),muoEta->at(iMuon),muoPhi->at(iMuon),MassMu);
        float qCone=0, ptCone=0;

        for(int ipf=0; ipf<nPF; ++ipf){
            float ptpfc = pfcPt->at(ipf);
            float etapfc = pfcEta->at(ipf);
            if( deltaR(etapfc, pfcPhi->at(ipf), muoEta->at(iMuon), muoPhi->at(iMuon)) > drCone) continue;
            if(ptpfc < 0.2) continue;
            if(fabs(etapfc) > 3.0) continue;  
            if(std::find(tkSsB.begin(), tkSsB.end(), pfcTrk->at(ipf)) != tkSsB.end()) continue;
            if(pfcTrk->at(ipf)>=0)
                if(fabs(dZ(pfcTrk->at(ipf), ssbPVT))>=1.0) 
                    continue;
      
            TLorentzVector a;
            a.SetPxPyPzE(pfcPx->at(ipf), pfcPy->at(ipf), pfcPz->at(ipf), pfcE->at(ipf));
            tCone += a;
            ++muoConeSize;

            qCone += pfcCharge->at(ipf) * pow(ptpfc, kappa);
            ptCone += pow(ptpfc, kappa);
            if(pfcCharge->at(ipf)==0) muoConeNF += pfcE->at(ipf);
            if(abs(pfcCharge->at(ipf))==1){
                muoConeNCH++;
                muoConeCF += pfcE->at(ipf);
            }
        }

        if(ptCone != 0) qCone /= ptCone;
        else qCone = 1;
        qCone *= trkCharge->at(itkmu);
        if(tCone.E()!=0){
            muoConeCF /= tCone.E();
            muoConeNF /= tCone.E();
        }

        muoConeQ = qCone;

        muoConePt = tCone.Pt();
        muoConeDr = deltaR(tCone.Eta(), tCone.Phi(), muoEta->at( iMuon ), muoPhi->at(iMuon));
        if(tCone.E() !=0 ) muoConeEnergyRatio = muoE->at(iMuon) / tCone.E();
        else muoConeEnergyRatio = 1;
        tCone -= tMu;
        muoConePtRel = muoPt->at( iMuon ) * (tMu.Vect().Unit() * tCone.Vect().Unit());


        //CONEClean variables
        float muoConeCleanPtRel = -1;
        float muoConeCleanDr = -1;
        float muoConeCleanEnergyRatio = -1;
        int   muoConeCleanSize = 0;
        float muoConeCleanQ = -1;
        float muoConeCleanPt = -1;
        float muoConeCleanNF = 0;
        float muoConeCleanCF = 0;
        int   muoConeCleanNCH = 0;

        TLorentzVector tConeClean(0.,0.,0.,0.);
        float qConeClean=0, ptConeClean=0;

        for(int ipf=0; ipf<nPF; ++ipf){
            float ptpfc = pfcPt->at(ipf);
            float etapfc = pfcEta->at(ipf);
            if( deltaR(etapfc, pfcPhi->at(ipf), muoEta->at(iMuon), muoPhi->at(iMuon)) > drCone) continue;
            if(ptpfc < 0.5) continue;
            if(fabs(etapfc) > 3.0) continue;
            if(pfcCharge->at(ipf) == 0) continue;
            if(std::find(tkSsB.begin(), tkSsB.end(), pfcTrk->at(ipf)) != tkSsB.end()) continue;
            if(pfcTrk->at(ipf)<0) continue;
            if(fabs(dZ(pfcTrk->at(ipf), ssbPVT))>=1.0) continue;
      
            TLorentzVector a;
            a.SetPxPyPzE(pfcPx->at(ipf), pfcPy->at(ipf), pfcPz->at(ipf), pfcE->at(ipf));
            tConeClean += a;
            ++muoConeCleanSize;

            qConeClean += pfcCharge->at(ipf) * pow(ptpfc, kappa);
            ptConeClean += pow(ptpfc, kappa);
            if(pfcCharge->at(ipf)==0) muoConeCleanNF += pfcE->at(ipf);
            if(abs(pfcCharge->at(ipf))==1){
                muoConeCleanNCH++;
                muoConeCleanCF += pfcE->at(ipf);
            }
        }

        if(ptConeClean != 0) qConeClean /= ptConeClean;
        else qConeClean = 1;
        qConeClean *= trkCharge->at(itkmu);
        if(tConeClean.E()!=0){
            muoConeCleanCF /= tConeClean.E();
            muoConeCleanNF /= tConeClean.E();
        }

        muoConeCleanQ = qConeClean;

        muoConeCleanPt = tConeClean.Pt();
        muoConeCleanDr = deltaR(tConeClean.Eta(), tConeClean.Phi(), muoEta->at( iMuon ), muoPhi->at(iMuon));
        if(tConeClean.E()!=0) muoConeCleanEnergyRatio = muoE->at(iMuon) / tConeClean.E();
        else muoConeCleanEnergyRatio = 1;
        tConeClean -= tMu;
        muoConeCleanPtRel = muoPt->at( iMuon ) * (tMu.Vect().Unit() * tConeClean.Vect().Unit());
        tConeClean += tMu; // for IP sign

        //------------------------------------------------FILLING------------------------------------------------
        (tWriter->muoPt) = muoPt->at( iMuon );
        (tWriter->muoEta) = muoEta->at( iMuon );
        (tWriter->muoPhi) = muoPhi->at(iMuon);
        (tWriter->muoCharge) = trkCharge->at(itkmu);

        (tWriter->muoDxy) = dSign(itkmu, tConeClean.Px(), tConeClean.Py())*abs(trkDxy->at(itkmu));
        (tWriter->muoDz) = dZ(itkmu, ssbPVT);
        (tWriter->muoExy) = trkExy->at(itkmu);
        (tWriter->muoEz) = trkEz->at(itkmu);

        (tWriter->muoSoftMvaValue) = computeMuonMva(iMuon);

        (tWriter->muoLund) = muoLund;
        (tWriter->muoAncestor) = muoAncestor; 

        (tWriter->muoDrB) = deltaR(tB.Eta(), tB.Phi(), muoEta->at(iMuon), muoPhi->at(iMuon));
        (tWriter->muoPFIso) = GetMuoPFiso(iMuon);

        (tWriter->osMuonJet) = osMuonJet;
        (tWriter->muoJetPt) = muoJetPt;
        (tWriter->muoJetPtRel) = muoJetPtRel;
        (tWriter->muoJetDr) = muoJetDr;
        (tWriter->muoJetEnergyRatio) = muoJetEnergyRatio;
        (tWriter->muoJetQ) = muoJetQ;
        (tWriter->muoJetSize) = muoJetSize;

        (tWriter->muoJetCSV) = muoJetCSV;
        (tWriter->muoJetDFprob) = muoJetDFprob;
        (tWriter->muoJetNDau) = muoJetNDau;
        (tWriter->muoJetNF) = muoJetNHF + muoJetNEF;
        (tWriter->muoJetCF) = muoJetCHF + muoJetCEF;
        (tWriter->muoJetNCH) = muoJetNCH;

        (tWriter->muoConePt) = muoConePt;
        (tWriter->muoConePtRel) = muoConePtRel;
        (tWriter->muoConeDr) = muoConeDr;
        (tWriter->muoConeEnergyRatio) = muoConeEnergyRatio;
        (tWriter->muoConeQ) = muoConeQ;
        (tWriter->muoConeSize) = muoConeSize;
        (tWriter->muoConeNF) = muoConeNF;
        (tWriter->muoConeCF) = muoConeCF;
        (tWriter->muoConeNCH) = muoConeNCH;

        (tWriter->muoConeCleanPt) = muoConeCleanPt;
        (tWriter->muoConeCleanPtRel) = muoConeCleanPtRel;
        (tWriter->muoConeCleanDr) = muoConeCleanDr;
        (tWriter->muoConeCleanEnergyRatio) = muoConeCleanEnergyRatio;
        (tWriter->muoConeCleanQ) = muoConeCleanQ;
        (tWriter->muoConeCleanSize) = muoConeCleanSize;
        (tWriter->muoConeCleanNF) = muoConeCleanNF;
        (tWriter->muoConeCleanCF) = muoConeCleanCF;
        (tWriter->muoConeCleanNCH) = muoConeCleanNCH;

        (tWriter->muoHowMany) = getNosMuons();

    }

    //------------------------------------------------TAG------------------------------------------------

    //CHARGE CORRELATION 
    if( muoAncestor >=0 ){
        if( TMath::Sign(1, ssBLund) == -1*trkCharge->at(itkmu) ){
            hmass_ssB_osCC->Fill(svtMass->at(ssbSVT), evtWeight);
            (tWriter->osMuonChargeInfo) = 1 ;
        }else{
            hmass_ssB_osWC->Fill(svtMass->at(ssbSVT), evtWeight);
            (tWriter->osMuonChargeInfo) = 0 ;
        }
    }else{
        hmass_ssB_osRC->Fill(svtMass->at(ssbSVT), evtWeight);
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

    if(use_gen){
        float eff = hmass_ssB_os->Integral() / hmass_ssB->Integral();
        float w  = hmass_ssB_osWT->Integral() / hmass_ssB_os->Integral();
        float power = eff*pow(1-2*w, 2);
        float tot = hmass_ssB_osCC->Integral() + hmass_ssB_osWC->Integral() + hmass_ssB_osRC->Integral();

        cout<<"CC   WC  RC"<<endl;
        cout<< hmass_ssB_osCC->Integral()/tot<<" "<< hmass_ssB_osWC->Integral()/tot<<" "<< hmass_ssB_osRC->Integral()/tot<<endl<<endl;

        cout<<"#B    eff%    w%    P%"<<endl;
        cout<< hmass_ssB->Integral()<<" "<<eff*100<<" "<<w*100<<" "<<power*100<<endl;
    }else{
        float eff = CountEventsWithFit(hmass_ssB_os, process) / CountEventsWithFit(hmass_ssB, process);
        float w  = CountEventsWithFit(hmass_ssB_osWT, process) / CountEventsWithFit(hmass_ssB_os, process);
        float power = eff*pow(1-2*w, 2);
        float tot = CountEventsWithFit(hmass_ssB_osCC, process ) + CountEventsWithFit(hmass_ssB_osWC, process) + CountEventsWithFit(hmass_ssB_osRC, process);

        cout<<"CC   WC  RC"<<endl;
        cout<< CountEventsWithFit(hmass_ssB_osCC, process)/tot<<" "<< CountEventsWithFit(hmass_ssB_osWC, process)/tot<<" "<< CountEventsWithFit(hmass_ssB_osRC, process)/tot<<endl<<endl;

        cout<<"#B    eff%    w%    P%"<<endl;
        cout<< CountEventsWithFit(hmass_ssB, process)<<" "<<eff*100<<" "<<w*100<<" "<<power*100<<endl;
    }



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
float PDAnalyzer::GetListPfcCharge(vector <int> *list, float kappa)
{
    float Q = 0;
    float pt = 0;
    for(int it:*list){
       Q += pfcCharge->at(it) * pow(pfcPt->at(it), kappa);
       pt += pow(pfcPt->at(it), kappa);
    }
    return Q/pt; 
}