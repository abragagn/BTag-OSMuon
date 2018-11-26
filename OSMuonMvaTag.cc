#include "PDSoftMuonMvaEstimator.h"
#include "AlbertoUtil.h"

OSMuonMvaTag::OSMuonMvaTag(TString weightsFile = "/lustre/cmswork/abragagn/weights/TMVAClassification_BDTOsMuon2016.weights.xml"):
    reader_("!Color:Silent"),
    ssIndex_(-1),
    osMuonIndex_(-1),
    osMuonTrackIndex_(-1)
{
    TMVA::PyMethodBase::PyInitialize();
    weightsFile_ = weightsFile,
    methodName_ = methodNameFromWeightName();
    setupReader();
}

OSMuonMvaTag::~OSMuonMvaTag() {}

// =====================================================================================
void OSMuonMvaTag::setupReader()
{

    reader_.AddVariable( "muoPt_", &muoPt_);
    reader_.AddVariable( "muoEta_", &absmuoEta_);
    reader_.AddVariable( "muoDxy_", &muoDxy_);
    reader_.AddVariable( "muoDz_", &muoDz_);
    reader_.AddVariable( "muoSoftMvaValue_", &muoSoftMvaValue_);
    reader_.AddVariable( "muoDrB_", &muoDrB_);
    reader_.AddVariable( "muoPFIso_", &muoPFIso_);

    if(weightsFile_.Contains("Jet")){
        reader_.AddVariable( "muoJetPt_", &muoJetPt_);
        reader_.AddVariable( "muoJetPtRel_", &muoJetPtRel_);
        reader_.AddVariable( "muoJetDr_", &muoJetDr_);
        reader_.AddVariable( "muoJetEnergyRatio_", &muoJetEnergyRatio_);
        reader_.AddVariable( "muoJetCSV_", &muoJetCSV_);
        if(!weightsFile_.Contains("2016")) reader_.AddVariable( "muoJetDFprob_", &muoJetDFprob_);
    }

    if(weightsFile_.Contains("Cone")){
        reader_.AddVariable( "muoConePt_", &muoConePt_);
        reader_.AddVariable( "muoConePtRel_", &muoConePtRel_);
        reader_.AddVariable( "muoConeDr_", &muoConeDr_);
        reader_.AddVariable( "muoConeEnergyRatio_", &mumuoConeEnergyRatio_oConeDr_);
        reader_.AddVariable( "muoConeSize_", &muoConeSize_);
    }

    reader_.AddVariable( "muoQCone_", &muoQCone_);

    reader_.BookMVA( methodName_, weightsFile_ );

    return;

}

TString OSMuonMvaTag::methodNameFromWeightName()
{
    TString prefix = "TMVAClassification_";
    int start = weightsName_.Index(prefix) + prefix.Length();
    int length = weightsName_.Index(".weights") - start;
    TString name( weightsName_(start, length) );
    return name;
}

int OSMuonMvaTag::getOsMuon(int iB = -999)
{
    if(iB == -999) iB = ssIndex_;

    vector <int> tkSsB = tracksFromSV(iB);

    int bestMuIndex = -1;
    float bestMuPt = minPtMuon;

    for(int iMuon = 0; iMuon < nMuons; ++iMuon ){

        int itkmu = muonTrack( iMuon, PDEnumString::muInner );
        if(itkmu<0) continue;

        if(std::find(tkSsB.begin(), tkSsB.end(), itkmu) != tkSsB.end()) continue;

        if(muoPt->at( iMuon ) < 2.) continue;
        if(abs(muoEta->at( iMuon )) > 2.4) continue;
       
        if(!isMvaMuon(iMuon, muonIdWpBarrel, muonIdWpEndcap)) continue;

        if(abs(dZ(itkmu, iSsPV)) > muoDzCut) continue;
        if(GetMuoPFiso(iMuon) > muoPFIsoCut)  continue;


        if(muoPt->at( iMuon ) > bestMuPt){
            bestMuPt = muoPt->at( iMuon );
            bestMuIndex = iMuon;
        }
    }

    osMuonIndex_ = bestMuIndex;
    osMuonTrackIndex_ = itkmu;
    return bestMuIndex;

}
void OSMuonMvaTag::computeVariables()
{
    int iB = ssIndex_;
    int iMuon = osMuonIndex_;

    int itkmu = muonTrack( iMuon, PDEnumString::muInner );
    TLorentzVector tB = GetTLorentzVecFromJpsiX(iSsB);

    DUMMY_ = -1;

    muoPt_ muoPt->at(iMuon);
    absMuoEta_ abs(muoEta->at(iMuon));
    muoDxy_= GetSignedDxy(iMuon, iSsPV);
    muoDz_= dZ(itkmu, iSsPV);
    muoSoftMvaValue_= computeMva(iMuon);
    muoDrB_= deltaR(tB.Eta(), tB.Phi(), muoEta->at(iMuon), muoPhi->at(iMuon));
    muoPFIso_= GetMuoPFiso(iMuon);

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
    qCone *= trkCharge->at(itkmu); 
    muoQCone_ = qCone;

    //JET VARIABLES

    int iJet = trkJet->at(itkmu);
    if(iJet<0 && trkPFC->at(itkmu)>=0) iJet=pfcJet->at(trkPFC->at(itkmu));  
    TVector3 vMu(muoPx->at(iMuon), muoPy->at(iMuon), muoPz->at(iMuon));

    muoJetPt_ = -1;
    muoJetPtRel_ = -1;
    muoJetDr_ = -1;
    muoJetEnergyRatio_ = -1;
    muoJetCSV_ = -1;
    muoJetDFprob_ = -1;
    muoJetSize_ = -1;


    if(iJet>=0){
        vector <int> jet_pfcs = pfCandFromJet( iJet );
        TVector3 vJet(jetPx->at(iJet), jetPy->at(iJet), jetPz->at(iJet));
        muoJetPt_ = jetPt->at(iJet);
        muoJetDr_ = deltaR(jetEta->at(iJet), jetPhi->at(iJet), muoEta->at( iMuon ), muoPhi->at(iMuon));
        muoJetEnergyRatio_ = muoE->at(iMuon) / jetE->at(iJet);
        vJet -= vMu;
        muoJetPtRel_ = muoPt->at( iMuon ) * (vMu.Unit() * vJet.Unit());
        muoJetSize_ = jet_pfcs.size();
        muoJetCSV_ = jetCSV->at(iJet);
        muoJetDFprob_ = GetJetProbb(iJet);
    }



    muoConePt_ = -1;
    muoConePtRel_ = -1;
    muoConeDr_ = -1;
    muoConeEnergyRatio_ = -1;
    muoConeSize_ = -1;

    TLorentzVector tCone, tMu;
    tCone.SetPtEtaPhiM(0.,0.,0.,0.);
    tMu.SetPtEtaPhiM(muoPt->at( iMuon ), muoEta->at( iMuon ), muoPhi->at( iMuon ), MassMu);

    for(int i=0; i<nPF; ++i){
        if( deltaR(pfcEta->at( i ), pfcPhi->at( i ), muoEta->at( iMuon ), muoPhi->at( iMuon )) > 0.4) continue;
        TLorentzVector a;
        a.SetPxPyPzE(pfcPx->at(i), pfcPy->at(i), pfcPz->at(i), pfcE->at(i));
        tCone += a;
        ++muoConeSize;
    }

    muoConePt_ = tCone.Pt();
    muoConeDr_ = deltaR(tCone.Eta(), tCone.Phi(), muoEta->at( iMuon ), muoPhi->at(iMuon));
    muoConeEnergyRatio_ = muoE->at(iMuon) / tCone.E();
    tCone -= tMu;
    muoConePtRel_ = muoPt->at( iMuon ) * (tMu.Vect().Unit() * tCone.Vect().Unit());

    return;
}

int OSMuonMvaTag::getOsMuonTag(int iB)
{
    ssIndex_ = iB;
    getOsMuon();
    if(osMuonIndex_ < 0) return 0;

    return -1*trkCharge->at(osMuonTrackIndex_);
    
}


float OSMuonMvaTag::getOsMuonMvaValue()
{
    computeVariables();
    return reader_.EvaluateMVA(methodName_);
}

float OSMuonMvaTag::getOsMuonMistagProb()
{
    //todo
}