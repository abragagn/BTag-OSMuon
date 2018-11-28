#include "PDSoftMuonMvaEstimator.h"
#include "AlbertoUtil.h"

OSMuonMvaTag::OSMuonMvaTag():
                osMuonTagReader_("!Color:Silent")
,               ssIndex_(-1)
,               osMuonIndex_(-1)
,               osMuonTrackIndex_(-1)
,               wpB_(0.)
,               wpE_(0.)
,               dzCut_(1.)
,               PFIsoCut_(5.)
{}

OSMuonMvaTag::~OSMuonMvaTag() {}

// =====================================================================================
void OSMuonMvaTag::setWeights(TString methodName, TString path)
{    
    TString year = "";
    if(methodName.Contains("2016")) year = "2016";
    if(methodName.Contains("2017")) year = "2017";
    if(methodName.Contains("2018")) year = "2018";

    TString weightsFile = path + year + "/TMVAClassification_" + methodName  + ".weights.xml";

    weightsFile_ = weightsFile;
    methodName_ = methodName;
}

void OSMuonMvaTag::setOsMuonCuts(float wpB, float wpE, float dzCut, float PFIsoCut)
{
    wpB_ = wpB;
    wpE_ = wpE;
    dzCut_ = dzCut;
    PFIsoCut_ = PFIsoCut;
}

void OSMuonMvaTag::inizializeOSMuonMvaReader(TString methodName, TString path = "/lustre/cmswork/abragagn/weights/" )
{

    TMVA::PyMethodBase::PyInitialize();
    setWeights(methodName, path);

    osMuonTagReader_.AddVariable( "muoPt", &muoPt_);
    osMuonTagReader_.AddVariable( "muoEta", &absmuoEta_);
    osMuonTagReader_.AddVariable( "muoDxy", &muoDxy_);
    osMuonTagReader_.AddVariable( "muoDz", &muoDz_);
    osMuonTagReader_.AddVariable( "muoSoftMvaValue", &muoSoftMvaValue_);
    osMuonTagReader_.AddVariable( "muoDrB", &muoDrB_);
    osMuonTagReader_.AddVariable( "muoPFIso", &muoPFIso_);

    if(weightsFile_.Contains("Jet")){
        osMuonTagReader_.AddVariable( "muoJetPt", &muoJetPt_);
        osMuonTagReader_.AddVariable( "muoJetPtRel", &muoJetPtRel_);
        osMuonTagReader_.AddVariable( "muoJetDr", &muoJetDr_);
        osMuonTagReader_.AddVariable( "muoJetEnergyRatio", &muoJetEnergyRatio_);
        osMuonTagReader_.AddVariable( "muoJetCSV", &muoJetCSV_);
        if(!weightsFile_.Contains("2016")) osMuonTagReader_.AddVariable( "muoJetDFprob", &muoJetDFprob_);
    }
    if(weightsFile_.Contains("Cone")){
        osMuonTagReader_.AddVariable( "muoConePt", &muoConePt_);
        osMuonTagReader_.AddVariable( "muoConePtRel", &muoConePtRel_);
        osMuonTagReader_.AddVariable( "muoConeDr", &muoConeDr_);
        osMuonTagReader_.AddVariable( "muoConeEnergyRatio", &muoConeEnergyRatio_);
        osMuonTagReader_.AddVariable( "muoConeSize", &muoConeSize_);
    }

    osMuonTagReader_.AddVariable( "muoQCone", &muoQCone_);

    osMuonTagReader_.BookMVA( methodName_, weightsFile_ );

}

TString OSMuonMvaTag::methodNameFromWeightName()
{
    TString prefix = "TMVAClassification_";
    int start = weightsFile_.Index(prefix) + prefix.Length();
    int length = weightsFile_.Index(".weights") - start;
    TString name( weightsFile_(start, length) );
    return name;
}

int OSMuonMvaTag::getOsMuon(int iB = -999)
{

    if(iB == -999) iB = ssIndex_;

    vector <int> tkSsB = tracksFromSV(iB);
    TLorentzVector tB = GetTLorentzVecFromJpsiX(iB);

    int bestMuIndex = -1;
    float bestMuPt = 2.;
    int itkmu = -1;

    for(int iMuon = 0; iMuon < nMuons; ++iMuon ){

        itkmu = muonTrack( iMuon, PDEnumString::muInner );
        if(itkmu<0) continue;

        if(std::find(tkSsB.begin(), tkSsB.end(), itkmu) != tkSsB.end()) continue;

        if(muoPt->at( iMuon ) < 2.) continue;
        if(abs(muoEta->at( iMuon )) > 2.4) continue;
       
        if(!isMvaMuon(iMuon, wpB_, wpE_)) continue;

        if(abs(dZ(itkmu, GetBestPV(iB, tB))) > dzCut_) continue;
        if(GetMuoPFiso(iMuon) > PFIsoCut_)  continue;

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
    TLorentzVector tB = GetTLorentzVecFromJpsiX(iB);
    vector <int> tkSsB = tracksFromSV(iB);

    DUMMY_ = -1;

    muoPt_ = muoPt->at(iMuon);
    absmuoEta_ = abs(muoEta->at(iMuon));
    muoDxy_= GetSignedDxy(iMuon, iB);
    muoDz_= dZ(itkmu, iB);
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
        ++muoConeSize_;
    }

    muoConePt_ = tCone.Pt();
    muoConeDr_ = deltaR(tCone.Eta(), tCone.Phi(), muoEta->at( iMuon ), muoPhi->at(iMuon));
    muoConeEnergyRatio_ = muoE->at(iMuon) / tCone.E();
    tCone -= tMu;
    muoConePtRel_ = muoPt->at( iMuon ) * (tMu.Vect().Unit() * tCone.Vect().Unit());

}

int OSMuonMvaTag::getOsMuonTag(int iB = -1)
{
    if(ssIndex_ < 0) ssIndex_ = iB;
    getOsMuon(ssIndex_);
    if(osMuonIndex_ < 0) return 0;

    return -1*trkCharge->at(osMuonTrackIndex_); 
}

float OSMuonMvaTag::getOsMuonMvaValue()
{
    computeVariables();
    return osMuonTagReader_.EvaluateMVA(methodName_);
}

/*float OSMuonMvaTag::getOsMuonMistagProb()
{
    //todo
}
*/