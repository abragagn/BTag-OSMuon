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
,               nMuonsSel_(0)
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

void OSMuonMvaTag::setOsMuonCuts(float wpB, float wpE, float dzCut)
{
    wpB_ = wpB;
    wpE_ = wpE;
    dzCut_ = dzCut;
}

void OSMuonMvaTag::inizializeOSMuonMvaTagReader(
    TString methodName, 
    TString path = "/lustre/cmswork/abragagn/mvaWeights/OsMuonTag/" )
{

    TMVA::PyMethodBase::PyInitialize();
    setWeights(methodName, path);

    osMuonTagReader_.AddVariable( "muoPt", &muoPt_);
    osMuonTagReader_.AddVariable( "abs_muoEta := fabs(muoEta)", &absmuoEta_);
    osMuonTagReader_.AddVariable( "muoDxy", &muoDxy_);
    osMuonTagReader_.AddVariable( "abs_muoDz := fabs(muoDz)", &absmuoDz_);
    osMuonTagReader_.AddVariable( "muoSoftMvaValue", &muoSoftMvaValue_);
    osMuonTagReader_.AddVariable( "muoDrB", &muoDrB_);
    osMuonTagReader_.AddVariable( "muoPFIso", &muoPFIso_);

    osMuonTagReader_.AddVariable( "muoJetConePt := muoJetPt != -1 ? muoJetPt : muoConePt", &muoJetConePt_);
    osMuonTagReader_.AddVariable( "muoJetConePtRel := muoJetPt != -1 ? muoJetPtRel : muoConePtRel", &muoJetConePtRel_);
    osMuonTagReader_.AddVariable( "muoJetConeDr", &muoJetConeDr_);
    osMuonTagReader_.AddVariable( "muoJetConeEnergyRatio := muoJetPt != -1 ? muoJetEnergyRatio : muoConeEnergyRatio", &muoJetConeEnergyRatio_);

    osMuonTagReader_.AddVariable( "muoJetCSV", &muoJetCSV_);
    if(!weightsFile_.Contains("2016")) osMuonTagReader_.AddVariable( "muoJetDFprob", &muoJetDFprob_);

    osMuonTagReader_.AddVariable( "muoJetConeSize := muoJetPt != -1 ? muoJetSize : muoConeSize", &muoJetConeSize_);
    osMuonTagReader_.AddVariable( "muoJetConeQ := muoJetPt != -1 ? muoJetQ : muoConeQ", &muoJetConeQ_);

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

int OSMuonMvaTag::getOsMuon()
{
    if(ssIndex_ < 0){ cout<<"SS NOT INITIALIZED"<<endl; return -2; }

    int iB = ssIndex_;
    int iPV = pvIndex_;

    vector <int> tkSsB = tracksFromSV(iB);
    TLorentzVector tB = GetTLorentzVecFromJpsiX(iB);
    if(pvIndex_ < 0) pvIndex_ = GetBestPV(iB, tB);

    int bestMuIndex = -1;
    float bestMuPt = 2.;
    int bestMuTrack = -1;
    nMuonsSel_ = 0;

    for(int iMuon = 0; iMuon < nMuons; ++iMuon ){

        int itkmu = muonTrack( iMuon, PDEnumString::muInner );
        if(itkmu<0) continue;

        if(std::find(tkSsB.begin(), tkSsB.end(), itkmu) != tkSsB.end()) continue;

        if(muoPt->at( iMuon ) < 2.) continue;
        if(fabs(muoEta->at( iMuon )) > 2.4) continue;
       
        if(!isMvaMuon(iMuon, wpB_, wpE_)) continue;

        if(fabs(dZ(itkmu, iPV)) > dzCut_) continue;
        if(deltaR(tB.Eta(), tB.Phi(), muoEta->at(iMuon), muoPhi->at(iMuon)) < 0.4) continue;
        //if(GetMuoPFiso(iMuon) > PFIsoCut_)  continue;

        nMuonsSel_++;

        if(muoPt->at( iMuon ) > bestMuPt){
            bestMuPt = muoPt->at( iMuon );
            bestMuIndex = iMuon;
            bestMuTrack = itkmu;
        }
    }

    osMuonIndex_ = bestMuIndex;
    osMuonTrackIndex_ = bestMuTrack;
    return bestMuIndex;

}

void OSMuonMvaTag::computeVariables()
{

    int iB = ssIndex_;
    int iMuon = osMuonIndex_;
    int iPV = pvIndex_;

    int itkmu = muonTrack( iMuon, PDEnumString::muInner );
    TLorentzVector tB = GetTLorentzVecFromJpsiX(iB);
    vector <int> tkSsB = tracksFromSV(iB);
    float kappa = 1;

    //JET VARIABLES

    int iJet = trkJet->at(itkmu);
    if(iJet<0 && trkPFC->at(itkmu)>=0) iJet=pfcJet->at(trkPFC->at(itkmu));  
    TVector3 vMu(muoPx->at(iMuon), muoPy->at(iMuon), muoPz->at(iMuon));

    float muoJetPt = -1;
    float muoJetPtRel = -1;
    float muoJetDr = -1;
    float muoJetEnergyRatio = -1;
    float muoJetCSV = -1;
    float muoJetDFprob = -1;
    float muoJetSize = -1;
    float muoJetQ = -1; 

    if(iJet>=0){
        vector <int> jet_pfcs = pfCandFromJet( iJet );
        TVector3 vJet(jetPx->at(iJet), jetPy->at(iJet), jetPz->at(iJet));
        muoJetPt = jetPt->at(iJet);
        muoJetDr = deltaR(jetEta->at(iJet), jetPhi->at(iJet), muoEta->at( iMuon ), muoPhi->at(iMuon));
        muoJetEnergyRatio = muoE->at(iMuon) / jetE->at(iJet);
        vJet -= vMu;
        muoJetPtRel = muoPt->at( iMuon ) * (vMu.Unit() * vJet.Unit());
        muoJetSize = jet_pfcs.size();
        muoJetCSV = jetCSV->at(iJet);
        muoJetDFprob = GetJetProbb(iJet);
        muoJetQ = GetJetCharge(iJet, kappa);
        muoJetQ *= trkCharge->at(itkmu); 
    }

    //CONE VARIABLES

    float muoConePt = -1;
    float muoConePtRel = -1;
    float muoConeDr = -1;
    float muoConeEnergyRatio = -1;
    float muoConeSize = 0;
    float muoConeQ = 0;

    TLorentzVector tCone, tMu;
    tCone.SetPtEtaPhiM(0.,0.,0.,0.);
    tMu.SetPtEtaPhiM(muoPt->at( iMuon ), muoEta->at( iMuon ), muoPhi->at( iMuon ), MassMu);

    float drCharge = 0.4;
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

    muoConePt = tCone.Pt();
    muoConeDr = deltaR(tCone.Eta(), tCone.Phi(), muoEta->at( iMuon ), muoPhi->at(iMuon));
    muoConeEnergyRatio = muoE->at(iMuon) / tCone.E();
    tCone -= tMu;
    muoConePtRel = muoPt->at( iMuon ) * (tMu.Vect().Unit() * tCone.Vect().Unit());
    if(ptCone != 0) qCone /= ptCone;
    else qCone = 1;
    qCone *= trkCharge->at(itkmu); 
    muoConeQ = qCone;

    DUMMY_ = -1;
    muoPt_ = muoPt->at(iMuon);
    absmuoEta_ = fabs(muoEta->at(iMuon));
    muoDxy_= GetSignedDxy(iMuon, iPV);
    absmuoDz_= fabs(dZ(itkmu, iPV));
    muoSoftMvaValue_= computeMuonMva(iMuon);
    muoDrB_= deltaR(tB.Eta(), tB.Phi(), muoEta->at(iMuon), muoPhi->at(iMuon));
    muoPFIso_= GetMuoPFiso(iMuon);
    muoJetConePt_ = muoJetPt != -1 ? muoJetPt : muoConePt;
    muoJetConePtRel_ = muoJetPt != -1 ? muoJetPtRel : muoConePtRel;
    muoJetConeDr_ = muoJetPt != -1 ? muoJetDr : muoConeDr;
    muoJetConeEnergyRatio_ = muoJetPt != -1 ? muoJetEnergyRatio : muoConeEnergyRatio;
    muoJetCSV_ = muoJetCSV;
    muoJetDFprob_ = muoJetDFprob;
    muoJetConeSize_ = muoJetPt != -1 ? muoJetSize : muoConeSize;
    muoJetConeQ_ = muoJetPt != -1 ? muoJetQ : muoConeQ;

}

int OSMuonMvaTag::getOsMuonTag()
{
    if(ssIndex_ < 0){ cout<<"SS NOT INITIALIZED"<<endl; return -999; }
    getOsMuon();
    if(osMuonIndex_ < 0) return 0;

    return -1*trkCharge->at(osMuonTrackIndex_); 
}

float OSMuonMvaTag::getOsMuonTagMvaValue()
{
    if(ssIndex_ < 0){ cout<<"SS NOT INITIALIZED"<<endl; return -999; }
    if(osMuonIndex_ < 0){ cout<<"WARNING: OS MU NOT INITIALIZED"<<endl; getOsMuon(); }
    
    computeVariables();
    return osMuonTagReader_.EvaluateMVA(methodName_);
}

/*float OSMuonMvaTag::getOsMuonTagMistagProb()
{
    //todo
}
*/
