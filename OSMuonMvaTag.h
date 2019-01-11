#ifndef OSMuonMvaTag_H
#define OSMuonMvaTag_H

#include <vector>
#include <set>
#include <map>
#include <string>

#include "PDSoftMuonMvaEstimator.h"
#include "PDAnalyzerUtil.h"
#include "AlbertoUtil.h"

#include "TString.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/PyMethodBase.h"

class OSMuonMvaTag : public virtual PDAnalyzerUtil
,                    public virtual AlbertoUtil
{

public:
    OSMuonMvaTag();
    ~OSMuonMvaTag();

    int     getOsMuon();
    int     getOsMuonTag();
    float   getOsMuonTagMvaValue();
    float   getOsMuonTagMistagProb();

    void    setSsForTag(int iB, int iPV) { ssIndex_ = iB; pvIndex_ = iPV;}
    void    setOsMuonCuts(float wpB, float wpE, float dzCut);
    void    inizializeOSMuonMvaTagReader(TString weightsFile, TString path);

    int     getNosMuons(){return nMuonsSel_;}
private:    
    TString methodNameFromWeightName();
    void    computeVariables();
    void    setWeights(TString methodName, TString path);

    TMVA::Reader osMuonTagReader_;
    TString weightsFile_;
    TString methodName_;

    int ssIndex_;
    int pvIndex_;
    int osMuonIndex_;
    int osMuonTrackIndex_;

    float wpB_;
    float wpE_;
    float dzCut_;
    float PFIsoCut_;

    int nMuonsSel_;

    //MVA Variables

    float muoPt_;
    float absmuoEta_;
    float muoDxy_;
    float absmuoDz_;
    float muoSoftMvaValue_;
    float muoDrB_;
    float muoPFIso_;

    float muoJetConePt_;
    float muoJetConePtRel_;
    float muoJetConeDr_;
    float muoJetConeEnergyRatio_;
    float muoJetCSV_;
    float muoJetDFprob_;
    float muoJetConeSize_;
    float muoJetConeQ_;

    float muoCharge_;

    float DUMMY_;

};

#endif
