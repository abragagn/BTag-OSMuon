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

class OSMuonMvaTag :    public virtual PDAnalyzerUtil
,                       public virtual AlbertoUtil
{

public:
    OSMuonMvaTag();
    ~OSMuonMvaTag();

    int     getOsMuon(int iB);
    int     getOsMuonTag(int iB);
    float   getOsMuonMvaValue();
    float   getOsMuonMistagProb();

    void    setSsForTag(int iB) { ssIndex_ = iB; }
    void    setOsMuonCuts(float wpB, float wpE, float dzCut, float PFIsoCut);
    void    inizializeOSMuonMvaReader(TString weightsFile, TString path);

private:    
    TString methodNameFromWeightName();
    void    computeVariables();
    void    setWeights(TString methodName, TString path);

    TMVA::Reader osMuonTagReader_;
    TString weightsFile_;
    TString methodName_;

    int ssIndex_;
    int osMuonIndex_;
    int osMuonTrackIndex_;

    float wpB_;
    float wpE_;
    float dzCut_;
    float PFIsoCut_;

    //MVA Variables

    float muoPt_;
    float absmuoEta_;
    float muoDxy_;
    float muoDz_;
    float muoSoftMvaValue_;
    float muoDrB_;
    float muoPFIso_;

    float muoJetPt_;
    float muoJetPtRel_;
    float muoJetDr_;
    float muoJetEnergyRatio_;
    float muoJetCSV_;
    float muoJetDFprob_;
    int   muoJetSize_;

    float muoQCone_;

    float muoConePt_;
    float muoConePtRel_;
    float muoConeDr_;
    float muoConeEnergyRatio_;
    int   muoConeSize_;

    float muoCharge_;

    float DUMMY_;

};

#endif
