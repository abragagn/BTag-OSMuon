#ifndef OSMuonMvaTag_H
#define OSMuonMvaTag_H

#include <vector>
#include <set>
#include <map>
#include <string>

#include "PDAnalyzerUtil.h"

#include "TString.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/PyMethodBase.h"

class OSMuonMvaTag :    public virtual PDAnalyzerUtil
{

public:
    OSMuonMvaTag();
    ~OSMuonMvaTag();

    float getEventMistagProb();
    int getEventTag();

private:
    TMVA::Reader reader_;
    TString weightsFile_;

    void OSMuonMvaTag::setupReader(TString);
    TString OSMuonMvaTag::methodNameFromWeightName();

    //MVA Variables

    float muoPt_;
    float muoEta_;
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
    int   muoJetSize_;

    float DUMMY_;

};

#endif
