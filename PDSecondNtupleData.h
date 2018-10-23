#ifndef PDSecondNtupleData_h
#define PDSecondNtupleData_h
#include <vector>
#include "NtuTool/Common/interface/TreeWrapper.h"
using namespace std;

class PDSecondNtupleData: public virtual TreeWrapper {

public:

void Reset() { autoReset(); }

PDSecondNtupleData() {


}
virtual ~PDSecondNtupleData() {
}

void initTree() {
    treeName = "PDsecondTree";

    setBranch( "ssbPt", &ssbPt, "ssbPt/F", &b_ssbPt );
    setBranch( "ssbEta", &ssbEta, "ssbEta/F", &b_ssbEta );
    setBranch( "ssbPhi", &ssbPhi, "ssbPhi/F", &b_ssbPhi );
    setBranch( "ssbDxy", &ssbDxy, "ssbDxy/F", &b_ssbDxy );
    setBranch( "ssbExy", &ssbExy, "ssbExy/F", &b_ssbExy );
    setBranch( "ssbDz", &ssbDz, "ssbDz/F", &b_ssbDz );
    setBranch( "ssbEz", &ssbEz, "ssbEz/F", &b_ssbEz );

    setBranch( "ssbMass", &ssbMass, "ssbMass/F", &b_ssbMass );
    setBranch( "JPsiMass", &JPsiMass, "JPsiMass/F", &b_JPsiMass );
    setBranch( "ssbSVT", &ssbSVT, "ssbSVT/I", &b_ssbSVT );
    setBranch( "ssbPVT", &ssbPVT, "ssbPVT/I", &b_ssbPVT );
    setBranch( "ssbLund", &ssbLund, "ssbLund/I", &b_ssbLund );
    setBranch( "ssbIsTight", &ssbIsTight, "ssbIsTight/I", &b_ssbIsTight );

    setBranch( "ssbLxy", &ssbLxy , "ssbLxy/F" , &b_ssbLxy );
    setBranch( "ssbCt2D", &ssbCt2D , "ssbCt2D/F" , &b_ssbCt2D );
    setBranch( "ssbCt2DErr", &ssbCt2DErr , "ssbCt2DErr/F" , &b_ssbCt2DErr );
    setBranch( "ssbCt2DSigmaUnit", &ssbCt2DSigmaUnit , "ssbCt2DSigmaUnit/F" , &b_ssbCt2DSigmaUnit );
    setBranch( "ssbCt3D", &ssbCt3D , "ssbCt3D/F" , &b_ssbCt3D );
    setBranch( "ssbCt3DErr", &ssbCt3DErr , "ssbCt3DErr/F" , &b_ssbCt3DErr );
    setBranch( "ssbCt3DSigmaUnit", &ssbCt3DSigmaUnit , "ssbCt3DSigmaUnit/F" , &b_ssbCt3DSigmaUnit );

    setBranch( "muoPt", &muoPt, "muoPt/F", &b_muoPt );
    setBranch( "muoEta", &muoEta, "muoEta/F", &b_muoEta );
    setBranch( "muoPhi", &muoPhi, "muoPhi/F", &b_muoPhi );
    setBranch( "muoCharge", &muoCharge, "muoCharge/I", &b_muoCharge );
    setBranch( "muoDxy", &muoDxy, "muoDxy/F", &b_muoDxy );
    setBranch( "muoExy", &muoExy, "muoExy/F", &b_muoExy );
    setBranch( "muoDz", &muoDz, "muoDz/F", &b_muoDz );
    setBranch( "muoEz", &muoEz, "muoEz/F", &b_muoEz );
    setBranch( "muoLund", &muoLund, "muoLund/I", &b_muoLund );
    setBranch( "muoAncestor", &muoAncestor, "muoAncestor/I", &b_muoAncestor );
    setBranch( "muoSoftMvaValue", &muoSoftMvaValue, "muoSoftMvaValue/F", &b_muoSoftMvaValue );


    setBranch( "osMuon", &osMuon, "osMuon/I", &b_osMuon );
    setBranch( "osMuonTag", &osMuonTag, "osMuonTag/I", &b_osMuonTag );
    setBranch( "osMuonChargeInfo", &osMuonChargeInfo, "osMuonChargeInfo/I", &b_osMuonChargeInfo );

    setBranch( "muoDR_B", &muoDR_B, "muoDR_B/F", &b_muoDR_B );
    setBranch( "muoDz_PV", &muoDz_PV, "muoDz_PV/F", &b_muoDz_PV );
    setBranch( "muoPFIso", &muoPFIso, "muoPFIso/F", &b_muoPFIso );
    setBranch( "muoJetPt", &muoJetPt, "muoJetPt/F", &b_muoJetPt );
    setBranch( "muoPtRel", &muoPtRel, "muoPtRel/F", &b_muoPtRel );
    setBranch( "muoDrJet", &muoDrJet, "muoDrJet/F", &b_muoDrJet );
    setBranch( "muoEnergyRatio", &muoEnergyRatio, "muoEnergyRatio/F", &b_muoEnergyRatio );
    setBranch( "muoQCone", &muoQCone, "muoQCone/F", &b_muoQCone );
    setBranch( "muoQJet", &muoQJet, "muoQJet/F", &b_muoQJet );
    setBranch( "muoJetCSV", &muoJetCSV, "muoJetCSV/F", &b_muoJetCSV );
    setBranch( "muoJetDFprob", &muoJetDFprob, "muoJetDFprob/F", &b_muoJetDFprob );
    setBranch( "muoJetSize", &muoJetSize, "muoJetSize/I", &b_muoJetSize );
    setBranch( "muoHowMany", &muoHowMany, "muoHowMany/I", &b_muoHowMany );

    setBranch( "evtNumber", &evtNumber, "evtNumber/I", &b_evtNumber );
    setBranch( "evtWeight", &evtWeight, "evtWeight/F", &b_evtWeight );
    setBranch( "hltJpsiMu", &hltJpsiMu , "hltJpsiMu/I" , &b_hltJpsiMu );
    setBranch( "hltJpsiTrkTrk", &hltJpsiTrkTrk , "hltJpsiTrkTrk/I" , &b_hltJpsiTrkTrk );
    setBranch( "hltJpsiTrk", &hltJpsiTrk , "hltJpsiTrk/I" , &b_hltJpsiTrk );
}

float ssbPt, ssbEta, ssbPhi, ssbMass, JPsiMass, ssbDxy, ssbExy, ssbDz, ssbEz, evtWeight;
float ssbLxy, ssbCt2D, ssbCt2DErr, ssbCt2DSigmaUnit, ssbCt3D, ssbCt3DErr, ssbCt3DSigmaUnit;
int ssbSVT, ssbPVT, ssbLund, evtNumber, hltJpsiMu, hltJpsiTrkTrk, hltJpsiTrk, ssbIsTight;

TBranch *b_ssbPt, *b_ssbEta, *b_ssbPhi, *b_ssbMass, *b_JPsiMass, *b_ssbDxy, *b_ssbExy, *b_ssbDz, *b_ssbEz;
TBranch *b_ssbLxy, *b_ssbCt2D, *b_ssbCt2DErr, *b_ssbCt2DSigmaUnit, *b_ssbCt3D, *b_ssbCt3DErr, *b_ssbCt3DSigmaUnit;
TBranch *b_ssbSVT, *b_ssbPVT, *b_ssbLund, *b_evtNumber, *b_hltJpsiMu, *b_hltJpsiTrkTrk, *b_hltJpsiTrk, *b_ssbIsTight, *b_evtWeight;


int muoLund, muoAncestor, muoCharge;
float muoSoftMvaValue;
float muoPt, muoEta, muoPhi, muoDxy, muoExy, muoDz, muoEz;

TBranch *b_muoPt, *b_muoEta, *b_muoPhi, *b_muoDxy, *b_muoExy, *b_muoDz, *b_muoEz, *b_muoCharge;
TBranch *b_muoLund, *b_muoAncestor, *b_muoSoftMvaValue;


int osMuon, osMuonTag, osMuonChargeInfo, muoJetSize, muoHowMany;
float muoDR_B, muoPFIso, muoPtRel, muoDrJet, muoEnergyRatio, muoQCone, muoDz_PV, muoQJet, muoJetCSV, muoJetDFprob, muoJetPt;

TBranch *b_osMuon, *b_osMuonTag, *b_osMuonChargeInfo, *b_muoJetSize, *b_muoHowMany;
TBranch *b_muoDR_B, *b_muoPFIso, *b_muoPtRel, *b_muoDrJet, *b_muoEnergyRatio, *b_muoQCone, *b_muoDz_PV, *b_muoQJet, *b_muoJetCSV, *b_muoJetDFprob, *b_muoJetPt;

private:

PDSecondNtupleData ( const PDSecondNtupleData& a );
PDSecondNtupleData& operator=( const PDSecondNtupleData& a );

};

#endif

