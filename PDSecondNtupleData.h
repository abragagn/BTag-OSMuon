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
    setBranch( "jpsiMass", &jpsiMass, "jpsiMass/F", &b_jpsiMass );
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


    setBranch( "osMuon", &osMuon, "osMuon/I", &b_osMuon );
    setBranch( "osMuonTag", &osMuonTag, "osMuonTag/I", &b_osMuonTag );
    setBranch( "osMuonChargeInfo", &osMuonChargeInfo, "osMuonChargeInfo/I", &b_osMuonChargeInfo );

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

    setBranch( "muoDrB", &muoDrB, "muoDrB/F", &b_muoDrB );
    setBranch( "muoDzPV", &muoDzPV, "muoDzPV/F", &b_muoDzPV );
    setBranch( "muoPFIso", &muoPFIso, "muoPFIso/F", &b_muoPFIso );

    setBranch( "muoJetPt", &muoJetPt, "muoJetPt/F", &b_muoJetPt );
    setBranch( "muoJetPtRel", &muoJetPtRel, "muoJetPtRel/F", &b_muoJetPtRel );
    setBranch( "muoJetDr", &muoJetDr, "muoJetDr/F", &b_muoJetDr );
    setBranch( "muoJetEnergyRatio", &muoJetEnergyRatio, "muoJetEnergyRatio/F", &b_muoJetEnergyRatio );
    setBranch( "muoJetQ", &muoJetQ, "muoJetQ/F", &b_muoJetQ );
    setBranch( "muoJetCSV", &muoJetCSV, "muoJetCSV/F", &b_muoJetCSV );
    setBranch( "muoJetDFprob", &muoJetDFprob, "muoJetDFprob/F", &b_muoJetDFprob );
    setBranch( "muoJetSize", &muoJetSize, "muoJetSize/I", &b_muoJetSize );

    setBranch( "muoSvtPt", &muoSvtPt, "muoSvtPt/F", &b_muoSvtPt );
    setBranch( "muoSvtPtRel", &muoSvtPtRel, "muoSvtPtRel/F", &b_muoSvtPtRel );
    setBranch( "muoSvtDr", &muoSvtDr, "muoSvtDr/F", &b_muoSvtDr );
    setBranch( "muoSvtEnergyRatio", &muoSvtEnergyRatio, "muoSvtEnergyRatio/F", &b_muoSvtEnergyRatio );
    setBranch( "muoSvtQ", &muoSvtQ, "muoSvtQ/F", &b_muoSvtQ );
    setBranch( "muoSvtCSV", &muoSvtCSV, "muoSvtCSV/F", &b_muoSvtCSV );
    setBranch( "muoSvtDFprob", &muoSvtDFprob, "muoSvtDFprob/F", &b_muoSvtDFprob );
    setBranch( "muoSvtSize", &muoSvtSize, "muoSvtSize/I", &b_muoSvtSize );

    setBranch( "muoQCone", &muoQCone, "muoQCone/F", &b_muoQCone );
    setBranch( "muoHowMany", &muoHowMany, "muoHowMany/I", &b_muoHowMany );

    setBranch( "evtNumber", &evtNumber, "evtNumber/I", &b_evtNumber );
    setBranch( "evtWeight", &evtWeight, "evtWeight/F", &b_evtWeight );
    setBranch( "hltJpsiMu", &hltJpsiMu , "hltJpsiMu/I" , &b_hltJpsiMu );
    setBranch( "hltJpsiTrkTrk", &hltJpsiTrkTrk , "hltJpsiTrkTrk/I" , &b_hltJpsiTrkTrk );
    setBranch( "hltJpsiTrk", &hltJpsiTrk , "hltJpsiTrk/I" , &b_hltJpsiTrk );
}

float ssbPt, ssbEta, ssbPhi, ssbMass, jpsiMass, ssbDxy, ssbExy, ssbDz, ssbEz, evtWeight;
float ssbLxy, ssbCt2D, ssbCt2DErr, ssbCt2DSigmaUnit, ssbCt3D, ssbCt3DErr, ssbCt3DSigmaUnit;
int ssbSVT, ssbPVT, ssbLund, evtNumber, hltJpsiMu, hltJpsiTrkTrk, hltJpsiTrk, ssbIsTight;

TBranch *b_ssbPt, *b_ssbEta, *b_ssbPhi, *b_ssbMass, *b_jpsiMass, *b_ssbDxy, *b_ssbExy, *b_ssbDz, *b_ssbEz;
TBranch *b_ssbLxy, *b_ssbCt2D, *b_ssbCt2DErr, *b_ssbCt2DSigmaUnit, *b_ssbCt3D, *b_ssbCt3DErr, *b_ssbCt3DSigmaUnit;
TBranch *b_ssbSVT, *b_ssbPVT, *b_ssbLund, *b_evtNumber, *b_hltJpsiMu, *b_hltJpsiTrkTrk, *b_hltJpsiTrk, *b_ssbIsTight, *b_evtWeight;


int muoLund, muoAncestor, muoCharge;
float muoSoftMvaValue;
float muoPt, muoEta, muoPhi, muoDxy, muoExy, muoDz, muoEz;

TBranch *b_muoPt, *b_muoEta, *b_muoPhi, *b_muoDxy, *b_muoExy, *b_muoDz, *b_muoEz, *b_muoCharge;
TBranch *b_muoLund, *b_muoAncestor, *b_muoSoftMvaValue;


int osMuon, osMuonTag, osMuonChargeInfo;
float muoDrB, muoDzPV, muoPFIso, muoQCone;

int muoJetSize, muoHowMany;
float muoJetPtRel, muoJetDr, muoJetEnergyRatio, muoJetQ, muoJetCSV, muoJetDFprob, muoJetPt;

int muoSvtSize;
float muoSvtPtRel, muoSvtDr, muoSvtEnergyRatio, muoSvtQ, muoSvtCSV, muoSvtDFprob, muoSvtPt;

TBranch *b_osMuon, *b_osMuonTag, *b_osMuonChargeInfo, *b_muoDrB, *b_muoPFIso, *b_muoQCone, *b_muoDzPV, *b_muoHowMany;
TBranch  *b_muoJetSize, *b_muoJetPtRel, *b_muoJetDr, *b_muoJetEnergyRatio, *b_muoJetQ, *b_muoJetCSV, *b_muoJetDFprob, *b_muoJetPt;
TBranch  *b_muoSvtSize, *b_muoSvtPtRel, *b_muoSvtDr, *b_muoSvtEnergyRatio, *b_muoSvtQ, *b_muoSvtCSV, *b_muoSvtDFprob, *b_muoSvtPt;

private:

PDSecondNtupleData ( const PDSecondNtupleData& a );
PDSecondNtupleData& operator=( const PDSecondNtupleData& a );

};

#endif

