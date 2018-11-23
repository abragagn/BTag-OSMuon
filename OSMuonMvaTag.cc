#include "PDSoftMuonMvaEstimator.h"
#include "AlbertoUtil.h"

OSMuonMvaTag::OSMuonMvaTag(TString weightsFile = "/lustre/cmswork/abragagn/weights/TMVAClassification_BDTOsMuon2016.weights.xml"):
    weightsFile_(weightsFile),
    reader_("!Color:Silent")
{
    TMVA::PyMethodBase::PyInitialize();
    setupReader();
}

OSMuonMvaTag::~OSMuonMvaTag() {}

// =====================================================================================
void OSMuonMvaTag::setupReader()
{

    reader_.AddVariable( "muoPt_", &muoPt_);
    reader_.AddVariable( "muoEta_", &muoEta_);
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
        reader_.AddVariable( "muoJetSize_", &muoJetSize_);
    }

    reader_.AddVariable( "muoQCone_", &muoQCone_);

    reader_.BookMVA( methodNameFromWeightName(), weightsFile_ );

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
