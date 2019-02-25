#ifndef PDAnalyzerUtil_H
#define PDAnalyzerUtil_H

#include "PDAnalysis/Ntu/interface/PDNtupleBranch.h"
#include "NtuAnalysis/Common/interface/NtuInterface.h"

#include "PDAnalysis/Ntu/interface/PDEnumString.h"

class PDAnalyzerUtil: public NtuInterface< PDNtupleBranch<WrapperBase> > {

 public:

  PDAnalyzerUtil();
  virtual ~PDAnalyzerUtil();

 protected:

  void setupNtuple();
  virtual void beginJob();

  // redefine "getEntry" to update HLT path lists when necessary
  // not fully working
//  virtual bool getEntry( int ientry );
  // pre-selection, with minimal data process 
  // before full ntuple entry reading
  virtual bool preSelect( int ientry );

  bool goodJet( unsigned int iJet );
  bool goodMu( unsigned int iMuon );

  bool hlt( PDEnumString::trigPath path );
  std::vector<number>* userInfo( PDEnumString::recoObject r,
                                 PDEnumString::infoType   t );

  const std::vector<int>& muonTracks( PDEnumString::trackType type );
  const std::vector<int>& pfCandFromJet( int iJet );
  const std::vector<int>& tracksFromJet( int iJet );
  const std::vector<int>& tracksFromPV( int iPV );
  const std::vector<int>& tracksFromSV( int iSV );
//  const std::vector<int>& tracksFullSV( int iSV );
  const std::vector<int>& subVtxFromSV( int iSV );

  int muonTrack( int iMu, PDEnumString::trackType type );
  int motherVtxSV( int iSV );

  int tagSVJet( int iJet );

  bool pfCandJetRel( int iPF, int iJet );
  bool tracksJetRel( int iTk, int iJet );
  bool tracksPVRel( int iTk, int iPV );
  bool tracksSVRel( int iTk, int iSV );
  bool subVtxSVRel( int jSV, int iSV );

  int findSVtx( int trkIndex, int prevSVtx = -1 );
  const std::vector<int>& sVtsWithTrack( int trkIndex );
  int ipIndex( int iTk, int iSV );
  int vpIndex( int iTk, int iSV );

  const std::string& svtName( int i );

  number dX ( int iTk );
  number dY ( int iTk );
  number dX ( int iTk, number vx, number vy );
  number dY ( int iTk, number vx, number vy );
  number dZ ( int iTk, number vx, number vy, number vz );
  number dXY( int iTk, number vx, number vy );
  int dSign( int iTk, number px, number py );
  int dSign( int iTk, number px, number py, number vx, number vy );
  int dSign( int iTk, int iJet ) { return dSign( iTk, jetPx->at( iJet ),
                                                      jetPy->at( iJet ) ); }
  int dSign( int iTk, int iJet, number vx, number vy )
                                 { return dSign( iTk, vx, vy, 
                                                      jetPx->at( iJet ),
                                                      jetPy->at( iJet ) ); }

  float jetCharge( int iJet, float k );
  float sVtxCharge( int iSV, float k );

 private:

  // cuts to select good jets
  int    jetNDauMin;
  int    jetNDauMax;
  double jetNHFMin;
  double jetNHFMax;
  double jetNEFMin;
  double jetNEFMax;
  double jetCHFMin;
  double jetCHFMax;
  double jetCEFMin;
  double jetCEFMax;
  double jetNCHMin;
  double jetNCHMax;
  double jetEtaCut;

  double muoChiMin;
  double muoPixHitMin;
  double muoMuoHitMin;
  double muoNTkLayMin;
  double muoNMatchMin;
  double muoDbMax;
  double muoDzMax;

  struct Id {
    unsigned int run;
    unsigned int event;
  };

  Id lastHLT;
  std::map<int,bool> hfMap;
  Id lastInfo;
  std::map< int,std::map< int,std::vector<number> > > uiMap;

  Id lastMuTrk;
  std::map< int,std::vector<int> > muTrkMap;

  Id lastPFCJet;
  std::map< int,std::vector<int> > pfjMap;
  Id lastTrkJet;
  std::map< int,std::vector<int> > tkjMap;
  Id lastTrkPV;
  std::map< int,std::vector<int> > tkpMap;
  Id lastTrkSV;
  std::map< int,std::vector<int> > tksMap;
//  Id lastTkfSV;
//  std::map< int,std::vector<int> > tkfMap;
  Id lastSubSV;
  std::map< int,std::vector<int> > subMap;
  Id lastMotSV;
  std::vector<int> msvList;
  Id lastVtsTk;
  std::map< int,std::vector<int> > vtkMap;

  Id lastJetSV;
  std::vector<int> jsvList;

  const std::vector<int>& listFromObject( int iObj, std::vector<int>* pObj,
                          std::map< int,std::vector<int> >& objMap,
                          Id& lastEID );

  Id lRelPFCJet;
  std::map< int,std::vector<bool> > pfjRelMap;
  Id lRelTrkJet;
  std::map< int,std::vector<bool> > tkjRelMap;
  Id lRelTrkPV;
  std::map< int,std::vector<bool> > tkpRelMap;
  Id lRelTrkSV;
  std::map< int,std::vector<bool> > tksRelMap;
//  Id lRelTkfSV;
//  std::map< int,std::vector<bool> > tkfRelMap;
  Id lRelSubSV;
  std::map< int,std::vector<bool> > subRelMap;

  Id lastIp;
  Id lastVp;
  std::map< int,std::map<int,int> > ipMap;
  std::map< int,std::map<int,int> > vpMap;

  bool comesFromObject( int iProd, int iObj, std::vector<int>* pObj,
                        std::map< int,std::vector<bool> >& objMap,
                        Id& lastEID );

  Id lastAssoc;
  std::vector<int> pfcGenP;
  std::vector<int> genReco;

  void grAssoc();

  bool chkNewEvent( const Id& eID );
  void updateEvent(       Id& eID );

  // dummy copy constructor and assignment
  PDAnalyzerUtil           ( const PDAnalyzerUtil& );
  PDAnalyzerUtil& operator=( const PDAnalyzerUtil& );

};


#endif

