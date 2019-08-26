// -*- C++ -*-
//
// Package:    PFAnalysis/PFAnalyzers
// Class:      PFTrackHFAnalyzer
//
/**\class PFTrackHFAnalyzer PFTrackHFAnalyzer.cc PFAnalysis/PFAnalyzers/plugins/PFTrackHFAnalyzer.cc

  Description: Analyzer of PFTracks/clusters in HF region
          The input step3 files will need: --outputCommand 'keep recoPFRecHits_particleFlow*_*_*','keep *_*pfTrack*_*_*'

  Implementation:
    [Notes on implementation]
*/
//
// Original Author:  Kenichi Hatakeyama
//         Created:  Wed, 17 Jul 2019 15:24:34 GMT
//
//

// system include files
#include <memory>
#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"  
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
//#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "boost/format.hpp"

#include "TH2.h"
#include "TH1.h"
#include "TProfile2D.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

bool mysort(const reco::PFCandidate& elem1, const reco::PFCandidate& elem2){ return elem1.energy() > elem2.energy() ; }

class PFTrackHFAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit PFTrackHFAnalyzer(const edm::ParameterSet&);
    ~PFTrackHFAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<reco::GenParticleCollection> genparToken_; 
    edm::EDGetTokenT<CaloParticleCollection> caloparToken_; 
    edm::EDGetTokenT<reco::VertexCollection> vertexToken_; 
    edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfcandToken_; 
    edm::EDGetTokenT<std::vector<reco::PFCluster>> pfclusterHFToken_; 
    edm::EDGetTokenT<std::vector<reco::PFRecHit>> pfrechitHFToken_; 
    edm::EDGetTokenT<std::vector<reco::PFRecTrack>> pftrackToken_; 
    edm::EDGetTokenT<reco::TrackCollection> trackToken_;  //used to select what tracks to read from configuration file
    edm::EDGetTokenT<HFRecHitCollection> hfrechitToken_;

    bool debug_;
    bool debugRecHit_;

    int nev = 0 ;
    std::vector<int> EventsToScan_;

    double ptlow_;
    double pthigh_;
    double etalow_;
    double etahigh_;

//    const std::vector<double> eventsToScan;
  
    TH1I * histo;
    map<TString, TH1*> m_Histos1D;
    map<TString, TH2*> m_Histos2D;
    map<TString, TProfile2D*> m_Profiles2D;
    void FillHist1D(const TString& histName, const Double_t& value, const double& weight);
    void FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight);
    void FillProfile2D(const TString& histName, const Double_t& value1, const Double_t& value2, const Double_t& value3, const double& weight);
    bool IsInCell(double eta, double phi, const CaloCellGeometry::CornersVec& CV);  
    int idphi(int iphiA, int iphiB);        



};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PFTrackHFAnalyzer::PFTrackHFAnalyzer(const edm::ParameterSet& iConfig)
  :
  genparToken_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_genpars"))),
  caloparToken_(consumes<CaloParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_calopars"))),
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_vertices"))),
  pfcandToken_(consumes<std::vector<reco::PFCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pfcands"))),
  pfclusterHFToken_(consumes<std::vector<reco::PFCluster>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pfclustersHF"))),
  pfrechitHFToken_(consumes<std::vector<reco::PFRecHit>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pfrechitsHF"))),
  pftrackToken_(consumes<std::vector<reco::PFRecTrack>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pftracks"))),
  trackToken_(consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_tracks"))),
  hfrechitToken_(consumes<HFRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_hfrechits"))),

  debug_(iConfig.getUntrackedParameter<bool>("debug")),
  debugRecHit_(iConfig.getUntrackedParameter<bool>("debugRecHit")),
  EventsToScan_(iConfig.getParameter<std::vector<int>>("EventsToScan")),
  ptlow_(iConfig.getParameter<double>("ptlow")),
  pthigh_(iConfig.getParameter<double>("pthigh")),
  etalow_(iConfig.getParameter<double>("etalow")),
  etahigh_(iConfig.getParameter<double>("etahigh"))

{
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  TString hname;

  hname = "gen_pt";   
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 250 , 0. , 250. );
  hname = "gen_eta";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 120 , 0. , 6. );


  hname = "pftrackHF_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 21 , -0.5 , 20.5 );
  hname = "pftrackHF_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 250 , 0. , 250. );
  hname = "pftrackHF_eta";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 120 , -6. , 6. );
  hname = "pftrackHF_pull";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , -50. , 50. );

  hname = "pfrechitHFEM_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , -0.5 , 199.5 );
  hname = "pfrechitHFHAD_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , -0.5 , 199.5 );

  hname = "pfrechitHFEM_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 200. );
  hname = "pfrechitHFHAD_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 200. );

  hname = "pfrechitHFEM_Emax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );
  hname = "pfrechitHFHAD_Emax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );
  hname = "pfrechitHF_Emax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );

  hname = "GProfileHFEMP_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "GProfileHFHADP_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "GProfileHFEMN_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "GProfileHFHADN_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "TProfileHFEMP_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "TProfileHFHADP_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "TProfileHFEMN_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "TProfileHFHADN_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );

  TFileDirectory subDir = fs->mkdir( "ScanEvents" );
  for (unsigned int i=0; i < EventsToScan_.size(); i++){
    hname = Form("GenEMP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("GenHADP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("ClustersEMP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("ClustersHADP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("EClustersEMP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("EClustersHADP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("ERClustersEMP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("ERClustersHADP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
  }
  
  hname = "pfclusHFEM_nhits";   
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 12 , -1.5 , 10.5 );
  hname = "pfclusHFHAD_nhits";  
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 12 , -1.5 , 10.5 );

  hname = "pfclusHFEM_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , -0.5 , 99.5 );
  hname = "pfclusHFHAD_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , -0.5 , 99.5 );

  hname = "pfclusHFEM_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 200. );
  hname = "pfclusHFHAD_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 200. );

  hname = "pfclusHFEM_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 200. );
  hname = "pfclusHFHAD_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 200. );

  hname = "pfclusHFEM_Emax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );
  hname = "pfclusHFHAD_Emax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. ); 
  hname = "pfclusHF_Emax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. ); 

  hname = "pfclusHFEM_ptmax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );
  hname = "pfclusHFHAD_ptmax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );
  hname = "pfclusHF_ptmax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );

  hname = "pfclusHF_nmatch";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 22 , -1.5 , 20.5, 22 , -1.5 , 20.5);
  hname = "pfclusHF_nmatch5";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 22 , -1.5 , 20.5, 22 , -1.5 , 20.5 );  
  hname = "pfclusHF_nmatch9";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 22 , -1.5 , 20.5, 22 , -1.5 , 20.5 );

  hname = "match_ptfrac";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100, 0. , 2.0 );   
  hname = "match_ptfrac5";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100, 0. , 2.0 );
  hname = "match_ptfrac9";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100, 0. , 2.0 );

  hname = "pfcandHFEM_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , -0.5 , 99.5 ); 
  hname = "pfcandHFHAD_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , -0.5 , 99.5 );
  hname = "pfcandHFCH_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , -0.5 , 99.5 );


  hname = "pfcandHFEM_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 100. );
  hname = "pfcandHFHAD_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 100. );
  hname = "pfcandHFCH_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 100. );

  hname = "pfcandHFEM_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 150 , 0. , 1500. );
  hname = "pfcandHFHAD_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 150 , 0. , 1500. );
  hname = "pfcandHFCH_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 150 , 0. , 1500. );


  hname = "pfcandHFEM_nelements" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 25 , -0.5, 24.5 );
  hname = "pfcandHFHAD_nelements" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 25 , -0.5, 24.5 );
  hname = "pfcandHFCH_nelements" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 25 , -0.5, 24.5 );

  hname = "pfcandHFEM_n9" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 16 , -0.5, 15.5 );
  hname = "pfcandHFHAD_n9" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 16 , -0.5, 15.5 );  
  hname = "pfcandHFCH_n9" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 16 , -0.5, 15.5 );

  hname = "Eratio_PFcand1toGen";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );
  hname = "Eratio_PFcand2toGen";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );
  hname = "Eratio_PFcand3toGen";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );

  hname = "Eratio_PFcandAlltoGen";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );
  hname = "Eratio_PFcandNonTtoGen";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );

  hname = "Eratio_PFcandEM1toGen";  
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );
  hname = "Eratio_PFcandHAD1toGen";  
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );
  hname = "Eratio_PFcandEMHAD1toGen";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );

  hname = "PovGenE_brem";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );
  hname = "PovGenE_nobrem";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );

  hname = "EovP_brem";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );
  hname = "EovP_nobrem";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );

}


PFTrackHFAnalyzer::~PFTrackHFAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PFTrackHFAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<reco::GenParticleCollection> genpars; iEvent.getByToken(genparToken_, genpars);
  edm::Handle<CaloParticleCollection> calopars; iEvent.getByToken(caloparToken_, calopars);
  edm::Handle<reco::VertexCollection> vertices; iEvent.getByToken(vertexToken_, vertices);
  edm::Handle<std::vector<reco::PFCandidate>> pfcands; iEvent.getByToken(pfcandToken_, pfcands);
  edm::Handle<std::vector<reco::PFCluster>> pfclustersHF; iEvent.getByToken(pfclusterHFToken_, pfclustersHF);
  edm::Handle<std::vector<reco::PFRecHit>> pfrechitsHF; iEvent.getByToken(pfrechitHFToken_, pfrechitsHF);
  edm::Handle<std::vector<reco::PFRecTrack>> pftracks; iEvent.getByToken(pftrackToken_, pftracks);
  edm::Handle<HFRecHitCollection> hfRecHits;   iEvent.getByToken(hfrechitToken_, hfRecHits);

  TString hname;

  nev++ ;
  std::vector<int>::iterator it;
  it = std::find (EventsToScan_.begin(), EventsToScan_.end(), nev);
  if (it != EventsToScan_.end())
  cout << " Will scan this event " << EventsToScan_[it-EventsToScan_.begin()] << endl ;


  if (debug_) 
  LogPrint("PFTrackHFAnalyzer") << "\n\n ======== genpars: ======== "  << genpars->size();
  double genP_pt = 0., genP_E = 0., genP_eta = 0., genP_phi = 0.,
         genN_pt = 0., genN_E = 0., genN_eta = 0., genN_phi = 0. ;
  int    genP_ieta1 = -1, genP_ieta2 = -1, genP_iphi1 = -1, genP_iphi2 = -1,
         genN_ieta1 = -1, genN_ieta2 = -1, genN_iphi1 = -1, genN_iphi2 = -1;
  for(const auto& genpar : *(genpars.product()) ){
    if (debug_)
    LogPrint("PFTrackHFAnalyzer") << boost::format("genpar (pt,eta,phi,E): (%6.2f, %6.2f, %6.2f, %6.2f)") % genpar.pt() % genpar.eta() % genpar.phi() % genpar.energy();
    if(genpar.eta()>0.){
      genP_pt  = genpar.pt();
      genP_E   = genpar.energy();
      genP_eta = genpar.eta();
      genP_phi = genpar.phi();
    } else {
      genN_pt  = genpar.pt();
      genN_E   = genpar.energy();
      genN_eta = genpar.eta();
      genN_phi = genpar.phi();
    }
  }

  if(genP_pt <ptlow_ || genP_pt > pthigh_)                 return ;
  if(fabs(genP_eta) <etalow_ || fabs(genP_eta) > etahigh_) return ;

  FillHist1D("gen_pt",  genP_pt,  1. );
  FillHist1D("gen_eta", genP_eta,  1. );

  if (debug_)
  LogPrint("PFTrackHFAnalyzer") << "\n =========== calopars: =========== "     << calopars->size();
  for(const auto& calopar : *(calopars.product()) ){
    if (debug_)
    LogPrint("PFTrackHFAnalyzer") << boost::format("calopar (pt,eta,phi): (%6.1f, %6.2f, %6.2f)") % calopar.pt() % calopar.eta() % calopar.phi();
  }


  if (debug_)
  LogPrint("PFTrackHFAnalyzer") << "\n =========== vertices: =========== "     << vertices->size();


  if (debug_)
  LogPrint("PFTrackHFAnalyzer") << "\n =========== pftracks: =========== "     << pftracks->size();
  int pftrackHF_n = 0;
  double trkP_P = 0., trkP_pt  = 0., trkP_pterror = 0., trkP_eta = 0., trkP_phi = 0.,
         trkN_P = 0., trkN_pt  = 0., trkN_pterror = 0., trkN_eta = 0., trkN_phi = 0.;
  int    trkP_ieta1 = -1, trkP_ieta2 = -1, trkP_iphi1 = -1, trkP_iphi2 = -1,
         trkN_ieta1 = -1, trkN_ieta2 = -1, trkN_iphi1 = -1, trkN_iphi2 = -1;
    
  for(const auto& pftrack : *(pftracks.product()) ){
    pftrackHF_n++;
    const reco::TrackRef trackref = pftrack.trackRef();
    if (debug_)
    LogPrint("PFTrackHFAnalyzer") << boost::format("pftrack (pt,eta,phi)@origin : (%6.1f +- %4.1f, %6.2f, %6.2f) ") % trackref->pt() % trackref->ptError() % trackref->eta() % trackref->phi() ;
    if(trackref->eta()>0.){   
      trkP_P   = trackref->p();
      trkP_pt  = trackref->pt(); 
      trkP_pterror = trackref->ptError();
      trkP_eta = trackref->eta();
      trkP_phi = trackref->phi();
      FillHist1D("pftrackHF_pt",  trkP_pt,  1. );
      FillHist1D("pftrackHF_eta", trkP_eta,  1. );
      FillHist1D("pftrackHF_pull",  ((trkP_pt-genP_pt)/trkP_pterror),  1. );
    } else {
      trkN_P   = trackref->p();
      trkN_pt  = trackref->pt();
      trkN_pterror = trackref->ptError();
      trkN_eta = trackref->eta();
      trkN_phi = trackref->phi();
      FillHist1D("pftrackHF_pt",  trkN_pt,  1. );
      FillHist1D("pftrackHF_eta", trkN_eta,  1. );
      FillHist1D("pftrackHF_pull",  ((trkN_pt-genN_pt)/trkN_pterror),  1. );
    }

/* ctmp      
    std::vector<reco::PFTrajectoryPoint> trajectoryPoints = pftrack.trajectoryPoints();
    
    constexpr reco::PFTrajectoryPoint::LayerType VFcalEntrance =
    reco::PFTrajectoryPoint::VFcalEntrance;
    
    const reco::PFTrajectoryPoint& tkAtHF =
    pftrack.extrapolatedPoint( VFcalEntrance );
    
    const double tracketa = tkAtHF.positionREP().Eta();
    const double trackphi = tkAtHF.positionREP().Phi();
    
    LogPrint("PFTrackHFAnalyzer") << boost::format("pftrack (pt,eta,phi)@origin (eta,phi)@HF: (%6.1f +- %4.1f, %6.2f, %6.2f) (%6.2f, %6.2f)")
    % trackref->pt() % trackref->ptError() % trackref->eta() % trackref->phi() % tracketa % trackphi;
*/

  }
  if(pftrackHF_n> 2) cout << " Warning!! More than two tracks. Scan this event # " << nev <<  endl ;
  FillHist1D("pftrackHF_n",  pftrackHF_n,  1. );

  if (debug_)
  LogPrint("PFTrackHFAnalyzer") << "\n  =========== HFRecHits: ===========  " <<  hfRecHits->size();;

  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  const CaloSubdetectorGeometry* HFGeom = geoHandle->getSubdetectorGeometry(DetId::Hcal, 4);

  for(const auto& rechit : *(hfRecHits.product()) ){
    int ieta = rechit.id().ieta() ;
    int iphi = rechit.id().iphi() ;
    int idep = rechit.id().depth() ;
    double E = rechit.energy() ;
    if (debugRecHit_)
    LogPrint("PFTrackHFAnalyzer") << boost::format(" hfrechit (ieta, iphi, depth, E): (%3d, %3d, %2d, %6.2f)") % ieta % iphi % idep % E ;

    auto thisCell = HFGeom->getGeometry(rechit.id().rawId());
    const CaloCellGeometry::CornersVec& cv = thisCell->getCorners();
    if (debugRecHit_)
    LogPrint("PFTrackHFAnalyzer") << boost::format(" corners (eta, eta, eta, eta, phi, phi, phi, phi): (%6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f)")
        % cv[0].eta() % cv[1].eta() % cv[2].eta() % cv[3].eta() % cv[0].phi() % cv[1].phi() % cv[2].phi() % cv[3].phi()  ;
    
    bool incell ;
    incell = IsInCell(genP_eta, genP_phi, cv) ;  
    if(incell) {
      if(idep == 1){
        genP_ieta1 = ieta ; genP_iphi1 = iphi ;
      } else {
        genP_ieta2 = ieta ; genP_iphi2 = iphi ;
      }
    } 
    incell = IsInCell(genN_eta, genN_phi, cv) ;
    if(incell) {
      if(idep == 1){ 
        genN_ieta1 = ieta ; genN_iphi1 = iphi ;
      } else {
        genN_ieta2 = ieta ; genN_iphi2 = iphi ;
      }
    }
    incell = IsInCell(trkP_eta, trkP_phi, cv) ;
    if(incell) {
      if(idep == 1){
        trkP_ieta1 = ieta ; trkP_iphi1 = iphi ;
      } else {
        trkP_ieta2 = ieta ; trkP_iphi2 = iphi ;
      }
    }
    incell = IsInCell(trkN_eta, trkN_phi, cv) ;
    if(incell) {
      if(idep == 1){
        trkN_ieta1 = ieta ; trkN_iphi1 = iphi ;
      } else {
        trkN_ieta2 = ieta ; trkN_iphi2 = iphi ;
      }
    }  
  } // end of loop over HFRecHits ...

  if (it != EventsToScan_.end()){ 
    hname = Form("GenEMP%i", nev);
    FillHist2D(hname, genP_ieta1,genP_iphi1,genP_E) ;
    hname = Form("GenHADP%i", nev);
    FillHist2D(hname, genP_ieta2,genP_iphi2,genP_E) ;
  }

  if (debug_){
    LogPrint("PFTrackHFAnalyzer") << boost::format("\nParticle in positive eta with (eta phi): (%6.2f, %6.2f)") % genP_eta % genP_phi;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth1 (ieta iphi): (%3d, %3d)") % genP_ieta1 % genP_iphi1  ;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth2 (ieta iphi): (%3d, %3d)") % genP_ieta2 % genP_iphi2  ;
    LogPrint("PFTrackHFAnalyzer") << boost::format(" Particle in negative eta with (eta phi): (%6.2f, %6.2f)") % genN_eta % genN_phi ;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth1 (ieta iphi): (%3d, %3d)") % genN_ieta1 % genN_iphi1  ;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth2 (ieta iphi): (%3d, %3d)") % genN_ieta2 % genN_iphi2  ;

    LogPrint("PFTrackHFAnalyzer") << boost::format("\nTrack in positive eta with (eta phi): (%6.2f, %6.2f)") % trkP_eta % trkP_phi;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth1 (ieta iphi): (%3d, %3d)") % trkP_ieta1 % trkP_iphi1  ;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth2 (ieta iphi): (%3d, %3d)") % trkP_ieta2 % trkP_iphi2  ;
    LogPrint("PFTrackHFAnalyzer") << boost::format(" Track in negative eta with (eta phi): (%6.2f, %6.2f)") % trkN_eta % trkN_phi ;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth1 (ieta iphi): (%3d, %3d)") % trkN_ieta1 % trkN_iphi1  ;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth2 (ieta iphi): (%3d, %3d)") % trkN_ieta2 % trkN_iphi2  ;
  }


  if (debug_)
  LogPrint("PFTrackHFAnalyzer") << "\n =========== pfrechitsHF: =========== "  << pfrechitsHF->size();
  int pfrechitHFEM_n = 0, pfrechitHFHAD_n = 0 ;
  double pfrechitHFEMP_Emax = -1.,  pfrechitHFHADP_Emax = -1.,
         pfrechitHFEMN_Emax = -1.,  pfrechitHFHADN_Emax = -1. ;

  for(const auto& pfrechit : *(pfrechitsHF.product()) ){
    int ieta = HcalDetId(pfrechit.detId()).ieta() ;
    int iphi = HcalDetId(pfrechit.detId()).iphi() ;
    int idep = HcalDetId(pfrechit.detId()).depth() ;
    double E = pfrechit.energy() ;
    if (debug_)
    LogPrint("PFTrackHFAnalyzer") << boost::format(" pfrechit (ieta, iphi, depth, E): (%3d, %3d, %2d, %6.2f)")  % ieta % iphi % idep % E;

    if(idep == 1) { // HFEM
      pfrechitHFEM_n++ ;
      FillHist1D("pfrechitHFEM_E", E, 1. );
      if(ieta > 0){
        if(E > pfrechitHFEMP_Emax ) pfrechitHFEMP_Emax = E ; 
        FillHist2D("GProfileHFEMP_E", ieta-genP_ieta1 , idphi(iphi, genP_iphi1), E/genP_E);
        FillHist2D("TProfileHFEMP_E", ieta-trkP_ieta1 , idphi(iphi, trkP_iphi1), E/trkP_P);
      } else {
        if(E > pfrechitHFEMN_Emax ) pfrechitHFEMN_Emax = E ;
        FillHist2D("GProfileHFEMN_E", ieta-genN_ieta1,  idphi(iphi, genN_iphi1), E/genN_E);
        FillHist2D("TProfileHFEMN_E", ieta-trkN_ieta1,  idphi(iphi, trkN_iphi1), E/trkN_P);
      }
    } else { // HFHAD
      pfrechitHFHAD_n++ ;
      FillHist1D("pfrechitHFHAD_E", E, 1. );
      if(ieta > 0){
        if(E > pfrechitHFHADP_Emax ) pfrechitHFHADP_Emax = E ;
        FillHist2D("GProfileHFHADP_E", ieta-genP_ieta2,  idphi(iphi, genP_iphi2), E/genP_E);
        FillHist2D("TProfileHFHADP_E", ieta-trkP_ieta2,  idphi(iphi, trkP_iphi2), E/trkP_P);
      } else {
        if(E > pfrechitHFHADN_Emax ) pfrechitHFHADN_Emax = E ;
        FillHist2D("GProfileHFHADN_E", ieta-genN_ieta2,  idphi(iphi, genN_iphi2), E/genN_E);
        FillHist2D("TProfileHFHADN_E", ieta-trkN_ieta2,  idphi(iphi, trkN_iphi2), E/trkN_P);
      }
    }
  } // end of loop obver pfrechitsHF
  if (debug_){
    LogPrint("PFTrackHFAnalyzer") << "  pfrechitsHFEM:  " << pfrechitHFEM_n;
    LogPrint("PFTrackHFAnalyzer") << "  pfrechitsHFHAD: " << pfrechitHFHAD_n;
  }
  FillHist1D("pfrechitHFEM_n",  pfrechitHFEM_n,  1. );
  FillHist1D("pfrechitHFHAD_n", pfrechitHFHAD_n, 1. );
  FillHist1D("pfrechitHFEM_Emax",  pfrechitHFEMP_Emax,  1. );
  FillHist1D("pfrechitHFEM_Emax",  pfrechitHFEMN_Emax,  1. );
  FillHist1D("pfrechitHFHAD_Emax", pfrechitHFHADP_Emax,  1. );
  FillHist1D("pfrechitHFHAD_Emax", pfrechitHFHADN_Emax,  1. );
  FillHist1D("pfrechitHF_Emax",  pfrechitHFEMP_Emax > pfrechitHFHADP_Emax ? pfrechitHFEMP_Emax : pfrechitHFHADP_Emax ,  1. );
  FillHist1D("pfrechitHF_Emax",  pfrechitHFEMN_Emax > pfrechitHFHADN_Emax ? pfrechitHFEMN_Emax : pfrechitHFHADN_Emax ,  1. );

  if (debug_)
  LogPrint("PFTrackHFAnalyzer") << "\n =========== pfclustersHF: =========== " << pfclustersHF->size();

  int pfclusHFEM_n  = 0, pfclusHFHAD_n  = 0,
      pfclusHFEMP_n = 0, pfclusHFHADP_n = 0 ;
  double pfclusHFEMP_Emax  = -1., pfclusHFEMN_Emax  = -1.,  pfclusHFHADP_Emax  = -1., pfclusHFHADN_Emax  = -1. ;
  double pfclusHFEMP_ptmax = -1., pfclusHFEMN_ptmax = -1.,  pfclusHFHADP_ptmax = -1., pfclusHFHADN_ptmax = -1. ;

  int pfclusHFEMP_nmatch = 0, pfclusHFHADP_nmatch = 0,  pfclusHFEMP_nmatch5 = 0, pfclusHFHADP_nmatch5 = 0,  pfclusHFEMP_nmatch9 = 0, pfclusHFHADP_nmatch9 = 0,
      pfclusHFEMN_nmatch = 0, pfclusHFHADN_nmatch = 0,  pfclusHFEMN_nmatch5 = 0, pfclusHFHADN_nmatch5 = 0,  pfclusHFEMN_nmatch9 = 0, pfclusHFHADN_nmatch9 = 0;

  double pfclusHFP_pttot = 0., pfclusHFP_pttot5 = 0., pfclusHFP_pttot9 = 0., 
         pfclusHFN_pttot = 0., pfclusHFN_pttot5 = 0., pfclusHFN_pttot9 = 0.;

  for(const auto& pfclus : *(pfclustersHF.product()) ){
    double eta = pfclus.eta() ;
    double phi = pfclus.phi() ;
    double pt  = pfclus.pt() ;
    double E   = pfclus.energy() ;
    int layer  = pfclus.layer() ; 
    int idep   = fabs(pfclus.depth() - 1.) < 0.001 ? 1 : 2 ;   
    const std::vector<reco::PFRecHitFraction> &fracs = pfclus.recHitFractions();   
    const std::vector<std::pair<DetId, float>> &hfracs = pfclus.hitsAndFractions();
    unsigned nhits = fracs.size();

    if (debug_)
    LogPrint("PFTrackHFAnalyzer") << boost::format("pfclus (pt,E,eta,phi,layer,depth,nhits): (%6.2f, %6.2f, %6.2f, %6.2f, %i, %i, %i)")
      % pt % E % eta % phi % layer % idep % nhits;

    bool matchP = false,  match4P = false,  match8P = false,
         matchN = false,  match4N = false,  match8N = false;

    for(unsigned i=0; i<nhits; i++) { 
//    const auto& id = hfracs[i].first.rawId();
      const reco::PFRecHitRef& pfRecHits = fracs[i].recHitRef();
      double rawenergy = pfRecHits->energy();
      double frac      = fracs[i].fraction();
      int ieta  = HcalDetId(pfRecHits->detId()).ieta() ;
      int iphi  = HcalDetId(pfRecHits->detId()).iphi() ;
    
      if (it != EventsToScan_.end()){ 
        if(idep == 1){
          hname = Form("ClustersEMP%i", nev);
          FillHist2D(hname, ieta,iphi,pfclusHFEMP_n+1) ; 
          hname = Form("EClustersEMP%i", nev);
          FillHist2D(hname, ieta,iphi,E) ;
          hname = Form("ERClustersEMP%i", nev);
          FillHist2D(hname, ieta,iphi,rawenergy*frac) ;
        } else {
          hname = Form("ClustersHADP%i", nev);
          FillHist2D(hname, ieta,iphi,pfclusHFHADP_n+1) ;
          hname = Form("EClustersHADP%i", nev);
          FillHist2D(hname, ieta,iphi,E) ;
          hname = Form("ERClustersHADP%i", nev);
          FillHist2D(hname, ieta,iphi,rawenergy*frac) ;
        } 
      } 

      if (debug_)
      LogPrint("PFTrackHFAnalyzer") << boost::format(" pfrechit (ieta, iphi, depth, E, frac): (%3d, %3d, %2d, %6.2f, %6.2f)")
      % HcalDetId(pfRecHits->detId()).ieta() % HcalDetId(pfRecHits->detId()).iphi()  % HcalDetId(pfRecHits->detId()).depth() % rawenergy % frac ;

      if(idep == 1){
        if( ieta == genP_ieta1 && iphi == genP_iphi1 ) matchP = true ;
        if( ieta == genN_ieta1 && iphi == genN_iphi1 ) matchN = true ;
        if( (abs(ieta - genP_ieta1) == 1 && iphi == genP_iphi1) || (abs(idphi(iphi, genP_iphi1)) == 1 && genP_ieta1 == ieta) ) match4P = true ;
        if( (abs(ieta - genN_ieta1) == 1 && iphi == genN_iphi1) || (abs(idphi(iphi, genN_iphi1)) == 1 && genN_ieta1 == ieta) ) match4N = true ;
        if( abs(ieta - genP_ieta1) == 1 && abs(idphi(iphi, genP_iphi1)) == 1 ) match8P = true ;
        if( abs(ieta - genN_ieta1) == 1 && abs(idphi(iphi, genN_iphi1)) == 1 ) match8N = true ;
      } else {
        if( ieta == genP_ieta2 && iphi == genP_iphi2 ) matchP = true ;
        if( ieta == genN_ieta2 && iphi == genN_iphi2 ) matchN = true ;
        if( (abs(ieta - genP_ieta2) == 1 && iphi == genP_iphi2) || (abs(idphi(iphi, genP_iphi2)) == 1 && genP_ieta2 == ieta) ) match4P = true ;
        if( (abs(ieta - genN_ieta2) == 1 && iphi == genN_iphi2) || (abs(idphi(iphi, genN_iphi2)) == 1 && genN_ieta2 == ieta) ) match4N = true ;
        if( abs(ieta - genP_ieta2) == 1 && abs(idphi(iphi, genP_iphi2)) == 1 ) match8P = true ;
        if( abs(ieta - genN_ieta2) == 1 && abs(idphi(iphi, genN_iphi2)) == 1 ) match8N = true ;
      }

      if (debug_ && (matchP || match4P || match8P) ){
        LogPrint("PFTrackHFAnalyzer") << boost::format(" => Cluster matched to genP particle  with (eta phi E) (ieta1 iphi1) (ieta2 iphi2):  (%6.2f, %6.2f, %6.2f), (%i, %i), (%i, %i)")
        % genP_eta % genP_phi % genP_E % genP_ieta1 % genP_iphi1 % genP_ieta2 % genP_iphi2 ;
      }
      if (debug_ && (matchN || match4N || match8N) ){
        LogPrint("PFTrackHFAnalyzer") << boost::format(" => Cluster matched to genN particle  with (eta phi E) (ieta1 iphi1) (ieta2 iphi2):  (%6.2f, %6.2f, %6.2f), (%i, %i), (%i, %i)")
        % genN_eta % genN_phi % genN_E % genN_ieta1 % genN_iphi1 % genN_ieta2 % genN_iphi2 ;
      }  
    }

    if(matchP) {
      if (idep == 1) pfclusHFEMP_nmatch++ ;
      else           pfclusHFHADP_nmatch++ ;
      pfclusHFP_pttot += pt;
    }
    if(matchN) {
      if (idep == 1) pfclusHFEMN_nmatch++ ;
      else           pfclusHFHADN_nmatch++ ;
      pfclusHFN_pttot += pt;
    }

    if(matchP || match4P ) {
      if (idep == 1) pfclusHFEMP_nmatch5++ ;
      else           pfclusHFHADP_nmatch5++ ;
      pfclusHFP_pttot5 += pt;
    }
      if(matchN || match4N ) {
      if (idep == 1) pfclusHFEMN_nmatch5++ ;
      else           pfclusHFHADN_nmatch5++ ;
      pfclusHFN_pttot5 += pt;
    }

    if(matchP || match4P || match8P ) {
      if (idep == 1) pfclusHFEMP_nmatch9++ ;
      else           pfclusHFHADP_nmatch9++ ;
      pfclusHFP_pttot9 += pt;
    } 
    if(matchN || match4N || match8N ) {
      if (idep == 1) pfclusHFEMN_nmatch9++ ;
      else           pfclusHFHADN_nmatch9++ ;
      pfclusHFN_pttot9 += pt;
    }

    if(idep == 1){
      pfclusHFEM_n++  ;
      if( eta > 0 )  pfclusHFEMP_n++  ;
      FillHist1D("pfclusHFEM_nhits",  nhits, 1.);
      FillHist1D("pfclusHFEM_E", E, 1.);
      FillHist1D("pfclusHFEM_pt",pt, 1.);
      if(eta > 0){
        if(E > pfclusHFEMP_Emax ) pfclusHFEMP_Emax  = E ;
        if(pt> pfclusHFEMP_ptmax) pfclusHFEMP_ptmax = pt ;
      } else {
        if(E > pfclusHFEMN_Emax ) pfclusHFEMN_Emax  = E ;
        if(pt> pfclusHFEMN_ptmax) pfclusHFEMN_ptmax = pt ;
      }
    } else {
      pfclusHFHAD_n++ ;
      if( eta > 0 )  pfclusHFHADP_n++  ;
      FillHist1D("pfclusHFHAD_nhits", nhits, 1.);
      FillHist1D("pfclusHFHAD_E", E, 1.);
      FillHist1D("pfclusHFHAD_pt",pt, 1.);
      if(eta > 0){
        if(E > pfclusHFHADP_Emax ) pfclusHFHADP_Emax  = E ;
        if(pt> pfclusHFHADP_ptmax) pfclusHFHADP_ptmax = pt ;
      } else {
        if(E > pfclusHFHADN_Emax ) pfclusHFHADN_Emax  = E ;
        if(pt> pfclusHFHADN_ptmax) pfclusHFHADN_ptmax = pt ;
      }
    } 
  } // end of loop over pfclusHF

  if (debug_){
    LogPrint("PFTrackHFAnalyzer") << "\n  pfclusHFEM_n:  " << pfclusHFEM_n;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHAD_n: " << pfclusHFHAD_n;

    LogPrint("PFTrackHFAnalyzer") << "\n  pfclusHFEMP_nmatch:  " << pfclusHFEMP_nmatch;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHADP_nmatch: " << pfclusHFHADP_nmatch;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFP_pttot:     " << pfclusHFP_pttot;

    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFEMN_nmatch:  " << pfclusHFEMN_nmatch;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHADN_nmatch: " << pfclusHFHADN_nmatch;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFN_pttot:     " << pfclusHFN_pttot;

    LogPrint("PFTrackHFAnalyzer") << "\n  pfclusHFEMP_nmatch5:  " << pfclusHFEMP_nmatch5;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHADP_nmatch5: " << pfclusHFHADP_nmatch5;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFP_pttot5:     " << pfclusHFP_pttot5;

    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFEMN_nmatch5:  " << pfclusHFEMN_nmatch5;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHADN_nmatch5: " << pfclusHFHADN_nmatch5;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFN_pttot5:     " << pfclusHFN_pttot5;

    LogPrint("PFTrackHFAnalyzer") << "\n  pfclusHFEMP_nmatch9:  " << pfclusHFEMP_nmatch9;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHADP_nmatch9: " << pfclusHFHADP_nmatch9;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFP_pttot9:     " << pfclusHFP_pttot9;
    
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFEMN_nmatch9:  " << pfclusHFEMN_nmatch9;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHADN_nmatch9: " << pfclusHFHADN_nmatch9;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFN_pttot9:     " << pfclusHFN_pttot9;

    cout << " pfclusHFEMP_ptmax  pfclusHFHADP_ptmax " << pfclusHFEMP_ptmax << "   " << pfclusHFHADP_ptmax << endl ;
    cout << " pfclusHFEMN_ptmax  pfclusHFHADN_ptmax " << pfclusHFEMN_ptmax << "   " << pfclusHFHADN_ptmax << endl ;
  }
  FillHist1D("pfclusHFEM_n",   pfclusHFEM_n, 1.);
  FillHist1D("pfclusHFHAD_n",  pfclusHFHAD_n, 1.);

  FillHist2D("pfclusHF_nmatch",  pfclusHFEMP_nmatch,  pfclusHFHADP_nmatch,  1.);
  FillHist2D("pfclusHF_nmatch",  pfclusHFEMN_nmatch,  pfclusHFHADN_nmatch,  1.);

  FillHist2D("pfclusHF_nmatch5", pfclusHFEMP_nmatch5, pfclusHFHADP_nmatch5, 1.);
  FillHist2D("pfclusHF_nmatch5", pfclusHFEMN_nmatch5, pfclusHFHADN_nmatch5, 1.);

  FillHist2D("pfclusHF_nmatch9", pfclusHFEMP_nmatch9, pfclusHFHADP_nmatch9, 1.);
  FillHist2D("pfclusHF_nmatch9", pfclusHFEMN_nmatch9, pfclusHFHADN_nmatch9, 1.);

  FillHist1D("match_ptfrac", pfclusHFP_pttot/genP_pt, 1.);
  FillHist1D("match_ptfrac", pfclusHFN_pttot/genN_pt, 1.);

  FillHist1D("match_ptfrac5", pfclusHFP_pttot5/genP_pt, 1.);
  FillHist1D("match_ptfrac5", pfclusHFN_pttot5/genN_pt, 1.);

  FillHist1D("match_ptfrac9", pfclusHFP_pttot9/genP_pt, 1.);
  FillHist1D("match_ptfrac9", pfclusHFN_pttot9/genN_pt, 1.);

  FillHist1D("pfclusHFEM_Emax",  pfclusHFEMP_Emax,  1. );
  FillHist1D("pfclusHFEM_Emax",  pfclusHFEMN_Emax,  1. );

  FillHist1D("pfclusHFHAD_Emax",  pfclusHFHADP_Emax,  1. );
  FillHist1D("pfclusHFHAD_Emax",  pfclusHFHADN_Emax,  1. );

  FillHist1D("pfclusHF_Emax",  pfclusHFEMP_Emax > pfclusHFHADP_Emax ? pfclusHFEMP_Emax : pfclusHFHADP_Emax ,  1. );
  FillHist1D("pfclusHF_Emax",  pfclusHFEMN_Emax > pfclusHFHADN_Emax ? pfclusHFEMN_Emax : pfclusHFHADN_Emax ,  1. );

  FillHist1D("pfclusHFEM_ptmax",  pfclusHFEMP_ptmax,  1. );
  FillHist1D("pfclusHFEM_ptmax",  pfclusHFEMN_ptmax,  1. );  
  
  FillHist1D("pfclusHFHAD_ptmax",  pfclusHFHADP_ptmax,  1. );
  FillHist1D("pfclusHFHAD_ptmax",  pfclusHFHADN_ptmax,  1. );
  
  FillHist1D("pfclusHF_ptmax",  pfclusHFEMP_ptmax > pfclusHFHADP_ptmax ? pfclusHFEMP_ptmax : pfclusHFHADP_ptmax ,  1. );
  FillHist1D("pfclusHF_ptmax",  pfclusHFEMN_ptmax > pfclusHFHADN_ptmax ? pfclusHFEMN_ptmax : pfclusHFHADN_ptmax ,  1. );


  if (debug_)
  LogPrint("PFTrackHFAnalyzer") << "\n =========== pfcands: =========== " << pfcands->size();
  int pfcandHFEM_n = 0, pfcandHFHAD_n = 0, pfcandHFCH_n = 0 ;
  vector<reco::PFCandidate> matchPFcandP ;
  vector<reco::PFCandidate> matchPFcandN ;

  for(const auto& pfcand : *(pfcands.product()) ){
//  if (pfcand.trackRef().isNonnull())  cout << "All PFCands: pfcand.trackRef().get() "  << pfcand.trackRef().get()  << endl ;
    double eta = pfcand.eta() ;
    double phi = pfcand.phi() ;
    double pt  = pfcand.pt() ;
    double E   = pfcand.energy() ;
    reco::PFCandidate::ParticleType id = pfcand.particleId();

    const reco::PFCandidate::ElementsInBlocks& theElements = pfcand.elementsInBlocks();
    int nblocks = theElements.size() ;
    if(nblocks < 1 ){
      LogPrint("PFTrackHFAnalyzer") << 
      boost::format(" Warning!! PFcandidate with no blocks:  pfcand (pt,E,eta,phi,id,nblocks): (%6.2f, %6.2f, %6.2f, %6.2f, %i, %i)") %
      pt % E % eta % phi % id  % nblocks  ;
      continue ;
    }

    const reco::PFBlockRef blockRef = theElements[0].first;
    const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
    unsigned int nele = elements.size() ;

    bool isHF = false, hasBREM = false ;
    for (unsigned int el = 0; el < nele; el++) {
      reco::PFBlockElement::Type type = elements[el].type();
      if(type == reco::PFBlockElement::HFEM || type == reco::PFBlockElement::HFHAD) isHF = true ; 
      if(type == reco::PFBlockElement::BREM) hasBREM = true ;
    }

    if (debug_){
      if(!isHF)
      LogPrint("PFTrackHFAnalyzer") << boost::format("[ pfcand (pt,E,eta,phi,id,nblocks, nelements): (%6.2f, %6.2f, %6.2f, %6.2f, %i, %i, %i)") %
      pt % E % eta % phi % id  % nblocks % nele ;
      else
      LogPrint("PFTrackHFAnalyzer") << boost::format("pfcand (pt,E,eta,phi,id,nblocks, nelements): (%6.2f, %6.2f, %6.2f, %6.2f, %i, %i, %i)") %
      pt % E % eta % phi % id  % nblocks % nele ;

      for (unsigned int el = 0; el < nele; el++) {
        reco::PFBlockElement::Type type = elements[el].type();
        if(!isHF) 
        LogPrint("PFTrackHFAnalyzer") << boost::format("  element type :  %i ]") % type ;
        else
        LogPrint("PFTrackHFAnalyzer") << boost::format("  element type :  %i ") % type ;
      }
    }

    // now only candidates with HF clusters ....
    if (!isHF) continue ;

    if (id == reco::PFCandidate::egamma_HF ){
      pfcandHFEM_n++ ;
      FillHist1D("pfcandHFEM_pt", pt, 1. );
      FillHist1D("pfcandHFEM_E",  E,  1. );
      FillHist1D("pfcandHFEM_nelements", nele, 1. );
    } else if (id == reco::PFCandidate::h_HF ) {
      pfcandHFHAD_n++ ;
      FillHist1D("pfcandHFHAD_pt", pt, 1. );
      FillHist1D("pfcandHFHAD_E",  E,  1. );
      FillHist1D("pfcandHFHAD_nelements", nele, 1. );
    } else {
      pfcandHFCH_n++ ;
      FillHist1D("pfcandHFCH_pt", pt, 1. );
      FillHist1D("pfcandHFCH_E",  E,  1. );
      FillHist1D("pfcandHFCH_nelements", nele, 1. );
    }

    bool matchCandP = false, matchCandN = false ;
    double ehf = 0., enonhf = 0.  ;
    double tp = -1. ;


    for (unsigned int el = 0; el < nele; el++) {
      reco::PFBlockElement::Type type = elements[el].type();
      if(type == reco::PFBlockElement::TRACK){
        if (debug_)
        LogPrint("PFTrackHFAnalyzer") << boost::format("   track in the block (pt, p, eta, phi, type): (%6.2f, %6.2f, %6.2f, %6.2f, %i)") %
        elements[el].trackRef()->pt() % elements[el].trackRef()->p() % elements[el].trackRef()->eta() % elements[el].trackRef()->phi() % type ;
        tp = elements[el].trackRef()->p() ;
      } else if(type == reco::PFBlockElement::HFEM || type == reco::PFBlockElement::HFHAD){
        if (debug_)
        LogPrint("PFTrackHFAnalyzer") << boost::format("   cluster in the block (pt, e, eta, phi, type): (%6.2f, %6.2f, %6.2f, %6.2f, %i)") %
        elements[el].clusterRef()->pt() % elements[el].clusterRef()->energy() % elements[el].clusterRef()->eta() % elements[el].clusterRef()->phi() % type ;
        const std::vector<reco::PFRecHitFraction> &fracs = elements[el].clusterRef()->recHitFractions();
        const std::vector<std::pair<DetId, float>> &hfracs = elements[el].clusterRef()->hitsAndFractions();
        ehf += elements[el].clusterRef()->energy() ;
        unsigned nhits = fracs.size();
        for(unsigned i=0; i<nhits; i++) {
          const reco::PFRecHitRef& pfRecHits = fracs[i].recHitRef();
          double rawenergy = pfRecHits->energy();
          double frac      = fracs[i].fraction();
          int ieta    = HcalDetId(pfRecHits->detId()).ieta() ;
          int iphi    = HcalDetId(pfRecHits->detId()).iphi() ;
          int idep    = HcalDetId(pfRecHits->detId()).depth() ; 
          if (debug_)
          LogPrint("PFTrackHFAnalyzer") << boost::format("      pfrechit (ieta, iphi, depth, E, frac): (%3d, %3d, %2d, %6.2f, %6.2f)")
          % ieta % iphi % idep % rawenergy % frac ; 

          if( idep == 1 && (abs(ieta - genP_ieta1)<= 1 && abs(idphi(iphi, genP_iphi1)) <=1 ) )   matchCandP = true ;
          if( idep == 2 && (abs(ieta - genP_ieta2)<= 1 && abs(idphi(iphi, genP_iphi2)) <=1 ) )   matchCandP = true ;
          if( idep == 1 && (abs(ieta - genN_ieta1)<= 1 && abs(idphi(iphi, genN_iphi1)) <=1 ) )   matchCandN = true ;
          if( idep == 2 && (abs(ieta - genN_ieta2)<= 1 && abs(idphi(iphi, genN_iphi2)) <=1 ) )   matchCandN = true ;
        }
      } else {
//        enonhf += elements[el].clusterRef()->energy() ;
        enonhf += 0. ;

      }
    }

    if(matchCandP)  matchPFcandP.push_back(pfcand);
    if(matchCandN)  matchPFcandN.push_back(pfcand);
    if( id == reco::PFCandidate::h ){
      if(hasBREM){
        FillHist1D("EovP_brem",     (ehf+enonhf)/tp, 1.);
      } else {
        FillHist1D("EovP_nobrem",   ehf/tp, 1.);
      }
      if(eta > 0.) {
	if(hasBREM) FillHist1D("PovGenE_brem",  tp/genP_E, 1.);
        else        FillHist1D("PovGenE_nobrem",tp/genP_E, 1.); 
      } else {
        if(hasBREM) FillHist1D("PovGenE_brem",  tp/genN_E, 1.);
        else        FillHist1D("PovGenE_nobrem",tp/genN_E, 1.);
      }
    }
  } // end of loop overy PFcandidates ...


  sort(matchPFcandP.begin(),matchPFcandP.end(), mysort);
  sort(matchPFcandN.begin(),matchPFcandN.end(), mysort);
  int nP = matchPFcandP.size() ;
  int nN = matchPFcandN.size() ;

  double etot = 0., etot_nontrk = 0.;
  int pfcandHFHAD_n9 = 0, pfcandHFEM_n9 = 0, pfcandHFCH_n9 = 0 ;
  for (int i=0; i<nP; i++){
    etot        += matchPFcandP[i].energy() ;
    reco::PFCandidate::ParticleType id = matchPFcandP[i].particleId();
    if (id == reco::PFCandidate::h_HF){
      pfcandHFHAD_n9++ ;
      etot_nontrk += matchPFcandP[i].energy() ;
    } else if (id == reco::PFCandidate::egamma_HF ) {
      pfcandHFEM_n9++ ;
      etot_nontrk += matchPFcandP[i].energy() ;
    } else if (id == reco::PFCandidate::h ) {
      pfcandHFCH_n9++ ;
    }
  }
  FillHist1D("pfcandHFEM_n9",  pfcandHFEM_n9, 1.);  
  FillHist1D("pfcandHFHAD_n9", pfcandHFHAD_n9, 1.);
  FillHist1D("pfcandHFCH_n9",  pfcandHFCH_n9, 1.);
  FillHist1D("Eratio_PFcandAlltoGen", etot/genP_E, 1.);
  FillHist1D("Eratio_PFcandNonTtoGen", etot_nontrk/genP_E, 1.);

  etot = 0., etot_nontrk = 0. ;
  pfcandHFHAD_n9 = 0; pfcandHFEM_n9 = 0, pfcandHFCH_n9 = 0 ;
  for (int i=0; i<nN; i++){
    etot += matchPFcandN[i].energy() ;
    reco::PFCandidate::ParticleType id = matchPFcandN[i].particleId();
    if (id == reco::PFCandidate::h_HF){
      pfcandHFHAD_n9++ ;
      etot_nontrk += matchPFcandN[i].energy() ;
    } else if (id == reco::PFCandidate::egamma_HF){
      pfcandHFEM_n9++ ;
      etot_nontrk += matchPFcandN[i].energy() ;
    } else if (id == reco::PFCandidate::h ) {
      pfcandHFCH_n9++ ;
    }
  }  
  FillHist1D("pfcandHFEM_n9",  pfcandHFEM_n9, 1.);
  FillHist1D("pfcandHFHAD_n9", pfcandHFHAD_n9, 1.);
  FillHist1D("pfcandHFCH_n9",  pfcandHFCH_n9, 1.);
  FillHist1D("Eratio_PFcandAlltoGen",  etot/genN_E, 1.);
  FillHist1D("Eratio_PFcandNonTtoGen", etot_nontrk/genN_E, 1.);

  if (nP > 0) {
    double r = matchPFcandP[0].energy()/genP_E ;
    FillHist1D("Eratio_PFcand1toGen", r, 1.);
    if(matchPFcandP[0].particleId() == reco::PFCandidate::egamma_HF) FillHist1D("Eratio_PFcandEM1toGen", r, 1.);
    const reco::PFCandidate::ElementsInBlocks& theElements = matchPFcandP[0].elementsInBlocks();
    const reco::PFBlockRef blockRef = theElements[0].first;
    const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
    if(matchPFcandP[0].particleId() == reco::PFCandidate::h_HF && elements.size() == 1) FillHist1D("Eratio_PFcandHAD1toGen",   r, 1.); 
    if(matchPFcandP[0].particleId() == reco::PFCandidate::h_HF && elements.size() == 2) FillHist1D("Eratio_PFcandEMHAD1toGen", r, 1.);

    if(r< 0.05) cout << " Scan this event # " << nev <<  endl ;
    
  } 

  if (nN > 0) {
    double r = matchPFcandN[0].energy()/genN_E ;
    FillHist1D("Eratio_PFcand1toGen", r, 1.);
    if(matchPFcandN[0].particleId() == reco::PFCandidate::egamma_HF) FillHist1D("Eratio_PFcandEM1toGen", r, 1. );
    const reco::PFCandidate::ElementsInBlocks& theElements = matchPFcandN[0].elementsInBlocks();
    const reco::PFBlockRef blockRef = theElements[0].first;
    const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
    if(matchPFcandN[0].particleId() == reco::PFCandidate::h_HF && elements.size() == 1) FillHist1D("Eratio_PFcandHAD1toGen",  r, 1.);
    if(matchPFcandN[0].particleId() == reco::PFCandidate::h_HF && elements.size() == 2) FillHist1D("Eratio_PFcandEMHAD1toGen",r, 1.); 
  }

  if (nP > 1) FillHist1D("Eratio_PFcand2toGen", (matchPFcandP[0].energy()+matchPFcandP[1].energy())/genP_E, 1. );
  if (nN > 1) FillHist1D("Eratio_PFcand2toGen", (matchPFcandN[0].energy()+matchPFcandN[1].energy())/genN_E, 1. );

  if (nP > 2) FillHist1D("Eratio_PFcand3toGen", (matchPFcandP[0].energy()+matchPFcandP[1].energy()+matchPFcandP[2].energy())/genP_E, 1. );
  if (nN > 2) FillHist1D("Eratio_PFcand3toGen", (matchPFcandN[0].energy()+matchPFcandN[1].energy()+matchPFcandN[2].energy())/genN_E, 1. );

  if (debug_){
    LogPrint("PFTrackHFAnalyzer") << "  pfcandsHF:    "  << pfcandHFEM_n+pfcandHFHAD_n+pfcandHFCH_n;
    LogPrint("PFTrackHFAnalyzer") << "  pfcandsHFEM:  "  << pfcandHFEM_n;
    LogPrint("PFTrackHFAnalyzer") << "  pfcandsHFHAD: "  << pfcandHFHAD_n;
    LogPrint("PFTrackHFAnalyzer") << "  pfcandsHFCH:  "  << pfcandHFCH_n;
  }

  FillHist1D("pfcandHFEM_n",  pfcandHFEM_n, 1.);
  FillHist1D("pfcandHFHAD_n", pfcandHFHAD_n, 1.);
  FillHist1D("pfcandHFCH_n",  pfcandHFCH_n, 1.);

  LogPrint("PFTrackHFAnalyzer") << "\n\n"  ; 

}


// ------------ method called once each job just before starting event loop  ------------
void
PFTrackHFAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
PFTrackHFAnalyzer::endJob()
{
  if (debug_) {
    // LogPrint("PFTrackHFAnalyzer") << "\n";

  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PFTrackHFAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

void PFTrackHFAnalyzer::FillHist1D(const TString& histName, const Double_t& value, const double& weight) 
{
  map<TString, TH1*>::iterator hid=m_Histos1D.find(histName);
  if (hid==m_Histos1D.end())
    cout << "%FillHist1D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value, weight);
}

void PFTrackHFAnalyzer::FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight) 
{
  map<TString, TH2*>::iterator hid=m_Histos2D.find(histName);
  if (hid==m_Histos2D.end())
    cout << "%FillHist2D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, weight);
}


void PFTrackHFAnalyzer::FillProfile2D(const TString& histName, const Double_t& value1, const Double_t& value2, const Double_t& value3, const double& weight) 
{
  map<TString, TProfile2D*>::iterator hid=m_Profiles2D.find(histName);
  if (hid==m_Profiles2D.end())
    cout << "%FillProfile2D -- Could not find profile with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, value3, weight);
}


bool PFTrackHFAnalyzer::IsInCell(double eta, double phi, const CaloCellGeometry::CornersVec& CV)
{

  if (eta < CV[0].eta() && eta > CV[2].eta()) {  
    if (phi < CV[0].phi() && phi > CV[2].phi()){
    return true ; 
    } else if (CV[0].phi() < CV[2].phi()) {
    if ( phi < CV[0].phi()) return true ;
    if ( phi > CV[2].phi()) return true ;
    }
  }

  return false ;
}


int PFTrackHFAnalyzer::idphi(int iphiA, int iphiB)
{ 
  int d = (iphiA - iphiB)/2 ;
  if(d >  18) d = 18 - d ;
  if(d < -18) d = 36 + d ;
  return d ; 
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFTrackHFAnalyzer);
