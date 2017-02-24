// -*- C++ -*-
//
// Package:    Validation/CTPPSFastValidation
// Class:      CTPPSFastValidation
//
/**\class CTPPSFastValidation CTPPSFastValidation.cc Validation/CTPPSFastValidation/plugins/CTPPSFastValidation.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Dilson De Jesus Damiao
//         Created:  Mon, 06 Jan 2017 11:54:27 GMT
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TMath.h>
#include <TStyle.h>
#include <cmath>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/CTPPSReco/interface/CTPPSFastRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSFastRecHitContainer.h"
#include "DataFormats/CTPPSReco/interface/CTPPSFastTrack.h"
#include "DataFormats/CTPPSReco/interface/CTPPSFastTrackContainer.h"
#include <vector>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.


class CTPPSFastValidation : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit CTPPSFastValidation(const edm::ParameterSet&);
  ~CTPPSFastValidation();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  void fillEventInfo(int);
  void store();
  void printcont();
  void format_h2(TH2F *h);
  // ----------member data ---------------------------
  edm::EDGetTokenT< edm::HepMCProduct > mcEventToken; // label of MC event
  edm::Handle< edm::HepMCProduct > EvtHandle ;
  edm::EDGetTokenT<edm::PSimHitContainer> _psimHitToken;
  edm::Handle<edm::PSimHitContainer> PSimHit;
  typedef std::vector<CTPPSFastRecHit> CTPPSFastRecHitContainer;
  edm::EDGetTokenT< CTPPSFastRecHitContainer > _recHitToken;
  edm::Handle<CTPPSFastRecHitContainer> recHits;
  typedef std::vector<CTPPSFastTrack> CTPPSFastTrackContainer;
  edm::EDGetTokenT< CTPPSFastTrackContainer > _tracksPPSToken;
  edm::Handle<CTPPSFastTrackContainer> tracksPPS;
  typedef edm::View<reco::PFJet> JetCollection;
  edm::Handle<JetCollection> jets;
  edm::EDGetTokenT<JetCollection> _jetsToken;
  //const edm::EDGetTokenT<edm::View<reco::PFJet> > _jetsToken;

  std::string _fPhysChannel;
  std::vector<const reco::PFJet*> JetsVector;
  std::vector<double> JetsVector_pt;
  std::vector<double> JetsVector_eta;
  std::vector<double> JetsVector_phi;
  double   Rjj = -999.0;
  double S = 13000.0;
  double MxPPS = -999.0;
  double Mjj = -999.0;
  static const int NMCPMAX = 100000;
  int count1 = 0;
  int count2 = 0;
  int EventKind, evNumber= 0;
  int contTracks = 0, contTracksP = 0, contTracksN = 0;
  int contrechits = 0, contsimhits = 0;
  float protonS1_p[NMCPMAX],protonS1_pt[NMCPMAX],protonS1_eta[NMCPMAX],protonS1_phi[NMCPMAX];
  double protonS1_t[NMCPMAX],protonS1_xi[NMCPMAX], vx[NMCPMAX],vy[NMCPMAX], vz[NMCPMAX];

  const double   ProtonMass = 0.93827;
  const double ProtonMassSQ = pow(ProtonMass,2);
  double fBeamEnergy;
  double fBeamMomentum;
  void Get_t_and_xi(const TLorentzVector* p,double& t, double& xi);
  void set_BeamEnergy(double e) {fBeamEnergy=e;fBeamMomentum = sqrt(fBeamEnergy*fBeamEnergy - ProtonMassSQ);};
  //TFileService
  edm::Service<TFileService> fs;
  TTree* AnalysisTree;

  // Proton Histograms
  TH2F *h_LHCT_xVsy_ARMFs1, *h_LHCT_xVsy_ARMFs2, *h_LHCT_xVsy_ARMBs1, *h_LHCT_xVsy_ARMBs2, *h_LHCT_xVsy_s2;
  TH2F *h_LHCT_xiVst_s2;
  TH2F *hxy_ARMB_simdet1, *hxy_ARMB_simdet2, *hxy_ARMB_simtof;
  TH2F *hxy_ARMF_simdet1, *hxy_ARMF_simdet2, *hxy_ARMF_simtof;
  TH2F *hxy_ARMB_recdet1, *hxy_ARMB_recdet2, *hxy_ARMB_rectof;
  TH2F *hxy_ARMF_recdet1, *hxy_ARMF_recdet2, *hxy_ARMF_rectof;
  TH2F *hxi_t_ARMF_trk, *hxi_t_ARMB_trk;
  TH2F *ht_ARMF_genreco, *ht_ARMB_genreco, *hxi_ARMF_genreco, *hxi_ARMB_genreco;
  TH2F *hMxPPSvsMjj;
  TH1F *h_LHCT_p_s2, *h_LHCT_pt_s2, *h_LHCT_pz_s2, *h_LHCT_py_s2, *h_LHCT_px_s2;
  TH1F *h_LHCT_vx_s2, *h_LHCT_vy_s2, *h_LHCT_vz_s2;
  TH1F *h_LHCT_eta_ARMF, *h_LHCT_eta_ARMB, *h_LHCT_phi_ARMF, *h_LHCT_phi_ARMB;
  TH1F *h_LHCT_px_ARMF, *h_LHCT_px_ARMB;
  TH1F *h_LHCT_py_ARMF, *h_LHCT_py_ARMB;
  TH1F *h_LHCT_pz_ARMF, *h_LHCT_pz_ARMB;
  TH1F *h_LHCT_pt_ARMF, *h_LHCT_pt_ARMB, *h_LHCT_p_ARMF, *h_LHCT_p_ARMB;
  TH1F *h_LHCT_xVertices, *h_LHCT_yVertices, *h_LHCT_zVertices;
  TH1F *h_LHCT_xi_s2, *h_LHCT_t_s2;
  TH1F *ht_ARMF_trk, *ht_ARMB_trk, *hxi_ARMF_trk, *hxi_ARMB_trk;
  TH1F *hx0_trk, *hy0_trk, *hip;
  TH1F * h_jet_eta, * h_jet_phi, * h_jet_pt;
  TH1F *hMjjdivMxPPS;
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
CTPPSFastValidation::CTPPSFastValidation(const edm::ParameterSet& iConfig)
//_jetsToken( consumes<edm::View<reco::PFJet> > (iConfig.getParameter<edm::InputTag>("jetsTag") ) )
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  mcEventToken    = mayConsume<edm::HepMCProduct>(iConfig.getUntrackedParameter<edm::InputTag>("MCEvent",std::string("")));
  _psimHitToken   = mayConsume<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("psimHitTag"));
  _recHitToken    = consumes<CTPPSFastRecHitContainer>(iConfig.getParameter<edm::InputTag>("recHitTag"));
  _tracksPPSToken = consumes<CTPPSFastTrackContainer>(iConfig.getParameter<edm::InputTag>("tracksPPSTag"));
  _jetsToken      = consumes<JetCollection>(iConfig.getParameter<edm::InputTag>("jetsTag"));
  _fPhysChannel   =iConfig.getParameter<std::string>("fPhysChannelTag");

  //evNumber=0;
}


CTPPSFastValidation::~CTPPSFastValidation()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

class GreaterPt{
public:
  bool operator()( const math::XYZTLorentzVector& a, const math::XYZTLorentzVector& b) {
    return a.pt() > b.pt();
  }
};

//
// member functions
//

// ------------ method called for each event  ------------

void CTPPSFastValidation::Get_t_and_xi(const TLorentzVector* proton,double& t,double& xi) {
  set_BeamEnergy(13000/2.);
  t = 0.;
  xi = -1.;
  if (!proton) return;
  double mom    = proton->P();
  if (mom>fBeamMomentum) mom=fBeamMomentum;
  double energy = proton->E();
  double theta  = (proton->Pz()>0)?proton->Theta():TMath::Pi()-proton->Theta();
  t      = -2.*(ProtonMassSQ-fBeamEnergy*energy+fBeamMomentum*mom*cos(theta));
  xi     = (1.0-energy/fBeamEnergy);
}

void CTPPSFastValidation::store(){
  //AnalysisTree->Fill();
}
void CTPPSFastValidation::printcont(){
  std::cout << "contTracks = " << contTracks << std::endl;
  std::cout << "contTracksP = " << contTracksP << std::endl;
  std::cout << "contTracksN = " << contTracksN << std::endl;
  std::cout << "contsimhits = " << contsimhits << std::endl;
  std::cout << "contrechits = " << contrechits << std::endl;
  std::cout << "store " << count2++ << std::endl;
}

void CTPPSFastValidation::fillEventInfo(int e){
  EventKind = e;
}

void CTPPSFastValidation::format_h2(TH2F* h){
    h->SetOption("COLZ");
    h->SetStats(kFALSE);
}

void
CTPPSFastValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  iEvent.getByToken( mcEventToken, EvtHandle );
  iEvent.getByToken( _psimHitToken, PSimHit);
  edm::PSimHitContainer::const_iterator itPSimHit;
  iEvent.getByToken( _recHitToken, recHits);
  edm::CTPPSFastRecHitContainer::const_iterator itRecHit;
  iEvent.getByToken( _tracksPPSToken, tracksPPS);
  edm::CTPPSFastTrackContainer::const_iterator itTracksPPS;

  Handle<edm::View<reco::PFJet> > jets;
  iEvent.getByToken(_jetsToken, jets);
  typename JetCollection::const_iterator jet;
  //std::cout<< jets.isValid() << std::endl;

  //math::XYZTLorentzVector p4jet[2];
  math::XYZTLorentzVector dijetSystem(0.,0.,0.,0.);
  JetsVector.clear();
  JetsVector_pt.clear();
  JetsVector_eta.clear();
  JetsVector_phi.clear();
  if(jets.isValid()) {
    //const reco::PFJet* jetAll = &(*jets);
    int ijet = 0;
    for(jet = jets->begin(); jet != jets->end(); jet++) {
      h_jet_eta->Fill(jet->eta());
      h_jet_phi->Fill(jet->phi());
      h_jet_pt->Fill(jet->pt());
      JetsVector_pt.push_back(jet->pt());
      JetsVector_eta.push_back(jet->eta());
      JetsVector_phi.push_back(jet->phi());
      if (ijet<2) dijetSystem += jet->p4();
      ijet++;
    }
    //   std::cout << dijetSystem.M() << std::endl;
    Mjj = dijetSystem.M();
  }

  const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
  EventKind = Evt->signal_process_id();
  std::vector<math::XYZTLorentzVector> protonCTPPS;
  protonCTPPS.clear();
  //-----------------------------------------------
  for(HepMC::GenEvent::particle_const_iterator i=Evt->particles_begin(); i != Evt->particles_end();++i) {
    int myId = (*i)->pdg_id();
    if(myId!=2212) continue;

    HepMC::GenVertex* pv = (*i)->production_vertex();
    HepMC::FourVector vertex = pv->position();
    HepMC::FourVector momentum=(*i)->momentum();
    const HepMC::FourVector p((*i)->momentum());
    protonCTPPS.push_back(math::XYZTLorentzVector(p.x(),p.y(),p.z(),p.t()));
    double t,xi;
    double px = momentum.x();
    double py = momentum.y();
    double pz = momentum.z();
    double e = sqrt(px*px+py*py+pz*pz+ProtonMassSQ);
    TLorentzVector* proton = new TLorentzVector(px,py,pz,e);
    Get_t_and_xi(proton,t,xi);
    //status() == 1 -> vertex at detector entrance
    if(vertex.eta() >= 8. && (*i)->status() == 1) { //Arm F
      h_LHCT_pt_ARMF->Fill(proton->Pt());
      h_LHCT_p_ARMF->Fill(proton->P());
      h_LHCT_px_ARMF->Fill(proton->Px());
      h_LHCT_py_ARMF->Fill(proton->Py());
      h_LHCT_pz_ARMF->Fill(proton->Pz());
      h_LHCT_eta_ARMF->Fill(proton->Eta());
      h_LHCT_phi_ARMF->Fill(proton->Phi());
      h_LHCT_xVsy_ARMFs1->Fill(vertex.x(),vertex.y());
    }
    if(vertex.eta() <= -8. && (*i)->status() == 1) { //Arm B
      h_LHCT_pt_ARMB->Fill(proton->Pt());
      h_LHCT_p_ARMB->Fill(proton->P());
      h_LHCT_px_ARMB->Fill(proton->Px());
      h_LHCT_py_ARMB->Fill(proton->Py());
      h_LHCT_pz_ARMB->Fill(proton->Pz());
      h_LHCT_eta_ARMB->Fill(proton->Eta());
      h_LHCT_phi_ARMB->Fill(proton->Phi());
      h_LHCT_xVsy_ARMBs1->Fill(vertex.x(),vertex.y());
    }
    //status() == 2 -> vertex at IP
    if((*i)->status() == 2) {
      Get_t_and_xi(proton,t,xi);
      h_LHCT_xi_s2->Fill(xi);
      h_LHCT_t_s2->Fill(t);
      h_LHCT_xiVst_s2->Fill(xi,t);
      h_LHCT_p_s2->Fill(proton->P());
      h_LHCT_pt_s2->Fill(proton->Pt());
      h_LHCT_pz_s2->Fill(proton->Pz());
      h_LHCT_py_s2->Fill(proton->Py());
      h_LHCT_px_s2->Fill(proton->Px());
      h_LHCT_vx_s2->Fill(vertex.x());
      h_LHCT_vy_s2->Fill(vertex.y());
      //std::cout << vertex.x() << "  "  << vertex.y() << std::endl;
      //std::cout << vertex.z() << std::endl;
      h_LHCT_vz_s2->Fill(vertex.z());
      double xiPPSArmF = -999., xiPPSArmB = -999.;
      if(tracksPPS.isValid()) {
        for (itTracksPPS = tracksPPS->begin(); itTracksPPS != tracksPPS->end(); ++itTracksPPS) {
          if(itTracksPPS->PZ()>0&&proton->Pz()>0) {
            contTracksP++;
            hx0_trk->Fill(itTracksPPS->X0());
            hy0_trk->Fill(itTracksPPS->Y0());
            hip->Fill(sqrt(itTracksPPS->X0()*itTracksPPS->X0()+itTracksPPS->Y0()*itTracksPPS->Y0()));
            ht_ARMF_trk->Fill(itTracksPPS->t());
            hxi_ARMF_trk->Fill(itTracksPPS->xi());
            xiPPSArmF = itTracksPPS->xi();
            hxi_t_ARMF_trk->Fill(itTracksPPS->xi(),itTracksPPS->t());
            ht_ARMF_genreco->Fill(t,itTracksPPS->t());
            hxi_ARMF_genreco->Fill(xi,itTracksPPS->xi());
          }
          if(itTracksPPS->PZ()<0&&proton->Pz()<0) {
            contTracksN++;
            hx0_trk->Fill(itTracksPPS->X0());
            hy0_trk->Fill(itTracksPPS->Y0());
            hip->Fill(sqrt(itTracksPPS->X0()*itTracksPPS->X0()+itTracksPPS->Y0()*itTracksPPS->Y0()));
            ht_ARMB_trk->Fill(itTracksPPS->t());
            hxi_ARMB_trk->Fill(itTracksPPS->xi());
            xiPPSArmB = itTracksPPS->xi();
            hxi_t_ARMB_trk->Fill(itTracksPPS->xi(),itTracksPPS->t());
            ht_ARMB_genreco->Fill(t,itTracksPPS->t());
            hxi_ARMB_genreco->Fill(xi,itTracksPPS->xi());
          }
        }
      }
      MxPPS = S*sqrt(xiPPSArmF*xiPPSArmB);
      if(xiPPSArmF > 0 && xiPPSArmB > 0 ){
        //std::cout  << JetsVector_pt.size() << " pt1 " << JetsVector_pt[0] << " pt0 " << JetsVector_pt[1] << " Mx  " << MxPPS << " Mjj " << Mjj << "Rjj  " << Mjj/MxPPS << " "  <<xiPPSArmF << " " << xiPPSArmB <<   std::endl;
        hMxPPSvsMjj->Fill(MxPPS,Mjj);
        hMjjdivMxPPS->Fill(Mjj/MxPPS);
      }
    }
  }
  if(PSimHit.isValid()) {
    for (itPSimHit = PSimHit->begin(); itPSimHit != PSimHit->end(); ++itPSimHit) {
      contsimhits++;
      if(itPSimHit->detUnitId() == 2031091712) hxy_ARMB_simdet1->Fill(itPSimHit->entryPoint().x(),itPSimHit->entryPoint().y());
      if(itPSimHit->detUnitId() == 2031616000) hxy_ARMB_simdet2->Fill(itPSimHit->entryPoint().x(),itPSimHit->entryPoint().y());
      if(itPSimHit->detUnitId() == 2063597568) hxy_ARMB_simtof->Fill(itPSimHit->entryPoint().x(),itPSimHit->entryPoint().y());
      if(itPSimHit->detUnitId() == 2014314496) hxy_ARMF_simdet1->Fill(itPSimHit->entryPoint().x(),itPSimHit->entryPoint().y());
      if(itPSimHit->detUnitId() == 2014838784) hxy_ARMF_simdet2->Fill(itPSimHit->entryPoint().x(),itPSimHit->entryPoint().y());
      if(itPSimHit->detUnitId() == 2046820352) hxy_ARMF_simtof->Fill(itPSimHit->entryPoint().x(),itPSimHit->entryPoint().y());
    }
  }
  if(recHits.isValid()) {
    for (itRecHit = recHits->begin(); itRecHit != recHits->end(); ++itRecHit) {
      contrechits++;
      if(itRecHit->detUnitId() == 2031091712) hxy_ARMB_recdet1->Fill(itRecHit->entryPoint().x(),itRecHit->entryPoint().y());
      if(itRecHit->detUnitId() == 2031616000) hxy_ARMB_recdet2->Fill(itRecHit->entryPoint().x(),itRecHit->entryPoint().y());
      if(itRecHit->detUnitId() == 2063597568) hxy_ARMB_rectof->Fill(itRecHit->entryPoint().x(),itRecHit->entryPoint().y());
      if(itRecHit->detUnitId() == 2014314496) hxy_ARMF_recdet1->Fill(itRecHit->entryPoint().x(),itRecHit->entryPoint().y());
      if(itRecHit->detUnitId() == 2014838784) hxy_ARMF_recdet2->Fill(itRecHit->entryPoint().x(),itRecHit->entryPoint().y());
      if(itRecHit->detUnitId() == 2046820352) hxy_ARMF_rectof->Fill(itRecHit->entryPoint().x(),itRecHit->entryPoint().y());
    }
  }  //store();
  //contTracks  = contTracksP+contTracksN;
  //printcont();
}


// ------------ method called once each job just before starting event loop  ------------
void
CTPPSFastValidation::beginJob()
{
  std::cout << _fPhysChannel << "  ---------- " << std::endl;
  // use TFileService for output to root file
  AnalysisTree = fs->make<TTree>("AnalysisTree","CTPPS Analysis Tree ");
  AnalysisTree->Branch("EventKind",&EventKind,"EventKind/I");
  // GenParticles at hadron level
  AnalysisTree->Branch("evNumber",&evNumber,"evNumber/I");
  AnalysisTree->Branch("protonS1_p",protonS1_p,"protonS1_p[evNumber]/F");
  AnalysisTree->Branch("protonS1_pt",protonS1_pt,"protonS1_pt[evNumber]/F");
  AnalysisTree->Branch("protonS1_eta",protonS1_eta,"protonS1_eta[evNumber]/F");
  AnalysisTree->Branch("protonS1_phi",protonS1_phi,"protonS1_phi[evNumber]/F");
  AnalysisTree->Branch("protonS1_t",protonS1_t,"protonS1_t[evNumber]/D");
  AnalysisTree->Branch("protonS1_xi",protonS1_xi,"protonS1_xi[evNumber]/D");
  AnalysisTree->Branch("vx",vx,"vx[evNumber]/D");
  AnalysisTree->Branch("vy",vy,"vx[evNumber]/D");
  AnalysisTree->Branch("vz",vz,"vx[evNumber]/D");

  h_LHCT_pt_ARMF = fs->make<TH1F>( "LHCT_pt_armF" , "LHCT_ptARMF; pt [GeV] ; nEvents", 100, .0, 2.0 );
  h_LHCT_px_ARMF = fs->make<TH1F>( "LHCT_px_armF" , "LHCT_pxARMF; px [GeV] ; nEvents", 100, -5.0, 5.0 );
  h_LHCT_py_ARMF = fs->make<TH1F>( "LHCT_py_armF" , "LHCT_pyARMF; py [GeV] ; nEvents", 100, -5.0, 5.0 );
  h_LHCT_pz_ARMF = fs->make<TH1F>( "LHCT_pz_armF" , "LHCT_pzARMF; pz [GeV] ; nEvents", 300, 5000.0, 6500.0 );
  h_LHCT_p_ARMF = fs->make<TH1F>( "LHCT_p_armF" , "LHCT_pARMF; p [GeV] ; nEvents", 300, 6500.0, 5000.0 );
  h_LHCT_eta_ARMF = fs->make<TH1F>( "LHCT_etaARMF" , "LHCT_etaARMF; #eta; nEvents" , 100, 8., 15.);
  h_LHCT_phi_ARMF = fs->make<TH1F>( "LHCT_phiARMF" , "LHCT_phiARMF; #phi; nEvents" , 100, -3.2, 3.2 );
  h_LHCT_xVsy_ARMFs1 = fs->make<TH2F>( "LHCT_xVsyARMFs1" , "LHCT_xVsyARMFs1; x [mm] ; y [mm]" , 100,-30,0.,100,-15,15);
  h_LHCT_pt_ARMB = fs->make<TH1F>( "LHCT_pt_armB" , "LHCT_ptARMB; pt [GeV] ; nEvents", 100, .0, 2.0 );
  h_LHCT_px_ARMB = fs->make<TH1F>( "LHCT_px_armB" , "LHCT_pxARMB; px [GeV] ; nEvents", 100, -5.0, 5.0 );
  h_LHCT_py_ARMB = fs->make<TH1F>( "LHCT_py_armB" , "LHCT_pyARMB; py [GeV] ; nEvents", 100, -5.0, 5.0 );
  h_LHCT_p_ARMB = fs->make<TH1F>( "LHCT_p_armB" , "LHCT_pARMB; p [GeV] ; nEvents", 300, 5000.0, 6500.0 );
  h_LHCT_pz_ARMB = fs->make<TH1F>( "LHCT_pz_armB" , "LHCT_pzARMB; p [GeV] ; nEvents", 300, -6500.0, -5000.0 );
  h_LHCT_eta_ARMB  = fs->make<TH1F>( "LHCT_etaARMB" , "LHCT_etaARMB; #eta; nEvents" , 100, -15.,-8.);
  h_LHCT_phi_ARMB = fs->make<TH1F>( "LHCT_phiARMB" , "LHCT_phiARMB; #phi; nEvents" , 100, -3.2, 3.2 );
  h_LHCT_xVsy_ARMBs1 = fs->make<TH2F>( "LHCT_xVsyARMBs1" , "LHCT_xVsyARMBs1; x [mm] ; y [mm]" ,100,-30,0.,100,-15,15);
  h_LHCT_xi_s2 = fs->make<TH1F>( "LHCT_xi_s2" , "LHCT_xis2; #xi ; nEvents", 100, 0., .21 );
  h_LHCT_t_s2  = fs->make<TH1F>( "LHCT_t_s2" , "LHCT_ts2; t [GeV^{2}] ; nEvents", 100, .0, 2.0 );
  h_LHCT_xiVst_s2 = fs->make<TH2F>( "LHCT_xiVst_s2" , "LHCT_xiVsts2; #xi ; t [GeV^{2}]" , 100, 0., .21 , 100, .0, 2.0 );
  h_LHCT_pt_s2 = fs->make<TH1F>( "LHCT_pt_s2" , "LHCT_pts2; pt [GeV] ; nEvents", 100, .0, 2.0);
  h_LHCT_pz_s2 = fs->make<TH1F>( "LHCT_pz_s2" , "LHCT_pzs2; p [GeV] ; nEvents", 1000, 6500.0, 6500.0 );
  h_LHCT_py_s2 = fs->make<TH1F>( "LHCT_py_s2" , "LHCT_pys2; p [GeV] ; nEvents", 1000, -2.0, 2.0 );
  h_LHCT_px_s2 = fs->make<TH1F>( "LHCT_px_s2" , "LHCT_pxs2; p [GeV] ; nEvents", 1000, -3.0, 3.0 );
  h_LHCT_p_s2 = fs->make<TH1F>( "LHCT_p_s2" , "LHCT_ps2; p [GeV] ; nEvents", 1000, 5200.0, 6500.0 );
  h_LHCT_vx_s2 = fs->make<TH1F>( "LHCT_x_s2" , "LHCT_xs2; x [mm] ; nEvents", 100, -1., 1. );
  h_LHCT_vy_s2 = fs->make<TH1F>( "LHCT_y_s2" , "LHCT_ys2; y [mm] ; nEvents", 100, -2., 2. );
  h_LHCT_vz_s2 = fs->make<TH1F>( "LHCT_z_s2" , "LHCT_zs2; z [mm] ; nEvents", 1000, -200.0, 200.0 );
  const int nXbins = 8;
  double xbin[nXbins+1] = {-18.195, -13.995, -9.795, -7.445, -5.695, -4.535, -3.515, -2.505, -1.695};
  double ybin = 4.2;
  float zdet1 = 203.8, zdet2 = 212.55, ztof = 215.7;
  hxy_ARMB_simdet1  = fs->make<TH2F>("hxy_ARMB_simdet1", Form("%s, SIM, Z = -%f, d = 15 #sigma;x [mm] ; y [mm]",_fPhysChannel.c_str(),zdet1), 100,-30,0.,100,-15,15);
  hxy_ARMB_simdet2  = fs->make<TH2F>("hxy_ARMB_simdet2", Form("%s, SIM, Z = -212.55 m, d = 15 #sigma;x [mm] ; y [mm]",_fPhysChannel.c_str()), 100,-30,0.,100,-15,15);
  hxy_ARMB_simtof  = fs->make<TH2F>("hxy_ARMB_simtof", Form("%s, SIM, Z = -215.7 m, d = 15 #sigma;x [mm] ; y [mm]",_fPhysChannel.c_str()),100,-30,0.,100,-15,15);
  hxy_ARMF_simdet1  = fs->make<TH2F>("hxy_ARMF_simdet1", Form("%s, SIM, Z = 203.8 m, d = 15 #sigma;x [mm] ; y [mm]",_fPhysChannel.c_str()), 100,-30,0.,100,-15,15);
  hxy_ARMF_simdet2  = fs->make<TH2F>("hxy_ARMF_simdet2", Form("%s, SIM, Z = 212.55 m, d = 15 #sigma;x [mm] ; y [mm]",_fPhysChannel.c_str()), 100,-30,0.,100,-15,15);
  hxy_ARMF_simtof  = fs->make<TH2F>("hxy_ARMF_simtof", Form("%s, SIM, Z = 215.7 m, d = 15 #sigma;x [mm] ; y [mm]",_fPhysChannel.c_str()),100,-30,0.,100,-15,15);
  hxy_ARMB_recdet1  = fs->make<TH2F>("hxy_ARMB_recdet1", Form("%s, RECO, Z = -203.8 m, d = 15 #sigma;x [mm] ; y [mm]",_fPhysChannel.c_str()) , 100,-30,0.,100,-15,15);
  hxy_ARMB_recdet2  = fs->make<TH2F>("hxy_ARMB_recdet2", Form("%s, RECO, Z = -212.55 m, d = 15 #sigma;x [mm] ; y [mm]",_fPhysChannel.c_str()) , 100,-30,0.,100,-15,15);
  hxy_ARMF_recdet1  = fs->make<TH2F>("hxy_ARMF_recdet1", Form("%s, RECO, Z = 203.8 m, d = 15 #sigma;x [mm] ; y [mm]",_fPhysChannel.c_str()), 100,-30,0.,100,-15,15);
  hxy_ARMF_recdet2  = fs->make<TH2F>("hxy_ARMF_recdet2",  Form("%s, RECO, Z = 212.55 m, d = 15 #sigma;x [mm] ; y [mm]",_fPhysChannel.c_str()), 100,-30,0.,100,-15,15);
  hxy_ARMF_rectof  = fs->make<TH2F>("hxy_ARMF_rectof",Form("%s, RECO, Z = 215.7 m, d = 15 #sigma;x [mm] ; y [mm]",_fPhysChannel.c_str()) ,nXbins+2, xbin,3,-1.5*ybin,1.5*ybin);
  hxy_ARMB_rectof  = fs->make<TH2F>("hxy_ARMB_rectof",Form("%s, RECO, Z = -215.7 m, d = 15 #sigma;x [mm] ; y [mm]",_fPhysChannel.c_str()) ,nXbins+2, xbin,3,-1.5*ybin,1.5*ybin);
  // Z > 0 -> ARMF && Z < 0 -> ARMB
  ht_ARMF_trk = fs->make<TH1F>( "ht_ARMF_trk" , Form("%s, CTPPS Tracks, Z > 0; |t_{reco}| [GeV^{2}] ; nEvents",_fPhysChannel.c_str()),100, .0, 2.0 );
  ht_ARMB_trk = fs->make<TH1F>( "ht_ARMB_trk" , Form("%s, CTPPS Tracks, Z < 0; |t_{reco}| [GeV^{2}] ; nEvents",_fPhysChannel.c_str()),100, .0, 2.0 );
  hxi_ARMF_trk= fs->make<TH1F>( "hxi_ARMF_trk" , Form("%s, CTPPS Tracks, Z > 0; #xi_{reco} ; nEvents",_fPhysChannel.c_str()),100, 0., .21 );
  hxi_ARMB_trk= fs->make<TH1F>( "hxi_ARMB_trk" , Form("%s, CTPPS Tracks, Z < 0; #xi_{reco} ; nEvents",_fPhysChannel.c_str()),100, 0., .21 );
  hxi_t_ARMF_trk = fs->make<TH2F>( "hxi_t_ARMF_trk" ,Form("%s, CTPPS Tracks, Z > 0;#xi_{reco} ; |t_{reco}| [GeV^{2}]",_fPhysChannel.c_str()),100, 0., .21 , 100, .0, 2.0 );
  hxi_t_ARMB_trk = fs->make<TH2F>( "hxi_t_ARMB_trk" ,Form("%s, CTPPS Tracks, Z < 0;#xi_{reco} ; |t_{reco}| [GeV^{2}]",_fPhysChannel.c_str()),100, 0., .21 , 100, .0, 2.0 );
  hx0_trk = fs->make<TH1F>( "hx0_trk" , "hx0_trk; X0 (mm) ; nEvents", 100, -1.,1.);
  hy0_trk = fs->make<TH1F>( "hy0_trk" , "hy0_trk; Y0 (mm) ; nEvents", 100, -1.,1.);
  hip = fs->make<TH1F>( "hip" , "hip; Impact parameter (mm) ; nEvents", 100, 0.,1.);
  hxi_ARMF_genreco  = fs->make<TH2F>("hxi_ARMF_genreco" , Form("%s, Proton Tracks, Z > 0; #xi_{gen} ; #xi_{reco} ",_fPhysChannel.c_str()), 100, 0., .21 ,  100, 0., .21 );
  hxi_ARMB_genreco  = fs->make<TH2F>("hxi_ARMB_genreco" , Form("%s, Proton Tracks, Z < 0; #xi_{gen} ; #xi_{reco} ",_fPhysChannel.c_str()), 100, 0., .21 ,  100, 0., .21 );
  ht_ARMF_genreco  = fs->make<TH2F>("ht_ARMF_genreco" , Form("%s, Proton Tracks, Z > 0; |t|_{gen} [GeV^{2}] ; |t|_{reco} [GeV^{2}] ",_fPhysChannel.c_str()), 100, .0, 5.0 ,100, .0, 5.0 );
  ht_ARMB_genreco  = fs->make<TH2F>("ht_ARMB_genreco" , Form("%s, Proton Tracks, Z < 0; |t|_{gen} [GeV^{2}] ; |t|_{reco} [GeV^{2}] ",_fPhysChannel.c_str()), 100, .0, 5.0 ,100, .0, 5.0 );
  h_jet_eta          = fs->make<TH1F>( "h_jet_eta",Form("%s, Jets #eta; #eta; nJets",_fPhysChannel.c_str()), 100,  -6.,    6.0 );
  h_jet_phi          = fs->make<TH1F>( "h_jet_phi", Form("%s, Jets #phi; #phi; nJets",_fPhysChannel.c_str()), 100, -3.2,    3.2 );
  h_jet_pt           = fs->make<TH1F>( "h_jet_pt", Form("%s, Jets p_{t}; p_{t} [GeV]; nJets",_fPhysChannel.c_str()), 1000,  0.0, 1000. );
  hMxPPSvsMjj = fs->make<TH2F>("hMxvsMjj" , "Mass_{X} vs M_{JJ}  distribution; M_{x}  [GeV];  M_{jj}  [GeV]" , 100, 0., 2000.,100, 0., 2000. );
  hMjjdivMxPPS = fs->make<TH1F>("hMjjdivMxPPS" , "R = M_{jj}/M_{X}  distribution; M_{jj}/M_{X}; Entries " , 100, 0.01, 2. );
  for (TH2F *hist : {h_LHCT_xiVst_s2, hMxPPSvsMjj,
    h_LHCT_xVsy_ARMBs1,  hxy_ARMB_simdet1, hxy_ARMB_simdet2, hxy_ARMB_simtof, hxy_ARMB_recdet1, hxy_ARMB_recdet2, hxy_ARMB_rectof,
      hxi_t_ARMB_trk, ht_ARMB_genreco, hxi_ARMB_genreco,
    h_LHCT_xVsy_ARMFs1,  hxy_ARMF_simdet1, hxy_ARMF_simdet2, hxy_ARMF_simtof, hxy_ARMF_recdet1, hxy_ARMF_recdet2, hxy_ARMF_rectof,
      hxi_t_ARMF_trk, ht_ARMF_genreco, hxi_ARMF_genreco})
    {format_h2(hist);}
    gStyle->SetPalette(1);
}

// ------------ method called once each job just after ending the event loop  ------------
void
CTPPSFastValidation::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CTPPSFastValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CTPPSFastValidation);
