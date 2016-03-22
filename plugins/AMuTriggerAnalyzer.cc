// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/AMuTriggerAnalyzer
// Class:      AMuTriggerAnalyzer
// 
/**\class AMuTriggerAnalyzer AMuTriggerAnalyzer.cc GGHAA2Mu2TauAnalysis/AMuTriggerAnalyzer/plugins/AMuTriggerAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mengyao Shi
//         Created:  Mon, 21 Mar 2016 12:00:57 GMT
//
//


// system include files
#include <memory>
#include <cmath>
// user include files
#include "TH1D.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "Tools/Common/interface/GenTauDecayID.h"
#include "Tools/Common/interface/Common.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class AMuTriggerAnalyzer : public edm::EDAnalyzer{
   public:
      explicit AMuTriggerAnalyzer(const edm::ParameterSet&);
      ~AMuTriggerAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> GenMatchedRecoMuonTag_;
      edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> MuPassTrigger_;
      edm::EDGetTokenT<reco::GenParticleCollection> genParticleTag_;
      std::map<std::string, TH1D*> histos1D_;
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
AMuTriggerAnalyzer::AMuTriggerAnalyzer(const edm::ParameterSet& iConfig):
  GenMatchedRecoMuonTag_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("GenMatchedRecoMuonTag"))),
  MuPassTrigger_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("MuPassTrigger"))),
  genParticleTag_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
  histos1D_()
{

}


AMuTriggerAnalyzer::~AMuTriggerAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
AMuTriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



   edm::Handle<edm::RefVector<std::vector<reco::Muon>>> MuonsPassTrigger;
   iEvent.getByToken(MuPassTrigger_,MuonsPassTrigger);

   edm::Handle<edm::RefVector<std::vector<reco::Muon>>> GenMatchedMuons;
   iEvent.getByToken(GenMatchedRecoMuonTag_,GenMatchedMuons);   

   edm::Handle<reco::GenParticleCollection> pGenParticles;
   iEvent.getByToken(genParticleTag_, pGenParticles);

   std::vector<reco::GenParticle*> genObjPtrs; 
   for(typename std::vector<reco::GenParticle>::const_iterator iGenObj=pGenParticles->begin(); iGenObj!=pGenParticles->end(); ++iGenObj)
   {
     const unsigned int absPDGID=fabs(iGenObj->pdgId());
     if(absPDGID==13)
       genObjPtrs.push_back(const_cast<reco::GenParticle*>(&(*iGenObj)));
   }
   for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iMuonPassTrigger=MuonsPassTrigger->begin();iMuonPassTrigger!=MuonsPassTrigger->end(); ++iMuonPassTrigger)
   {
     int nearestGenObjKey0=-1;
     const reco::GenParticle* nearestGenObj0=
       Common::nearestObject(*iMuonPassTrigger, genObjPtrs, nearestGenObjKey0);
     histos1D_["numerator"]->Fill(reco::deltaR((*iMuonPassTrigger)->eta(), (*iMuonPassTrigger)->phi(), nearestGenObj0->eta(), nearestGenObj0->phi())); 
   }
   for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iGenMatchedMuon=GenMatchedMuons->begin(); iGenMatchedMuon!=GenMatchedMuons->end(); ++iGenMatchedMuon)
   {
     int nearestGenObjKey=-1;
     const reco::GenParticle* nearestGenObj=
      Common::nearestObject(*iGenMatchedMuon, genObjPtrs,nearestGenObjKey);
     histos1D_["denominator"]->Fill(reco::deltaR((*iGenMatchedMuon)->eta(),(*iGenMatchedMuon)->phi(), nearestGenObj->eta(), nearestGenObj->phi()));
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
AMuTriggerAnalyzer::beginJob()
{
  edm::Service<TFileService> fileService;
  histos1D_["denominator"]=fileService->make<TH1D>("denominator", "dR Of RecoMuonWithLooseID45GeV2.1etaCutGenMatchedToAMu vs. nearestGenMuon ", 100, 0.0, 0.05);
  histos1D_["numerator"]=fileService->make<TH1D>("numerator","Denominator*FiredHLT*TriggerMatch(PT match)",100, 0.0, 0.05);
  histos1D_["Efficiency"]=fileService->make<TH1D>("Efficiency","Efficiency of hlt+match", 100,0.0,0.05);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AMuTriggerAnalyzer::endJob() 
{
  histos1D_["Efficiency"]->Divide(histos1D_["numerator"], histos1D_["denominator"]);
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AMuTriggerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AMuTriggerAnalyzer);
