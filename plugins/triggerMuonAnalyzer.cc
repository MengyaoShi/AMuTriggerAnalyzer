// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/triggerMuonAnalyzer
// Class:      triggerMuonAnalyzer
// 
/**\class triggerMuonAnalyzer triggerMuonAnalyzer.cc GGHAA2Mu2TauAnalysis/triggerMuonAnalyzer/plugins/triggerMuonAnalyzer.cc

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

class triggerMuonAnalyzer : public edm::EDAnalyzer{
   public:
      explicit triggerMuonAnalyzer(const edm::ParameterSet&);
      ~triggerMuonAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> triggerMuon_;
      edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> Mu1Mu2_;
      edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> Mu1Mu2NoVeto_;
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
triggerMuonAnalyzer::triggerMuonAnalyzer(const edm::ParameterSet& iConfig):
  triggerMuon_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("triggerMuon"))),
  Mu1Mu2_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("Mu1Mu2"))),
  Mu1Mu2NoVeto_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("Mu1Mu2NoVeto"))),
  histos1D_()
{
}


triggerMuonAnalyzer::~triggerMuonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
triggerMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle<edm::RefVector<std::vector<reco::Muon>>> ptriggerMuon;
   iEvent.getByToken(triggerMuon_,ptriggerMuon);

   edm::Handle<edm::RefVector<std::vector<reco::Muon>>> pMu1Mu2;
   iEvent.getByToken(Mu1Mu2_, pMu1Mu2);

   edm::Handle<edm::RefVector<std::vector<reco::Muon>>> pMu1Mu2NoVeto;
   iEvent.getByToken(Mu1Mu2NoVeto_, pMu1Mu2NoVeto);

  std::vector<reco::Muon*> Mu1Mu2Ptrs;
  for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iMu1Mu2=pMu1Mu2->begin(); iMu1Mu2!=pMu1Mu2->end();++iMu1Mu2)
  {
     Mu1Mu2Ptrs.push_back(const_cast<reco::Muon*>(&(**iMu1Mu2)));
  }

 
   for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iMu1Mu2NoVeto=pMu1Mu2NoVeto->begin(); iMu1Mu2NoVeto!=pMu1Mu2NoVeto->end();++iMu1Mu2NoVeto)
   {
      histos1D_["pt_Mu1Mu2"]->Fill((*iMu1Mu2NoVeto)->pt());
   }
   for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator itriggerMuon=ptriggerMuon->begin(); itriggerMuon!=ptriggerMuon->end(); ++itriggerMuon)
   {

     int nearestRecoObjKey=-1;
     const reco::Muon* nearestMu1Mu2=
       Common::nearestObject(*itriggerMuon, Mu1Mu2Ptrs, nearestRecoObjKey);
     double deltaR=0.0;
     double invMass=0.0;
     deltaR=reco::deltaR(**itriggerMuon, *nearestMu1Mu2);
     invMass=((*itriggerMuon)->p4()+nearestMu1Mu2->p4()).M();
     histos1D_["deltaR"]->Fill(deltaR);
     histos1D_["pt_trigger"]->Fill((*itriggerMuon)->pt()); 
     histos1D_["pt_Mu2"]->Fill(nearestMu1Mu2->pt());
     histos1D_["invMass"]->Fill(invMass);
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
triggerMuonAnalyzer::beginJob()
{
  edm::Service<TFileService> fileService;
  histos1D_["deltaR"]=fileService->make<TH1D>("deltaR of trigger muon and nearest Mu1 or Mu2","deltaR of trigger muon and nearest Mu1 or Mu2", 10, 0.0, 0.1);
  histos1D_["pt_trigger"]=fileService->make<TH1D>("pt of trigger", "pt of trigger muon(H750a09)", 10, 00.0, 300);
  histos1D_["pt_Mu2"]=fileService->make<TH1D>("pt of reco muon","pt of trigger muon's nearest muon choosing from Mu1 or Mu2 (H750a09)",10, 0.0, 300);
  histos1D_["pt_Mu1Mu2"]=fileService->make<TH1D>("pt of all Mu1Mu2","pt of all Mu1Mu2 (H750a09)",10, 0.0, 300);
  histos1D_["invMass"]=fileService->make<TH1D>("invMass", "Invariant Mass of Mu1 Mu2", 100, 0, 10);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
triggerMuonAnalyzer::endJob() 
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
triggerMuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(triggerMuonAnalyzer);
