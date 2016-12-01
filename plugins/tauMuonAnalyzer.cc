// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/tauMuonAnalyzer
// Class:      tauMuonAnalyzer
// 
/**\class tauMuonAnalyzer tauMuonAnalyzer.cc GGHAA2Mu2TauAnalysis/tauMuonAnalyzer/plugins/tauMuonAnalyzer.cc

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

class tauMuonAnalyzer : public edm::EDAnalyzer{
   public:
      explicit tauMuonAnalyzer(const edm::ParameterSet&);
      ~tauMuonAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::GenParticleCollection>  genParticleTag_;
      edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> thirdHighestPtMuon_;
      edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> Mu1Mu2_;
      std::map<std::string, TH1D*> histos1D_;
      std::vector<double> Mu3PtBins_;
      std::vector<double> Mu3dRBins_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//`
tauMuonAnalyzer::tauMuonAnalyzer(const edm::ParameterSet& iConfig):
  genParticleTag_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleTag"))),  
  thirdHighestPtMuon_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("thirdHighestPtMuon"))),
  Mu1Mu2_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("Mu1Mu2"))),
  histos1D_(),
  Mu3PtBins_(iConfig.getParameter<std::vector<double> >("Mu3PtBins")),
  Mu3dRBins_(iConfig.getParameter<std::vector<double>>("Mu3dRBins"))
{
}


tauMuonAnalyzer::~tauMuonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
tauMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle<edm::RefVector<std::vector<reco::Muon>>> thirdHighestPtMuon;
   iEvent.getByToken(thirdHighestPtMuon_,thirdHighestPtMuon);

   edm::Handle<edm::RefVector<std::vector<reco::Muon>>> pMu1Mu2;
   iEvent.getByToken(Mu1Mu2_, pMu1Mu2);

  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByToken(genParticleTag_, pGenParticles);

  std::vector<reco::GenParticle*> genObjPtrs;
  for (typename reco::GenParticleCollection::const_iterator iGenObj = pGenParticles->begin(); 
       iGenObj != pGenParticles->end(); ++iGenObj) {
    genObjPtrs.push_back(const_cast<reco::GenParticle*>(&(*iGenObj)));
  }

  std::vector<reco::Muon*> Mu1Mu2Ptrs;
  for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iMu1Mu2=pMu1Mu2->begin(); iMu1Mu2!=pMu1Mu2->end();++iMu1Mu2)
  {
     Mu1Mu2Ptrs.push_back(const_cast<reco::Muon*>(&(**iMu1Mu2)));
  }

 
   for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator ithirdHighestPtMuon=thirdHighestPtMuon->begin();ithirdHighestPtMuon!=thirdHighestPtMuon->end(); ++ithirdHighestPtMuon)
   {
     int nearestGenObjKey=-1;
     const reco::GenParticle* nearestGenObj=
       Common::nearestObject(*ithirdHighestPtMuon, genObjPtrs, nearestGenObjKey);
       while(abs(nearestGenObj->motherRef()->pdgId())==13)
           nearestGenObj=(&(*(nearestGenObj->motherRef())));

     int nearestRecoObjKey=-1;
     const reco::Muon* nearestMu1Mu2=
       Common::nearestObject(*ithirdHighestPtMuon, Mu1Mu2Ptrs, nearestRecoObjKey);
     double deltaR=0.0;
     deltaR=reco::deltaR(**ithirdHighestPtMuon, *nearestMu1Mu2);
     histos1D_["deltaR_reco"]->Fill(deltaR);
     std::cout<< " Mu3's Nearest Gen Muon's mother=="<<nearestGenObj->motherRef()->pdgId()<<std::endl;
     if(abs(nearestGenObj->motherRef()->pdgId())==15)
       histos1D_["pt_Mu3_if"]->Fill((*ithirdHighestPtMuon)->pt());
     else
   	{
          std::cout<< "pt of Mu3=="<< (*ithirdHighestPtMuon)->pt()<<std::endl;
	  std::cout<< "eta of Mu3=="<< (*ithirdHighestPtMuon)->eta()<<std::endl;
          for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iMu1Mu2=pMu1Mu2->begin(); iMu1Mu2!=pMu1Mu2->end();++iMu1Mu2)
          {
            std::cout<< "pt of Mu1Mu2=="<< (*iMu1Mu2)->pt()<<std::endl;
	  }
          histos1D_["pt_Mu3_else"]->Fill((*ithirdHighestPtMuon)->pt());
	}
     histos1D_["pt_gen"]->Fill(nearestGenObj->pt()); 
     histos1D_["pt_reco"]->Fill((*ithirdHighestPtMuon)->pt());
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
tauMuonAnalyzer::beginJob()
{
  edm::Service<TFileService> fileService;
  histos1D_["deltaR_reco"]=fileService->make<TH1D>("deltaR of Mu3 and nearest Mu1 or Mu2","deltaR of Mu3 and nearest Mu1 or Mu2", Mu3dRBins_.size()-1, &Mu3dRBins_[0]);
  histos1D_["deltaR_reco"]->Sumw2();
  histos1D_["pt_gen"]=fileService->make<TH1D>("pt of gen", "pt of Mu3's nearest generation-level-muon (H750a09)", Mu3PtBins_.size()-1,&Mu3PtBins_[0]);
  histos1D_["pt_gen"]->Sumw2();
  histos1D_["pt_reco"]=fileService->make<TH1D>("pt of reco muon","pt of Mu3 (H750a09)",Mu3PtBins_.size()-1,&Mu3PtBins_[0]);
  histos1D_["pt_reco"]->Sumw2();
  histos1D_["pt_Mu3_else"]=fileService->make<TH1D>("pt of Mu3 else", "pt of Mu3 if Mu3 is not from tau", Mu3PtBins_.size()-1, &Mu3PtBins_[0]);
  histos1D_["pt_Mu3_else"]->Sumw2();
  histos1D_["pt_Mu3_if"]=fileService->make<TH1D>("pt of Mu3 if", "pt of Mu3 if Mu3 is from tau", Mu3PtBins_.size()-1, &Mu3PtBins_[0]);
  histos1D_["pt_Mu3_if"]->Sumw2();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
tauMuonAnalyzer::endJob() 
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
tauMuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(tauMuonAnalyzer);
