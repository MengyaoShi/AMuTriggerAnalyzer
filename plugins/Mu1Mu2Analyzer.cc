// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/Mu1Mu2Analyzer
// Class:      Mu1Mu2Analyzer
// 
/**\class Mu1Mu2Analyzer Mu1Mu2Analyzer.cc GGHAA2Mu2TauAnalysis/Mu1Mu2Analyzer/plugins/Mu1Mu2Analyzer.cc

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

class Mu1Mu2Analyzer : public edm::EDAnalyzer{
   public:
      explicit Mu1Mu2Analyzer(const edm::ParameterSet&);
      ~Mu1Mu2Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::GenParticleCollection>  genParticleTag_;
      edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> Mu1Mu2_;
      edm::EDGetTokenT<reco::PFCandidateCollection> particleFlow_;
      std::map<std::string, TH1D*> histos1D_;
      std::map<std::string, TH2D*> histos2D_;
      int countbothfromW=0;
      std::vector<double> Mu2PtBins_;
      std::vector<double> invMassBins_; 
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
Mu1Mu2Analyzer::Mu1Mu2Analyzer(const edm::ParameterSet& iConfig):
  genParticleTag_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleTag"))),  
  Mu1Mu2_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("Mu1Mu2"))),
  particleFlow_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("particleFlow"))),
  histos1D_(),
  histos2D_(),
  Mu2PtBins_(iConfig.getParameter<std::vector<double> >("Mu2PtBins")),
  invMassBins_(iConfig.getParameter<std::vector<double>>("invMassBins"))
{
}


Mu1Mu2Analyzer::~Mu1Mu2Analyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Mu1Mu2Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<edm::RefVector<std::vector<reco::Muon>>> pMu1Mu2;
   iEvent.getByToken(Mu1Mu2_, pMu1Mu2);

  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByToken(genParticleTag_, pGenParticles);

  edm::Handle<reco::PFCandidateCollection> pPFCandidates;
  iEvent.getByToken(particleFlow_, pPFCandidates);

  std::vector<reco::GenParticle*> genObjPtrs;
  for (typename reco::GenParticleCollection::const_iterator iGenObj = pGenParticles->begin(); 
       iGenObj != pGenParticles->end(); ++iGenObj) {
    genObjPtrs.push_back(const_cast<reco::GenParticle*>(&(*iGenObj)));
  }
 
  std::vector<reco::Muon*> Mu1Mu2Ptrs;
  double invMass=0;
  reco::Muon* HighestPtMu1Mu2;
  reco::Muon* LowestPtMu1Mu2;
  invMass=((*pMu1Mu2)[0]->p4()+(*pMu1Mu2)[1]->p4()).M();
  if((*pMu1Mu2)[0]->pt()> (*pMu1Mu2)[1]->pt()){
    HighestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[0])));
    LowestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[1])));
    }
  else
    {
    HighestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[1])));
    LowestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[0])));
    } 
  double Mu2Pt=0;
  double dR=0.0;
  dR=deltaR(*LowestPtMu1Mu2, *HighestPtMu1Mu2);
  Mu2Pt=LowestPtMu1Mu2->pt();
  histos2D_["dRVsMu2Pt"]->Fill(dR, Mu2Pt); 
  histos2D_["Mu1PtMu2Pt"]->Fill(HighestPtMu1Mu2->pt(), Mu2Pt);
  histos1D_["Mu1Pt"]->Fill(HighestPtMu1Mu2->pt());
  histos1D_["Mu2Pt"]->Fill(Mu2Pt);
  histos1D_["dRMu1Mu2"]->Fill(dR);
  histos1D_["dRMu1Mu2Wider"]->Fill(dR);
  std::vector<reco::PFCandidate*> PFCandidatePtrs;
  double PtSum=0;
  double Iso=0;
  for(typename reco::PFCandidateCollection::const_iterator iPFCandidate = pPFCandidates->begin(); iPFCandidate != pPFCandidates->end(); ++iPFCandidate){
    PFCandidatePtrs.push_back(const_cast<reco::PFCandidate*> (&(*iPFCandidate)));
    if(deltaR(*HighestPtMu1Mu2, *iPFCandidate )<0.4&&(deltaR(*HighestPtMu1Mu2, *iPFCandidate)>0.0001&&((*HighestPtMu1Mu2).pt()-(*iPFCandidate).pt())/((*HighestPtMu1Mu2).pt())>0.0001)&&(deltaR(*LowestPtMu1Mu2, *iPFCandidate)>0.0001&&((*LowestPtMu1Mu2).pt()-(*iPFCandidate).pt())/((*LowestPtMu1Mu2).pt())>0.0001))
      PtSum+=(*iPFCandidate).pt();
    else continue;
  }
 
  Iso=PtSum/(HighestPtMu1Mu2->pt());
  for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iMu1Mu2=pMu1Mu2->begin(); iMu1Mu2!=pMu1Mu2->end();++iMu1Mu2)
  {
     Mu1Mu2Ptrs.push_back(const_cast<reco::Muon*>(&(**iMu1Mu2)));
  }
  
   std::cout<<"*******************New Event************************"<<std::endl;
   int countMufromW=0; 
   for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iMu1Mu2=pMu1Mu2->begin();iMu1Mu2!=pMu1Mu2->end(); ++iMu1Mu2)
   {
    
     int nearestGenObjKey=-1;
     const reco::GenParticle* nearestGenObj=
       Common::nearestObject(*iMu1Mu2, genObjPtrs, nearestGenObjKey);
       while(abs(nearestGenObj->motherRef()->pdgId())==13)
           nearestGenObj=(&(*(nearestGenObj->motherRef())));

     std::cout<< " Mu1 Mu2's Nearest Gen Muon's mother=="<<nearestGenObj->motherRef()->pdgId()<<std::endl;
     if(abs(nearestGenObj->motherRef()->pdgId())==24)
       countMufromW++;
          std::cout<< "pt of Mu1 Mu2=="<< (*iMu1Mu2)->pt()<<std::endl;
	  std::cout<< "eta of Mu1 Mu2=="<< (*iMu1Mu2)->eta()<<std::endl;
     histos1D_["pt_gen"]->Fill(nearestGenObj->pt()); 
     histos1D_["pt_reco"]->Fill((*iMu1Mu2)->pt());
   }
   if(countMufromW==2)
   {
      countbothfromW++;
      histos1D_["pt_ifMu1Mu2arebothfromW"]->Fill(LowestPtMu1Mu2->pt());
   }
      
   histos1D_["invMass"]->Fill(invMass);
   histos2D_["Iso"]->Fill(Iso,LowestPtMu1Mu2->pt());
   histos1D_["Iso1D"]->Fill(Iso);
}


// ------------ method called once each job just before starting event loop  ------------
void 
Mu1Mu2Analyzer::beginJob()
{
  edm::Service<TFileService> fileService;
  histos1D_["pt_gen"]=fileService->make<TH1D>("pt of gen", "pt of Mu1 Mu2's nearest generation-level-muon (H750a09)", 100, 0.0, 100);
  histos1D_["pt_reco"]=fileService->make<TH1D>("pt of reco muon","pt of Mu1 Mu2 (H750a09)",100, 0.0, 100);
  histos1D_["invMass"]=fileService->make<TH1D>("invMass of Mu1 Mu2","invMass of Mu1 Mu2 (H750a09)",invMassBins_.size()-1, &invMassBins_[0]);
//  histos1D_["invMass"]->Sumw2();
  histos2D_["Iso"]=fileService->make<TH2D>("Iso", "Iso", 50, 0.0, 5.0, 300, 0, 300);
  histos2D_["Iso"]->SetXTitle("relative Isolation within cone 0.4");
  histos2D_["Iso"]->SetYTitle("Pt of Mu2");
  histos2D_["Mu1PtMu2Pt"]=fileService->make<TH2D>("Mu1PtMu2Pt","Mu1Pt vs Mu2Pt", 100, 0, 500, 40, 0, 200);
  histos1D_["Mu2Pt"]=fileService->make<TH1D>("pt of Mu2", "Pt of RecoMu2",Mu2PtBins_.size()-1,&Mu2PtBins_[0]);
  histos1D_["Mu2Pt"]->Sumw2();
  histos1D_["Mu1Pt"]=fileService->make<TH1D>("pt of Mu1", "pt of Mu1", 100, 0, 500);
  histos1D_["Iso1D"]=fileService->make<TH1D>("Iso1D", "Iso1D", 50, 0, 5.0);
  histos1D_["dRMu1Mu2"]=fileService->make<TH1D>("dRMu1Mu2", "dRMu1Mu2", 50, 0, 5.0);
  histos1D_["dRMu1Mu2Wider"]=fileService->make<TH1D>("dRMu1Mu2Wider", "dRMu1Mu2Wider", 50, 0, 5.0);
  histos2D_["dRVsMu2Pt"]=fileService->make<TH2D>("dRVsMu2Pt", "dRVsMu2Pt", 50, 0, 5.0, 50, 0, 50.0);
  histos1D_["pt_ifMu1Mu2arebothfromW"]=fileService->make<TH1D>("ptMu2_if","pt of Mu2 if both Mu1 Mu2 are from W", 50, 0, 150);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Mu1Mu2Analyzer::endJob() 
{
   std::cout<<"CountBothfromW="<<countbothfromW<<std::endl;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Mu1Mu2Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Mu1Mu2Analyzer);
