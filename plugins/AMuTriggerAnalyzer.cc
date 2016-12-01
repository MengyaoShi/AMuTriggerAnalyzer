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
#include "TEfficiency.h"
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
  edm::EDGetTokenT<reco::GenParticleRefVector> selectedGenParticleTag_; 
   std::vector<double> PtBins_;
   std::vector<double> EtaBins_;
   std::vector<double> dRBins_;
      std::map<std::string, TH1D*> histos1D_;
      std::map<std::string, TH2D*> histos2D_;
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
  selectedGenParticleTag_(consumes<reco::GenParticleRefVector>(iConfig.getParameter<edm::InputTag>("selectedGenParticleTag"))),
  PtBins_(iConfig.getParameter<std::vector<double>>("PtBins")),
  EtaBins_(iConfig.getParameter<std::vector<double>>("EtaBins")),
  dRBins_(iConfig.getParameter<std::vector<double>>("dRBins")),
  histos1D_(),
  histos2D_()
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
   
   std::cout<<"In Analyzer, MuonsPassTrigger->size()="<<MuonsPassTrigger->size()<<std::endl;

   edm::Handle<edm::RefVector<std::vector<reco::Muon>>> GenMatchedMuons;
   iEvent.getByToken(GenMatchedRecoMuonTag_,GenMatchedMuons);   
 
   edm::Handle<reco::GenParticleRefVector> pSelectedGenParticles;
   iEvent.getByToken(selectedGenParticleTag_, pSelectedGenParticles);
 
   std::vector<reco::GenParticle*> genObjPtrs;
  for (typename reco::GenParticleRefVector::const_iterator iGenObj = pSelectedGenParticles->begin(); 
       iGenObj != pSelectedGenParticles->end(); ++iGenObj) {
    genObjPtrs.push_back(const_cast<reco::GenParticle*>((*iGenObj).get()));
  }
  
   double deltaRGen=0.0;
   deltaRGen=reco::deltaR(*genObjPtrs[0], *genObjPtrs[1]);


   if(GenMatchedMuons->size()!=0)
   {
     histos1D_["denominator_dR"]->Fill(deltaRGen);
     double maxGenMatchedMuonPt=0;
     reco::MuonRef maxPtGenMatchedMuon;
     for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iGenMatchedMuon=GenMatchedMuons->begin(); iGenMatchedMuon!=GenMatchedMuons->end(); ++iGenMatchedMuon)
     {
       if((*iGenMatchedMuon)->pt() > maxGenMatchedMuonPt)
       {
         maxGenMatchedMuonPt=(*iGenMatchedMuon)->pt();
         maxPtGenMatchedMuon=(*iGenMatchedMuon);
       }
       else
         continue;
     }
     histos1D_["denominator_pt"]->Fill(maxPtGenMatchedMuon->pt());
     if(MuonsPassTrigger->size()!=0 )
     {
       histos1D_["numerator_dR"]->Fill(deltaRGen);
     /*double maxTriggerMuonPt=0;
 *      reco::MuonRef maxPtTriggerMuon;
 *           for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iMuonPassTrigger=MuonsPassTrigger->begin();iMuonPassTrigger!=MuonsPassTrigger->end(); ++iMuonPassTrigger)
 *                {
 *                       if((*iMuonPassTrigger)->pt()> maxTriggerMuonPt)
 *                              {
 *                                       maxTriggerMuonPt=(*iMuonPassTrigger)->pt();
 *                                                maxPtTriggerMuon=(*iMuonPassTrigger);
 *                                                       }
 *                                                              else
 *                                                                       continue;
 *                                                                            }*/
       histos1D_["numerator_pt"]->Fill(maxPtGenMatchedMuon->pt());
     }
   }

   if(MuonsPassTrigger->size() <( GenMatchedMuons->size()) && (MuonsPassTrigger->size()==0))
   {
     std::cout<<"Inefficiency source! Pt of GenMatchedMuon"<<std::endl;
     std::cout<<"MuonsPassTrigger->size()="<<MuonsPassTrigger->size()<<"; GenMatchedMuons->size()="<<GenMatchedMuons->size()<<std::endl;
     for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iGenMatchedMuon=GenMatchedMuons->begin(); iGenMatchedMuon!=GenMatchedMuons->end(); ++iGenMatchedMuon)
     {
       std::cout<<"deltaRGen="<< deltaRGen<<std::endl;
       std::cout<<"Eta="<< (*iGenMatchedMuon)->eta()<<std::endl;
     }
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
AMuTriggerAnalyzer::beginJob()
{
  std::cout<<"in beginJob"<<std::endl;
  edm::Service<TFileService> fileService;

  histos1D_["denominator_dR"]=fileService->make<TH1D>("denominator_dR", "dR ", dRBins_.size()-1, &dRBins_[0] );
  histos1D_["denominator_dR"]->SetXTitle("dR");
  histos1D_["denominator_dR"]->Sumw2();

  histos1D_["numerator_dR"]=fileService->make<TH1D>("numerator_dR","dR", dRBins_.size()-1, &dRBins_[0]);
  histos1D_["numerator_dR"]->SetXTitle("dR");
  histos1D_["numerator_dR"]->Sumw2();

  histos1D_["denominator_pt"]=fileService->make<TH1D>("denominator_pt", "pt of RecoMuonWithSoftID45GeV2.1etaCutGenMatchedToAMu H750a09 ", PtBins_.size()-1, &PtBins_[0] );
  histos1D_["denominator_pt"]->SetXTitle("pt (see title description)(GeV)");
  histos1D_["denominator_pt"]->Sumw2();  

  histos1D_["numerator_pt"]=fileService->make<TH1D>("numerator_pt","pt of RecoMuonWithSoftID45GeV2.1etaCutGenMatchedToAMu H750a09", PtBins_.size()-1, &PtBins_[0]);
  histos1D_["numerator_pt"]->SetXTitle("pt");
  histos1D_["numerator_pt"]->Sumw2();

histos1D_["denominator_eta"]=fileService->make<TH1D>("denominator_eta", "eta of RecoMuonWithSoftID45GeV2.1etaCutGenMatchedToAMu H750a09", EtaBins_.size()-1, &EtaBins_[0] );
  histos1D_["denominator_eta"]->SetXTitle("eta (see title description)");
  histos1D_["denominator_eta"]->Sumw2();

  histos1D_["numerator_eta"]=fileService->make<TH1D>("numerator_eta","eta of RecoMuonWithSoftID45GeV2.1etaCutGenMatchedToAMu H750a09", EtaBins_.size()-1, &EtaBins_[0]);
  histos1D_["numerator_eta"]->SetXTitle("eta");
  histos1D_["numerator_eta"]->Sumw2();

  histos1D_["Efficiency_dR"]=fileService->make<TH1D>("Efficiency_dR","Efficiency_dR)", dRBins_.size()-1, &dRBins_[0]);

  histos1D_["Efficiency_pt"]=fileService->make<TH1D>("Efficiency_pt","Efficiency of pass hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q(no HLT fire required) with change of pt of object pass hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q (H750a09)", PtBins_.size()-1, &PtBins_[0]);
  histos1D_["Efficiency_pt"]->SetXTitle("Pt");
  histos1D_["Efficiency_pt"]->SetYTitle("Efficiency (see title description)");

  histos1D_["Efficiency_eta"]=fileService->make<TH1D>("Efficiency_eta","Efficiency of pass hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q(no HLT fire required) with change of eta of object pass hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q (H750a09)", EtaBins_.size()-1, &EtaBins_[0]);
  histos1D_["Efficiency_eta"]->SetXTitle("Eta");
  histos1D_["Efficiency_eta"]->SetYTitle("Efficiency (see title description)");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
AMuTriggerAnalyzer::endJob() 
{
 histos1D_["Efficiency_dR"]->Divide(histos1D_["numerator_dR"], histos1D_["denominator_dR"], 1, 1, "B");
  histos1D_["Efficiency_pt"]->Divide(histos1D_["numerator_pt"], histos1D_["denominator_pt"], 1, 1, "B");
  histos1D_["Efficiency_eta"]->Divide(histos1D_["numerator_eta"], histos1D_["denominator_eta"], 1, 1, "B");

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
