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
       void reset(const bool);   

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> GenMatchedRecoMuonTag_;
      edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> MuPassTrigger_;
      edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> GenMatchedRecoMuonTag2_;
   std::vector<double> EffBins_;
   std::vector<double> Chi2Bins_;
   std::vector<double> PtBins_;
   std::vector<double> EtaBins_;
      std::map<std::string, TH1D*> histos1D_;
      std::map<std::string, TH2D*> histos2D_;
  TFile* out_;
 std::string outFileName_;
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
  GenMatchedRecoMuonTag2_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("GenMatchedRecoMuonTag2"))),
  EffBins_(iConfig.getParameter<std::vector<double> >("DRBins")),
  Chi2Bins_(iConfig.getParameter<std::vector<double>>("Chi2Bins")),
  PtBins_(iConfig.getParameter<std::vector<double>>("PtBins")),
  EtaBins_(iConfig.getParameter<std::vector<double>>("EtaBins")),
  histos1D_(),
  histos2D_(),
  outFileName_(iConfig.getParameter<std::string>("outFileName"))
{
   reset(false);
}


AMuTriggerAnalyzer::~AMuTriggerAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  reset(true);
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
  
   edm::Handle<edm::RefVector<std::vector<reco::Muon>>> GenMatchedMuons2;
   iEvent.getByToken(GenMatchedRecoMuonTag2_,GenMatchedMuons2);

   std::vector<reco::Muon*> recoMuon2Ptrs;
   for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iGenMatchedMuon2=GenMatchedMuons2->begin(); iGenMatchedMuon2!= GenMatchedMuons2->end(); ++iGenMatchedMuon2)
   {
      recoMuon2Ptrs.push_back(const_cast<reco::Muon*>(iGenMatchedMuon2->get()));
   }
   

   for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iMuonPassTrigger=MuonsPassTrigger->begin();iMuonPassTrigger!=MuonsPassTrigger->end(); ++iMuonPassTrigger)
   {
     int nearestRecoObjKey1=-1;
     const reco::Muon* nearestRecoObj1=
       Common::nearestObject(*iMuonPassTrigger, recoMuon2Ptrs, nearestRecoObjKey1);
     if(nearestRecoObj1!=NULL)
     histos1D_["numerator"]->Fill(reco::deltaR((*iMuonPassTrigger)->eta(), (*iMuonPassTrigger)->phi(), nearestRecoObj1->eta(), nearestRecoObj1->phi()));

     histos1D_["numerator_pt"]->Fill((*iMuonPassTrigger)->pt());
     histos1D_["numerator_eta"]->Fill((*iMuonPassTrigger)->eta());
     histos2D_["twoDnumerator"]->Fill(reco::deltaR((*iMuonPassTrigger)->eta(), (*iMuonPassTrigger)->phi(), nearestRecoObj1->eta(), nearestRecoObj1->phi()), (*iMuonPassTrigger)->eta());
   }
   for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iGenMatchedMuon=GenMatchedMuons->begin(); iGenMatchedMuon!=GenMatchedMuons->end(); ++iGenMatchedMuon)
   {
     int nearestRecoObjKey2=-1;
     const reco::Muon* nearestRecoObj2=
      Common::nearestObject(*iGenMatchedMuon, recoMuon2Ptrs,nearestRecoObjKey2);
     if(nearestRecoObj2!=NULL)
     {
     histos1D_["denominator"]->Fill(reco::deltaR((*iGenMatchedMuon)->eta(),(*iGenMatchedMuon)->phi(), nearestRecoObj2->eta(), nearestRecoObj2->phi()));
     histos2D_["2dplot"]->Fill(reco::deltaR((*iGenMatchedMuon)->eta(),(*iGenMatchedMuon)->phi(), nearestRecoObj2->eta(), nearestRecoObj2->phi()),nearestRecoObj2->pt()); 
     }
//     histos1D_["deno_Chi2"]->Fill((*iGenMatchedMuon)->globalTrack()->normalizedChi2());
     histos1D_["deno_Chi2"]->Fill((*iGenMatchedMuon)->globalTrack()->normalizedChi2());
     histos1D_["denominator_pt"]->Fill((*iGenMatchedMuon)->pt());
     histos1D_["denominator_eta"]->Fill((*iGenMatchedMuon)->eta());
     histos2D_["twoDdenominator"]->Fill((reco::deltaR((*iGenMatchedMuon)->eta(),(*iGenMatchedMuon)->phi(), nearestRecoObj2->eta(), nearestRecoObj2->phi())), (*iGenMatchedMuon)->eta());
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
AMuTriggerAnalyzer::beginJob()
{
out_ = new TFile(outFileName_.c_str(), "RECREATE");
  edm::Service<TFileService> fileService;
  histos1D_["denominator"]=fileService->make<TH1D>("denominator", "dR Of RecoMuonWithSoftID45GeV2.1etaCutGenMatchedToAMu vs. SecondHighestPTGenMatchedRecoMuonPassSoftIDOppositeSign H750a09",  EffBins_.size()-1, &EffBins_[0]);
  histos1D_["denominator"]->Sumw2();
  histos1D_[ "denominator" ]->SetXTitle( "dR (see title description)" );

  histos1D_["numerator"]=fileService->make<TH1D>("numerator","Denominator*FiredHLT*TriggerMatch(PT match) H750a09", EffBins_.size()-1, &EffBins_[0]);
  histos1D_["numerator"]->SetXTitle("dR");
  histos1D_["numerator"]->Sumw2();

  histos1D_["deno_Chi2"]=fileService->make<TH1D>("deno_Chi2","Chi2 of RecoMuonGenMatchedtoAMu H750a09", Chi2Bins_.size()-1, &Chi2Bins_[0]);
  histos1D_["deno_Chi2"]->SetXTitle("Chi2");
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

  histos1D_["Efficiency"]=fileService->make<TH1D>("Efficiency","Efficiency of pass hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q(no HLT fire required) with change of DeltaR H750a09", EffBins_.size()-1, &EffBins_[0]);
  histos1D_["Efficiency"]->SetXTitle("dR");
  histos1D_["Efficiency"]->SetYTitle("Efficiency (see title description)");
 
  histos1D_["Efficiency_pt"]=fileService->make<TH1D>("Efficiency_pt","Efficiency of pass hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q(no HLT fire required) with change of pt of object pass hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q (H750a09)", PtBins_.size()-1, &PtBins_[0]);
  histos1D_["Efficiency_pt"]->SetXTitle("Pt");
  histos1D_["Efficiency_pt"]->SetYTitle("Efficiency (see title description)");

  histos1D_["Efficiency_eta"]=fileService->make<TH1D>("Efficiency_eta","Efficiency of pass hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q(no HLT fire required) with change of eta of object pass hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q (H750a09)", EtaBins_.size()-1, &EtaBins_[0]);
  histos1D_["Efficiency_eta"]->SetXTitle("Eta");
  histos1D_["Efficiency_eta"]->SetYTitle("Efficiency (see title description)");

  histos2D_["2dplot"]=fileService->make<TH2D>("2Dplot","2Dplot of deltaR mentioned before vs pt of mu2 H750a09",100,0.0,0.5, 100, 0, 200);
  histos2D_["2dplot"]->SetXTitle("dR");
  histos2D_["2dplot"]->SetYTitle("pt of Mu2(Second HighestPTGenMatchedRecoMuonPassSoftIDOppositeSign)");

  histos2D_["twoDdenominator"]=fileService->make<TH2D>("twoDdenominator","Delta R between Mu1 Mu2 vs Eta of Mu1 ",  EffBins_.size()-1, &EffBins_[0],EtaBins_.size()-1, &EtaBins_[0] );
  histos2D_["twoDdenominator"]->SetXTitle("dR between Mu1 Mu2");
  histos2D_["twoDdenominator"]->SetYTitle("Eta of Mu1");

histos2D_["twoDnumerator"]=fileService->make<TH2D>("twoDnumerator","Delta R between Mu1 Mu2 with Mu1 pass Subfilter VS Eta of Mu1 pass Subfilter", EffBins_.size()-1, &EffBins_[0],EtaBins_.size()-1, &EtaBins_[0] );
  histos2D_["twoDnumerator"]->SetXTitle("dR between Mu1 Mu2 with Mu1 pass Subfilter");
  histos2D_["twoDnumerator"]->SetYTitle("Eta of Mu1 pass subfilter ");
 
  histos2D_["Efficiency_2D"]=fileService->make<TH2D>("Efficiency_2D", "Efficiency of pass hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q(no HLT fire required)(H750a09)",EffBins_.size()-1, &EffBins_[0],EtaBins_.size()-1, &EtaBins_[0] );
  histos2D_["Efficiency_2D"]->SetXTitle("dR of mu1 mu2");
  histos2D_["Efficiency_2D"]->SetYTitle("Eta of Mu1");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AMuTriggerAnalyzer::endJob() 
{
  histos1D_["Efficiency"]->Divide(histos1D_["numerator"], histos1D_["denominator"],1, 1,"B");
  histos1D_["Efficiency_pt"]->Divide(histos1D_["numerator_pt"], histos1D_["denominator_pt"], 1, 1, "B");
  histos1D_["Efficiency_eta"]->Divide(histos1D_["numerator_eta"], histos1D_["denominator_eta"], 1, 1, "B");
  histos2D_["Efficiency_2D"]->Divide(histos2D_["twoDnumerator"],histos2D_["twoDdenominator"], 1, 1,"B");

  TEfficiency *Eff2D= new TEfficiency(*histos2D_["twoDnumerator"],*histos2D_["twoDdenominator"]);
out_->cd(); 
 TCanvas c1("example", "",600,600);
   c1.SetFillStyle(1001);
   c1.SetFillColor(kWhite);
   Eff2D->Draw("AP");
   c1.Write();
  out_->Write();
  out_->Close();
}
void 
AMuTriggerAnalyzer::reset(const bool doDelete)
{
  if (doDelete && (out_ != NULL)) delete out_;
  out_ = NULL;
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
