// -*- C++ -*-
// 
// Package:    seedAnalyzer/seedAnalyzer
// Class:      seedAnalyzer
// 
/**\class seedAnalyzer seedAnalyzer.cc seedAnalyzer/seedAnalyzer/plugins/seedAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Leonardo GIANNINI
//         Created:  Thu, 01 Dec 2016 15:40:12 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"


#include "DataFormats/JetReco/interface/PFJet.h"
#include <DataFormats/TrackReco/interface/Track.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"

#include <RecoVertex/VertexPrimitives/interface/BasicVertexState.h>
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/SharedTracks.h"

#include <DataFormats/Candidate/interface/Candidate.h>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/CandidatePtrTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "seedAnalyzer/seedAnalyzer/interface/trackVars2.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

/*#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"


#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
*/

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.


class AnalyzerSignedIP_MINIAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit AnalyzerSignedIP_MINIAOD(const edm::ParameterSet&);
      ~AnalyzerSignedIP_MINIAOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      // ----------member data ---------------------------
      
      int funzionePippo (int vivaLaVita);
      int seedToJetMatching(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi);
      std::vector<int> seedToAllJets30Matching(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi, std::vector<double> jet_pt);
      int seedToJet30Matching(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi, std::vector<double> jet_pt);
      std::pair<int,float> ClosestJet(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi);
      int genLevel(std::vector<pat::PackedGenParticle> genP, reco::TransientTrack seed, std::vector<int>& MomFlav, std::vector<int>& BChain, std::vector<double>& allChi);      
      int genMap(std::vector<pat::PackedGenParticle> genP, std::vector<reco::TransientTrack> tracks, std::vector<int>& genPnumbers, std::vector<double>& allChi);      
      double trackEff(std::vector<pat::PackedGenParticle> genP, std::vector<reco::TransientTrack> tracks, double &BChain_eff, double &DChain_eff, double &B_eff,
      double &low_pt_eff, double &high_pt_eff, std::vector<double> &minGenChiSquare);
      
      //genMap(genParticles, selectedTracks, genP_index, Match_chisquare );
      bool genLevelMM(std::vector<pat::PackedGenParticle> genP, std::vector<reco::TransientTrack> tracks, std::vector<int>& genPnumbers, std::vector<double>& assigned_chisquares);
      
      std::vector<float> firstHadronInChain(pat::PackedGenParticle gParticle, const reco::Candidate* mother);
      
      edm::Service<TFileService> file;
//       edm::RefVector<std::vector<pat::PackedGenParticle> > mothers;
//       reco::Candidate mother;

      TTree *tree;
      
      double m;
      float n;
      
      double match_eff=0;
      double tracking_eff=-10;
      
      double BChain_eff=0;
      double DChain_eff=0;
      double B_eff=0;
      double low_pt_eff=0;
      double high_pt_eff=0;
      
      int b_found=0, d_found=0;
      int id_mom;
      int bd_i=0;
      
      double gen_pv_x=-10;
      double gen_pv_y=-10;
      double gen_pv_z=-10;
      
      double pv_x=-10;
      double pv_y=-10;
      double pv_z=-10;
      
      int evt=0;
      int lumi=0;
      int run=0;
      
      int n_seed=0;
      
      std::vector<int> genP_index;
//      std::vector<int> genP_id_MM;
      
      std::vector<double> jet_pt;
      std::vector<double> jet_px;
      std::vector<double> jet_py;
      std::vector<double> jet_pz;
      std::vector<double> jet_eta;
      std::vector<double> jet_phi;
      std::vector<double> jet_mass;
      std::vector<int> jet_flavour;
      std::vector<int> jet_seedIndex;
      
      std::vector<double> seed_pt;
      std::vector<double> seed_eta;
      std::vector<double> seed_phi;
      std::vector<double> seed_mass;
      
      std::vector<double> seed_dz;
      std::vector<double> seed_dxy;
      std::vector<double> seed_3D_ip;
      std::vector<double> seed_3D_sip;
      std::vector<double> seed_2D_ip;
      std::vector<double> seed_2D_sip;
      std::vector<double> seed_3D_signedIp;
      std::vector<double> seed_3D_signedSip;
      std::vector<double> seed_2D_signedIp;
      std::vector<double> seed_2D_signedSip;
      std::vector<int> seed_JetMatch;
      
      std::vector<double> seed_chi2reduced;
      std::vector<double> seed_nPixelHits;
      std::vector<double> seed_nHits;
      std::vector<double> seed_jetAxisDistance;
      std::vector<double> seed_jetAxisDlength;
      std::vector<float> seed_ClosestJet_dR;
      
      
      
      //vars from matching genParticles if any
      std::vector<double> seed_MC_pt;
      std::vector<double> seed_MC_eta;
      std::vector<double> seed_MC_phi;
      std::vector<double> seed_MC_mass;
      std::vector<double> seed_MC_dz;
      std::vector<double> seed_MC_dxy;
      std::vector<int> seed_MC_MomPdgId;
      std::vector<int> seed_MC_MomFlavour;
      std::vector<int> seed_MC_BChain;
      std::vector<int> seed_MC_DChain; //no B in this case
      std::vector<double> seed_MC_vx;
      std::vector<double> seed_MC_vy;
      std::vector<double> seed_MC_vz;
      std::vector<double> seed_MC_pvd;
      //end vars from matching genParticles if any
      
      std::vector<int> nearTracks_Nvtx;
      std::vector<int> nearTracks_nTracks;//max 20
      std::vector<double> nearTracks_pt;
      std::vector<double> nearTracks_eta;
      std::vector<double> nearTracks_phi;
      std::vector<double> nearTracks_mass;
      std::vector<double> nearTracks_dz;
      std::vector<double> nearTracks_dxy;
      std::vector<double> nearTracks_3D_ip;
      std::vector<double> nearTracks_3D_sip;
      std::vector<double> nearTracks_2D_ip;
      std::vector<double> nearTracks_2D_sip;
//       std::vector<double> nearTracks_3D_signedIp;
//       std::vector<double> nearTracks_3D_signedSip;
//       std::vector<double> nearTracks_2D_signedIp;
//       std::vector<double> nearTracks_2D_signedSip;
      
      std::vector<double> nearTracks_PCAdist;
      std::vector<double> nearTracks_PCAdsig;
      
      std::vector<double> nearTracks_PCAonSeed_x;
      std::vector<double> nearTracks_PCAonSeed_y;
      std::vector<double> nearTracks_PCAonSeed_z;
      
      std::vector<double> nearTracks_PCAonSeed_xerr;
      std::vector<double> nearTracks_PCAonSeed_yerr;
      std::vector<double> nearTracks_PCAonSeed_zerr;
      
      std::vector<double> nearTracks_PCAonTrack_x;
      std::vector<double> nearTracks_PCAonTrack_y;
      std::vector<double> nearTracks_PCAonTrack_z;
      
      std::vector<double> nearTracks_PCAonTrack_xerr;
      std::vector<double> nearTracks_PCAonTrack_yerr;
      std::vector<double> nearTracks_PCAonTrack_zerr; 

      std::vector<double> nearTracks_dotprodTrack;
      std::vector<double> nearTracks_dotprodSeed;
      std::vector<double> nearTracks_dotprodTrackSeed2D;
      std::vector<double> nearTracks_dotprodTrackSeed3D;
      std::vector<double> nearTracks_dotprodTrackSeedVectors2D;
      std::vector<double> nearTracks_dotprodTrackSeedVectors3D;
      
      std::vector<double> nearTracks_PCAonSeed_pvd;
      std::vector<double> nearTracks_PCAonTrack_pvd;
      
      
      //vars from matching genParticles if any
      std::vector<double> nearTracks_MC_pt;
      std::vector<double> nearTracks_MC_eta;
      std::vector<double> nearTracks_MC_phi;
      std::vector<double> nearTracks_MC_dz;
      std::vector<double> nearTracks_MC_dxy;
      
      std::vector<int> nearTracks_MC_MomPdgId;
      std::vector<int> nearTracks_MC_MomFlavour;
      std::vector<int> nearTracks_MC_BChain;
      std::vector<int> nearTracks_MC_DChain; //no B in this case

      std::vector<double> nearTracks_MC_Track_vx;
      std::vector<double> nearTracks_MC_Track_vy;
      std::vector<double> nearTracks_MC_Track_vz;
      
      std::vector<double> nearTracks_MC_fromSeedVtx;
      std::vector<double> nearTracks_MC_fromSeedChain;
      std::vector<double> nearTracks_MC_pvd;
      //end vars from matching genParticles if any

      std::vector<double> Match_chisquare;
//      std::vector<double> Match_chisquare2;
      std::vector<double> minChiSquare_genP;
      std::vector<trackVars2> nearTracks;
      trackVars2 myTrack;

      float min3DIPValue=0.005;
      float min3DIPSignificance=1.2;
      int max3DIPValue=9999.;
      int max3DIPSignificance=9999.;
      
      /*
       * IVF configuration for seeding
       * seedMax3DIPSignificance = cms.double(9999.),
           seedMax3DIPValue = cms.double(9999.),
           seedMin3DIPSignificance = cms.double(1.2),
           seedMin3DIPValue = cms.double(0.005),*/
      
      std::pair<int, float> pair;
      
      
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> > CandidateToken;
      edm::EDGetTokenT<std::vector< reco::GenParticle > > prunedGenParticleToken;
      edm::EDGetTokenT<std::vector< pat::PackedGenParticle > > genParticleToken;
      edm::EDGetTokenT<std::vector<pat::Jet> > jetsToken;

      edm::EDGetTokenT<reco::VertexCollection> token_primaryVertex;
      edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;
//       edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetToken;
     
    
      
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



/*
•pT > 1 GeV 
•chi2/ndof < 5 
•|dz| < 17 cm 
•|dxy| < 2 cm 
•Nhits(pixel)≥ 1 
•Nhits(total)≥ 1
•Distance(jet axis, track) < 0.07 cm a.k.a. “distance”•
Distance(PV, closest point of the track to the jet axis) < 5 cm a.k.a. “length” 
•Then tracks are sorted w.r.t. IP2D significance and are used in the CSV 
computation
*/



AnalyzerSignedIP_MINIAOD::AnalyzerSignedIP_MINIAOD(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");
    tree=file->make<TTree>("tree","tree");
    
    tree->Branch("lumi",&lumi, "lumi/I");
    tree->Branch("evt",&evt, "evt/I");
    tree->Branch("run",&run, "run/I");
    
    tree->Branch("gen_pv_x",&gen_pv_x, "gen_pv_x/D");
    tree->Branch("gen_pv_y",&gen_pv_y, "gen_pv_y/D");
    tree->Branch("gen_pv_z",&gen_pv_z, "gen_pv_z/D");
    tree->Branch("pv_x",&pv_x, "pv_x/D");
    tree->Branch("pv_y",&pv_y, "pv_y/D");
    tree->Branch("pv_z",&pv_z, "pv_z/D");
    
    tree->Branch("match_eff", &match_eff, "match_eff/D");
    tree->Branch("tracking_eff", &tracking_eff, "tracking_eff/D");
    
    tree->Branch("BChain_eff", &BChain_eff, "BChain_eff/D");
    tree->Branch("DChain_eff", &DChain_eff, "DChain_eff/D");
    tree->Branch("B_eff", &B_eff, "B_eff/D");
    tree->Branch("low_pt_eff", &low_pt_eff, "low_pt_eff/D");
    tree->Branch("high_pt_eff", &high_pt_eff, "high_pt_eff/D");
    
    
    tree->Branch("jet_pt",&jet_pt);
    tree->Branch("jet_eta",&jet_eta);
    tree->Branch("jet_phi",&jet_phi);
    tree->Branch("jet_mass",&jet_mass);
    tree->Branch("jet_flavour",&jet_flavour);
    tree->Branch("jet_seedIndex",&jet_seedIndex);
    
    tree->Branch("n_seed",&n_seed, "n_seed/I");

    tree->Branch("seed_pt",&seed_pt);
    tree->Branch("seed_eta",&seed_eta);
    tree->Branch("seed_phi",&seed_phi);
    tree->Branch("seed_mass",&seed_mass);
    
    tree->Branch("seed_dz", &seed_dz);
    tree->Branch("seed_dxy", &seed_dxy);
    tree->Branch("seed_3D_ip", &seed_3D_ip);
    tree->Branch("seed_3D_sip", &seed_3D_sip);
    tree->Branch("seed_2D_ip", &seed_2D_ip);
    tree->Branch("seed_2D_sip", &seed_2D_sip);
    
    tree->Branch("seed_3D_signedIp", &seed_3D_signedIp);
    tree->Branch("seed_3D_signedSip", &seed_3D_signedSip);
    tree->Branch("seed_2D_signedIp", &seed_2D_signedIp);
    tree->Branch("seed_2D_signedSip", &seed_2D_signedSip);
    
    tree->Branch("seed_JetMatch", &seed_JetMatch);
    tree->Branch("seed_chi2reduced",&seed_chi2reduced);
    tree->Branch("seed_nPixelHits",&seed_nPixelHits);
    tree->Branch("seed_nHits",&seed_nHits);
    tree->Branch("seed_jetAxisDistance",&seed_jetAxisDistance);
    tree->Branch("seed_jetAxisDlength",&seed_jetAxisDlength);
    tree->Branch("seed_ClosestJet_distance", &seed_ClosestJet_dR);
    
    tree->Branch("seed_MC_pt",&seed_MC_pt);
    tree->Branch("seed_MC_eta",&seed_MC_eta);
    tree->Branch("seed_MC_phi",&seed_MC_phi);
    tree->Branch("seed_MC_mass",&seed_MC_mass);
    tree->Branch("seed_MC_dz",&seed_MC_dz);
    tree->Branch("seed_MC_dxy",&seed_MC_dxy);
    tree->Branch("seed_MC_MomFlavour",&seed_MC_MomFlavour);
    tree->Branch("seed_MC_MomPdgId",&seed_MC_MomPdgId);
    tree->Branch("seed_MC_BChain",&seed_MC_BChain);
    tree->Branch("seed_MC_DChain",&seed_MC_DChain);
    tree->Branch("seed_MC_vx", &seed_MC_vx);
    tree->Branch("seed_MC_vy", &seed_MC_vy);
    tree->Branch("seed_MC_vz", &seed_MC_vz);
    tree->Branch("seed_MC_pvd", &seed_MC_pvd);  
    
    
    tree->Branch("nearTracks_Nvtx", &nearTracks_Nvtx);
    tree->Branch("nearTracks_nTracks", &nearTracks_nTracks);
    tree->Branch("nearTracks_pt", &nearTracks_pt);
    tree->Branch("nearTracks_eta", &nearTracks_eta);
    tree->Branch("nearTracks_phi", &nearTracks_phi);
    tree->Branch("nearTracks_mass", &nearTracks_mass);
    tree->Branch("nearTracks_dz", &nearTracks_dz);
    tree->Branch("nearTracks_dxy", &nearTracks_dxy);
    tree->Branch("nearTracks_3D_ip", &nearTracks_3D_ip);
    tree->Branch("nearTracks_3D_sip", &nearTracks_3D_sip);
    tree->Branch("nearTracks_2D_ip", &nearTracks_2D_ip);
    tree->Branch("nearTracks_2D_sip", &nearTracks_2D_sip);
//     tree->Branch("nearTracks_3D_signedIp", &nearTracks_3D_signedIp);
//     tree->Branch("nearTracks_3D_signedSip", &nearTracks_3D_signedSip);
//     tree->Branch("nearTracks_2D_signedIp", &nearTracks_2D_signedIp);
//     tree->Branch("nearTracks_2D_signedSip", &nearTracks_2D_signedSip);

    tree->Branch("nearTracks_PCAdist", &nearTracks_PCAdist);
    tree->Branch("nearTracks_PCAdsig", &nearTracks_PCAdsig);
    
    tree->Branch("nearTracks_PCAonSeed_x", &nearTracks_PCAonSeed_x);
    tree->Branch("nearTracks_PCAonSeed_y", &nearTracks_PCAonSeed_y);
    tree->Branch("nearTracks_PCAonSeed_z", &nearTracks_PCAonSeed_z);

    tree->Branch("nearTracks_PCAonSeed_xerr", &nearTracks_PCAonSeed_xerr);
    tree->Branch("nearTracks_PCAonSeed_yerr", &nearTracks_PCAonSeed_yerr);
    tree->Branch("nearTracks_PCAonSeed_zerr", &nearTracks_PCAonSeed_zerr);

    tree->Branch("nearTracks_PCAonTrack_x", &nearTracks_PCAonTrack_x);
    tree->Branch("nearTracks_PCAonTrack_y", &nearTracks_PCAonTrack_y);
    tree->Branch("nearTracks_PCAonTrack_z", &nearTracks_PCAonTrack_z);

    tree->Branch("nearTracks_PCAonTrack_xerr", &nearTracks_PCAonTrack_xerr);
    tree->Branch("nearTracks_PCAonTrack_yerr", &nearTracks_PCAonTrack_yerr);
    tree->Branch("nearTracks_PCAonTrack_zerr", &nearTracks_PCAonTrack_zerr); 

    tree->Branch("nearTracks_dotprodTrack", &nearTracks_dotprodTrack);
    tree->Branch("nearTracks_dotprodSeed", &nearTracks_dotprodSeed);
    tree->Branch("nearTracks_dotprodTrackSeed2D", &nearTracks_dotprodTrackSeed2D);
    tree->Branch("nearTracks_dotprodTrackSeed3D", &nearTracks_dotprodTrackSeed3D);
    tree->Branch("nearTracks_dotprodTrackSeedVectors2D", &nearTracks_dotprodTrackSeedVectors2D);
    tree->Branch("nearTracks_dotprodTrackSeedVectors3D", &nearTracks_dotprodTrackSeedVectors3D);
    
    tree->Branch("nearTracks_PCAonSeed_pvd", &nearTracks_PCAonSeed_pvd);
    tree->Branch("nearTracks_PCAonTrack_pvd", &nearTracks_PCAonTrack_pvd);
    
    tree->Branch("nearTracks_MC_pt", &nearTracks_MC_pt);
    tree->Branch("nearTracks_MC_eta", &nearTracks_MC_eta);
    tree->Branch("nearTracks_MC_phi", &nearTracks_MC_phi);
    tree->Branch("nearTracks_MC_dz", &nearTracks_MC_dz);
    tree->Branch("nearTracks_MC_dxy", &nearTracks_MC_dxy);
    tree->Branch("nearTracks_MC_MomFlavour", &nearTracks_MC_MomFlavour);
    tree->Branch("nearTracks_MC_MomPdgId", &nearTracks_MC_MomPdgId);
    tree->Branch("nearTracks_MC_BChain", &nearTracks_MC_BChain);
    tree->Branch("nearTracks_MC_DChain", &nearTracks_MC_DChain);

    tree->Branch("nearTracks_MC_Track_vx", &nearTracks_MC_Track_vx);
    tree->Branch("nearTracks_MC_Track_vy", &nearTracks_MC_Track_vy);
    tree->Branch("nearTracks_MC_Track_vz", &nearTracks_MC_Track_vz);
    
    tree->Branch("nearTracks_MC_fromSeedVtx", &nearTracks_MC_fromSeedVtx);
    tree->Branch("nearTracks_MC_fromSeedChain", &nearTracks_MC_fromSeedChain);
    tree->Branch("nearTracks_MC_pvd", &nearTracks_MC_pvd);
    tree->Branch("Match_chisquare", &Match_chisquare);
    
    tree->Branch("minChiSquare_genP", &minChiSquare_genP);

    CandidateToken = consumes<edm::View<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    token_primaryVertex = consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    token_beamSpot = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
    prunedGenParticleToken = consumes<std::vector< reco::GenParticle > >(edm::InputTag("prunedGenParticles"));
    genParticleToken = consumes<std::vector< pat::PackedGenParticle > >(edm::InputTag("packedGenParticles"));
    jetsToken = consumes<std::vector<pat::Jet> >(edm::InputTag("slimmedJets"));
//     jetToken = consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("jetMCSrc"));
    

}


AnalyzerSignedIP_MINIAOD::~AnalyzerSignedIP_MINIAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
AnalyzerSignedIP_MINIAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    evt=iEvent.id().event();
    lumi=iEvent.id().luminosityBlock();
    run=iEvent.id().run();
    
    jet_pt.clear();
    jet_px.clear();
    jet_py.clear();
    jet_pz.clear();
    jet_eta.clear();
    jet_phi.clear();
    jet_mass.clear();
    jet_flavour.clear();
    jet_seedIndex.clear();

    seed_pt.clear();
    seed_eta.clear();
    seed_phi.clear();
    seed_mass.clear();
    seed_dz.clear();
    seed_dxy.clear();
    seed_3D_ip.clear();
    seed_3D_sip.clear();
    seed_2D_ip.clear();
    seed_2D_sip.clear();
    seed_3D_signedIp.clear();
    seed_3D_signedSip.clear();
    seed_2D_signedIp.clear();
    seed_2D_signedSip.clear();
    seed_JetMatch.clear();
    seed_chi2reduced.clear();
    seed_nPixelHits.clear();
    seed_nHits.clear();
    seed_jetAxisDistance.clear();
    seed_jetAxisDlength.clear();    
    seed_ClosestJet_dR.clear();
    
    seed_MC_MomFlavour.clear();
    seed_MC_MomPdgId.clear();
    seed_MC_BChain.clear();
    seed_MC_DChain.clear();
    
    nearTracks_nTracks.clear();
    
    nearTracks_Nvtx.clear();
    nearTracks_pt.clear();
    nearTracks_eta.clear();
    nearTracks_phi.clear();
    nearTracks_dz.clear();
    nearTracks_dxy.clear();
    nearTracks_mass.clear();
    nearTracks_3D_ip.clear();
    nearTracks_3D_sip.clear();
    nearTracks_2D_ip.clear();
    nearTracks_2D_sip.clear();
//     nearTracks_3D_signedIp.clear(); //signed only in jet
//     nearTracks_3D_signedSip.clear();
//     nearTracks_2D_signedIp.clear();
//     nearTracks_2D_signedSip.clear();
    nearTracks_PCAdist.clear();
    nearTracks_PCAdsig.clear();
    nearTracks_PCAonSeed_x.clear();
    nearTracks_PCAonSeed_y.clear();
    nearTracks_PCAonSeed_z.clear();
    nearTracks_PCAonSeed_xerr.clear();
    nearTracks_PCAonSeed_yerr.clear();
    nearTracks_PCAonSeed_zerr.clear();
    nearTracks_PCAonTrack_x.clear();
    nearTracks_PCAonTrack_y.clear();
    nearTracks_PCAonTrack_z.clear();
    nearTracks_PCAonTrack_xerr.clear();
    nearTracks_PCAonTrack_yerr.clear();
    nearTracks_PCAonTrack_zerr.clear();
    nearTracks_dotprodTrack.clear();
    nearTracks_dotprodSeed.clear();
    nearTracks_dotprodTrackSeed2D.clear();
    nearTracks_dotprodTrackSeed3D.clear();
    nearTracks_dotprodTrackSeedVectors2D.clear();
    nearTracks_dotprodTrackSeedVectors3D.clear();
    
    
    nearTracks_PCAonSeed_pvd.clear();
    nearTracks_PCAonTrack_pvd.clear();
    
    //vars from matching genParticles if any
    nearTracks_MC_pt.clear();
    nearTracks_MC_eta.clear();
    nearTracks_MC_phi.clear();
    nearTracks_MC_dz.clear();
    nearTracks_MC_dxy.clear();

    nearTracks_MC_MomFlavour.clear();
    nearTracks_MC_MomPdgId.clear();
    nearTracks_MC_BChain.clear();
    nearTracks_MC_DChain.clear();

    seed_MC_vx.clear();
    seed_MC_vy.clear();
    seed_MC_vz.clear();
    seed_MC_pvd.clear();
    nearTracks_MC_Track_vx.clear();
    nearTracks_MC_Track_vy.clear();
    nearTracks_MC_Track_vz.clear();
    
    nearTracks_MC_fromSeedVtx.clear();
    nearTracks_MC_fromSeedChain.clear();
    nearTracks_MC_pvd.clear();
    
    seed_MC_pt.clear();
    seed_MC_eta.clear();
    seed_MC_phi.clear();
    seed_MC_mass.clear();
    seed_MC_dz.clear();
    seed_MC_dxy.clear();
    //end vars from matching genParticles if any
    
    Match_chisquare.clear();


    using namespace edm;
    using namespace reco;
    
    
       //1st way to get objects
    Handle<std::vector<pat::PackedGenParticle> > genParticlesCollection;
    iEvent.getByToken(genParticleToken, genParticlesCollection);
    std::vector<pat::PackedGenParticle>  genParticles = *genParticlesCollection.product();
    
    Handle<std::vector<reco::GenParticle> > prunedGenParticlesCollection;
    iEvent.getByToken(prunedGenParticleToken, prunedGenParticlesCollection);
    std::vector<reco::GenParticle>  prunedGenParticles = *prunedGenParticlesCollection.product();
    
    Handle<std::vector<pat::Jet> > jetsCollection;
    iEvent.getByToken(jetsToken, jetsCollection);
    std::vector<pat::Jet>  ak4jets = *jetsCollection.product();
    
    
    //2nd way to get objects
    Handle<edm::View<pat::PackedCandidate> > tracks;
    iEvent.getByToken(CandidateToken, tracks);

    edm::Handle<reco::VertexCollection > primaryVertices;
    iEvent.getByToken(token_primaryVertex, primaryVertices);
    
//     edm::Handle<reco::JetFlavourInfoMatchingCollection> jetMC;    
//     iEvent.getByToken(jetToken, jetMC);
    
    
   
//   std::cout << tracks->size() << std::endl;
   
   if(primaryVertices->size()!=0){
        const reco::Vertex &pv = (*primaryVertices)[0];
        GlobalPoint pvp(pv.x(),pv.y(),pv.z());
//        std::cout << pv.x() << std::endl;
        
        gen_pv_x=genParticles[2].vx();
        gen_pv_y=genParticles[2].vy();
        gen_pv_z=genParticles[2].vz();
        
        GlobalPoint gpvp(gen_pv_x,gen_pv_y,gen_pv_z);
        
        pv_x=pv.x();
        pv_y=pv.y();
        pv_z=pv.z();
        
        std::cout << "check the primary  " << std::endl;
        std::cout <<" gen -->" << gen_pv_x << ", "<< gen_pv_y << ", "<< gen_pv_z << "   " ;
        std::cout <<" reco -->" << pv_x << ", "<< pv_y << ", "<< pv_z <<std::endl;
        
        
        //before track etc: do jetsCollection
        for (std::vector<pat::Jet>::const_iterator iter = ak4jets.begin(); iter != ak4jets.end(); ++iter) {
            unsigned int fl = std::abs(iter->partonFlavour());

            jet_flavour.push_back(fl);
            jet_pt.push_back(iter->pt());
            jet_eta.push_back(iter->eta());
            jet_phi.push_back(iter->phi());
            jet_mass.push_back(iter->mass());
            jet_px.push_back(iter->px());
            jet_py.push_back(iter->py());
            jet_pz.push_back(iter->pz());
            
            std::cout<< "pt: " <<iter->pt()<< "eta: " << iter->eta()<< "phi: " << iter->phi()<< "m: " << iter->mass() << std::endl;
            
            
            
        }
        
        
        edm::ESHandle<TransientTrackBuilder> trackBuilder;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);

        edm::Handle<BeamSpot> beamSpot;
        iEvent.getByToken(token_beamSpot,beamSpot);
        
        
        std::cout << "Transient tracks" << std::endl;
        std::vector<TransientTrack> selectedTracks;
        
        std::cout << "Masses" << std::endl;
        std::vector<float> masses;
        
        //if primary: build transient tracks form packedCandidates
        
        for(typename edm::View<pat::PackedCandidate>::const_iterator track = tracks->begin(); track != tracks->end(); ++track) {
        
            unsigned int k=track - tracks->begin();
             if ((*tracks)[k].pt()>10) std::cout << "all tracks Candidate pt"  << (*tracks)[k].pt() << " eta  " 
            <<(*tracks)[k].eta() << "  phi " <<(*tracks)[k].phi() << " distance DZ " << std::fabs(pv_z-(*tracks)[k].vz()) << " bestTrack  " <<  (*tracks)[k].bestTrack() << std::endl;
            if ((*tracks)[k].pt()>10) std::cout << "all tracks Candidate vtx"  <<  "  z  " << pv_z << "  " << (*tracks)[k].vz() 
            << "  x  "  << pv_x << "  " << (*tracks)[k].vx() << "  y  " << pv_y << "  " << (*tracks)[k].vy() <<std::endl;

            if ((*tracks)[k].pt()>10 and (*tracks)[k].bestTrack() != 0){
            std::cout << " " << trackBuilder->build(tracks->ptrAt(k)).track().pt() << " " << trackBuilder->build(tracks->ptrAt(k)).track().eta() <<  " " << trackBuilder->build(tracks->ptrAt(k)).track().phi()   << std::endl;
         std::cout <<" vertex for the track  " <<  trackBuilder->build(tracks->ptrAt(k)).track().vx() <<  " "<< trackBuilder->build(tracks->ptrAt(k)).track().vy() <<  " " << trackBuilder->build(tracks->ptrAt(k)).track().vz() << std::endl;
         std::cout <<" reco primary -->" << pv_x << ", "<< pv_y << ", "<< pv_z <<std::endl;
         std::cout<< " dxy "<< trackBuilder->build(tracks->ptrAt(k)).track().dxy() << " dz " << trackBuilder->build(tracks->ptrAt(k)).track().dz() << " wrt to PV " << trackBuilder->build(tracks->ptrAt(k)).track().dxy(pv.position())<< "  " << trackBuilder->build(tracks->ptrAt(k)).track().dz(pv.position())<<std::endl;
            }
            
            
            /*Questi tagli ip tag*/
            
//             if (track.pt() > m_cutMinPt &&
//             track.hitPattern().numberOfValidHits() >= m_cutTotalHits &&         // min num tracker hits
//             track.hitPattern().numberOfValidPixelHits() >= m_cutPixelHits &&
//             track.normalizedChi2() < m_cutMaxChiSquared &&
//             std::abs(track.dxy(pv->position())) < m_cutMaxTIP &&
//             std::abs(track.dz(pv->position())) < m_cutMaxLIP)
            
//                 if( fabs(impactParameters[i].distanceToJetAxis.value()) < m_cutMaxDistToAxis  &&        // distance to JetAxis 0.07
// 	  (impactParameters[i].closestToJetAxis - pv).mag() < m_cutMaxDecayLen  &&      // max decay len 5

// 	  ) 
            if((*tracks)[k].bestTrack() != 0 &&  (*tracks)[k].pt()>0.5) {
            if (std::fabs(pv_z-trackBuilder->build(tracks->ptrAt(k)).track().vz())<0.5){
            selectedTracks.push_back(trackBuilder->build(tracks->ptrAt(k)));
            masses.push_back(tracks->ptrAt(k)->mass());}}
            }
        
        std::cout << "function entering .." << std::endl;
//         tracking_eff=trackEff(genParticles, selectedTracks, BChain_eff, DChain_eff, B_eff, low_pt_eff, high_pt_eff, minChiSquare_genP);
        
        //selected tracks are all the good tracks, now match to the gen particles for later reference
//        genMap(genParticles, selectedTracks, genP_index, Match_chisquare );
//         genLevelMM(genParticles, selectedTracks, genP_index, Match_chisquare );
        
        match_eff=0;
        
        for(std::vector<TransientTrack>::const_iterator it = selectedTracks.begin(); it != selectedTracks.end(); it++){
        
        int index= it - selectedTracks.begin();
        
        std::cout << " see how the matching went  : new particle" << "  eff:  " << match_eff << "/" << selectedTracks.size() <<  "  at index " << index << std::endl;
        
//        std::cout << " " << it->track().pt() << " " << it->track().eta() <<  " " << it->track().phi()   << std::endl;
        std::cout << " " << it->track().pt() << " " << it->track().eta() <<  " " << it->track().phi()   << std::endl;
        std::cout <<" vertex for the track  " <<  it->track().vx() <<  " "<< it->track().vy() <<  " " << it->track().vz() << std::endl;
        std::cout <<" reco primary -->" << pv_x << ", "<< pv_y << ", "<< pv_z <<std::endl;
//         std::cout << " " << index << " == " << it - selectedTracks.begin() <<  " " << genP_index[index]   << "  " << genP_index.size() <<std::endl;
        
        
        
        int g_index=0;//genP_index[index];
        if (g_index>0){
        match_eff=match_eff+1;
        
        std::cout << " " << genParticles[g_index].pt() << " " << genParticles[g_index].eta() <<  " " << genParticles[g_index].phi()   << std::endl;
        }
        else
        {
        std::cout << " no particle match" << std::endl;
        }
        
        }
        match_eff=match_eff/selectedTracks.size();
        
        
        
        std::vector<TransientTrack> seeds;        
        //if primary: build transient tracks form packedCandidates
        std::cout << "seeds" << std::endl;

        for(std::vector<TransientTrack>::const_iterator it = selectedTracks.begin(); it != selectedTracks.end(); it++){
        std::vector<int> jet_index_array = seedToAllJets30Matching(*it, jet_eta,jet_phi, jet_pt); 
        
        for (unsigned j_i=0; j_i < jet_index_array.size(); j_i++) {
        
            std::cout << "array of indexes " << jet_index_array[j_i] << std::endl;
        int seed_jet_match_i=jet_index_array[j_i];

//         int seed_jet_match_i=seedToJetMatching(*it, jet_eta,jet_phi); 
        std::pair<bool,Measurement1D> ip = IPTools::absoluteImpactParameter3D(*it, pv);        
        std::pair<bool,Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(*it, pv);
        
            std::cout << ip.second.significance() << " " <<  ( it->track().normalizedChi2()<5.) <<  (seed_jet_match_i>=0 )<< (std::fabs(it->track().dxy(pv.position())) < 2 )<<  (std::fabs(it->track().dz(pv.position())) < 17) << std::endl;
        
            if(ip.first && ip.second.value() >= 0.0 && ip.second.significance() >= 1.0 &&
            ip.second.value() <= max3DIPValue && ip.second.significance() <= max3DIPSignificance &&
            it->track().normalizedChi2()<5. &&  seed_jet_match_i>=0 && std::fabs(it->track().dxy(pv.position())) < 2 &&
            std::fabs(it->track().dz(pv.position())) < 17
            /*&&
            
            std::abs(track.dxy(pv->position())) < m_cutMaxTIP &&
           std::abs(track.dz(pv->position())) < m_cutMaxLIP)
            ip2d.second.value() >= min3DIPValue && ip2d.second.significance() >= min3DIPSignificance &&
            ip2d.second.value() <= max3DIPValue && ip2d.second.significance() <= max3DIPSignificance*/
            
            )
            { 
                
            GlobalVector direction(jet_px[seed_jet_match_i], jet_py[seed_jet_match_i], jet_pz[seed_jet_match_i]);

            std::pair<bool,Measurement1D> ipSigned = IPTools::signedImpactParameter3D(*it,direction, pv);        
            std::pair<bool,Measurement1D> ip2dSigned = IPTools::signedTransverseImpactParameter(*it,direction, pv);
            
            std::cout << "new seed : 2d+3d sign." << ip2d.second.significance() << "  " << ip.second.significance() <<"  "<< std::endl;
            std::cout << "new seed : in jet : " << jet_pt[seed_jet_match_i] <<"  "<<jet_eta[seed_jet_match_i]<<"  "<<jet_phi[seed_jet_match_i]<<"  "<< std::endl;
//              std::cout << "new seed " <<  it-selectedTracks.begin() << " ref " << "qualcosa key"  << " " << ip.second.value() <<std::endl;
//              std::cout << ip.second.significance() << " " << it->track().hitPattern().trackerLayersWithMeasurement() <<std::endl;
//              std::cout << " " << it->track().pt() << " " << it->track().eta() << std::endl;
//             int index= it - selectedTracks.begin();
            int g_index=0;//genP_index[index];
            if (g_index>0){
//                match_eff=match_eff+1;
                std::cout << " " << it->track().pt() << " " << it->track().eta() <<  " " << it->track().phi()   << std::endl;
                std::cout << " " << genParticles[g_index].pt() << " " << genParticles[g_index].eta() <<  " " << genParticles[g_index].phi()   << std::endl;
                }
            else
            std::cout << "no match for seed " << it->track().pt() << " " << it->track().eta() <<  " " << it->track().phi()   << std::endl;
//            std::cout << "new seed " <<  it-selectedTracks.begin() << " ref " << it->trackBaseRef().key()  << " " << ip.second.value() << " " << ip.second.significance() << " " << it->track().hitPattern().trackerLayersWithMeasurement() << " " << it->track().pt() << " " << it->track().eta() << std::endl;
            seeds.push_back(*it);
            
            seed_pt.push_back(it->track().pt());
            seed_eta.push_back(it->track().eta());
            seed_phi.push_back(it->track().phi());
            seed_mass.push_back(masses[it-selectedTracks.begin()]);
            
            seed_dz.push_back(it->track().dz());
            seed_dxy.push_back(it->track().dxy());
            seed_3D_ip.push_back(ip.second.value());
            seed_3D_sip.push_back(ip.second.significance());
            seed_2D_ip.push_back(ip2d.second.value());
            seed_2D_sip.push_back(ip2d.second.significance());
            seed_3D_signedIp.push_back(ipSigned.second.value());
            seed_3D_signedSip.push_back(ipSigned.second.significance());
            seed_2D_signedIp.push_back(ip2dSigned.second.value());
            seed_2D_signedSip.push_back(ip2dSigned.second.significance());
            
                       
            seed_JetMatch.push_back(seed_jet_match_i);
            
            pair = ClosestJet(*it, jet_eta,jet_phi);
            seed_ClosestJet_dR.push_back(pair.second);
            
            seed_chi2reduced.push_back(it->track().normalizedChi2());
            seed_nPixelHits.push_back(it->track().hitPattern().numberOfValidPixelHits());
            seed_nHits.push_back(it->track().hitPattern().numberOfValidHits());  
            
            if(seed_jet_match_i>=0){
                
            //GlobalVector direction(jet_px[seed_jet_match_i], jet_py[seed_jet_match_i], jet_pz[seed_jet_match_i]); 
            std::pair<double, Measurement1D> jet_distance =IPTools::jetTrackDistance(*it, direction, pv);
            seed_jetAxisDistance.push_back(std::fabs(jet_distance.second.value())); 
            std::cout<<jet_distance.second.value()<<"      distance"<<jet_distance.second.significance()<<std::endl;
            
            TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(it->impactPointState(),pv, direction,it->field());
            if (closest.isValid()) seed_jetAxisDlength.push_back((closest.globalPosition() - pvp).mag()); 
            else seed_jetAxisDlength.push_back(-99); 
            
            }
            
           else{seed_jetAxisDistance.push_back(-100.); seed_jetAxisDlength.push_back(-100.);  }
          
//             jet_seedIndex.push_back();

            // MC Truth for seeds only in RECO
            
            /*if (genP_index[it - selectedTracks.begin()]>0){
            int g_i=genP_index[it - selectedTracks.begin()];
            seed_MC_pt.push_back(genParticles[g_i].pt());
            seed_MC_eta.push_back(genParticles[g_i].eta());
            seed_MC_phi.push_back(genParticles[g_i].phi());
            seed_MC_mass.push_back(genParticles[g_i].mass());
            seed_MC_dxy.push_back((-genParticles[g_i].vx() * genParticles[g_i].py() + genParticles[g_i].vy() * genParticles[g_i].px()) / genParticles[g_i].pt());
            seed_MC_dz.push_back(genParticles[g_i].vz() - (genParticles[g_i].vx() * genParticles[g_i].px() + genParticles[g_i].vy() * genParticles[g_i].py()) / genParticles[g_i].pt() * (genParticles[g_i].pz() / genParticles[g_i].pt()));
            seed_MC_vx.push_back(genParticles[g_i].vx());
            seed_MC_vy.push_back(genParticles[g_i].vy());
            seed_MC_vz.push_back(genParticles[g_i].vz());
            GlobalPoint seed_vtx(genParticles[g_i].vx(), genParticles[g_i].vy(), genParticles[g_i].vz());
            seed_MC_pvd.push_back((seed_vtx-gpvp).mag()); //"write lambda var for gen particle");
            

            //mother id
            const reco::Candidate* mother=genParticles[g_i].mother(0);

            if (genParticles[g_i].numberOfMothers()>0 )
            {             
            id_mom = abs(mother->pdgId());
            id_mom=std::max((id_mom/1000) % 10,(id_mom/100) % 10);}
            else
            id_mom = -1;

            seed_MC_MomFlavour.push_back(id_mom);
            
            if (genParticles[g_i].numberOfMothers()>0 )
            seed_MC_MomPdgId.push_back(mother->pdgId());
            else
            seed_MC_MomPdgId.push_back(39);
            //bchain
            
            b_found=0;
            d_found=0;
            while (mother!=NULL && b_found==0 && mother->status()<3)
            {
            id_mom = abs(mother->pdgId());
            std::cout<<"seed chain verify "<< mother->pdgId() <<"  status   "<< mother->status() <<std::endl;
            id_mom=std::max((id_mom/1000) % 10,(id_mom/100) % 10);
            if (id_mom==5)
            b_found=1;
            if (id_mom==4)
            d_found=1;
            mother=mother->mother(0);
            }

            seed_MC_BChain.push_back(b_found);
            seed_MC_DChain.push_back(d_found*(!b_found));
            
            bd_i = std::max(b_found,d_found);
            
            }
            else{
            
            seed_MC_pt.push_back(-100);
            seed_MC_eta.push_back(-100);
            seed_MC_phi.push_back(-100);
            seed_MC_mass.push_back(-100);
            seed_MC_dxy.push_back(-100);
            seed_MC_dz.push_back(-100);
            seed_MC_vx.push_back(-100);
            seed_MC_vy.push_back(-100);
            seed_MC_vz.push_back(-100);
            seed_MC_pvd.push_back(-100);
            seed_MC_MomFlavour.push_back(-100);
            seed_MC_BChain.push_back(-100);
            seed_MC_DChain.push_back(-100);
            seed_MC_MomPdgId.push_back(39);
            bd_i = 0;
            }*/

            
//            int myseed = genLevel(genParticles, *it, seed_MC_MomFlavour, seed_MC_BChain, Match_chisquare);
//             int s_i = genP_index[it - selectedTracks.begin()];
            nearTracks.clear();
            std::cout<< "nearTracks" << std::endl;

            for(std::vector<reco::TransientTrack>::const_iterator tt = selectedTracks.begin();tt!=selectedTracks.end(); ++tt )   {
                VertexDistance3D distanceComputer;
                TwoTrackMinimumDistance dist;
                if(*tt==*it) continue;
                if(std::fabs(pv_z-tt->track().vz())>0.1) continue;
                if(dist.calculate(tt->impactPointState(),it->impactPointState())) {
                    GlobalPoint ttPoint          = dist.points().first;
//                    std::cout << ttPoint << std::endl;
                    GlobalError ttPointErr       = tt->impactPointState().cartesianError().position();
                    GlobalPoint seedPosition     = dist.points().second;
                    GlobalError seedPositionErr  = it->impactPointState().cartesianError().position();
                    Measurement1D m = distanceComputer.distance(VertexState(seedPosition,seedPositionErr), VertexState(ttPoint, ttPointErr));
                    GlobalPoint cp(dist.crossingPoint()); 
//                    std::cout << ttPointErr.cxx() << "  errore come fatto  " << std::endl;

                    float PCAseedFromPV =  (dist.points().second-pvp).mag();
                    float PCAtrackFromPV =  (dist.points().first-pvp).mag();
//                    std::cout << distanceFromPV << std::endl;
                    
                    float distance = dist.distance();
//                    float significance = m.significance();
//                    std::cout << distance << std::endl;
//                    std::cout << m.value() << "   " << m.significance() << " m 1d measurement " << std::endl;
                    
                    
                    GlobalVector trackDir2D(tt->impactPointState().globalDirection().x(),tt->impactPointState().globalDirection().y(),0.); 
                    GlobalVector seedDir2D(it->impactPointState().globalDirection().x(),it->impactPointState().globalDirection().y(),0.); 
                    GlobalVector trackPCADir2D(dist.points().first.x()-pvp.x(),dist.points().first.y()-pvp.y(),0.); 
                    GlobalVector seedPCADir2D(dist.points().second.x()-pvp.x(),dist.points().second.y()-pvp.y(),0.); 
//                    //SK:UNUSED//    float dotprodTrackSeed2D = trackDir2D.unit().dot(seedDir2D.unit());

                    float dotprodTrack = (dist.points().first-pvp).unit().dot(tt->impactPointState().globalDirection().unit());
                    float dotprodSeed = (dist.points().second-pvp).unit().dot(it->impactPointState().globalDirection().unit());
                    
                    float dotprodTrackSeed2D = trackDir2D.unit().dot(seedDir2D.unit());
                    float dotprodTrackSeed3D = it->impactPointState().globalDirection().unit().dot(tt->impactPointState().globalDirection().unit());
                    float dotprodTrackSeed2DV = trackPCADir2D.unit().dot(seedPCADir2D.unit());
                    float dotprodTrackSeed3DV = (dist.points().second-pvp).unit().dot((dist.points().first-pvp).unit());
                    
                    std::pair<bool,Measurement1D> t_ip = IPTools::absoluteImpactParameter3D(*tt,pv);        
                    std::pair<bool,Measurement1D> t_ip2d = IPTools::absoluteTransverseImpactParameter(*tt,pv);
                   
                    
                    myTrack.set_values(tt->track().pt(), tt->track().eta(), tt->track().phi(),  tt->track().dz(), tt->track().dxy(), distance, 
                    m.significance(), /*put significance etc*/ 
                    
                    seedPosition.x(), seedPosition.y(), seedPosition.z(), seedPositionErr.cxx(), seedPositionErr.cyy(), seedPositionErr.czz(),
                     ttPoint.x(),  ttPoint.y(),  ttPoint.z(),  ttPointErr.cxx(),  ttPointErr.cyy(),  ttPointErr.czz(),
                    dotprodTrack, dotprodSeed
                    
                    );
                    myTrack.set_index(-1);
                    myTrack.set_distances(PCAseedFromPV, PCAtrackFromPV);
                    myTrack.set_vars(masses[tt-selectedTracks.begin()],t_ip2d.second.value() , t_ip2d.second.significance(),
                    t_ip.second.value() , t_ip.second.significance(), dotprodTrackSeed2D, dotprodTrackSeed3D, dotprodTrackSeed2DV, dotprodTrackSeed3DV );
                    nearTracks.push_back(myTrack);
                   
//                    float w = distanceFromPV*distanceFromPV/(pvDistance*distance);
//                    bool selected = (m.significance() < clusterMaxSignificance && 
//                    dotprodSeed > clusterMinAngleCosine && //Angles between PV-PCAonSeed vectors and seed directions
//                    dotprodTrack > clusterMinAngleCosine && //Angles between PV-PCAonTrack vectors and track directions
//                    //                    dotprodTrackSeed2D > clusterMinAngleCosine && //Angle between track and seed
//                    //      distance*clusterScale*tracks.size() < (distanceFromPV+pvDistance)*(distanceFromPV+pvDistance)/pvDistance && // cut scaling with track density
//                    distance*distanceRatio < distanceFromPV && // cut scaling with track density
//                    distance < clusterMaxDistance);  // absolute distance cut

//                    #ifdef VTXDEBUG
//                    std::cout << tt->trackBaseRef().key() << " :  " << (selected?"+":" ")<< " " << m.significance() << " < " << clusterMaxSignificance <<  " &&  " << 
//                    dotprodSeed  << " > " <<  clusterMinAngleCosine << "  && " << 
//                    dotprodTrack  << " > " <<  clusterMinAngleCosine << "  && " << 
//                    dotprodTrackSeed2D  << " > " <<  clusterMinAngleCosine << "  &&  "  << 
//                    distance*distanceRatio  << " < " <<  distanceFromPV << "  crossingtoPV: " << distanceFromPV << " dis*scal " <<  distance*distanceRatio << "  <  " << distanceFromPV << " dist: " << distance << " < " << clusterMaxDistance <<  std::endl; // cut scaling with track density
//                    #endif           
//                    if(selected)
//                    {
//                    result.push_back(*tt);
//                    seedingPoint = GlobalPoint(cp.x()*w+seedingPoint.x(),cp.y()*w+seedingPoint.y(),cp.z()*w+seedingPoint.z());  
//                    sumWeights+=w; 
//                    }
                }
            }
//            
           for(unsigned int i=0; i< nearTracks.size(); i++) {
           std::cout << nearTracks[i].dist << "  "; 
           }
           std::cout <<  "  " << std::endl; 
           
           for(unsigned int i=0; i< nearTracks.size(); i++) {
           std::cout << nearTracks[i].dsig << "  "; 
           }
           std::cout <<  "  " << std::endl; 
//            
           std::sort (nearTracks.begin(), nearTracks.end(), sortfunction2());
//            
//            for(unsigned int i=0; i< nearTracks.size(); i++) {
//            std::cout << nearTracks[i].dist << "  "; 
//            }
            std::cout <<  "  " << std::endl; 
            std::cout <<  "  " << std::endl; 
            nearTracks.resize(20);
            for(unsigned int i=0; i< nearTracks.size(); i++) {
            std::cout << nearTracks[i].dist << "  "; 
            }
            std::cout <<  "  " << std::endl; 
            for(unsigned int i=0; i< nearTracks.size(); i++) {
            std::cout << nearTracks[i].pt << "  ";}
            std::cout <<  "  " << std::endl; 
            
            
            for(unsigned int i=0; i< nearTracks.size(); i++) {
            nearTracks_nTracks.push_back(nearTracks.size());
            nearTracks_Nvtx.push_back(seeds.size()-1);
            nearTracks_pt.push_back(nearTracks[i].pt);
            nearTracks_eta.push_back(nearTracks[i].eta);
            nearTracks_phi.push_back(nearTracks[i].phi);
            nearTracks_dz.push_back(nearTracks[i].dz);
            nearTracks_dxy.push_back(nearTracks[i].dxy);
            nearTracks_mass.push_back(nearTracks[i].mass);
            nearTracks_3D_ip.push_back(nearTracks[i].t3Dip);
            nearTracks_3D_sip.push_back(nearTracks[i].t3Dsip);
            nearTracks_2D_ip.push_back(nearTracks[i].t2Dip);
            nearTracks_2D_sip.push_back(nearTracks[i].t2Dsip);
            nearTracks_PCAdist.push_back(nearTracks[i].dist);
            nearTracks_PCAdsig.push_back(nearTracks[i].dsig);
            nearTracks_PCAonSeed_x.push_back(nearTracks[i].PCA_sx);
            nearTracks_PCAonSeed_y.push_back(nearTracks[i].PCA_sy);
            nearTracks_PCAonSeed_z.push_back(nearTracks[i].PCA_sz);
            nearTracks_PCAonSeed_xerr.push_back(nearTracks[i].PCA_sxerr);
            nearTracks_PCAonSeed_yerr.push_back(nearTracks[i].PCA_syerr);
            nearTracks_PCAonSeed_zerr.push_back(nearTracks[i].PCA_szerr);
            nearTracks_PCAonTrack_x.push_back(nearTracks[i].PCA_tx);
            nearTracks_PCAonTrack_y.push_back(nearTracks[i].PCA_ty);
            nearTracks_PCAonTrack_z.push_back(nearTracks[i].PCA_tz);
            nearTracks_PCAonTrack_xerr.push_back(nearTracks[i].PCA_txerr);
            nearTracks_PCAonTrack_yerr.push_back(nearTracks[i].PCA_tyerr);
            nearTracks_PCAonTrack_zerr.push_back(nearTracks[i].PCA_tzerr);
            nearTracks_dotprodTrack.push_back(nearTracks[i].dotprodTrack);
            nearTracks_dotprodSeed.push_back(nearTracks[i].dotprodSeed);
            nearTracks_dotprodTrackSeed2D.push_back(nearTracks[i].dotprodTrackSeed2D);
            nearTracks_dotprodTrackSeed3D.push_back(nearTracks[i].dotprodTrackSeed3D);
            nearTracks_dotprodTrackSeedVectors2D.push_back(nearTracks[i].dotprodTrackSeedVectors2D);
            nearTracks_dotprodTrackSeedVectors3D.push_back(nearTracks[i].dotprodTrackSeedVectors3D);
            nearTracks_PCAonSeed_pvd.push_back(nearTracks[i].seedPCA_pv);
            nearTracks_PCAonTrack_pvd.push_back(nearTracks[i].trackPCA_pv);

            //MC Truth for Tracks (only in RECO)
                
            /*if (nearTracks[i].index>0){
                std::cout<<"track "<<i<< std::endl;
            int g_i=nearTracks[i].index;
            nearTracks_MC_pt.push_back(genParticles[g_i].pt());
            nearTracks_MC_eta.push_back(genParticles[g_i].eta());
            nearTracks_MC_phi.push_back(genParticles[g_i].phi());
            nearTracks_MC_dxy.push_back((-genParticles[g_i].vx() * genParticles[g_i].py() + genParticles[g_i].vy() * genParticles[g_i].px()) / genParticles[g_i].pt());
            nearTracks_MC_dz.push_back(genParticles[g_i].vz() - (genParticles[g_i].vx() * genParticles[g_i].px() + genParticles[g_i].vy() * genParticles[g_i].py()) / genParticles[g_i].pt() * (genParticles[g_i].pz() / genParticles[g_i].pt()));
            nearTracks_MC_Track_vx.push_back(genParticles[g_i].vx());
            nearTracks_MC_Track_vy.push_back(genParticles[g_i].vy());
            nearTracks_MC_Track_vz.push_back(genParticles[g_i].vz());
            nearTracks_MC_fromSeedVtx.push_back((genParticles[g_i].vx()==seed_MC_vx[seed_MC_vx.size() - 1])&&(genParticles[g_i].vy()==seed_MC_vy[seed_MC_vy.size() - 1])
            &&(genParticles[g_i].vz()==seed_MC_vz[seed_MC_vz.size() - 1]));            
                       
            GlobalPoint track_vtx(genParticles[g_i].vx(), genParticles[g_i].vy(), genParticles[g_i].vz());
            nearTracks_MC_pvd.push_back((track_vtx-gpvp).mag());
            const reco::Candidate* mother =genParticles[g_i].mother(0);

            //mother id
            std::cout<<"mothers size "<<mother << "  and statusss  " <<mother->status()<< "genpv   "<< gpvp << std::endl;
            if (genParticles[g_i].numberOfMothers()>0 )
            {id_mom = abs(mother->pdgId());
            id_mom=std::max((id_mom/1000) % 10,(id_mom/100) % 10);}
            else
            id_mom = -1;
            
            nearTracks_MC_MomFlavour.push_back(id_mom);
            
            if (genParticles[g_i].numberOfMothers()>0 )
            nearTracks_MC_MomPdgId.push_back(mother->pdgId());
            else
            nearTracks_MC_MomPdgId.push_back(39);
            
            //bchain
            
            b_found=0;
            d_found=0;
            while (mother!=NULL && b_found==0 && mother->status()<3)
            {
            id_mom = abs(mother->pdgId());
            std::cout<<"track "<<i<<" chain verify "<< mother->pdgId() <<"  status   "<< mother->status() <<std::endl;
            id_mom=std::max((id_mom/1000) % 10,(id_mom/100) % 10);
            if (id_mom==5)
            b_found=1;
            if (id_mom==4)
            d_found=1;
            std::cout<<"chain verify "<< b_found << "  " <<d_found <<"  "<<std::endl;
            mother=mother->mother(0);
            }

            nearTracks_MC_BChain.push_back(b_found);
            nearTracks_MC_DChain.push_back(d_found*(!b_found));
            
            
            
            //same chain
            int same_chain=0; 
            
            if(s_i>0&&g_i>0){
            if ((b_found || d_found)&&bd_i)
                
            {   
                //change to stop at first displaced decay instead o first hadron?
                std::vector<float> v1=firstHadronInChain(genParticles[g_i], mother);
                std::vector<float> v2=firstHadronInChain(genParticles[s_i], mother);
                same_chain=std::equal( v1.begin(), v1.end(), v2.begin() );
                std::cout << "verify bfound  "<< b_found << " dfound  "<< d_found << " same chain  "<< same_chain << "  bdi  "<< bd_i <<std::endl;
            }
            nearTracks_MC_fromSeedChain.push_back(same_chain);
            }
            else nearTracks_MC_fromSeedChain.push_back(-100); 
            
            
            
            
            
            }
            else{
            nearTracks_MC_pt.push_back(-100);
            nearTracks_MC_eta.push_back(-100);
            nearTracks_MC_phi.push_back(-100);
            nearTracks_MC_dxy.push_back(-100);
            nearTracks_MC_dz.push_back(-100);
            nearTracks_MC_Track_vx.push_back(-100);
            nearTracks_MC_Track_vy.push_back(-100);
            nearTracks_MC_Track_vz.push_back(-100);
            nearTracks_MC_MomFlavour.push_back(-100);
            nearTracks_MC_BChain.push_back(-100);
            nearTracks_MC_DChain.push_back(-100);
            nearTracks_MC_fromSeedVtx.push_back(-100);
            nearTracks_MC_fromSeedChain.push_back(-100);
            nearTracks_MC_MomPdgId.push_back(39);
            nearTracks_MC_pvd.push_back(-100);
            }*/
            }
            
            
            std::cout << "check the seed again before " << it->track().phi() << " " << it->track().pt() << " " << it->track().eta() << std::endl;

            
            }
            
            
            
            
            
        n_seed=seeds.size();
        }
        }
        }
   tree->Fill();



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}



int AnalyzerSignedIP_MINIAOD::funzionePippo (int vivaLaVita) {
return 1;
}


std::vector<float> AnalyzerSignedIP_MINIAOD::firstHadronInChain(pat::PackedGenParticle gParticle, const reco::Candidate* mother){
    std::vector<float> vtx;
    vtx.push_back(0);
    vtx.push_back(0);
    vtx.push_back(0);
    vtx.push_back(0);
    vtx.push_back(0);
    vtx.push_back(0);
    vtx.push_back(0);
    std::cout<< "init func" <<std::endl;
//     edm::RefVector<std::vector<pat::PackedGenParticle> > mothers;
    mother=gParticle.mother(0);
//     int b_found=0;
//     int d_found=0;
    int first_h=0;
    int id_mom;
    int id_before;
    std::cout<< "mothers" <<std::endl;
    while (gParticle.numberOfMothers()>0 && first_h==0)
    {
         std::cout<< "motherspdid  " << mother->pdgId()<< "  status  " << mother->status() << std::endl;
//          std::cout << "decay vertex" <<" vx  "<<mother->daughter(0)->vx()<<" vy  "<<mother->daughter(0)->vy()<<" vz  "<<mother->daughter(0)->vz()<<std::endl;
         std::cout << "first h" <<" p  "<<mother->pt()<<" eta  "<<mother->eta()<<" phi  "<<mother->phi()<<std::endl;
    id_mom = abs(mother->pdgId());
    id_mom=std::max((id_mom/1000) % 10,(id_mom/100) % 10);
    if (id_mom==5){
//     b_found=1;
         std::cout<< "mothersVpdid  " << mother->mother(0)->pdgId()<<std::endl;
    id_before=std::abs(mother->mother(0)->pdgId());
    id_before=std::max((id_before/1000) % 10,(id_before/100) % 10);
    if (id_before!=5){  
        first_h=1;
        std::cout << "first h" <<std::endl;
//         std::cout << "first h" <<" vx  "<<mother->daughter(0)->vx()<<" vy  "<<mother->daughter(0)->vy()<<" vz  "<<mother->daughter(0)->vz()<<std::endl;
        std::cout << "first h" <<" p  "<<mother->pt()<<" eta  "<<mother->eta()<<" phi  "<<mother->phi()<<std::endl;
        for (unsigned int o=0; o<mother->numberOfDaughters(); o++){
            std::cout<<"dau___"<<mother->daughter(o)->pdgId()<<"---"<<mother->daughter(o)->pt()<<" --  "<<mother->daughter(o)->eta()<<" --  "<<mother->daughter(o)->phi()<<std::endl;
        }
        vtx[0]=(mother->pt());
        vtx[1]=(mother->eta());
        vtx[2]=(mother->phi());
        vtx[3]=(mother->daughter(0)->vx());
        vtx[4]=(mother->daughter(0)->vy());
        vtx[5]=(mother->daughter(0)->vz());
        vtx[6]=(mother->pdgId());
    }
    }
    
    if (id_mom==4){
//     d_found=1;
        std::cout<< "mothersVpdid  " << mother->mother(0)->pdgId()<<std::endl;
    id_before=std::abs(mother->mother(0)->pdgId());
    id_before=std::max((id_before/1000) % 10,(id_before/100) % 10);
    if (id_before!=5 && id_before!=4)  {
        first_h=1;
        std::cout << "first h" <<std::endl;
        std::cout << "first h" <<" vx  "<<mother->daughter(0)->vx()<<" vy  "<<mother->daughter(0)->vy()<<" vz  "<<mother->daughter(0)->vz()<<std::endl;
        std::cout << "first h" <<" p  "<<mother->pt()<<" eta  "<<mother->eta()<<" phi  "<<mother->phi()<<std::endl;
        for (unsigned int o=0; o<mother->numberOfDaughters(); o++){
            std::cout<<"dau___"<<mother->daughter(o)->pdgId()<<"---"<<mother->daughter(o)->pt()<<" --  "<<mother->daughter(o)->eta()<<" --  "<<mother->daughter(o)->phi()<<std::endl;
        }
        vtx[0]=(mother->pt());
        vtx[1]=(mother->eta());
        vtx[2]=(mother->phi());
        vtx[3]=(mother->daughter(0)->vx());
        vtx[4]=(mother->daughter(0)->vy());
        vtx[5]=(mother->daughter(0)->vz());
        vtx[6]=(mother->pdgId());
    }
    }
    std::cout << "torefvect" <<std::endl;
    mother=mother->mother(0);
    std::cout << "refvect" << mother->numberOfMothers()<<std::endl;
    }
    
    
return vtx;   
}


std::vector<int> AnalyzerSignedIP_MINIAOD::seedToAllJets30Matching(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi, std::vector<double> jet_pt)
{
  float angular_distance =  0.4;
  float dist_temp;
  std::vector<int> index_array;
  for (unsigned int i = 0; i <jet_eta.size(); i++) {
     if (jet_pt[i]>30) {
     dist_temp=std::sqrt(std::pow(jet_eta[i]-seed.track().eta(),2) + std::pow(jet_phi[i]-seed.track().phi(),2) );
     if (dist_temp<=angular_distance) {
         std::cout<< "jet macthibg function"<< jet_eta[i]<<"  " <<seed.track().eta()<<"  " << jet_phi[i]<<"  " <<seed.track().phi()<<std::endl;
//          angular_distance=dist_temp;
         index_array.push_back(i);
    }}
      
}
return index_array;
}



int AnalyzerSignedIP_MINIAOD::seedToJetMatching(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi) {
  float angular_distance =  0.4;
  float dist_temp;
  int index = -1;
  for (unsigned int i = 0; i <jet_eta.size(); i++) {
     if (index==-1) {
     dist_temp=std::sqrt(std::pow(jet_eta[i]-seed.track().eta(),2) + std::pow(jet_phi[i]-seed.track().phi(),2) );
     if (dist_temp<angular_distance) {
         std::cout<< "jet macthibg function"<< jet_eta[i]<<"  " <<seed.track().eta()<<"  " << jet_phi[i]<<"  " <<seed.track().phi()<<std::endl;
//          angular_distance=dist_temp;
         index=i;
    }}
      
}
std::cout<<angular_distance<<"  index  "<<index<<std::endl;
return index;
}

int AnalyzerSignedIP_MINIAOD::seedToJet30Matching(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi, std::vector<double> jet_pt) {
  float angular_distance =  0.4;
  float dist_temp;
  int index = -1;
  for (unsigned int i = 0; i <jet_eta.size(); i++) {
     if (jet_pt[i]>30){
     dist_temp=std::sqrt(std::pow(jet_eta[i]-seed.track().eta(),2) + std::pow(jet_phi[i]-seed.track().phi(),2) );
     if (dist_temp<angular_distance) {
         
         angular_distance=dist_temp;
         index=i;
    }
    }
      
}
std::cout<<angular_distance<<"  index 30 "<<index<<std::endl;
return index;
}

std::pair<int,float> AnalyzerSignedIP_MINIAOD::ClosestJet(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi) {
  float angular_distance =  15;
  float dist_temp;
  int index = -1;
  std::pair<int, float> values;
  for (unsigned int i = 0; i <jet_eta.size(); i++) {

     dist_temp=std::sqrt(std::pow(jet_eta[i]-seed.track().eta(),2) + std::pow(jet_phi[i]-seed.track().phi(),2) );
     if (dist_temp<angular_distance) {

         angular_distance=dist_temp;
         index=i;
    }

}
values.first=index;
values.second=angular_distance;
std::cout<<angular_distance<<"  index  "<<index<<std::endl;
return values;
}





double AnalyzerSignedIP_MINIAOD::trackEff(std::vector<pat::PackedGenParticle> genP, std::vector<reco::TransientTrack> tracks, double &BChain_eff, double &DChain_eff, double &B_eff,
      double &low_pt_eff, double &high_pt_eff, std::vector<double> &minGenChiSquare){

const reco::Candidate* mother; 
minGenChiSquare.clear();
int id_mom = -1;
int b_found=0;
int d_found=0;
double eff=0;
int sample_size=0;
int B_sample_size=0;
int BChain_sample_size=0;
int DChain_sample_size=0;
int highpt_sample_size=0;
int lowpt_sample_size=0;
B_eff=0;
BChain_eff=0;
DChain_eff=0;
high_pt_eff=0;
low_pt_eff=0;

if (genP.size()==0)
return 0;

std::vector<trackGenMatch2> chi_squares; 
edm::RefVector<std::vector<pat::PackedGenParticle> > mothers; 
trackGenMatch2 chi_square; 

std::cout << "function entered .." << std::endl;
for (unsigned int g = 0; g < genP.size(); g++){
    chi_squares.clear();
    if (std::fabs(genP[g].eta())<2.4&&genP[g].status()==1&&genP[g].charge()!=0){

        std::cout<< genP[g].pdgId() << "P pdgid " <<std::endl;        
        std::cout<< genP[g].mother(0)->pdgId() << " mother pdgid " <<std::endl;
        std::cout<< genP[g].numberOfMothers() << " Nmothers " <<std::endl;
        if (genP[g].mother(0)==NULL) continue;        
        if (genP[g].mother(0)->numberOfDaughters()<1)  continue; 
        if (genP[g].mother(0)->pdgId()==2212)  continue;
        auto newgenP=genP[g].mother(0)->daughter(0);
        
        for(unsigned int ig=0; ig <genP[g].mother(0)->numberOfDaughters(); ig++){
            newgenP=genP[g].mother(0)->daughter(ig);
            std::cout<< newgenP->pdgId() << " newP pdgid in cycle  " << genP[g].mother(0)->numberOfDaughters()<< std::endl;         
            if (newgenP->pdgId()==genP[g].pdgId() && std::fabs(newgenP->pt()-genP[g].pt())<0.01) {
                
             
             break;   
            }
                
                
            
        }
        
        std::cout<< newgenP->pdgId() << " newP pdgid " <<std::endl;       
        
        double dxy = (-newgenP->vx() * newgenP->py() + newgenP->vy() * newgenP->px()) / newgenP->pt();
        double dsz =  newgenP->vz() * newgenP->pt() / newgenP->p() - (newgenP->vx() * newgenP->px() + newgenP->vy() * newgenP->py()) / newgenP->pt() * newgenP->pz() / newgenP->p();

        
        sample_size=sample_size+1;
//        std::cout<<"genp checled"<<std::endl;
//         double dxy = (-genP[g].vx() * genP[g].py() + genP[g].vy() * genP[g].px()) / genP[g].pt();
//         double dsz =  genP[g].vz() * genP[g].pt() / genP[g].p() - (genP[g].vx() * genP[g].px() + genP[g].vy() * genP[g].py()) / genP[g].pt() * genP[g].pz() / genP[g].p();
        
       std::cout<<"genp checled  eta "<<genP[g].eta()<< " pt" <<genP[g].pt()<< "  phi "<< genP[g].phi()<< std::endl; 
       std::cout<<"genp checled  eta "<<newgenP->eta()<< " pt" <<newgenP->pt()<< "  phi "<< newgenP->phi()<< std::endl; 
        std::cout<<"genp checled  vx "<<newgenP->vx()<< " vy " <<newgenP->vy()<< " vz "<< newgenP->vz()<< std::endl;
//         std::cout<<"genp checled  dxy "<<newgenP->dxy()<< " dz " <<newgenP->dz()<< std::endl;
        for(std::vector<reco::TransientTrack>::const_iterator it = tracks.begin(); it != tracks.end(); it++){
        double chi, tot=0;

        chi= (it->track().eta()-newgenP->eta())/(it->track().etaError());
        chi=chi*chi;
        tot=tot+chi;
        chi= (it->track().pt()-newgenP->pt())/(it->track().ptError());
        chi=chi*chi;
        tot=tot+chi;
        chi= (it->track().phi()-newgenP->phi())/(it->track().phiError());
        chi=chi*chi;
        tot=tot+chi;
        chi= (it->track().dsz()-dsz)/(it->track().dszError());
        chi=chi*chi;
        tot=tot+chi;
        chi= (it->track().dxy()-dxy)/(it->track().dxyError());
        chi=chi*chi;
        tot=tot+chi;
        chi_square.set_chi(tot);
//        chi_square.set_GenIndex(g);
        chi_squares.push_back(chi_square);
//         std::cout<<tot<<std::endl;
        
        }
    std::cout<<"sorting"<<std::endl;
    if (chi_squares.size()==0)
        return 0;
    std::sort (chi_squares.begin(), chi_squares.end(), sortgen2());
    std::cout<<"sorted ok"<<std::endl;

    minGenChiSquare.push_back(chi_squares[0].chi_square);
    mother=genP[g].mother(0);
    
   //bchain

    b_found=0;
    d_found=0;
    while (mother!=NULL && b_found==0 && mother->status()<3)
    {
    id_mom = abs(mother->pdgId());
    id_mom=std::max((id_mom/1000) % 10,(id_mom/100) % 10);
    if (id_mom==5)
    b_found=1;
    if (id_mom==4)
    d_found=1;
    mother=mother->mother(0);
    }


    //mother id
    mother=genP[g].mother(0);
        
    if (mother!=NULL )
    {id_mom = abs(mother->pdgId());
    id_mom=std::max((id_mom/1000) % 10,(id_mom/100) % 10);}
    else
    id_mom = -1;
            
    
    if (genP[g].pt()<0.5){
    lowpt_sample_size++;
    if (chi_squares[0].chi_square<50)
    low_pt_eff++;
    }
    else{
    highpt_sample_size++;
    if (chi_squares[0].chi_square<50)
    high_pt_eff++; }
    
    if (id_mom==5){
    B_sample_size++;
    if (chi_squares[0].chi_square<50)
    B_eff++;}
    
    if (b_found){
    BChain_sample_size++;
    if (chi_squares[0].chi_square<50)
    BChain_eff++;}
    
    if (d_found){
    DChain_sample_size++;
    if (chi_squares[0].chi_square<50)
    DChain_eff++;}

   
    
    if (chi_squares[0].chi_square<50){
    eff=eff+1;
    }
    
    }
    }


std::cout<<"division"<<std::endl;
B_eff=B_eff/B_sample_size;
BChain_eff=BChain_eff/BChain_sample_size;
DChain_eff=DChain_eff/DChain_sample_size;
high_pt_eff=high_pt_eff/highpt_sample_size;
low_pt_eff=low_pt_eff/lowpt_sample_size;


eff=eff/sample_size;
return eff;
}



bool AnalyzerSignedIP_MINIAOD::genLevelMM(std::vector<pat::PackedGenParticle> genP, std::vector<reco::TransientTrack> tracks, std::vector<int>& genPnumbers, std::vector<double>& assigned_chisquares){ 
    std::multimap<float,std::pair<int,int>> chisquares_map;
    genPnumbers.clear();
    assigned_chisquares.clear();
    for(std::vector<reco::TransientTrack>::const_iterator it = tracks.begin(); it != tracks.end(); it++){
        genPnumbers.push_back(-1);
        //baseline == no assignment
    for (unsigned int g = 0; g < genP.size(); g++){
        
        if (std::fabs(genP[g].eta())<5&&genP[g].status()==1){
            
        const reco::Candidate* mother=genP[g].mother(0);        
        if (mother==NULL) break;
        
        auto newgenP=mother->daughter(0);
        if (newgenP==NULL) break;
        
        double dxy = (-newgenP->vx() * newgenP->py() + newgenP->vy() * newgenP->px()) / newgenP->pt();
        double dsz =  newgenP->vz() * newgenP->pt() / newgenP->p() - (newgenP->vx() * newgenP->px() + newgenP->vy() * newgenP->py()) / newgenP->pt() * newgenP->pz() / newgenP->p();

        double chi, tot=0;

        chi= (it->track().eta()-newgenP->eta())/(it->track().etaError());
        chi=chi*chi;
        tot=tot+chi;
        chi= (it->track().pt()-newgenP->pt())/(it->track().ptError());
        chi=chi*chi;
        tot=tot+chi;
        chi= (it->track().phi()-newgenP->phi())/(it->track().phiError());
        chi=chi*chi;
        tot=tot+chi;
        chi= (it->track().dsz()-dsz)/(it->track().dszError());
        chi=chi*chi;
        tot=tot+chi;
        chi= (it->track().dxy()-dxy)/(it->track().dxyError());
        chi=chi*chi;

        tot=tot+chi;
        
        if (tot<1000){chisquares_map.insert(std::make_pair(tot, std::make_pair(it-tracks.begin(),g)));}

        }
    }
    }
    
    for(std::multimap<float,std::pair<int,int>>::const_iterator im = chisquares_map.begin(); im != chisquares_map.end(); im++){
        
              
        std::cout<<"@@  chi  "<<im->first<<"  track_ind  "<< im->second.first <<"  gen_index  "<< im->second.second<<std::endl;
        //check if assigned--otherwise assigned
        if (genPnumbers[im->second.first]==-1){
           if (im->first<1000) {
               
               genPnumbers[im->second.first]=im->second.second;
               assigned_chisquares.push_back(im->first);
	       std::cout<< "matching happened @@ chi2 "<<im->first<<std::endl;
        }
        }

    }
    
    
    
return 1;
}


int AnalyzerSignedIP_MINIAOD::genMap(std::vector<pat::PackedGenParticle> genP, std::vector<reco::TransientTrack> tracks, std::vector<int>& genPnumbers, std::vector<double>& allChi){

std::vector<trackGenMatch2> chi_squares;  
trackGenMatch2 chi_square; 
edm::RefVector<std::vector<pat::PackedGenParticle> > mothers;
edm::RefVector<std::vector<pat::PackedGenParticle> > daughters;
genPnumbers.clear();

for(std::vector<reco::TransientTrack>::const_iterator it = tracks.begin(); it != tracks.end(); it++){
    chi_squares.clear();
    std::cout << "vertex for the track" <<  it->track().vx() <<  " "<< it->track().vy() <<  " " << it->track().vz() << std::endl;
    int j=-1;
    for (unsigned int g = 0; g < genP.size(); g++){

        if (std::fabs(genP[g].eta())<5){
        double dxy = (-genP[g].vx() * genP[g].py() + genP[g].vy() * genP[g].px()) / genP[g].pt();
        double dsz =  genP[g].vz() * genP[g].pt() / genP[g].p() - (genP[g].vx() * genP[g].px() + genP[g].vy() * genP[g].py()) / genP[g].pt() * genP[g].pz() / genP[g].p();

        double chi, tot=0;

        chi= (it->track().eta()-genP[g].eta())/(it->track().etaError());
        chi=chi*chi;
        tot=tot+chi;
        chi= (it->track().pt()-genP[g].pt())/(it->track().ptError());
        chi=chi*chi;
        tot=tot+chi;
        chi= (it->track().phi()-genP[g].phi())/(it->track().phiError());
        chi=chi*chi;
        tot=tot+chi;
        chi= (it->track().dsz()-dsz)/(it->track().dszError());
        chi=chi*chi;
        tot=tot+chi;
        chi= (it->track().dxy()-dxy)/(it->track().dxyError());
        chi=chi*chi;
        tot=tot+chi;
        
        

        chi_square.set_chi(tot);
        chi_square.set_GenIndex(g);
        chi_square.set_numberOfDaughters(genP[g].numberOfDaughters());
        chi_square.set_Status(genP[g].status());
        chi_squares.push_back(chi_square);

        }
    }
    
    std::sort (chi_squares.begin(), chi_squares.end(), sortgen2());
    bool to_match=1;
    bool match=0;
    unsigned int i=0;
    while(match==0 and i < chi_squares.size() and to_match)
    {
    if (chi_squares[i].Status==1 && chi_squares[i].chi_square<50) /*detector-stable particle*/{
//    if (chi_squares[i].chi_square<50) /*detector-stable particle*/{
    match=1;
    }
    else if (chi_squares[i].chi_square>50){to_match=0;}
    i++;
    }
    if (match){
    allChi.push_back(chi_squares[i-1].chi_square);
    j=chi_squares[i-1].GenIndex;
    genPnumbers.push_back(j);
    std::cout << " see how the matching went: index  " << chi_squares[i-1].GenIndex << ", chisqaure  " <<  chi_squares[i-1].chi_square<< std::endl;
    std::cout << " a bit of a check if sorting is ok " << chi_squares[i-1].GenIndex << ", chisqaure  " <<  chi_squares[i-1].chi_square
    << ", chisqaure 2 " <<  chi_squares[i].chi_square<< ", chisqaure 3 " <<  chi_squares[i+1].chi_square<< std::endl;}
    
    else{allChi.push_back(chi_squares[0].chi_square);
//    j=chi_squares[i-1].GenIndex;
    genPnumbers.push_back(-1);
    std::cout << " no matching here: index  " << chi_squares[0].GenIndex << ", chisqaure  " <<  chi_squares[0].chi_square<< std::endl;
    std::cout << " a bit of a check if sorting is ok " << chi_squares[0].GenIndex << ", chisqaure  " <<  chi_squares[0].chi_square
    << ", chisqaure 2 " <<  chi_squares[1].chi_square<< ", chisqaure 3 " <<  chi_squares[2].chi_square<< std::endl;
    }
}



return 1;

}






int AnalyzerSignedIP_MINIAOD::genLevel(std::vector<pat::PackedGenParticle> genP, reco::TransientTrack seed, std::vector<int>& MomFlav, std::vector<int>& BChain, std::vector<double>& allChi){

std::vector<trackGenMatch2> chi_squares;  
trackGenMatch2 chi_square; 

edm::RefVector<std::vector<pat::PackedGenParticle> > mothers;
edm::RefVector<std::vector<pat::PackedGenParticle> > daughters;

int id, id_mom=-111111, id_dau=-111111;

for (unsigned int i = 0; i < genP.size(); i++)
{

if (fabs(genP[i].eta())<5){
////dxy = (-vx() * py() + vy() * px()) / pt();
////dsz =  vz() * pt() / p() - (vx() * px() + vy() * py()) / pt() * pz() / p();

double dxy = (-genP[i].vx() * genP[i].py() + genP[i].vy() * genP[i].px()) / genP[i].pt();
double dsz =  genP[i].vz() * genP[i].pt() / genP[i].p() - (genP[i].vx() * genP[i].px() + genP[i].vy() * genP[i].py()) / genP[i].pt() * genP[i].pz() / genP[i].p();

//std::cout << dxy << " " << dsz << " " << genP[i].phi() << " " << "PRIMA RIGA" << std::endl;
//std::cout << genP[i].pt() << " " << genP[i].eta() << " " << genP[i].phi() << " " << "2 RIGA" << std::endl;
//std::cout << seed.track().phi() << " " << seed.track().pt() << " " << seed.track().eta() << " " << "3 RIGA" << std::endl;
//std::cout << seed.track().dsz() << " " << seed.track().dxy() << " " << seed.track().eta() << " " << std::endl;
//std::cout << seed.track().phiError() << " " << seed.track().ptError() << " " << seed.track().etaError() << " " << std::endl;
//std::cout << seed.track().dszError() << " " << seed.track().dxyError() << " " << seed.track().etaError() << " " << std::endl;

double chi, tot=0;

chi= (seed.track().eta()-genP[i].eta())/(seed.track().etaError());
chi=chi*chi;
tot=tot+chi;
chi= (seed.track().pt()-genP[i].pt())/(seed.track().ptError());
chi=chi*chi;
tot=tot+chi;
chi= (seed.track().phi()-genP[i].phi())/(seed.track().phiError());
chi=chi*chi;
tot=tot+chi;
chi= (seed.track().dsz()-dsz)/(seed.track().dszError());
chi=chi*chi;
tot=tot+chi;
chi= (seed.track().dxy()-dxy)/(seed.track().dxyError());
chi=chi*chi;
tot=tot+chi;



//std::cout << tot << " chisqaure  "<< std::endl;


chi_square.set_chi(tot);
allChi.push_back(tot);
chi_square.set_GenIndex(i);
chi_square.set_numberOfDaughters(genP[i].numberOfDaughters());
chi_square.set_Status(genP[i].status());


chi_squares.push_back(chi_square);


}
}


std::sort (chi_squares.begin(), chi_squares.end(), sortgen2());

bool to_match=1;
bool match=0;
unsigned int i=0; 

while(match==0 and i < chi_squares.size() and to_match)
{
if (chi_squares[i].Status==1) /*detector-stable particle*/{
match=1;

}

std::cout << i << " while here" << match << "  " << chi_squares[i].Status << " "  << chi_squares[i].GenIndex   << std::endl;
i++;
}

std::cout << i-1 << " end of while" << std::endl;
i=chi_squares[i-1].GenIndex;
std::cout << i-1 << " reset while" << std::endl;

const reco::Candidate* daughter=genP[i-1].daughter(0);
const reco::Candidate* mother=genP[i-1].mother(0);

id = genP[i-1].pdgId();

//mother id

if (genP[i-1].numberOfMothers()>0 )
{id_mom = mother->pdgId();
id_mom=std::max((id_mom/1000) % 10,(id_mom/100) % 10);}
else
id_mom = -1;

MomFlav.push_back(id_mom);
//mother id



if (genP[i-1].numberOfDaughters()>0 )
id_dau = daughter->pdgId();
else
id_dau = -111111111;

std::cout << id << " " << id_mom << " " << id_dau << " " << std::endl;
std::cout << genP[i-1].numberOfMothers() << " " << genP[i-1].numberOfDaughters() << " "  << " " << std::endl;
std::cout << genP[i].status() << " " << std::endl;

//b in chain
int b_found=0;
while (mother!=NULL && b_found==0 && mother->status()<3)
{
id_mom = mother->pdgId();
id_mom=std::max((id_mom/1000) % 10,(id_mom/100) % 10);
if (id_mom==5)
b_found=1;
mother=mother->mother(0);
}

BChain.push_back(b_found);

if (match)
return i;
else
return -1;

}



// ------------ method called once each job just before starting event loop  ------------
void 
AnalyzerSignedIP_MINIAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AnalyzerSignedIP_MINIAOD::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AnalyzerSignedIP_MINIAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzerSignedIP_MINIAOD);
