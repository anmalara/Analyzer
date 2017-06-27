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
#include "DataFormats/GeometrySurface/interface/Line.h"

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


#include "RecoBTag/TrackProbability/interface/HistogramProbabilityEstimator.h"

class HistogramProbabilityEstimator;
#include <typeinfo>


#include "CondFormats/BTauObjects/interface/TrackProbabilityCalibration.h"
#include "CondFormats/DataRecord/interface/BTagTrackProbability2DRcd.h"
#include "CondFormats/DataRecord/interface/BTagTrackProbability3DRcd.h"
#include "FWCore/Framework/interface/EventSetupRecord.h"
#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "FWCore/Framework/interface/EventSetupRecordKey.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.


class AnalyzerSignedIP_MINIAOD_wArr : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit AnalyzerSignedIP_MINIAOD_wArr(const edm::ParameterSet&);
      ~AnalyzerSignedIP_MINIAOD_wArr();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      // ----------member data ---------------------------
      
      int funzionePippo (int vivaLaVita);
      int seedToJetMatching(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi);
      std::vector<int> seedToAllJets20Matching(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi, std::vector<double> jet_pt);
      int seedToJet30Matching(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi, std::vector<double> jet_pt);
      std::pair<int,float> ClosestJet(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi);
      int genLevel(std::vector<pat::PackedGenParticle> genP, reco::TransientTrack seed, std::vector<int>& MomFlav, std::vector<int>& BChain, std::vector<double>& allChi);      
      int genMap(std::vector<pat::PackedGenParticle> genP, std::vector<reco::TransientTrack> tracks, std::vector<int>& genPnumbers, std::vector<double>& allChi);  
      template<class T> void insert_inArray(T (&tArray)[10], unsigned int size, unsigned int index, double element);
	  template<class T> void insert_inArray200(T (&tArray)[200], unsigned int size, unsigned int index, double element);
      int jet_flavour_function(const pat::Jet& jet, bool usePhysForLightAndUndefined);
      double trackEff(std::vector<pat::PackedGenParticle> genP, std::vector<reco::TransientTrack> tracks, double &BChain_eff, double &DChain_eff, double &B_eff,
      double &low_pt_eff, double &high_pt_eff, std::vector<double> &minGenChiSquare);
      
      //genMap(genParticles, selectedTracks, genP_index, Match_chisquare );
      bool genLevelMM(std::vector<pat::PackedGenParticle> genP, std::vector<reco::TransientTrack> tracks, std::vector<int>& genPnumbers, std::vector<double>& assigned_chisquares);
      
      std::vector<float> firstHadronInChain(pat::PackedGenParticle gParticle, const reco::Candidate* mother);
      
      void checkEventSetup(const edm::EventSetup & iSetup);
      
      edm::Service<TFileService> file;
//       edm::RefVector<std::vector<pat::PackedGenParticle> > mothers;
//       reco::Candidate mother;

      TTree *tree;
      std::auto_ptr<HistogramProbabilityEstimator> m_probabilityEstimator;
      bool m_computeProbabilities=1;
      unsigned long long  m_calibrationCacheId2D; 
      unsigned long long m_calibrationCacheId3D;
      
      int count_seeds=0;
      
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
      std::vector<float> jet_CSVv2;
      std::vector<float> jet_CMVA;
      std::vector<float> jet_TCHE;
      std::vector<float> jet_SSVHE;
      std::vector<reco::Jet> Jets;
     
      double jetpt;      
      double jeteta;
      double jetphi;
      double jetmass;
      int jetflavour;
      int jetNseeds;
      float jetCSVv2;
      float jetCMVA;
      float jetTCHE;
      float jetSSVHE;
      
      double seed_pt[10];
      double seed_eta[10];
      double seed_phi[10];
      double seed_mass[10];
      
      double seed_dz[10];
      double seed_dxy[10];
      double seed_3D_ip[10];
      double seed_3D_sip[10];
      double seed_2D_ip[10];
      double seed_2D_sip[10];
      double seed_3D_signedIp[10];
      double seed_3D_signedSip[10];
      double seed_2D_signedIp[10];
      double seed_2D_signedSip[10];
      double seed_3D_TrackProbability[10];
      double seed_2D_TrackProbability[10];
      int seed_JetMatch[10];
      
                
      double seed_chi2reduced[10];
      double seed_nPixelHits[10];
      double seed_nHits[10];
      double seed_jetAxisDistance[10];
      double seed_jetAxisDlength[10];
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
      double nearTracks_pt[200];
      double nearTracks_eta[200];
      double nearTracks_phi[200];
      double nearTracks_mass[200];
      double nearTracks_dz[200];
      double nearTracks_dxy[200];
      double nearTracks_3D_ip[200];
      double nearTracks_3D_sip[200];
      double nearTracks_2D_ip[200];
      double nearTracks_2D_sip[200];
      double nearTracks_PCAdist[200];
      double nearTracks_PCAdsig[200];      
      double nearTracks_PCAonSeed_x[200];
      double nearTracks_PCAonSeed_y[200];
      double nearTracks_PCAonSeed_z[200];      
      double nearTracks_PCAonSeed_xerr[200];
      double nearTracks_PCAonSeed_yerr[200];
      double nearTracks_PCAonSeed_zerr[200];      
      double nearTracks_PCAonTrack_x[200];
      double nearTracks_PCAonTrack_y[200];
      double nearTracks_PCAonTrack_z[200];      
      double nearTracks_PCAonTrack_xerr[200];
      double nearTracks_PCAonTrack_yerr[200];
      double nearTracks_PCAonTrack_zerr[200]; 
      double nearTracks_dotprodTrack[200];
      double nearTracks_dotprodSeed[200];
      double nearTracks_dotprodTrackSeed2D[200];
      double nearTracks_dotprodTrackSeed3D[200];
      double nearTracks_dotprodTrackSeedVectors2D[200];
      double nearTracks_dotprodTrackSeedVectors3D[200];      
      double nearTracks_PCAonSeed_pvd[200];
      double nearTracks_PCAonTrack_pvd[200];
      double nearTracks_PCAjetAxis_dist[200];
      double nearTracks_PCAjetMomenta_dotprod[200];
      double nearTracks_PCAjetDirs_DEta[200];
      double nearTracks_PCAjetDirs_DPhi[200];

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


      float min3DIPValue=0.005;
      float min3DIPSignificance=1.2;
      int max3DIPValue=9999.;
      int max3DIPSignificance=9999.;
	  
	  trackVars2 myTrack;	  
      std::vector<trackVars2> nearTracks; 
      
      /*
       * IVF configuration for seeding
       * seedMax3DIPSignificance = cms.double(9999.),
           seedMax3DIPValue = cms.double(9999.),
           seedMin3DIPSignificance = cms.double(1.2),
           seedMin3DIPValue = cms.double(0.005),*/
      
      std::pair<int, float> pair;
      std::multimap<double,std::pair<const reco::TransientTrack*,const std::vector<trackVars2>>> SortedSeedsMap; /*,std::vector<trackVars2>>> SortedSeedsMap;
      std::multimap<double,std::pair<&reco::TransientTrack ,std::vector<trackVars2>>> SortedSeedsMap;*/
      
      
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

/*Questi tagli ip tag*/

// if (track.pt() > m_cutMinPt &&
// track.hitPattern().numberOfValidHits() >= m_cutTotalHits &&         // min num tracker hits
// track.hitPattern().numberOfValidPixelHits() >= m_cutPixelHits &&
// track.normalizedChi2() < m_cutMaxChiSquared &&
// std::abs(track.dxy(pv->position())) < m_cutMaxTIP &&
// std::abs(track.dz(pv->position())) < m_cutMaxLIP)
// 
//     if( fabs(impactParameters[i].distanceToJetAxis.value()) < m_cutMaxDistToAxis  &&        // distance to JetAxis 0.07
// (impactParameters[i].closestToJetAxis - pv).mag() < m_cutMaxDecayLen  &&      // max decay len 5
// 
// )



AnalyzerSignedIP_MINIAOD_wArr::AnalyzerSignedIP_MINIAOD_wArr(const edm::ParameterSet& iConfig)

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
   
    
    tree->Branch("jet_pt",&jetpt, "jet_pt/D");
    tree->Branch("jet_eta",&jeteta, "jet_eta/D");
    tree->Branch("jet_phi",&jetphi, "jet_phi/D");
    tree->Branch("jet_mass",&jetmass, "jet_mass/D");
    tree->Branch("jet_flavour",&jetflavour, "jet_flavour/I");
    tree->Branch("jet_Nseeds",&jetNseeds, "jet_Nseeds/I"); 
    
    tree->Branch("jet_CSVv2", &jetCSVv2,  "jet_CSVv2/F");
    tree->Branch("jet_CMVA", &jetCMVA,  "jet_CMVA/F");
    tree->Branch("jet_TCHE", &jetTCHE,  "jet_TCHE/F");
    tree->Branch("jet_SSVHE", &jetSSVHE,  "jet_SSVHE/F");
    
    
    tree->Branch("n_seed",&n_seed, "n_seed/I");

    tree->Branch("seed_pt",&seed_pt, "seed_pt[10]/D");
    tree->Branch("seed_eta",&seed_eta, "seed_eta[10]/D");
    tree->Branch("seed_phi",&seed_phi, "seed_phi[10]/D");
    tree->Branch("seed_mass",&seed_mass, "seed_mass[10]/D");
    
    tree->Branch("seed_dz", &seed_dz, "seed_dz[10]/D");
    tree->Branch("seed_dxy", &seed_dxy, "seed_dxy[10]/D");
    tree->Branch("seed_3D_ip", &seed_3D_ip, "seed_3D_ip[10]/D");
    tree->Branch("seed_3D_sip", &seed_3D_sip, "seed_3D_sip[10]/D");
    tree->Branch("seed_2D_ip", &seed_2D_ip, "seed_2D_ip[10]/D");
    tree->Branch("seed_2D_sip", &seed_2D_sip, "seed_2D_sip[10]/D");
    
    tree->Branch("seed_3D_signedIp", &seed_3D_signedIp, "seed_3D_signedIp[10]/D");
    tree->Branch("seed_3D_signedSip", &seed_3D_signedSip, "seed_3D_signedSip[10]/D");
    tree->Branch("seed_2D_signedIp", &seed_2D_signedIp, "seed_2D_signedIp[10]/D");
    tree->Branch("seed_2D_signedSip", &seed_2D_signedSip, "seed_2D_signedSip[10]/D");
    tree->Branch("seed_3D_TrackProbability", &seed_3D_TrackProbability, "seed_3D_TrackProbability[10]/D");
    tree->Branch("seed_2D_TrackProbability", &seed_2D_TrackProbability, "seed_2D_TrackProbability[10]/D");
    
    tree->Branch("seed_JetMatch", &seed_JetMatch, "seed_JetMatch[10]/D");
    tree->Branch("seed_chi2reduced",&seed_chi2reduced, "seed_chi2reduced[10]/D");
    tree->Branch("seed_nPixelHits",&seed_nPixelHits, "seed_nPixelHits[10]/D");
    tree->Branch("seed_nHits",&seed_nHits, "seed_nHits[10]/D");
    tree->Branch("seed_jetAxisDistance",&seed_jetAxisDistance, "seed_jetAxisDistance[10]/D");
    tree->Branch("seed_jetAxisDlength",&seed_jetAxisDlength, "seed_jetAxisDlength[10]/D");
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
    tree->Branch("nearTracks_pt", &nearTracks_pt, "nearTracks_pt[200]/D");
    tree->Branch("nearTracks_eta", &nearTracks_eta, "nearTracks_eta[200]/D");
    tree->Branch("nearTracks_phi", &nearTracks_phi, "nearTracks_phi[200]/D");
    tree->Branch("nearTracks_mass", &nearTracks_mass, "nearTracks_mass[200]/D");
    tree->Branch("nearTracks_dz", &nearTracks_dz, "nearTracks_dz[200]/D");
    tree->Branch("nearTracks_dxy", &nearTracks_dxy, "nearTracks_dxy[200]/D");
    tree->Branch("nearTracks_3D_ip", &nearTracks_3D_ip, "nearTracks_3D_ip[200]/D");
    tree->Branch("nearTracks_3D_sip", &nearTracks_3D_sip, "nearTracks_3D_sip[200]/D");
    tree->Branch("nearTracks_2D_ip", &nearTracks_2D_ip, "nearTracks_2D_ip[200]/D");
    tree->Branch("nearTracks_2D_sip", &nearTracks_2D_sip, "nearTracks_2D_sip[200]/D");

    tree->Branch("nearTracks_PCAdist", &nearTracks_PCAdist, "nearTracks_PCAdist[200]/D");
    tree->Branch("nearTracks_PCAdsig", &nearTracks_PCAdsig, "nearTracks_PCAdsig[200]/D");
    
    tree->Branch("nearTracks_PCAonSeed_x", &nearTracks_PCAonSeed_x, "nearTracks_PCAonSeed_x[200]/D");
    tree->Branch("nearTracks_PCAonSeed_y", &nearTracks_PCAonSeed_y, "nearTracks_PCAonSeed_y[200]/D");
    tree->Branch("nearTracks_PCAonSeed_z", &nearTracks_PCAonSeed_z, "nearTracks_PCAonSeed_z[200]/D");

    tree->Branch("nearTracks_PCAonSeed_xerr", &nearTracks_PCAonSeed_xerr, "nearTracks_PCAonSeed_xerr[200]/D");
    tree->Branch("nearTracks_PCAonSeed_yerr", &nearTracks_PCAonSeed_yerr, "nearTracks_PCAonSeed_yerr[200]/D");
    tree->Branch("nearTracks_PCAonSeed_zerr", &nearTracks_PCAonSeed_zerr, "nearTracks_PCAonSeed_zerr[200]/D");

    tree->Branch("nearTracks_PCAonTrack_x", &nearTracks_PCAonTrack_x, "nearTracks_PCAonTrack_x[200]/D");
    tree->Branch("nearTracks_PCAonTrack_y", &nearTracks_PCAonTrack_y, "nearTracks_PCAonTrack_y[200]/D");
    tree->Branch("nearTracks_PCAonTrack_z", &nearTracks_PCAonTrack_z, "nearTracks_PCAonTrack_z[200]/D");

    tree->Branch("nearTracks_PCAonTrack_xerr", &nearTracks_PCAonTrack_xerr, "nearTracks_PCAonTrack_xerr[200]/D");
    tree->Branch("nearTracks_PCAonTrack_yerr", &nearTracks_PCAonTrack_yerr, "nearTracks_PCAonTrack_yerr[200]/D");
    tree->Branch("nearTracks_PCAonTrack_zerr", &nearTracks_PCAonTrack_zerr, "nearTracks_PCAonTrack_zerr[200]/D"); 

    tree->Branch("nearTracks_dotprodTrack", &nearTracks_dotprodTrack, "nearTracks_dotprodTrack[200]/D");
    tree->Branch("nearTracks_dotprodSeed", &nearTracks_dotprodSeed, "nearTracks_dotprodSeed[200]/D");
    tree->Branch("nearTracks_dotprodTrackSeed2D", &nearTracks_dotprodTrackSeed2D, "nearTracks_dotprodTrackSeed2D[200]/D");
    tree->Branch("nearTracks_dotprodTrackSeed3D", &nearTracks_dotprodTrackSeed3D, "nearTracks_dotprodTrackSeed3D[200]/D");
    tree->Branch("nearTracks_dotprodTrackSeedVectors2D", &nearTracks_dotprodTrackSeedVectors2D, "nearTracks_dotprodTrackSeedVectors2D[200]/D");
    tree->Branch("nearTracks_dotprodTrackSeedVectors3D", &nearTracks_dotprodTrackSeedVectors3D, "nearTracks_dotprodTrackSeedVectors3D[200]/D");
    
    tree->Branch("nearTracks_PCAonSeed_pvd", &nearTracks_PCAonSeed_pvd, "nearTracks_PCAonSeed_pvd[200]/D");
    tree->Branch("nearTracks_PCAonTrack_pvd", &nearTracks_PCAonTrack_pvd, "nearTracks_PCAonTrack_pvd[200]/D");
    tree->Branch("nearTracks_PCAjetAxis_dist",&nearTracks_PCAjetAxis_dist,"nearTracks_PCAjetAxis_dist[200]/D");
    tree->Branch("nearTracks_PCAjetMomenta_dotprod",&nearTracks_PCAjetMomenta_dotprod,"nearTracks_PCAjetMomenta_dotprod[200]/D");
    tree->Branch("nearTracks_PCAjetDirs_DEta",&nearTracks_PCAjetDirs_DEta,"nearTracks_PCAjetDirs_DEta[200]/D");
    tree->Branch("nearTracks_PCAjetDirs_DPhi",&nearTracks_PCAjetDirs_DPhi,"nearTracks_PCAjetDirs_DPhi[200]/D");
    
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


AnalyzerSignedIP_MINIAOD_wArr::~AnalyzerSignedIP_MINIAOD_wArr()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
AnalyzerSignedIP_MINIAOD_wArr::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
std::cout <<  " inizio evento  " << std::endl;
    if (m_computeProbabilities) checkEventSetup(iSetup);
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
    jet_CSVv2.clear();   
    jet_CMVA.clear();   
    jet_TCHE.clear();   
    jet_SSVHE.clear();  
    Jets.clear();


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
            int fl = jet_flavour_function(*iter, 1);
//             unsigned int fl = std::abs(iter->partonFlavour());
            
            if (iter->pt()>20){

            jet_flavour.push_back(fl);
            jet_pt.push_back(iter->pt());
            jet_eta.push_back(iter->eta());
            jet_phi.push_back(iter->phi());
            jet_mass.push_back(iter->mass());
            jet_px.push_back(iter->px());
            jet_py.push_back(iter->py());
            jet_pz.push_back(iter->pz());
            
            jet_CSVv2.push_back(iter->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
            jet_CMVA.push_back(iter->bDiscriminator("pfCombinedMVAV2BJetTags"));
            jet_TCHE.push_back(iter->bDiscriminator("pfTrackCountingHighEffBJetTags"));
            jet_SSVHE.push_back(iter->bDiscriminator("pfSimpleInclusiveSecondaryVertexHighEffBJetTags"));
            Jets.push_back(*(iter));
            
            }
            
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
        
     
        for (unsigned j_i=0; j_i < jet_pt.size(); j_i++) {
            
            //clear vars to fill
            
            
            std::fill_n(seed_pt, 10, 0.);//seed_pt.clear();
            std::fill_n(seed_eta, 10, 0.);//seed_eta.clear();
            std::fill_n(seed_phi, 10, 0.);//seed_phi.clear();
            std::fill_n(seed_mass, 10, 0.);//seed_mass.clear();
            std::fill_n(seed_dz, 10, 0.);//seed_dz.clear();
            std::fill_n(seed_dxy, 10, 0.);//seed_dxy.clear();
            std::fill_n(seed_3D_ip, 10, 0.);//seed_3D_ip.clear();
            std::fill_n(seed_3D_sip, 10, 0.);//seed_3D_sip.clear();
            std::fill_n(seed_2D_ip, 10, 0.);//seed_2D_ip.clear();
            std::fill_n(seed_2D_sip, 10, 0.);//seed_2D_sip.clear();
            std::fill_n(seed_3D_signedIp, 10, 0.);//seed_3D_signedIp.clear();
            std::fill_n(seed_3D_signedSip, 10, 0.);//seed_3D_signedSip.clear();
            std::fill_n(seed_2D_signedIp, 10, 0.);//seed_2D_signedIp.clear();
            std::fill_n(seed_2D_signedSip, 10, 0.);//seed_2D_signedSip.clear();
            std::fill_n(seed_JetMatch, 10, 0.);//seed_JetMatch.clear();
            std::fill_n(seed_chi2reduced, 10, 0.);//seed_chi2reduced.clear();
            std::fill_n(seed_nPixelHits, 10, 0.);//seed_nPixelHits.clear();
            std::fill_n(seed_nHits, 10, 0.);//seed_nHits.clear();
            std::fill_n(seed_jetAxisDistance, 10, 0.);//seed_jetAxisDistance.clear();
            std::fill_n(seed_jetAxisDlength, 10, 0.);//seed_jetAxisDlength.clear();    
            std::fill_n(seed_3D_TrackProbability, 10, 0.);
            std::fill_n(seed_2D_TrackProbability, 10, 0.);
            seed_ClosestJet_dR.clear();            
            seed_MC_MomFlavour.clear();
            seed_MC_MomPdgId.clear();
            seed_MC_BChain.clear();
            seed_MC_DChain.clear();            
            nearTracks_nTracks.clear();            
            nearTracks_Nvtx.clear();
            std::fill_n(nearTracks_pt, 200, 0.);//nearTracks_pt.clear();
            std::fill_n(nearTracks_eta, 200, 0.);//nearTracks_eta.clear();
            std::fill_n(nearTracks_phi, 200, 0.);//nearTracks_phi.clear();
            std::fill_n(nearTracks_dz, 200, 0.);//nearTracks_dz.clear();
            std::fill_n(nearTracks_dxy, 200, 0.);//nearTracks_dxy.clear();
            std::fill_n(nearTracks_mass, 200, 0.);//nearTracks_mass.clear();
            std::fill_n(nearTracks_3D_ip, 200, 0.);//nearTracks_3D_ip.clear();
            std::fill_n(nearTracks_3D_sip, 200, 0.);//nearTracks_3D_sip.clear();
            std::fill_n(nearTracks_2D_ip, 200, 0.);//nearTracks_2D_ip.clear();
            std::fill_n(nearTracks_2D_sip, 200, 0.);//nearTracks_2D_sip.clear();
            std::fill_n(nearTracks_PCAdist, 200, 0.);//nearTracks_PCAdist.clear();
            std::fill_n(nearTracks_PCAdsig, 200, 0.);//nearTracks_PCAdsig.clear();
            std::fill_n(nearTracks_PCAonSeed_x, 200, 0.);//nearTracks_PCAonSeed_x.clear();
            std::fill_n(nearTracks_PCAonSeed_y, 200, 0.);//nearTracks_PCAonSeed_y.clear();
            std::fill_n(nearTracks_PCAonSeed_z, 200, 0.);//nearTracks_PCAonSeed_z.clear();
            std::fill_n(nearTracks_PCAonSeed_xerr, 200, 0.);//nearTracks_PCAonSeed_xerr.clear();
            std::fill_n(nearTracks_PCAonSeed_yerr, 200, 0.);//nearTracks_PCAonSeed_yerr.clear();
            std::fill_n(nearTracks_PCAonSeed_zerr, 200, 0.);//nearTracks_PCAonSeed_zerr.clear();
            std::fill_n(nearTracks_PCAonTrack_x, 200, 0.);//nearTracks_PCAonTrack_x.clear();
            std::fill_n(nearTracks_PCAonTrack_y, 200, 0.);//nearTracks_PCAonTrack_y.clear();
            std::fill_n(nearTracks_PCAonTrack_z, 200, 0.);//nearTracks_PCAonTrack_z.clear();
            std::fill_n(nearTracks_PCAonTrack_xerr, 200, 0.);//nearTracks_PCAonTrack_xerr.clear();
            std::fill_n(nearTracks_PCAonTrack_yerr, 200, 0.);//nearTracks_PCAonTrack_yerr.clear();
            std::fill_n(nearTracks_PCAonTrack_zerr, 200, 0.);//nearTracks_PCAonTrack_zerr.clear();
            std::fill_n(nearTracks_dotprodTrack, 200, 0.);//nearTracks_dotprodTrack.clear();
            std::fill_n(nearTracks_dotprodSeed, 200, 0.);//nearTracks_dotprodSeed.clear();
            std::fill_n(nearTracks_dotprodTrackSeed2D, 200, 0.);//nearTracks_dotprodTrackSeed2D.clear();
            std::fill_n(nearTracks_dotprodTrackSeed3D, 200, 0.);//nearTracks_dotprodTrackSeed3D.clear();
            std::fill_n(nearTracks_dotprodTrackSeedVectors2D, 200, 0.);//nearTracks_dotprodTrackSeedVectors2D.clear();
            std::fill_n(nearTracks_dotprodTrackSeedVectors3D, 200, 0.);//nearTracks_dotprodTrackSeedVectors3D.clear();           
            std::fill_n(nearTracks_PCAonSeed_pvd, 200, 0.);//nearTracks_PCAonSeed_pvd.clear();
            std::fill_n(nearTracks_PCAonTrack_pvd, 200, 0.);//nearTracks_PCAonTrack_pvd.clear();
            std::fill_n(nearTracks_PCAjetAxis_dist, 200, 0.);
            std::fill_n(nearTracks_PCAjetMomenta_dotprod, 200, 0.);
            std::fill_n(nearTracks_PCAjetDirs_DEta, 200, 0.);
            std::fill_n(nearTracks_PCAjetDirs_DPhi, 200, 0.);
            
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
       
            jetpt=jet_pt[j_i];
            jeteta=jet_eta[j_i];
            jetphi=jet_phi[j_i];
            jetmass=jet_mass[j_i];
            jetflavour=jet_flavour[j_i];
            jetCSVv2=jet_CSVv2[j_i];
            jetCMVA=jet_CMVA[j_i];
            jetTCHE=jet_TCHE[j_i];
            jetSSVHE=jet_SSVHE[j_i];

            SortedSeedsMap.clear();
            
                       
            
        count_seeds=0;
        int seed_jet_match_i=j_i;
        GlobalVector direction(jet_px[seed_jet_match_i], jet_py[seed_jet_match_i], jet_pz[seed_jet_match_i]);
        
        for(std::vector<TransientTrack>::const_iterator it = selectedTracks.begin(); it != selectedTracks.end(); it++){
            
        float angular_distance=std::sqrt(std::pow(jet_eta[j_i]-it->track().eta(),2) + std::pow(jet_phi[j_i]-it->track().phi(),2) );
        
        if (angular_distance<0.4){  
            
        bool match_ok= (seed_jet_match_i>=0 );   
        
        
        std::pair<bool,Measurement1D> ip = IPTools::absoluteImpactParameter3D(*it, pv);        
        std::pair<bool,Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(*it, pv);
		std::pair<double, Measurement1D> jet_dist =IPTools::jetTrackDistance(*it, direction, pv);                   
        TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(it->impactPointState(),pv, direction,it->field());
		float length=777.77;
        if (closest.isValid()) length=(closest.globalPosition() - pvp).mag(); 

        
            std::cout << ip.second.significance() << " " <<  ( it->track().normalizedChi2()<5.) <<  match_ok<< (std::fabs(it->track().dxy(pv.position())) < 2 )<<  (std::fabs(it->track().dz(pv.position())) < 17) << std::endl;
        
            if(ip.first && ip.second.value() >= 0.0 && ip.second.significance() >= 1.0 &&
            ip.second.value() <= max3DIPValue && ip.second.significance() <= max3DIPSignificance &&
            it->track().normalizedChi2()<5. && std::fabs(it->track().dxy(pv.position())) < 2 &&
            std::fabs(it->track().dz(pv.position())) < 17  && jet_dist.second.value()<0.07 && length<5. )

            {
            
            std::pair<bool,Measurement1D> ipSigned = IPTools::signedImpactParameter3D(*it,direction, pv); 
            count_seeds++;
            nearTracks.clear();
            std::cout<< "nearTracks" << nearTracks.size() << "after clear" << std::endl;

            for(std::vector<reco::TransientTrack>::const_iterator tt = selectedTracks.begin();tt!=selectedTracks.end(); ++tt )   {
                VertexDistance3D distanceComputer;
                TwoTrackMinimumDistance dist;
                if(*tt==*it) continue;
                if(std::fabs(pv_z-tt->track().vz())>0.1) continue;
                if(dist.calculate(tt->impactPointState(),it->impactPointState())) {
                    GlobalPoint ttPoint          = dist.points().first;
                    GlobalError ttPointErr       = tt->impactPointState().cartesianError().position();
                    GlobalPoint seedPosition     = dist.points().second;
                    GlobalError seedPositionErr  = it->impactPointState().cartesianError().position();
                    Measurement1D m = distanceComputer.distance(VertexState(seedPosition,seedPositionErr), VertexState(ttPoint, ttPointErr));
                    GlobalPoint cp(dist.crossingPoint()); 

                    GlobalVector PairMomentum(it->track().px()+tt->track().px(), it->track().py()+tt->track().py(), it->track().pz()+tt->track().pz());
                    GlobalVector  PCA_pv(cp-pvp);

                    float PCAseedFromPV =  (dist.points().second-pvp).mag();
                    float PCAtrackFromPV =  (dist.points().first-pvp).mag();               
                    float distance = dist.distance();


                    GlobalVector trackDir2D(tt->impactPointState().globalDirection().x(),tt->impactPointState().globalDirection().y(),0.); 
                    GlobalVector seedDir2D(it->impactPointState().globalDirection().x(),it->impactPointState().globalDirection().y(),0.); 
                    GlobalVector trackPCADir2D(dist.points().first.x()-pvp.x(),dist.points().first.y()-pvp.y(),0.); 
                    GlobalVector seedPCADir2D(dist.points().second.x()-pvp.x(),dist.points().second.y()-pvp.y(),0.); 


                    float dotprodTrack = (dist.points().first-pvp).unit().dot(tt->impactPointState().globalDirection().unit());
                    float dotprodSeed = (dist.points().second-pvp).unit().dot(it->impactPointState().globalDirection().unit());                    
                    float dotprodTrackSeed2D = trackDir2D.unit().dot(seedDir2D.unit());
                    float dotprodTrackSeed3D = it->impactPointState().globalDirection().unit().dot(tt->impactPointState().globalDirection().unit());
                    float dotprodTrackSeed2DV = trackPCADir2D.unit().dot(seedPCADir2D.unit());
                    float dotprodTrackSeed3DV = (dist.points().second-pvp).unit().dot((dist.points().first-pvp).unit());

                    std::pair<bool,Measurement1D> t_ip = IPTools::absoluteImpactParameter3D(*tt,pv);        
                    std::pair<bool,Measurement1D> t_ip2d = IPTools::absoluteTransverseImpactParameter(*tt,pv);

                    myTrack.set_values(tt->track().pt(), tt->track().eta(), tt->track().phi(),  tt->track().dz(pv.position()), tt->track().dxy(pv.position()), distance, 
                    m.significance(), seedPosition.x(), seedPosition.y(), seedPosition.z(), seedPositionErr.cxx(), seedPositionErr.cyy(), seedPositionErr.czz(), 
                    ttPoint.x(),  ttPoint.y(),  ttPoint.z(),  ttPointErr.cxx(),  ttPointErr.cyy(),  ttPointErr.czz(), dotprodTrack, dotprodSeed );
                    myTrack.set_index(-1);
                    myTrack.set_distances(PCAseedFromPV, PCAtrackFromPV);
                    myTrack.set_vars(masses[tt-selectedTracks.begin()],t_ip2d.second.value() , t_ip2d.second.significance(),
                    t_ip.second.value() , t_ip.second.significance(), dotprodTrackSeed2D, dotprodTrackSeed3D, dotprodTrackSeed2DV, dotprodTrackSeed3DV ); 
                    
                    Line::PositionType pos(pvp);
                    Line::DirectionType dir(direction);
                    Line::DirectionType pairMomentumDir(PairMomentum);
                    Line jetLine(pos,dir);   
                    Line PCAMomentumLine(cp,pairMomentumDir);
                    float PCA_JetAxis_dist=jetLine.distance(cp).mag();
                    float dotprodMomenta=PairMomentum.unit().dot(direction.unit());
                    float dEta=std::fabs(PCA_pv.eta()-jeteta);
                    float dPhi=std::fabs(PCA_pv.phi()-jetphi);
                    
                    myTrack.setSeedMass(masses[it-selectedTracks.begin()]);                    
                    myTrack.set_JetAxisVars(PCA_JetAxis_dist,dotprodMomenta,dEta,dPhi);
                    nearTracks.push_back(myTrack);
                    /* std::cout <<  nearTracks.back().dist << " dist primann " << myTrack.dist << " " << nearTracks.size() << std::endl; */
                    
                    
                    
                }
                
                
                
                
            }
            
/*             std::cout <<nearTracks[0].dist << "  dist0  "<<std::endl; 
            for(unsigned int i=0; i< nearTracks.size(); i++) {
            
            std::cout << nearTracks[i].dist << "  dist  "; 
            std::cout << nearTracks[i].dsig << "  dsig ";
            }
            std::cout <<  "  " << std::endl; 
            
            std::cout<< "sorting" << std::endl */;            
            std::sort (nearTracks.begin(), nearTracks.end(), sortfunction2());

            nearTracks.resize(20);

            for(unsigned int i=0; i< nearTracks.size(); i++) {
            std::cout << nearTracks[i].dist << "  disat  "; 
            std::cout << nearTracks[i].dsig << "  dsiag ";
            }
            std::cout <<  "  " << std::endl;  
            
            std::cout << "check the seed again before " << it->track().phi() << " " << it->track().pt() << " " << it->track().eta() << std::endl;
            
            std::cout << "check the seed again before " << ipSigned.second.significance()<< " inserted   " << std::endl;
            SortedSeedsMap.insert(std::make_pair(-ipSigned.second.significance(), std::make_pair(&(*it), nearTracks)));
            
            }
            
            
            
        }
        }
        
		unsigned int seeds_max_counter=0;
        for(std::multimap<double,std::pair<const reco::TransientTrack*,const std::vector<trackVars2>>>::const_iterator im = SortedSeedsMap.begin(); im != SortedSeedsMap.end(); im++){
            
            if(seeds_max_counter>=10) break;

            std::cout<<-(im->first)<< "  sorted sips" <<std::endl; 			
            std::cout << "check the seed again before filling" << im->second.first->track().phi() << " " << im->second.first->track().pt() << " " << im->second.first->track().eta() << std::endl;
            std::pair<bool,Measurement1D> ipSigned = IPTools::signedImpactParameter3D(*im->second.first,direction, pv);        
            std::pair<bool,Measurement1D> ip2dSigned = IPTools::signedTransverseImpactParameter(*im->second.first,direction, pv);  
            std::pair<bool,Measurement1D> ip = IPTools::absoluteImpactParameter3D(*im->second.first, pv);        
            std::pair<bool,Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(*im->second.first, pv);	
            std::cout << "new seed : 2d+3d sign." << ip2d.second.significance() << "  " << ip.second.significance() <<"  "<< std::endl;
            std::cout << "new seed : in jet : " << jet_pt[seed_jet_match_i] <<"  "<<jet_eta[seed_jet_match_i]<<"  "<<jet_phi[seed_jet_match_i]<<"  "<< std::endl;

            if (m_computeProbabilities) {
                    
                    
                    const reco::Track & track = im->second.first->track();
                    
                    //probability with 3D ip
                    
                    std::pair<bool,double> probability = m_probabilityEstimator->probability(0, 0,ip.second.significance(),track,Jets[seed_jet_match_i],pv);
                    double prob3D=(probability.first ? probability.second : -1.);
                    std::cout<<prob3D<<std::endl;

                    //probability with 2D ip
                    
                    probability = m_probabilityEstimator->probability(0, 1,ip2d.second.significance(),im->second.first->track(),Jets[seed_jet_match_i],pv);
                    double prob2D=(probability.first ? probability.second : -1.);
                    std::cout<<prob2D<<std::endl;
                    
                    
                    seed_3D_TrackProbability[seeds_max_counter]=prob3D;
                    seed_2D_TrackProbability[seeds_max_counter]=prob2D;
                 
            } 
            
            seed_pt[seeds_max_counter]=im->second.first->track().pt();
            seed_eta[seeds_max_counter]=im->second.first->track().eta();
            seed_phi[seeds_max_counter]=im->second.first->track().phi();
            seed_mass[seeds_max_counter]=im->second.second.at(0).seedMass;
            seed_dz[seeds_max_counter]=im->second.first->track().dz(pv.position());
            seed_dxy[seeds_max_counter]=im->second.first->track().dxy(pv.position());
            seed_3D_ip[seeds_max_counter]=ip.second.value();
            seed_3D_sip[seeds_max_counter]=ip.second.significance();
            seed_2D_ip[seeds_max_counter]=ip2d.second.value();
            seed_2D_sip[seeds_max_counter]=ip2d.second.significance();
            seed_3D_signedIp[seeds_max_counter]=ipSigned.second.value();
            seed_3D_signedSip[seeds_max_counter]=ipSigned.second.significance();
            seed_2D_signedIp[seeds_max_counter]=ip2dSigned.second.value();
            seed_2D_signedSip[seeds_max_counter]=ip2dSigned.second.significance();
            seed_JetMatch[seeds_max_counter]=seed_jet_match_i;			
            seed_chi2reduced[seeds_max_counter]=im->second.first->track().normalizedChi2();
            seed_nPixelHits[seeds_max_counter]=im->second.first->track().hitPattern().numberOfValidPixelHits();
            seed_nHits[seeds_max_counter]=im->second.first->track().hitPattern().numberOfValidHits();

            std::pair<double, Measurement1D> jet_distance =IPTools::jetTrackDistance(*im->second.first, direction, pv);
            seed_jetAxisDistance[seeds_max_counter]=std::fabs(jet_distance.second.value());

            std::cout<<jet_distance.second.value()<<"      distance"<<jet_distance.second.significance()<<std::endl;            
            TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(im->second.first->impactPointState(),pv, direction,im->second.first->field());
            if (closest.isValid()) seed_jetAxisDlength[seeds_max_counter]=(closest.globalPosition() - pvp).mag(); 
            else seed_jetAxisDlength[seeds_max_counter]= -99; 


            std::cout <<  im->second.second.size() << std::endl; 	
            for(unsigned int i=0; i< im->second.second.size(); i++) {

            if((seeds_max_counter*20+i)>=200) break;
			nearTracks_nTracks.push_back(20);
            nearTracks_Nvtx.push_back(count_seeds);
            std::cout << im->second.second.at(i).dist << "  dist  "; 
            std::cout << im->second.second.at(i).dsig << "  dsig ";
            nearTracks_pt[seeds_max_counter*20+i]=im->second.second.at(i).pt;
            nearTracks_eta[seeds_max_counter*20+i]=im->second.second.at(i).eta;
            nearTracks_phi[seeds_max_counter*20+i]=im->second.second.at(i).phi;
            nearTracks_dz[seeds_max_counter*20+i]=im->second.second.at(i).dz;
            nearTracks_dxy[seeds_max_counter*20+i]=im->second.second.at(i).dxy;
            nearTracks_mass[seeds_max_counter*20+i]=im->second.second.at(i).mass;
            nearTracks_3D_ip[seeds_max_counter*20+i]=im->second.second.at(i).t3Dip;
            nearTracks_3D_sip[seeds_max_counter*20+i]=im->second.second.at(i).t3Dsip;
            nearTracks_2D_ip[seeds_max_counter*20+i]=im->second.second.at(i).t2Dip;
            nearTracks_2D_sip[seeds_max_counter*20+i]=im->second.second.at(i).t2Dsip;
            nearTracks_PCAdist[seeds_max_counter*20+i]=im->second.second.at(i).dist;
            nearTracks_PCAdsig[seeds_max_counter*20+i]=im->second.second.at(i).dsig;
            nearTracks_PCAonSeed_x[seeds_max_counter*20+i]=im->second.second.at(i).PCA_sx;
            nearTracks_PCAonSeed_y[seeds_max_counter*20+i]=im->second.second.at(i).PCA_sy;
            nearTracks_PCAonSeed_z[seeds_max_counter*20+i]=im->second.second.at(i).PCA_sz;
            nearTracks_PCAonSeed_xerr[seeds_max_counter*20+i]=im->second.second.at(i).PCA_sxerr;
            nearTracks_PCAonSeed_yerr[seeds_max_counter*20+i]=im->second.second.at(i).PCA_syerr;
            nearTracks_PCAonSeed_zerr[seeds_max_counter*20+i]=im->second.second.at(i).PCA_szerr;
            nearTracks_PCAonTrack_x[seeds_max_counter*20+i]=im->second.second.at(i).PCA_tx;
            nearTracks_PCAonTrack_y[seeds_max_counter*20+i]=im->second.second.at(i).PCA_ty;
            nearTracks_PCAonTrack_z[seeds_max_counter*20+i]=im->second.second.at(i).PCA_tz;
            nearTracks_PCAonTrack_xerr[seeds_max_counter*20+i]=im->second.second.at(i).PCA_txerr;
            nearTracks_PCAonTrack_yerr[seeds_max_counter*20+i]=im->second.second.at(i).PCA_tyerr;
            nearTracks_PCAonTrack_zerr[seeds_max_counter*20+i]=im->second.second.at(i).PCA_tzerr;
            nearTracks_dotprodTrack[seeds_max_counter*20+i]=im->second.second.at(i).dotprodTrack;
            nearTracks_dotprodSeed[seeds_max_counter*20+i]=im->second.second.at(i).dotprodSeed;
            nearTracks_dotprodTrackSeed2D[seeds_max_counter*20+i]=im->second.second.at(i).dotprodTrackSeed2D;
            nearTracks_dotprodTrackSeed3D[seeds_max_counter*20+i]=im->second.second.at(i).dotprodTrackSeed3D;
            nearTracks_dotprodTrackSeedVectors2D[seeds_max_counter*20+i]=im->second.second.at(i).dotprodTrackSeedVectors2D;
            nearTracks_dotprodTrackSeedVectors3D[seeds_max_counter*20+i]=im->second.second.at(i).dotprodTrackSeedVectors3D;

            nearTracks_PCAonSeed_pvd[seeds_max_counter*20+i]=im->second.second.at(i).seedPCA_pv;
            nearTracks_PCAonTrack_pvd[seeds_max_counter*20+i]=im->second.second.at(i).trackPCA_pv;
            std::cout <<  "after lambda  " << seeds_max_counter*20+i<< std::endl;
            std::cout << im->second.second.at(i).PCA_JetAxis_distance << "  " << im->second.second.at(i).PCAPair_Jet_dotprod << " "
            << im->second.second.at(i).PCAAxis_JetAxis_DEta <<" " <<im->second.second.at(i).PCAAxis_JetAxis_DPhi << std::endl;

            nearTracks_PCAjetAxis_dist[seeds_max_counter*20+i]=im->second.second.at(i).PCA_JetAxis_distance;
            nearTracks_PCAjetMomenta_dotprod[seeds_max_counter*20+i]=im->second.second.at(i).PCAPair_Jet_dotprod;


            nearTracks_PCAjetDirs_DEta[seeds_max_counter*20+i]=im->second.second.at(i).PCAAxis_JetAxis_DEta;
            nearTracks_PCAjetDirs_DPhi[seeds_max_counter*20+i]=im->second.second.at(i).PCAAxis_JetAxis_DPhi;


            }
            seeds_max_counter++;
            std::cout <<  im->second.second.size()<<  "fine traccia seed  " << std::endl;
            std::cout <<  "  " << std::endl;  

        }
        std::cout <<  " fine jet  " << std::endl;
        jetNseeds=count_seeds;
        SortedSeedsMap.clear();
        std::cout <<  " mmclear " << std::endl;
        tree->Fill();
        std::cout <<  " filled " << std::endl;
        }
        }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
std::cout <<  " fine evento  " << std::endl;


}



int AnalyzerSignedIP_MINIAOD_wArr::funzionePippo (int vivaLaVita) {
return 1;
}


std::vector<float> AnalyzerSignedIP_MINIAOD_wArr::firstHadronInChain(pat::PackedGenParticle gParticle, const reco::Candidate* mother){
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


std::vector<int> AnalyzerSignedIP_MINIAOD_wArr::seedToAllJets20Matching(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi, std::vector<double> jet_pt)
{
  float angular_distance =  0.4;
  float dist_temp;
  std::vector<int> index_array;
  for (unsigned int i = 0; i <jet_eta.size(); i++) {
     if (jet_pt[i]>20) {
     dist_temp=std::sqrt(std::pow(jet_eta[i]-seed.track().eta(),2) + std::pow(jet_phi[i]-seed.track().phi(),2) );
     if (dist_temp<=angular_distance) {
         std::cout<< "jet macthibg function"<< jet_eta[i]<<"  " <<seed.track().eta()<<"  " << jet_phi[i]<<"  " <<seed.track().phi()<<std::endl;
//          angular_distance=dist_temp;
         index_array.push_back(i);
    }}
      
}
return index_array;
}



int AnalyzerSignedIP_MINIAOD_wArr::seedToJetMatching(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi) {
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

int AnalyzerSignedIP_MINIAOD_wArr::seedToJet30Matching(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi, std::vector<double> jet_pt) {
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

std::pair<int,float> AnalyzerSignedIP_MINIAOD_wArr::ClosestJet(reco::TransientTrack seed, std::vector<double> jet_eta, std::vector<double> jet_phi) {
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





double AnalyzerSignedIP_MINIAOD_wArr::trackEff(std::vector<pat::PackedGenParticle> genP, std::vector<reco::TransientTrack> tracks, double &BChain_eff, double &DChain_eff, double &B_eff,
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



bool AnalyzerSignedIP_MINIAOD_wArr::genLevelMM(std::vector<pat::PackedGenParticle> genP, std::vector<reco::TransientTrack> tracks, std::vector<int>& genPnumbers, std::vector<double>& assigned_chisquares){ 
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


int AnalyzerSignedIP_MINIAOD_wArr::genMap(std::vector<pat::PackedGenParticle> genP, std::vector<reco::TransientTrack> tracks, std::vector<int>& genPnumbers, std::vector<double>& allChi){

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






int AnalyzerSignedIP_MINIAOD_wArr::genLevel(std::vector<pat::PackedGenParticle> genP, reco::TransientTrack seed, std::vector<int>& MomFlav, std::vector<int>& BChain, std::vector<double>& allChi){

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

template<class T> void AnalyzerSignedIP_MINIAOD_wArr::insert_inArray(T (&tArray)[10], unsigned int size, unsigned int index, double element)
{
    if  (index>(size-1)) return;
    
    if (index==(size-1))
    {
        tArray[size-1]=element;
    }
    else{
        tArray[size-1]=tArray[size-2];
        insert_inArray(tArray, size-1, index, element);
    }
    return;
}

template<class T> void AnalyzerSignedIP_MINIAOD_wArr::insert_inArray200(T (&tArray)[200], unsigned int size, unsigned int index, double element)
{
	std::cout<< "insert_inArray200"  << std::endl;
	std::cout<< size << "  " << index << std::endl;
    
    if  (index>(size-1)) return;
    
    unsigned int block_index=index-(index%20);
    
    if (tArray[block_index]==0 && tArray[block_index+1]==0 &&  tArray[block_index+2]==0) tArray[index]=element;
    
    else{
        for (unsigned int j = 180+(index%20) ; j>index; j=j-20)
        {
            std::cout<< j << std::endl;
            tArray[j]=tArray[j-20];
         }
        tArray[index]=element;
    }
    return;
}


// Jet Flvour
// BB 5
// B  5
// C  4
// UD 1
// S  3
// G  21
// NAN 0


int AnalyzerSignedIP_MINIAOD_wArr::jet_flavour_function(const pat::Jet& jet, bool usePhysForLightAndUndefined) {
    int hflav = abs(jet.hadronFlavour());
    int pflav = abs(jet.partonFlavour());
    int physflav = 0;
    if(jet.genParton()) physflav=abs(jet.genParton()->pdgId());
    std::size_t nbs = jet.jetFlavourInfo().getbHadrons().size();
    std::size_t ncs = jet.jetFlavourInfo().getcHadrons().size();

    if(hflav == 5) { //B jet
        if(nbs > 1) return 5;
        else if(nbs == 1) {
            return 5;
        }
        else {
            if(usePhysForLightAndUndefined){
                if(physflav == 21) return 21;
                else if(physflav == 3) return 3;
                else if(physflav == 2 || physflav ==1) return 1;
                else return 0;
            }
            else return 0;
        }
    }
    else if(hflav == 4) { //C jet
        return 4;
    }
    else { //not a heavy jet
        if(std::abs(pflav) == 4 || std::abs(pflav) == 5 || nbs || ncs) {
            if(usePhysForLightAndUndefined){
                if(physflav == 21) return 21;
                else if(physflav == 3) return 3;
                else if(physflav == 2 || physflav ==1) return 1;
                else return 0;
            }
            else return 0;
        }
        else if(usePhysForLightAndUndefined){
            if(physflav == 21) return 21;
            else if(physflav == 3) return 3;
            else if(physflav == 2 || physflav ==1) return 1;
            else return 0;
        }
        else {
            if(pflav == 21) return 21;
            else if(pflav == 3) return 3;
            else if(pflav == 2 || pflav ==1) return 1;
            else return 0;
        }
    }
}

void AnalyzerSignedIP_MINIAOD_wArr::checkEventSetup(const edm::EventSetup & iSetup) {
  using namespace edm;
  using namespace edm::eventsetup;

   const EventSetupRecord & re2D= iSetup.get<BTagTrackProbability2DRcd>();
   const EventSetupRecord & re3D= iSetup.get<BTagTrackProbability3DRcd>();
   unsigned long long cacheId2D= re2D.cacheIdentifier();
   unsigned long long cacheId3D= re3D.cacheIdentifier();

   if(cacheId2D!=m_calibrationCacheId2D || cacheId3D!=m_calibrationCacheId3D  )  //Calibration changed
   {
     //iSetup.get<BTagTrackProbabilityRcd>().get(calib);
     ESHandle<TrackProbabilityCalibration> calib2DHandle;
     iSetup.get<BTagTrackProbability2DRcd>().get(calib2DHandle);
     ESHandle<TrackProbabilityCalibration> calib3DHandle;
     iSetup.get<BTagTrackProbability3DRcd>().get(calib3DHandle);

     const TrackProbabilityCalibration *  ca2D= calib2DHandle.product();
     const TrackProbabilityCalibration *  ca3D= calib3DHandle.product();

     m_probabilityEstimator.reset(new HistogramProbabilityEstimator(ca3D,ca2D));

   }
   m_calibrationCacheId3D=cacheId3D;
   m_calibrationCacheId2D=cacheId2D;
}

// ------------ method called once each job just before starting event loop  ------------
void 
AnalyzerSignedIP_MINIAOD_wArr::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AnalyzerSignedIP_MINIAOD_wArr::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AnalyzerSignedIP_MINIAOD_wArr::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzerSignedIP_MINIAOD_wArr);
