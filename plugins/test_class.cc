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
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include <stdint.h>
#include <stdio.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs. 
class test_class : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit test_class(const edm::ParameterSet&);
      ~test_class();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      // ----------member data ---------------------------
      
//       int funzionePippo (int vivaLaVita);
      
      edm::Service<TFileService> file;
     
      TTree *tree;
      
      int evt=0;
      int lumi=0;
      int run=0;
      int min3DIPValue=0;
      int min3DIPSignificance=2;
      int max3DIPValue=100;
      int max3DIPSignificance=10;
//       std::vector<double> seed_3D_ip;
//       std::vector<double> seed_3D_sip;
      std::vector<float> res;


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
test_class::test_class(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
    usesResource("TFileService");
    tree=file->make<TTree>("tree","tree");
    
    tree->Branch("lumi",&lumi, "lumi/I");
    tree->Branch("evt",&evt, "evt/I");
    tree->Branch("run",&run, "run/I");  
/*    
    tree->Branch("seed_3D_ip",&seed_3D_ip);
    tree->Branch("seed_3D_sip",&seed_3D_sip);*/
    tree->Branch("seed_Test",&res);
    
    
    
}


test_class::~test_class()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}



//
// member functions
//

// ------------ method called for each event  ------------
void
test_class::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    
          std::vector<trackVars2> nearTracks;
      trackVars2 myTrack;

    evt=iEvent.id().event();
    lumi=iEvent.id().luminosityBlock();
    run=iEvent.id().run();    
    std::cout<< "nearTracks" <<  "   "  << nearTracks.size() << "before fill" << std::endl;
    nearTracks.push_back(myTrack);
    std::cout <<  nearTracks.back().dist << " dist prima " << myTrack.dist << " " << nearTracks.size() << std::endl;
    printf("%zu", SIZE_MAX);
   tree->Fill();
}
            
            
            
void 
test_class::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
test_class::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
test_class::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(test_class);
