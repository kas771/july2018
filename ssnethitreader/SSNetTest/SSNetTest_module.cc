/**
 * @file   SSNetTest_module.cc
 * @brief  Takes an SSNetHits.root file 
 * @author ksutton
 * 
*/


// LArSoft includes
 #include "lardataobj/Simulation/SimChannel.h"
 #include "larsim/Simulation/LArG4Parameters.h"
 #include "lardataobj/RecoBase/Hit.h"
 #include "lardataobj/RecoBase/Cluster.h"
 #include "larcore/Geometry/Geometry.h"
 #include "larcore/Geometry/GeometryCore.h"
 #include "nusimdata/SimulationBase/MCParticle.h"
 #include "nusimdata/SimulationBase/MCTruth.h"
 #include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
 #include "lardataobj/RecoBase/Vertex.h"
 #include "lardataobj/RecoBase/PFParticle.h"
 #include "lardataobj/RecoBase/Track.h"
 #include "lardataobj/RecoBase/Shower.h"
 
// Framework includes
 #include "canvas/Utilities/Exception.h"
 #include "art/Framework/Core/EDAnalyzer.h"
 #include "art/Framework/Principal/Event.h"
 #include "art/Framework/Principal/Handle.h"
 #include "art/Framework/Services/Registry/ServiceHandle.h"
 #include "art/Framework/Services/Optional/TFileService.h"
 #include "art/Framework/Core/ModuleMacros.h"
 #include "canvas/Persistency/Common/FindManyP.h"

 // utility libraries
 #include "messagefacility/MessageLogger/MessageLogger.h"
 #include "fhiclcpp/ParameterSet.h"
 #include "fhiclcpp/types/Table.h"
 #include "fhiclcpp/types/Atom.h"
 #include "cetlib/pow.h" // cet::sum_of_squares()

// ROOT includes. Note: To look up the properties of the ROOT classes,
 // use the ROOT web site; e.g.,// <http://root.cern.ch/root/html532/ClassIndex.html>
 #include "TH1.h"
 #include "TH2.h"
 #include "TTree.h"
 #include "TLorentzVector.h"
 #include "TVector3.h"

 // C++ Includes
 #include <map>
 #include <vector>
 #include <string>
 #include <cmath>
 #include <memory>

namespace{
double DetectorDiagonal(geo::GeometryCore const& geom);
//bool inROI(double* vertex_pos, double* hit_pos);//checks if hit is in ROI
}

namespace lar {
namespace example {
  class SSNetTest : public art::EDAnalyzer
  {
  public:
	struct Config {
      
      // save some typing:
      	using Name = fhicl::Name;
     	using Comment = fhicl::Comment;
      // one Atom for each parameter
        fhicl::Atom<art::InputTag> SimulationLabel {Name("SimulationLabel"),
      		Comment("tag of the input data product with the detector simulation information")
        };
      
      fhicl::Atom<int> PxType {
      Name("PxType"),
      Comment("")
      };
                                                                 
     fhicl::Atom<art::InputTag> SSHitLabel {
      Name("SSHitLabel"),
      Comment("tag of the input SSnet data product with reconstructed hits")
      };
                                          
     fhicl::Atom<art::InputTag> HitLabel {
      Name("HitLabel"),
      Comment("tag of the input data product with reconstructed hits")
      };
   
     /*fhicl::Atom<art::InputTag> HitLabel[1] {
      Name("ShrHitLabel"),
      Comment("tag of the input data product with reconstructed hits")
      };*/
    
       
      fhicl::Atom<art::InputTag> ClusterLabel {
      Name("ClusterLabel"),
      Comment("tag of the input data product with reconstructed clusters")
      };
                                                                                                                                          
      fhicl::Atom<int> PDGcode {
      Name("PDGcode"),
      Comment("particle type (PDG ID) of the primary particle to be selected")
      };     //                                                                                                                                                                         
     fhicl::Atom<double> BinSize {
	Name("BinSize"),                                                                                                                                                                    Comment("dx [cm] used for the dE/dx calculation")                                                                                                                                                                                       };
     fhicl::Atom<art::InputTag> VertexLabel {
	Name("Vertex"),
	Comment("currently using Pandora vertex")                                                                                                                                                                                       };
     fhicl::Atom<art::InputTag> PFPLabel {
	Name("PFP"),
	Comment("tag of the PFParticle producer")                                                                                                                                                                                       };
     fhicl::Atom<art::InputTag> TrackLabel {
	Name("TrackLabel"),
	Comment("tag of the track producer")                                                                                                                                                                                       };
     fhicl::Atom<art::InputTag> ShowerLabel {
	Name("ShowerLabel"),
	Comment("tag of the shower producer")                                                                                                                                                                                       };
    
};//struct
    
     using Parameters = art::EDAnalyzer::Table<Config>;

    explicit SSNetTest(Parameters const& config); 
  
    virtual void beginJob() override;
    virtual void beginRun(const art::Run& run) override;
    virtual void analyze (const art::Event& event) override;

  private:

    // The parameters we'll read from the .fcl file.
    art::InputTag fSimulationProducerLabel; ///< The name of the producer that tracked simulated particles through the detector
    art::InputTag fHitProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fSSHitProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fClusterProducerLabel;    ///< The name of the producer that created clusters
    int fSelectedPDG;                     ///< PDG code of particle we'll focus on
    double fBinSize;                      ///< For dE/dx work: the value of dx. 
    int fPxType;
    art::InputTag fVertexProducerLabel; ///<name of the producer that created the vertices
    art::InputTag fPFPLabel; ///name of PFParticle producer
    art::InputTag fTrackProducerLabel; ///name of track producer
    art::InputTag fShowerProducerLabel; ///name of track producer


// vector of shower-like hit indices
//   std::vector<size_t> _shrhits;

 // Pointers to the histograms we'll create. 
    TH1D* fPDGCodeHist;     ///< PDG code of all particles
    TH1D* fMomentumHist;    ///< momentum [GeV] of all selected particles
    TH1D* fTrackLengthHist; ///< true length [cm] of all selected particles

    // The n-tuples we'll create.
    TTree* fmytree;     ///< tuple with simulated data

    // The variables that will go into the n-tuple.
    int fEvent;     ///< number of the event being processed
    int fRun;       ///< number of the run being processed
    int fSubRun;    ///< number of the sub-run being processed
    int fSimPDG;       ///< PDG ID of the particle begin processed
    int fSimTrackID;   ///< GEANT ID of the particle begin processed

    Float_t _xpos, _ypos, _zpos; // xyz of vertex

    double fStartXYZT[4]; ///< (x,y,z,t) of the true start of the particle
    double fEndXYZT[4];   ///< (x,y,z,t) of the true end of the particle
    double fStartPE[4];   ///< (Px,Py,Pz,E) at the true start of the particle
    double fEndPE[4];     ///< (Px,Py,Pz,E) at the true end of the particle

    geo::GeometryCore const* fGeometry;    
     }; // class SSNetTest

SSNetTest::SSNetTest(Parameters const& config) // Initialize member data here.
	: EDAnalyzer(config)
 	, fSimulationProducerLabel(config().SimulationLabel())
    	, fHitProducerLabel       (config().HitLabel())
    	, fSSHitProducerLabel       (config().SSHitLabel())
    	, fClusterProducerLabel   (config().ClusterLabel())
    	, fSelectedPDG            (config().PDGcode())
    	, fBinSize                (config().BinSize())
//	, fVertexProducerLabel (config().Vertex())
       // , fLArCVLocation  	  (config().LArCVLocation())
        , fPxType          	  (config().PxType())
        , fVertexProducerLabel    (config().VertexLabel()) 
	, fPFPLabel		  (config().PFPLabel())  
	, fTrackProducerLabel		  (config().TrackLabel())  	
	, fShowerProducerLabel		  (config().ShowerLabel())  	
	// fInHitProducer   = p.get<std::string>("InHitProducer","gaushit");
        //fPxThresholdHigh = p.get<double>     ("PxThresholdHigh"        );
        //fPxThresholdLow  = p.get<double>     ("PxThresholdLow"         );
       // produces< std::vector< recob::Hit > >(Form("shrhit%zu",(size_t)(fPxThresholdHigh*100.)));
        //produces< std::vector< recob::Hit > >(Form("shrhit%zu",(size_t)(fPxThresholdLow*100.) ));

  {
    // get a pointer to the geometry service provider
     fGeometry = lar::providerFrom<geo::Geometry>();
  }

  
  //-----------------------------------------------------------------------
  void SSNetTest::beginJob()
  {
    // Get the detector length, to determine the maximum bin edge of one
    // of the histograms.
    const double detectorLength = DetectorDiagonal(*fGeometry);

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;
  
    fPDGCodeHist     = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
    fMomentumHist    = tfs->make<TH1D>("mom",     ";particle Momentum (GeV);",    100, 0.,    10.);
    fTrackLengthHist = tfs->make<TH1D>("length",  ";particle track length (cm);", 200, 0, detectorLength);
 
    // Define our n-tuples, which are limited forms of ROOT
    // TTrees. Start with the TTree itself.
    fmytree     = tfs->make<TTree>("SSNetTestSimulation",    "SSNetTestSimulation");

    // Define the branches (columns) of our simulation n-tuple. When
    // we write a variable, we give the address of the variable to
    // TTree::Branch.
    fmytree->Branch("Event",       &fEvent,          "Event/I");
    fmytree->Branch("SubRun",      &fSubRun,         "SubRun/I");
    fmytree->Branch("Run",         &fRun,            "Run/I");
    fmytree->Branch("TrackID",     &fSimTrackID,     "TrackID/I");
    fmytree->Branch("PDG",         &fSimPDG,         "PDG/I");
    // When we write arrays, we give the address of the array to
    // TTree::Branch; in C++ this is simply the array name.
    fmytree->Branch("StartXYZT",   fStartXYZT,       "StartXYZT[4]/D");
    fmytree->Branch("EndXYZT",     fEndXYZT,         "EndXYZT[4]/D");
    fmytree->Branch("StartPE",     fStartPE,         "StartPE[4]/D");
    fmytree->Branch("EndPE",       fEndPE,           "EndPE[4]/D");
    // For a variable-length array: include the number of bins.
    //fmytree->Branch("NdEdx",       &fSimNdEdxBins,   "NdEdx/I");
    // ROOT can understand fairly well vectors of numbers (and little more)
    //fmytree->Branch("dEdx",        &fSimdEdxBins);

}

   
  //-----------------------------------------------------------------------
  void SSNetTest::beginRun(const art::Run& /*run*/)
  {
 //   std::cout<<"flag 1"<<std::endl;
    art::ServiceHandle<sim::LArG4Parameters> larParameters;
   // fElectronsToGeV = 1./larParameters->GeVToElectrons();
  //std::cout<<"flag 2"<<std::endl;
}

  //-----------------------------------------------------------------------
  void SSNetTest::analyze(const art::Event& event) 
  {
    std::cout<<"------- starting event --------"<<std::endl;
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

   //read in PFParticles
   art::ValidHandle< std::vector<recob::PFParticle> > PFPHandle
   = event.getValidHandle<std::vector<recob::PFParticle>>
        (fPFPLabel);

  //read in clusters
  art::ValidHandle< std::vector<recob::Cluster> > clusterHandle
      = event.getValidHandle<std::vector<recob::Cluster>>
        (fClusterProducerLabel);

  //read in tracks
  art::ValidHandle< std::vector<recob::Track> > trackHandle
      = event.getValidHandle<std::vector<recob::Track>>
        (fTrackProducerLabel);

 //read in showers
  art::ValidHandle< std::vector<recob::Shower> > showerHandle
      = event.getValidHandle<std::vector<recob::Shower>>
        (fShowerProducerLabel);
/*
 //read in vertices
  art::ValidHandle< std::vector<recob::Vertex> > vertexHandle
      = event.getValidHandle<std::vector<recob::Vertex>>
        (fVertexProducerLabel);
*/
 // grab vertices, tracks, and showers associated to PFParticle
  const art::FindManyP<recob::Vertex> pfp_vtx_assn_v(PFPHandle, event, fPFPLabel);  
  const art::FindManyP<recob::Track  > pfp_trk_assn_v(PFPHandle, event, fPFPLabel);
  const art::FindManyP<recob::Shower> pfp_shr_assn_v(PFPHandle, event, fPFPLabel);
/* 
// grab showers  associated to vertices
   art::FindManyP<recob::Shower> vtx_shr_assn_v(vertexHandle, event, fVertexProducerLabel);
*/ 
// grab hits associated to track
  art::FindManyP<recob::Hit> trk_hit_assn_v(trackHandle, event, fTrackProducerLabel);
 
  // grab hits associated to shower
  art::FindManyP<recob::Hit> shr_hit_assn_v(showerHandle, event, fShowerProducerLabel);


/*
 *Perform hit-matching between hits and sshits
 */

//load hits
 art::Handle< std::vector<recob::Hit> > hitHandle;
 if (!event.getByLabel(fHitProducerLabel, hitHandle)) return;
 
//load sshits  
 art::Handle< std::vector<recob::Hit> > sshitHandle;
 if (!event.getByLabel(fSSHitProducerLabel, sshitHandle)) return;

// map associated peak time of a hit to a pointer for the sshit
std::map<int, const recob::Hit* > _hitmap;
 
 //int m= 0; //number of matched sshits
int n = sshitHandle->size(); //the total number of sshits
// For every Hit:
    for ( auto const& hit : (*hitHandle) )
      {
	//int plane_h = hit.View();
	
	//get the time and channel
	float t = hit.PeakTime();
	auto hitChannelNumber = hit.Channel();
//	std::cout<<"the channel number = "<<hitChannelNumber<<std::endl;
     	
	//int idx = 0; //the sshit index
	//for every SSHit
	for ( auto const& sshit : (*sshitHandle) ){
		//if the sshit isn't already matched
		auto search = _hitmap.find(sshit.PeakTime());
		if ( search != _hitmap.end() ) continue;
		//if the channel matches
		if (sshit.Channel() == hitChannelNumber){	
			//if the peak time matches
			if (sshit.PeakTime() == t){
				//save the peak time as the key and a pointer to the sshit
				_hitmap[t] = &sshit;
				//m++; //increment number of matched hits
			}//peak time
		}//channel
	}//for each SSHits
 } // for each Hit

//std::cout<<"the number of items in the map = "<<_hitmap.size()<<std::endl;
std::cout<<"number of sshits = "<<n<<", number of matches = "<<_hitmap.size()<<std::endl;
 
/*
 * loop over PFParticles and find 1 shower 1 track topologies
 *
 */

  std::vector<art::Ptr<recob::Track>> my_trks; //vector of pointers to tracks in case of 1 shower 1 track events
  std::vector<art::Ptr<recob::Shower>> my_shrs; //vector of pointers to tracks in case of 1 shower 1 track events
  std::vector<double*> my_vtxs; //vector of pointers to the first element in a double for the vertex XYZ position (3D)
 
 // std::cout<<"number of PFParticles: "<<PFPHandle->size()<<std::endl;
  //for each PFP
  for ( size_t pfp_index = 0; pfp_index != PFPHandle->size(); ++pfp_index ){
	auto const& pfp = PFPHandle->at(pfp_index);
	int pdg = pfp.PdgCode();
	
	//only want primary particles
	if (pfp.IsPrimary() == false){continue;}	
	
	//check for proton
	if (pdg == fSelectedPDG){
		//std::cout<<"found proton"<<std::endl;
		std::cout<<"PDG  = "<<pdg<<std::endl;
	}
	
	//check parent of particle
	//std::cout<<"The parent particle index is : "<<pfp.Parent()<<std::endl;

	//check number of daughters
	auto const& daughters = pfp.Daughters(); //get collection of daughter particles
	if (daughters.size() != 2){continue;} //want precisely 1 track 1 shower topology

	bool has_shower = false; //is true if at least one daughter is a shower
	bool has_track = false;	//is true if at least one daughter is a track
	//size_t pfp_daughter_trk_ind = -1; //index of PFP in PFPHandle corresponding to daughter track
	//size_t pfp_daughter_shr_ind = -1; //index of PFP in PFPHandle corresponding to daughter shower
	art::Ptr<recob::Track> this_track; //pointer to track
	art::Ptr<recob::Shower> this_shower; //pointer to shower

	std::cout<<"number of daughters = "<<daughters.size()<<std::endl;
	for (size_t d = 0; d != daughters.size(); ++d ){//for each daughter
		size_t daughter_index = daughters.at(d);
		//std::cout<<"the daughter particle(s) ID :"<<daughter_index<<std::endl;
		auto const& daughterPFP = PFPHandle->at(daughter_index);
		
		//check if the daughter is a track
		auto const& pfp_trk_v = pfp_trk_assn_v.at(daughter_index);
		std::cout<<"the number of tracks associated to this daughter = "<<pfp_trk_v.size()<<std::endl;
		if(pfp_trk_v.size() == 1){
			has_track = true;
			//get the track
			this_track = pfp_trk_v.at(0);
		}
		
		//check if the daughter is a shower
		auto const& pfp_shr_v = pfp_shr_assn_v.at(daughter_index);
		std::cout<<"the number of showers associated to this daughter = "<<pfp_shr_v.size()<<std::endl;
		if(pfp_shr_v.size() == 1){
			has_shower = true;
			this_shower = pfp_shr_v.at(0);
			//pfp_daughter_shr_ind = daughter_index;
	}
	
		std::cout<<"the daughter PDG code is :"<<daughterPFP.PdgCode()<<std::endl;
	}//for each daughter

	//skip if there isn't exactly one track and one shower
	if (has_shower == false || has_track == false){
		//std::cout<<"this event doesn't have one track and one shower"<<std::endl;
		continue;
	}else{
		my_trks.push_back(this_track);//add track pointer to vector
		my_shrs.push_back(this_shower);  //add shower pointer to vector
	}
	

	//if(pfp_daughter_shr_ind <= (size_t)-1 || pfp_daughter_trk_ind <= (size_t)-1){continue;}//check that track+shower vars were filled
 
	//get PFP corresponding to shower and track
	// auto const& shrPFP = PFPHandle->at(pfp_daughter_shr_ind);
	//auto const& trkPFP = PFPHandle->at(pfp_daughter_trk_ind);


	//for this event, get the vertex
	auto const& pfp_vtx_v = pfp_vtx_assn_v.at(pfp_index);

	//for each vertex
	for (size_t vtx_index = 0; vtx_index != pfp_vtx_v.size(); ++vtx_index ){
		auto vtx = *(pfp_vtx_v.at(vtx_index));
		std::cout<<"vertex ID  = "<<vtx.ID()<<std::endl;

		//get the vertex position
		double xyz[3] = {} ;
        	vtx.XYZ(&xyz[0]);
        	_xpos = xyz[0];
       		_ypos = xyz[1];
        	_zpos = xyz[2];
        	
		//store the vertex postion
		my_vtxs.push_back(xyz);
        //	std::cout<<"vertex pos xyz = "<<_xpos<<", "<<_ypos<<", "<< _zpos<<std::endl;

	} //for each vertex
}//for each PFP

std::cout<<"the number of stored vertices = "<<my_vtxs.size()<<std::endl;

/* 
 * Remove shrhits matched with hits in the shower and look at shrhits matched to the track
 */

//std::cout<<"checking shower"<<std::endl;
//if(showerHandle->size() != 0) { //skip events with no shower
  for ( size_t shr_index = 0; shr_index != showerHandle->size(); ++shr_index ){
	auto const& shr = showerHandle->at(shr_index);
	auto const& start  = shr.ShowerStart();
	for (size_t s = 0; s != my_shrs.size(); ++s ){
		auto const& current_shr = *(my_shrs.at(s));
		auto const& current_start = current_shr.ShowerStart();
		if (start == current_start){
			std::cout<<"matched showers"<<std::endl;
			std::cout<<"the shower length is  = "<<shr.Length()<<", and the opening angle is = "<<shr.OpenAngle()<<std::endl;
			
			 //get the associated hits for the shower
			 const std::vector<art::Ptr<recob::Hit> > shr_hit_v = shr_hit_assn_v.at(shr_index);
			std::cout<<"the number of hits in the shower = "<<shr_hit_v.size()<<std::endl;		 
			
			//for each hit
			 for (size_t h=0; h < shr_hit_v.size(); h++){
			 	auto const& this_hit = *(shr_hit_v.at(h));
				float this_time = this_hit.PeakTime();
				
				//if the peak time is in the hit map, remove the shrhits from the map
				if(_hitmap.count(this_time) >= 1){
					_hitmap.erase(this_time);
				}

			}//for each hit
		std::cout<<"number of remaining matched shr hits = "<<_hitmap.size()<<std::endl;

		}//if the showers match
	}//for each shower from a 1 shower 1 track topology 

	//auto const& length = shr.Length();
	//std::cout<<"shower length = "<<length<<std::endl;
	
//	auto const& start  = shr.ShowerStart();
//	std::cout<<"shower start xyz = "<<start.X()<<", "<<start.Y()<<", "<< start.Z()<<std::endl;
}//for each shower

//for each track
  for ( size_t trk_index = 0; trk_index != trackHandle->size(); ++trk_index ){
	auto const& trk = trackHandle->at(trk_index);
	auto const& start  = trk.Start();
	for (size_t tr = 0; tr != my_trks.size(); ++tr ){
		auto const& current_trk = *(my_trks.at(tr));
		auto const& current_start = current_trk.Start();
		int number_matched_trk_hits = 0;
		if (start == current_start){
			std::cout<<"matched tracks"<<std::endl;
			//std::cout<<"the shower length is  = "<<shr.Length()<<", and the opening angle is = "<<shr.OpenAngle()<<std::endl;
			
			 //get the associated hits for the track
			 const std::vector<art::Ptr<recob::Hit> > trk_hit_v = trk_hit_assn_v.at(trk_index);
			std::cout<<"the number of hits in the track = "<<trk_hit_v.size()<<std::endl;		 
			
			//for each hit
			 for (size_t h=0; h < trk_hit_v.size(); h++){
			 	auto const& this_hit = *(trk_hit_v.at(h));
				float this_time = this_hit.PeakTime();
				
				//if the peak time is in the hit map, remove the shrhits from the map
				if(_hitmap.count(this_time) >= 1){
					//_hitmap.erase(this_time);
					number_matched_trk_hits++;
				}

			}//for each hit
		//std::cout<<"number of remaining matched shr hits = "<<_hitmap.size()<<std::endl;
		std::cout<<"the number of matched shr hits in the track = "<<number_matched_trk_hits<<std::endl;
		}//if the tracks match
	}//for each track from a 1 shower 1 track topology 
}//for each track

/*
 * Look at remaining shrhits in ROI
 *
 */
/*
//for each vertex position in the vector
for ( size_t vtx_index = 0; vtx_index != my_vtxs.size(); ++vtx_index ){
	double xyz_pos[3] = {};
	xyz_pos =  *(my_vtxs.at(vtx_index));
	//double X = xyz_pos;
	//double Y = ++xyz_pos;
	
	double X, Y; 
	X = xyz_pos[0];
	Y = xyz_pos[1];
	std::cout<<"the vertex XYZ = "<<X<<Y<<std::endl;
}//loop over vertices
*/

 auto particleHandle
      = event.getValidHandle<std::vector<simb::MCParticle>>
      (fSimulationProducerLabel);
  
 if (!particleHandle){std::cout<<"missing particle handle"<<std::endl;}
    
 
std::map< int, const simb::MCParticle* > particleMap;
//std::cout<<"flag 3.1"<<std::endl;
    for ( auto const& particle : (*particleHandle) )
      {
	// For the methods you can call to get particle information,
	// see ${NUTOOLS_INC}/SimulationBase/MCParticle.h.
	fSimTrackID = particle.TrackId();
	//std::cout<<"TrackID: "<<fSimTrackID<<std::endl;
//std::cout<<"flag 3.1.1"<<std::endl;
	// Add the address of the MCParticle to the map, with the track ID as the key.
	particleMap[fSimTrackID] = &particle; 
//std::cout<<"flag 3.1.2"<<std::endl;
	// Histogram the PDG code of every particle in the event.
	fSimPDG = particle.PdgCode();
	//if(fSimPDG == fSelectedPDG){std::cout<<"Found particle with PDG: "<<fSimPDG<<std::endl;}
	fPDGCodeHist->Fill( fSimPDG );
//std::cout<<"flag 3.2"<<std::endl;


	//if ( particle.Process() == "primary"  &&  fSimPDG == fSelectedPDG )
	if (fSimPDG == fSelectedPDG ) { //for all events containing a proton
	  //std::cout<<"flag 3.2.1"<<std::endl;
	    // A particle has a trajectory, consisting of a set of
	    // 4-positions and 4-mommenta.
	    const size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();

	    // For trajectories, as for vectors and arrays, the
	    // first point is #0, not #1.
	    const int last = numberTrajectoryPoints - 1;
	    const TLorentzVector& positionStart = particle.Position(0);
	    const TLorentzVector& positionEnd   = particle.Position(last);
	    const TLorentzVector& momentumStart = particle.Momentum(0);
	    const TLorentzVector& momentumEnd   = particle.Momentum(last);
//std::cout<<"flag 3.2.2"<<std::endl;

	    // Make a histogram of the starting momentum.
	   // std::cout<<"P: "<<momentumStart.P()<<std::endl;
	  //  fMomentumHist->Fill( momentumStart.P() );
//std::cout<<"flag 3.2.2.1"<<std::endl;

	    // Fill arrays with the 4-values. (Don't be fooled by
	    // the name of the method; it just puts the numbers from
	    // the 4-vector into the array.)
	    positionStart.GetXYZT( fStartXYZT );
//std::cout<<"flag 3.2.2.2"<<std::endl;

	    positionEnd.GetXYZT( fEndXYZT );
	//std::cout<<"flag 3.2.2.3"<<std::endl;

	    momentumStart.GetXYZT( fStartPE );
	//std::cout<<"flag 3.2.2.4"<<std::endl;

	    momentumEnd.GetXYZT( fEndPE );
//std::cout<<"flag 3.2.3"<<std::endl;

	    // Use a polar-coordinate view of the 4-vectors to
	    // get the track length.
	    const double trackLength = ( positionEnd - positionStart ).Rho();
	    LOG_DEBUG("SSNetTest")
	      << "Track length: " << trackLength << " cm";
	    
	    // Fill a histogram of the track length.
	    fTrackLengthHist->Fill( trackLength ); 
//std::cout<<"flag 3.2.4"<<std::endl

            LOG_DEBUG("SSNetTest")
	      << "track ID=" << fSimTrackID << " (PDG ID: " << fSimPDG << ") "
	      << trackLength << " cm long, momentum " << momentumStart.P();
	   //std::cout<<"filling tree"<<std::endl;
	    fmytree->Fill();

	  } // if primary and PDG selected by user
      } // loop over all particles in the event. 
 //   std::cout<<"flag 4"<<std::endl;
    // std::cout<<"flag 6"<<std::endl;
    // We have a map of dE/dx vectors. Write each one of them to the
    // reconstruction n-tuple. 
    //for ( const auto& dEdxEntry : dEdxMap )
      //{
	// At this point, we've filled in all the reconstruction
	// n-tuple's variables. Write it.
	//fReconstructionNtuple->Fill();
      //}
 
    // First, read in the clusters.
    // Again, here we require clusters to be available by getting a ValidHandle.
    // Otherwise, the following code would not work.
      // const art::FindManyP<recob::PFParticle> findManyPFP(clusterHandle, event, fClusterProducerLabel);

    const art::FindManyP<recob::Hit> findManyHits(clusterHandle, event, fClusterProducerLabel);

    if ( findManyHits.isValid() )
      {
        for ( size_t cluster_index = 0; cluster_index != clusterHandle->size(); ++cluster_index )
	  {
            auto const& hits = findManyHits.at( cluster_index );
	    mf::LogInfo("SSNetTest")  
	      << "Cluster ID=" << clusterHandle->at( cluster_index ).ID()
	      << " has " << hits.size() << " hits";
	  }
      } // findManyHits valid
    else
      {
	mf::LogError("SSNetTest")  
	  << "findManyHits recob::Hit for recob::Cluster failed;"
	  << " cluster label='" << fClusterProducerLabel << "'";
      }
 } // SSNetTest::analyze()
 
  
  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see SSNetTest.fcl for more information.
  DEFINE_ART_MODULE(SSNetTest)

} // namespace example
} // namespace lar 

namespace {
 // time to define that function...
 double DetectorDiagonal(geo::GeometryCore const& geom) {
	const double length = geom.DetLength();
	const double width = 2. * geom.DetHalfWidth();
	const double height = 2. * geom.DetHalfHeight();
	return std::sqrt(cet::sum_of_squares(length, width, height));
	} // DetectorDiagonal()

//takes the 3D vertex position and 3D hit position
//returns true if the hit is within the radius of the vertex, false otherwise 
/* bool inROI(double* vertex_pos, double* hit_pos){
	return false;
	}*/	
} // local namespace
