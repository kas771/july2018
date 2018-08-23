////////////////////////////////////////////////////////////////////////
//
// EventFilter class
//
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "TMath.h"
/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

///filters for events, etc
namespace filter {

  class EventFilter : public art::EDFilter  {

  public:
    
    explicit EventFilter(fhicl::ParameterSet const& ); 
    virtual ~EventFilter();      
    
    bool filter(art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

    
    std::vector < unsigned int > SetOfBadEvents()   const { return fBadEvents;}
    std::vector < unsigned int > SetOfBadRuns()     const { return fBadRuns;  }
    
  private:
    
    std::vector < unsigned int >            fBadEvents; ///< list of bad events
    std::vector < unsigned int >            fBadRuns;   ///< list of bad runs

    std::vector < unsigned int >            fSelEvents; ///< list of selected events
    std::vector < unsigned int >            fSelRuns;   ///< list of selected runs
    std::vector < unsigned int >            fSelSubRuns;   ///< list of selected runs
    std::string fEventList;
    int         fSelection; //0: reject events based on input
                            //>0: accept events based on txt file
                            //<0: reject events based on txt file
                                    

  }; //class EventFilter
}

void filter::EventFilter::beginJob()
{
}

//-------------------------------------------------
void filter::EventFilter::endJob()
{
}

///////////////////////////////////////////////////////

void filter::EventFilter::reconfigure(fhicl::ParameterSet const& p)
{
  fBadEvents  = p.get < std::vector <unsigned int> >("BadEvents");	       
  fBadRuns    = p.get < std::vector <unsigned int> >("BadRuns");	
  
  fSelection = p.get< int >("Selection");
  fEventList = p.get< std::string >("EventList");
  fSelEvents.clear();
  fSelRuns.clear();
  std::ifstream in;
  in.open(fEventList.c_str());
  char line[1024];
  while(1){
    in.getline(line,1024);
    if (!in.good()) break;
    unsigned int n0, n1, n2;
    sscanf(line,"%u %u %u",&n0,&n1,&n2);
    fSelRuns.push_back(n0);
    fSelSubRuns.push_back(n1);
    fSelEvents.push_back(n2);
  }
  in.close();
       
}


filter::EventFilter::EventFilter(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
}

filter::EventFilter::~EventFilter()
{
}

bool filter::EventFilter::filter(art::Event &evt)
{   
  unsigned int evtNo = (unsigned int) evt.id().event();
  unsigned int runNo = (unsigned int) evt.run();
  unsigned int subrunNo = (unsigned int) evt.subRun();
  if (fSelection==0){
    std::vector <unsigned int> sobe = SetOfBadEvents();
    std::vector <unsigned int> sobr = SetOfBadRuns();
    if (sobe.size() != sobr.size()) {
      throw cet::exception("EventFilter.cxx: ") << " BadEvent and BadRun list must be same length. Line " <<__LINE__ << ", " << __FILE__ << "\n";
    }
    
    for (unsigned int ii=0; ii<sobe.size(); ++ii){
      if(sobe.at(ii)==evtNo && sobr.at(ii)==runNo) 
	{
	  mf::LogInfo("EventFilter: ") << "\t\n Skipping run/event " << runNo <<"/"<< evtNo << " by request.\n";
	  return false;
	}
    }
    return true;  
  }
  else{
    for (unsigned int ii = 0; ii<fSelRuns.size(); ii++){
      if (fSelRuns[ii] == runNo && fSelSubRuns[ii] == subrunNo && fSelEvents[ii] == evtNo){
	//std::cout<<"true"<<std::endl;
	if (fSelection>0){
	  return true;
	}
	else{
	  return false;
	}
      }
    }
    if (fSelection>0){
      return false;
    }
    else {
      return true;
    }
  }
}

namespace filter {

   DEFINE_ART_MODULE(EventFilter)

} //namespace filt
