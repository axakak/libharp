// traceDate.cpp
//
//
#ifndef TRACE_DATA_H
#define TRACE_DATA_H

//Standared lib
#include <vector>
#include <fstream>

#include <yaml-cpp/yaml.h>


using namespace std;

struct Event
{
  double x,y,z,time;
};

ostream& operator << (ostream& stream, Event e)
{
  stream << "x: " << e.x << " "
         << "y: " << e.y << " "
         << "z: " << e.z << " "
         << "time: " << e.time;

  return stream;
}


class TraceData 
{
public:

  //default constructor
  TraceData();

  void loadYamlFile(char* traceFile);

  //get functions
  int getPaitentID() {return patientID;}
  double getTotalTime() {return totalTime;}
  Event getEvent(int i) {return events[i];}

  size_t size() {return events.size();}
private:

  int patientID;
  double totalTime;
  vector<Event> events;
};


TraceData::TraceData(): patientID(0), totalTime(0)
{
}


void TraceData::loadYamlFile(char* traceFile)
{
  YAML::Node doc = YAML::LoadFile(traceFile);

  //TODO: check for data

  patientID = doc["patient-id"].as<int>();
  totalTime = doc["total-time"].as<double>();
  
  //load additional nodes

  
  //load events
  Event event;
  const YAML::Node& eventsNode = doc["events"];

  for(size_t i=1; i < eventsNode.size(); i++)
  {
    const YAML::Node& node = eventsNode[i];

    event.x = node[0].as<double>();
    event.y = node[1].as<double>();
    event.z = node[2].as<double>();
    event.time = node[3].as<double>();
 
    events.push_back(event);
  }


}

#endif //TRACE_DATA_H
