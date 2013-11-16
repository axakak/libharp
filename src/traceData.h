// traceDate.cpp
//
//
#ifndef TRACE_DATA_H
#define TRACE_DATA_H

//Standared libs
#include <vector>
#include <fstream>

//3rd party libs
#include <yaml-cpp/yaml.h>


#define PI 3.1415926535897932

using namespace std;

struct Event
{
  double x,y,z,time;
};


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

  void normalizeEvents();

private:
  int patientID;
  double totalTime;
  vector<Event> events;
};

#endif //TRACE_DATA_H
