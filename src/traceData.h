// traceDate.cpp
//
//
#ifndef TRACE_DATA_H
#define TRACE_DATA_H

#include <fstream>
#include <yaml-cpp/yaml.h>


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

private:

  int patientID;
  double totalTime;

  Event *dataSet;
};


TraceData::TraceData(): patientID(0), totalTime(0), dataSet(NULL)
{
}


void TraceData::loadYamlFile(char* traceFile)
{
  YAML::Node doc = YAML::LoadFile(traceFile);

  patientID = doc["patient-id"].as<int>();
  totalTime = doc["total-time"].as<double>();
}

#endif //TRACE_DATA_H
