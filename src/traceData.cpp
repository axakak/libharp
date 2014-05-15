#include "traceData.h"

#include <iostream>
#include <sstream>
#include <iomanip>


TraceData::TraceData(): totalTime(0)
{
}


void TraceData::loadYamlFile(char* traceFile)
{
  YAML::Node doc = YAML::LoadFile(traceFile);

  events.clear();

  //load metadata
  //TODO: implement data validation
  patientID = doc["patient-id"].as<std::string>();
  date = doc["date"].as<std::string>();
  location = doc["location"].as<std::string>();
  patternType = doc["pattern"]["type"].as<std::string>();
  patternLevel = doc["pattern"]["level"].as<int>();
  workspace = doc["workspace"].as<std::string>();
  coordinateSpace = doc["coordinate-space"].as<std::string>();
  totalTime = doc["total-time"].as<double>();

  //load events
  const YAML::Node& eventsNode = doc["events"];

  for(size_t i = 1; i < eventsNode.size(); i++)
    events.push_back(eventsNode[i].as<Event>());
}


void TraceData::exportYamlFile(char* traceFile) const
{
  YAML::Emitter out;

  out << YAML::BeginDoc
      << YAML::BeginMap
      << YAML::Key << "patient-id" << YAML::Value << patientID
      << YAML::Key << "date" << YAML::Value << date
      << YAML::Key << "location" << YAML::Value << location
      << YAML::Key << "pattern"
      << YAML::BeginMap
        << YAML::Key << "type" << YAML::Value << patternType
        << YAML::Key << "level" << YAML::Value << patternLevel
      << YAML::EndMap
      << YAML::Key << "workspace" << YAML::Value << workspace
      << YAML::Key << "coordinate-space" << YAML::Value << coordinateSpace
      << YAML::Key << "total-time" << YAML::Value << totalTime
      << YAML::Key << "events"
      << YAML::Value << YAML::BeginSeq
      << YAML::Flow
      << YAML::BeginSeq << "x" << "y" << "z" << "time" << YAML::EndSeq;

  for(auto& event : events)
    out << event;

  out << YAML::EndSeq
      << YAML::EndMap;


  ofstream file(traceFile);
  file << "%YAML 1.2\n" << out.c_str();
  file.close();
}


void TraceData::normalizeEvents()
{
  Event eventMax, eventMin, eventScale;
  vector<Event>::iterator eventVecItr;

  eventMax = eventMin = events.front();
  eventMax.time = events.back().time;

  //find min and max for each event spacial dimention
  for(eventVecItr = events.begin(); eventVecItr != events.end(); eventVecItr++)
  {
    if(eventVecItr->x > eventMax.x)
      eventMax.x = eventVecItr->x;
    else if(eventVecItr->x < eventMin.x) 
      eventMin.x = eventVecItr->x;

    if(eventVecItr->y > eventMax.y)
      eventMax.y = eventVecItr->y;
    else if(eventVecItr->y < eventMin.y) 
      eventMin.y = eventVecItr->y;

    if(eventVecItr->z > eventMax.z)
      eventMax.z = eventVecItr->z;
    else if(eventVecItr->z < eventMin.z) 
      eventMin.z = eventVecItr->z;
  }

  eventScale.x = 1 / (eventMax.x - eventMin.x);
  eventScale.y = 1 / (eventMax.y - eventMin.y);
  eventScale.z = 1 / (eventMax.z - eventMin.z);
  eventScale.time = (2*PI) / eventMax.time;

  for(eventVecItr = events.begin(); eventVecItr != events.end(); eventVecItr++)
  {
    eventVecItr->x = (eventVecItr->x - eventMin.x) * eventScale.x;
    eventVecItr->y = (eventVecItr->y - eventMin.y) * eventScale.y;
    eventVecItr->z = (eventVecItr->z - eventMin.z) * eventScale.z;
    eventVecItr->time = eventVecItr->time * eventScale.time;
  }
}


ostream& operator<< (ostream& stream, Event& e)
{
  stream << "x: " << fixed << setprecision(4) << setw(8) << e.x << " "
         << "y: " << setw(8) << e.y << " "
         << "z: " << setw(8) << e.z << " "
         << "time: " << setprecision(1) << setw(8) << e.time;

  return stream;
}


YAML::Emitter& operator<< (YAML::Emitter& out, const Event& v)
{
  //yaml-cpp does not support fixed floating-point formatting,
  //osrtingstream is only used for formatting. 
  //TODO: replace with faster formatting procedure. -ak

  ostringstream s;
  s.precision(4);
  s.setf(ios_base::fixed, ios_base::floatfield);

  out << YAML::Flow << YAML::BeginSeq;
  s << v.x;
  out << s.str();
  s.str("");
  s << v.y;
  out << s.str();
  s.str("");
  s << v.z;
  out << s.str();
  s.str("");
  s << setprecision(1) << v.time;
  out << s.str();
  out << YAML::EndSeq;

  return out;
}


