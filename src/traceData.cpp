#include "traceData.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <random>
#include <algorithm>

TraceData::TraceData(): totalTime(0)
{
}

TraceData::TraceData(const string& traceFile)
{
  loadYamlFile(traceFile);
}

void TraceData::loadYamlFile(const string& traceFile)
{
  YAML::Node doc = YAML::LoadFile(traceFile);

  fileName = traceFile;

  events.clear();

  //load metadata
  //TODO: implement data validation
  patientID = doc["patient-id"].as<std::string>();
  date = doc["date"].as<std::string>();
  location = doc["location"].as<std::string>();

  if(doc["pattern"])
    patternType = doc["pattern"]["type"].as<std::string>();
    patternLevel = doc["pattern"]["level"].as<int>();

  if(doc["classification"])
    classificationGroup = doc["classification"]["group"].as<int>();
  else
    classificationGroup = 0;

  workspace = doc["workspace"].as<std::string>();
  coordinateSpace = doc["coordinate-space"].as<std::string>();
  totalTime = doc["total-time"].as<double>();

  //load events
  //FIXME: some old files label "events" as "data" and will not load
  const YAML::Node& eventsNode = doc["events"];

  for(size_t i = 1; i < eventsNode.size(); i++)
    events.push_back(eventsNode[i].as<Event>());

  findMinMaxBounds();
}


void TraceData::exportYamlFile(const string& traceFile) const
{
  ofstream file(traceFile);
  file << "%YAML 1.2\n" << exportYamlString();
  file.close();
}


string TraceData::exportYamlString() const
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
      << YAML::Key << "classification"
      << YAML::BeginMap
        << YAML::Key << "group" << YAML::Value << classificationGroup
      << YAML::EndMap
      << YAML::Key << "workspace" << YAML::Value << workspace
      << YAML::Key << "workspace" << YAML::Value << workspace
      << YAML::Key << "coordinate-space" << YAML::Value << coordinateSpace
      << YAML::Key << "total-time" << YAML::Value << totalTime
      << YAML::Key << "events"
      << YAML::Value << YAML::BeginSeq
      << YAML::Flow
      << YAML::BeginSeq << "x" << "y" << "z" << "time" << YAML::EndSeq;

  for(auto &event : events)
    out << event;

  out << YAML::EndSeq
      << YAML::EndMap;

  return out.c_str();
}


void TraceData::exportCsvFile(const string& traceFile) const
{
  ofstream file(traceFile);

  file.precision(4);
  file.setf(ios_base::fixed, ios_base::floatfield);

  file << "x,y,z,time" << endl;

  for(auto &event : events)
  {
    file << event.x << "," << event.y << ","
         << event.z << "," << event.time << endl;
  }

  file.close();
}


string TraceData::exportCsvString() const
{
  stringstream csvStringStream;

  csvStringStream.precision(4);
  csvStringStream.setf(ios_base::fixed, ios_base::floatfield);

  for(auto &event : events)
  {
    csvStringStream << event.x << "," << event.y << ","
                    << event.z << "," << event.time << endl;
  }

  return csvStringStream.str();
}


void TraceData::findMinMaxBounds()
{
  maxBound = minBound = events.front();
  maxBound.time = events.back().time;

  //find min and max for each event spacial dimension
  for(auto &event : events)
  {
    if(event.x > maxBound.x)
      maxBound.x = event.x;
    else if(event.x < minBound.x)
      minBound.x = event.x;

    if(event.y > maxBound.y)
      maxBound.y = event.y;
    else if(event.y < minBound.y)
      minBound.y = event.y;

    if(event.z > maxBound.z)
      maxBound.z = event.z;
    else if(event.z < minBound.z)
      minBound.z = event.z;
  }
}


void TraceData::normalizeEvents()
{
  Event eventScale;

  //TODO: find better way to normalize z, values are very scattered -ak
  eventScale.x = 1 / (maxBound.x - minBound.x);
  eventScale.y = 1 / (maxBound.y - minBound.y);
  eventScale.z = 1 / (maxBound.z - minBound.z);
  eventScale.time = (2*g_pi) / maxBound.time;

  for(auto &event : events)
  {
    event.x = (event.x - minBound.x) * eventScale.x;
    event.y = (event.y - minBound.y) * eventScale.y;
    event.z = (event.z - minBound.z) * eventScale.z;
    event.time = (event.time - minBound.time) * eventScale.time;
  }
}


void TraceData::shuffleEvents()
{
  random_device rd;
  mt19937 g(rd());
  shuffle(events.begin(), events.end(), g);
}


void TraceData::insertEvents(TraceData& td)
{
  events.insert(events.begin(), td.events.begin(), td.events.end());
}


YAML::Emitter& operator<< (YAML::Emitter& out, const Event& v)
{
  //yaml-cpp doesn't support fixed floating-point output(i.e no trailing zeros),
  //so output is jagged //TODO: implement fixed point formatting when avaiable.

  out << YAML::Flow
      << YAML::BeginSeq
      << YAML::Precision(8) << v.x
      << YAML::Precision(8) << v.y
      << YAML::Precision(8) << v.z
      << YAML::Precision(8) << v.time
      << YAML::EndSeq;

  return out;
}


ostream& operator<< (ostream& stream, Event& e)
{
  stream << "x: " << fixed << setprecision(4) << setw(8) << e.x << " "
         << "y: " << setw(8) << e.y << " "
         << "z: " << setw(8) << e.z << " "
         << "time: " << setw(8) << e.time;

  return stream;
}


void TraceData::loadTraceDataList(const string& tdFileList, vector<TraceData>& tdv)
{
  cout << "\x1B[36m==>\x1B[0m "
       <<  "Loading traces from " << tdFileList << endl;

  ifstream fileList(tdFileList);
  string traceFile;
  int eventCount = 0;

  // for each file in list
  while(getline(fileList, traceFile))
  {
    cout << "Loading trace " << traceFile << endl;
    tdv.emplace_back(traceFile);
    eventCount += tdv.back().size();
  }

  cout << "Loading complete... " << tdv.size() << " traces containing "
       << eventCount << " total events" << endl;

  fileList.close();
}

void TraceData::exportTracesYamlFile(const string& fileName, vector<TraceData>& tdv)
{
  cout << "Exporting " << tdv.size() << " traces to " << fileName;

  ofstream file(fileName);
  file << "%YAML 1.2";

  for(auto &trace : tdv)
    file << endl << trace.exportYamlString();

  file.close();

  cout << "\x1B[32m\tComplete\x1B[0m" << endl;
}
