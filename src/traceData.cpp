#include "traceData.h"

TraceData::TraceData(): patientID(0), totalTime(0)
{
}


void TraceData::loadYamlFile(char* traceFile)
{
  YAML::Node doc = YAML::LoadFile(traceFile);

  //load metadata
  //TODO: check for data
  patientID = doc["patient-id"].as<int>();
  totalTime = doc["total-time"].as<double>();
  
  //TODO:load additional nodes

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


void TraceData::normalizeEvents()
{
  Event eMax, eMin, eScale;
  vector<Event>::iterator eit;

  eMax = eMin = events.front();
  eMax.time = events.back().time;

  //find min and max for each event spacial dimention
  for(eit = events.begin(); eit != events.end(); eit++)
  {
    if(eit->x > eMax.x)
      eMax.x = eit->x;
    else if(eit->x < eMin.x) 
      eMin.x = eit->x;

    if(eit->y > eMax.y)
      eMax.y = eit->y;
    else if(eit->y < eMin.y) 
      eMin.y = eit->y;

    if(eit->z > eMax.z)
      eMax.z = eit->z;
    else if(eit->z < eMin.z) 
      eMin.z = eit->z;
  }

  eScale.x = 1 / (eMax.x - eMin.x);
  eScale.y = 1 / (eMax.y - eMin.y);
  eScale.z = 1 / (eMax.z - eMin.z);
  eScale.time = (2*PI) / eMax.time;

  for(eit = events.begin(); eit != events.end(); eit++)
  {
    eit->x = (eit->x - eMin.x) * eScale.x;
    eit->y = (eit->y - eMin.y) * eScale.y;
    eit->z = (eit->z - eMin.z) * eScale.z;
    eit->time = eit->time * eScale.time;
  }
}


/******************************************************************************
* Non-member Functions
******************************************************************************/

ostream& operator << (ostream& stream, Event e)
{
  stream << "x: " << e.x << " "
         << "y: " << e.y << " "
         << "z: " << e.z << " "
         << "time: " << e.time;

  return stream;
}

