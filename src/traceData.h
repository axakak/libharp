// traceDate.h
#ifndef TRACE_DATA_H
#define TRACE_DATA_H

//Standard libs
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
  void exportYamlFile(char* traceFile);


  //mutator functions
  int  getPaitentID();
  void setPaitentID(int pID);

  double getTotalTime();
  void setTotalTime(double time);

  Event& getEvent(int i);

  size_t size();

  void normalizeEvents();

private:
  int patientID;
  string date;
  string location;
  string patternType;
  int patternLevel;
  string workspace;
  string coordinateSpace;

  double totalTime;
  vector<Event> events;
};

ostream& operator<< (ostream& stream, Event& e);

YAML::Emitter& operator << (YAML::Emitter& out, const Event& v);


namespace YAML
{
  template<>
  struct convert<Event>
  {
    static Node encode(const Event& rhs)
    {
      Node node;
      node.push_back(rhs.x);
      node.push_back(rhs.y);
      node.push_back(rhs.z);
      node.push_back(rhs.time);
      return node;
    }

    static bool decode(const Node& node, Event& rhs)
    {
      if(!node.IsSequence() || node.size() != 4)
        return false;

      rhs.x = node[0].as<double>();
      rhs.y = node[1].as<double>();
      rhs.z = node[2].as<double>();
      rhs.time = node[3].as<double>();
      return true;
    }
  };
}


#endif //TRACE_DATA_H
