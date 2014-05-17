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
  void exportYamlFile(char* traceFile) const;

  // access methods
  string getPaitentID() const;
  string getDate() const; //TODO: select standared date format
  string getLocation() const;
  string getPatternType() const;
  int getPatternLevel() const;
  string getWorkspace() const;
  string getCoordinateSpace() const;
  double getTotalTime() const;

  // modifier methods
  void setPaitentID(string pID);
  void setDate(string);
  void setLocation(string);
  void setPatternType(string);
  void setPatternLevel(int);
  void setWorkspace(string);
  void setCoordinateSpace(string);
  void setTotalTime(double time);

  void clear();

  // std::vector methods
  Event& at(size_t pos);
  const Event& at(size_t pos) const;
  Event& operator[](size_t pos);
  void push_back(const Event&);
  void pop_back();
  bool empty() const;
  size_t size() const;

  void normalizeEvents();

private:
  string patientID;
  string date;
  string location;
  string patternType;
  int patternLevel;
  string workspace;
  string coordinateSpace;
  double totalTime;
  vector<Event> events;
};

// inline Methods

inline string TraceData::getPaitentID() const
{
  return patientID;
}


inline string TraceData::getDate() const
{
  return date;
}

inline string TraceData::getLocation() const
{
  return location;
}


inline string TraceData::getPatternType() const
{
  return patternType;
}


inline int TraceData::getPatternLevel() const
{
  return patternLevel;
}


inline string TraceData::getWorkspace() const
{
  return workspace;
}


inline string TraceData::getCoordinateSpace() const
{
  return coordinateSpace;
}


inline double TraceData::getTotalTime() const
{
  return totalTime;
}


inline void TraceData::setPaitentID(string str)
{
  patientID = str;
}

inline void TraceData::setDate(string str)
{
  date = str;
}


inline void TraceData::setLocation(string str)
{
  location = str;
}


inline void TraceData::setPatternType(string str)
{
  patternType = str;
}


inline void TraceData::setPatternLevel(int pLevel)
{
  patternLevel = pLevel;
}


inline void TraceData::setWorkspace(string str)
{
  workspace = str;
}


inline void TraceData::setCoordinateSpace(string str)
{
  coordinateSpace = str;
}


inline void TraceData::setTotalTime(double time)
{
 totalTime = time;
}


inline Event& TraceData::at(size_t pos)
{
  return events.at(pos);
}


inline const Event& TraceData::at(size_t pos) const
{
  return events.at(pos);
}


inline Event& TraceData::operator[](size_t pos)
{
  return events[pos];
}


inline void TraceData::push_back(const Event& value)
{
  events.push_back(value);
}


inline void TraceData::pop_back()
{
  events.pop_back();
}


inline bool TraceData::empty() const
{
  return events.empty();
}


inline size_t TraceData::size() const
{
  return events.size();
}


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
