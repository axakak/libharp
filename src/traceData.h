// traceDate.h
#ifndef TRACE_DATA_H
#define TRACE_DATA_H

//Standard libs
#include <vector>
#include <fstream>
#include <string>

//3rd party libs
#include <yaml-cpp/yaml.h>

#if _MSC_VER >= 1700 && _MSC_VER < 1800 //for Visual Studio 2012
  const double g_pi = 3.14159265358979323846;
#else
  constexpr double g_pi = 3.14159265358979323846;
#endif

using namespace std;


struct Event
{
  Event(): x(0), y(0), z(0), time(0){}
  Event(double X, double Y, double Z, double Time): x(X), y(Y), z(Z), time(Time){}

  Event operator-(const Event& rhs) const
  {
    return Event(x-rhs.x, y-rhs.y, z-rhs.z, time-rhs.time);
  }

  Event operator+(const Event& rhs) const
  {
    return Event(x+rhs.x, y+rhs.y, z+rhs.z, time+rhs.time);
  }

  Event operator*(const double rhs) const
  {
    return Event(x*rhs,y*rhs,z*rhs,time*rhs);
  }

  Event& operator+=(const Event& rhs)
  {
    x = x + rhs.x;
    y = y + rhs.y;
    z = z + rhs.z;
    time = time + rhs.time;

    return *this;
  }

  double x,y,z,time;
};


class TraceData
{
public:
  //default constructor
  TraceData();
  TraceData(const string& traceFile);

  // I/O functions
  void loadYamlFile(const string& traceFile);
  void exportYamlFile(const string& traceFile) const;
  string exportYamlString() const;
  void exportCsvFile(const string& traceFile) const;
  string exportCsvString() const;

  // access methods
  string getPaitentID() const;
  string getDate() const; //TODO: select standared date format
  string getLocation() const;
  string getPatternType() const;
  int getPatternLevel() const;
  int getClassificationGroup() const;
  string getWorkspace() const;
  string getCoordinateSpace() const;
  double getTotalTime() const;
  Event getMaxBound() const;
  Event getMinBound() const;

  // modifier methods
  void setPaitentID(string pID);
  void setDate(string);
  void setLocation(string);
  void setPatternType(string);
  void setPatternLevel(int);
  void setClassificationGroup(int);
  void setWorkspace(string);
  void setCoordinateSpace(string);
  void setTotalTime(double time);

  void clear();

  // std::vector methods
  Event& at(size_t pos);
  const Event& at(size_t pos) const;
  Event& operator[](size_t pos);
  const Event& operator[](size_t pos) const;
  void push_back(const Event&);
  void pop_back();
  bool empty() const;
  size_t size() const;

  void normalizeEvents();
  void shuffleEvents();
  void insertEvents(TraceData& td);

private:
  void findMinMaxBounds();

  string fileName;

  string patientID;
  string date;
  string location;
  string patternType;
  int patternLevel;
  int classificationGroup;
  string workspace;
  string coordinateSpace;
  double totalTime;

  vector<Event> events;
  Event maxBound;
  Event minBound;
};

/************************************************************
 *  Inline Methods
 ***********************************************************/

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


inline int TraceData::getClassificationGroup() const
{
  return classificationGroup;
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


inline Event TraceData::getMaxBound() const
{
  return maxBound;
}


inline Event TraceData::getMinBound() const
{
  return minBound;
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


inline void TraceData::setClassificationGroup(int classGroup)
{
  classificationGroup = classGroup;
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


inline const Event& TraceData::operator[](size_t pos) const
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


/************************************************************
 * YAML parser/emitter
 ***********************************************************/

YAML::Emitter& operator<< (YAML::Emitter& out, const Event& v);


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
