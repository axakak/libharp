//c-rbf.h
#ifndef C_RBF_H
#define C_RBF_H

#include "traceData.h"

#include <cmath>
#include <vector>
#include <list>
#include <unordered_map>

using namespace std;

struct Complex
{
  Complex(): amplitude(0), phase(0){}

  double amplitude, phase;
};

/************************************************************
 * Spatio-temporal Objects
 ***********************************************************/

class SpatioTemporalNeuron
{
public:
  SpatioTemporalNeuron(double X, double Y, double Z, double Time)
    :weight(X, Y, Z, Time), error(0){}

  SpatioTemporalNeuron(const Event& w)
    :weight(w), error(0){}

  //accessor/mutator methods
  void setWeight(const Event& event);
  void setGain(double amp, double ph=0);
  void setError(double err);
  void setIndex(int i);
  const Event& getWeight() const;
  const Complex& getGain() const;
  double getError() const;

  //evaluation methods
  void computeGain(const Event& event);

  //training methods
  double computeDistance(const Event& event) const;
  void accumulateError(double err);
  void scaleError(double factor);
  void incramentEdgeAges();
  void adapt(const Event& event);
  void connect(SpatioTemporalNeuron* stNeuron);
  void disconnect(SpatioTemporalNeuron* stNeuron);
  void disconnectOld(int maxAge);
  void neighbourWithLargestError(const SpatioTemporalNeuron* stNeuron);
  bool noEdges() const;
  void assignEvent(int cGroup, const Event& event);

  void exportConnectionsYaml(YAML::Emitter& e);

private:
  Event weight;
  Complex gain;
  double error;
  int index;
  unordered_map<SpatioTemporalNeuron*, int> edges;
  unordered_multimap<int, const Event*> eventCluster;
};


class SpatioTemporalLayer
{
public:
  string exportNeuronsYamlString();
  void exportNeuronsYamlFile(const string& fileName);

  void evaluate(const Event& event);
  void train(TraceData& td);
  void initRandomNeurons(int count);
  void assignEventsToNeurons(const vector<TraceData>& traces);

  SpatioTemporalNeuron* findNeuronNearestToEvent(const Event& event);
  std::pair<SpatioTemporalNeuron*, SpatioTemporalNeuron*>
   findNeuronsNearestToEvent(const Event& event);

  /*
  SpatioTemporalNeuron& operator[](size_t pos);
  const SpatioTemporalNeuron& operator[](size_t pos) const;
  */

private:
  //list container was selected for constant insert and delete times,
  //which is only important during training
  list<SpatioTemporalNeuron> neurons;
};


/************************************************************
 * Class Objects
 ***********************************************************/

class ClassNeuron
{
public:
  void computeGain(const SpatioTemporalLayer& stl);
  const Complex& getGain() const;

private:
  double omega(double x, double w);

  //training methods
  void computeWeights(SpatioTemporalLayer& stl);

  /* Vector of N weights, one for each ST neuron */
  vector<Complex> weights;
  Complex gain;
  //int classGroup;
};


class ClassLayer
{
public:
  void train(SpatioTemporalLayer& stl, const vector<TraceData>& tdv);
  void evaluate(const SpatioTemporalLayer& stl);

  ClassNeuron& operator[](size_t pos);

private:
  vector<ClassNeuron> neurons;//TODO: change to map?
};


/************************************************************
 *  Network Objects
 ***********************************************************/

class CRBFNeuralNetwork
{
public:
  CRBFNeuralNetwork(): trace(){}

  CRBFNeuralNetwork(const CRBFNeuralNetwork& orig);

  void evaluateTrace();

  void train(const string& traceFileList);
  void exportYamlFile(const string& traceFile);
  void exportCsvFile(const string& tracesFile) const;
  void loadTraceFileList(const string& traceFileList);
  void exportTracesYamlFile(const string& tracesFile) const;
  void loadCRBFNeuralNetworkFile(const string& crbfFile);

private:
  double computeErrorIncrement(const Complex& gain);

  /* Input */
  TraceData trace;

  vector<TraceData> traces;

  /* Neural Network Layers */
  SpatioTemporalLayer stLayer;
  ClassLayer cLayer;

  /* Output */
  vector<double> cumulativeErrors;
};


/************************************************************
 *  Inline Methods
 ***********************************************************/

inline void SpatioTemporalNeuron::setGain(double amp, double ph)
{
  gain.amplitude = amp;
  gain.phase = ph;
}


inline void SpatioTemporalNeuron::setWeight(const Event& event)
{
  weight = event;
}


inline void SpatioTemporalNeuron::setError(double err)
{
  error = err;
}


inline void SpatioTemporalNeuron::setIndex(int i)
{
  index = i;
}


inline void SpatioTemporalNeuron::accumulateError(double err)
{
  error += err;
}


inline void SpatioTemporalNeuron::scaleError(double factor)
{
  error *= factor;
}


inline const Complex& SpatioTemporalNeuron::getGain() const
{
  return gain;
}


inline const Event& SpatioTemporalNeuron::getWeight() const
{
  return weight;
}


inline double SpatioTemporalNeuron::getError() const
{
  return error;
}

inline bool SpatioTemporalNeuron::noEdges() const
{
  return edges.empty();
}

/*
inline SpatioTemporalNeuron& SpatioTemporalLayer::operator[](size_t pos)
{
  return neurons[pos];
}


inline const SpatioTemporalNeuron& SpatioTemporalLayer::operator[](size_t pos) const
{
  return neurons[pos];
}
*/

inline const Complex& ClassNeuron::getGain() const
{
  return gain;
}


inline ClassNeuron& ClassLayer::operator[](size_t pos)
{
  return neurons[pos];
}


/************************************************************
 * YAML parser/emitter
 ***********************************************************/

YAML::Emitter& operator<< (YAML::Emitter& out, const SpatioTemporalNeuron& v);

namespace YAML
{
  template<>
  struct convert<SpatioTemporalNeuron>
  {
    static Node encode(const SpatioTemporalNeuron& rhs)
    {
      Node node;
      node.push_back(rhs.getWeight());
      return node;
    }

    static bool decode(const Node& node, SpatioTemporalNeuron& rhs)
    {
      rhs.setWeight(node.as<Event>());
      return true;
    }
  };
}

#endif //C_RBF_H
