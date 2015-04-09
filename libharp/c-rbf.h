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
  Complex(double amp, double ph): amplitude(amp), phase(ph){}

  double amplitude, phase;
};


/*******************************************************************************
 * Spatio-temporal Objects
 ******************************************************************************/

class SpatioTemporalNeuron
{
public:
  SpatioTemporalNeuron(double X, double Y, double Z, double Time)
    :weight(X, Y, Z, Time), error(0), index(0){}

  SpatioTemporalNeuron(const Event& w)
    :weight(w), error(0), index(0){}

  //accessor/mutator methods
  void setWeight(const Event& event);
  void setGain(double amp, double ph=0);
  void setError(double err);
  void setIndex(int i);
  const Event& getWeight() const;
  double getError() const;
  int getIndex() const {return index;}

  //evaluation methods
  Complex computeGain(const Event& event) const;

  //training methods
  double computeDistance(const Event& event) const;
  void accumulateError(double err);
  void scaleError(double factor);
  void incramentEdgeAges();
  void adapt(const Event& event);
  void connect(SpatioTemporalNeuron* stNeuron);
  void disconnect(SpatioTemporalNeuron* stNeuron);
  bool disconnectOld(int aMax);
  void neighbourWithLargestError(const SpatioTemporalNeuron* stNeuron);
  bool noEdges() const;

  void exportConnectionsYaml(YAML::Emitter& e);

private:
  Event weight;
  double error;
  int index;
  unordered_map<SpatioTemporalNeuron*, int> edges;
};


class SpatioTemporalLayer
{
public:
  string exportYamlString();
  void exportNeuronsYamlFile(const string& filename);

  vector<Complex> evaluate(const Event& event) const;
  void train(vector<TraceData>& tdv, int ageMax, int insertionInterval, int reportCount);
  void initRandomNeurons(int count);
  void loadFile(const string& filename);

  SpatioTemporalNeuron* findNeuronWithLeastGain(const Event& event);
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


/*******************************************************************************
 * Class Objects
 ******************************************************************************/

class ClassNeuron
{
public:
  ClassNeuron(int cGroup):classGroup(cGroup){}


  //evaluation methods
  Complex computeGain(const vector<Complex>& stlGains) const;

  //training methods
  void computeWeight(const SpatioTemporalNeuron* stn, unordered_multimap<int, Complex> cgm);

  //expot methods
  void exportWeightsYaml(YAML::Emitter& e) const;

  //accessor methods
  int getClassGroup() const {return classGroup;}

private:
  double omega(double x, double w) const;

  int classGroup;

  /* Vector of N weights, one for each ST neuron */
  vector<Complex> weights;
};


class ClassLayer
{
public:
  void train(SpatioTemporalLayer& stl, const vector<TraceData>& tdv);
  vector<Complex> evaluate(const vector<Complex>& stLayerGains);
  void loadFile(const string& filename);
  string exportYamlString();

  ClassNeuron& operator[](size_t pos);
  size_t size() const {return neurons.size();}

private:
  vector<ClassNeuron> neurons;
};


/*******************************************************************************
 *  Network Objects
 ******************************************************************************/

class CRBFNeuralNetwork
{
public:
  CRBFNeuralNetwork(){};

  vector<double> evaluateTrace(const string& traceFile);

  void train(const string& traceFileList);
  void exportYamlFile(const string& traceFile);
  void loadCRBFNeuralNetworkFile(const string& crbfFile);

private:
  double computeErrorIncrement(const Complex& gain, int tdSize) const;

  /* Neural Network Layers */
  SpatioTemporalLayer stLayer;
  ClassLayer cLayer;
};


/*******************************************************************************
 *  Inline Methods
 ******************************************************************************/

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



inline ClassNeuron& ClassLayer::operator[](size_t pos)
{
  return neurons[pos];
}


/*******************************************************************************
 * YAML parser/emitter
 ******************************************************************************/

YAML::Emitter& operator<< (YAML::Emitter& out, const SpatioTemporalNeuron& v);

#endif //C_RBF_H
