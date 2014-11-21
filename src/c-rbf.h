//c-rbf.h
#ifndef C_RBF_H
#define C_RBF_H

#include "traceData.h"

#include <cmath>

using namespace std;

class vector;

struct Complex
{
  double amplitude, phase;
};

/************************************************************
 * Spatio-temporal Objects
 ***********************************************************/

class SpatioTemporalNeuron
{
public:
  SpatioTemporalNeuron(double X, double Y, double Z, double Time)
      :weight(X, Y, Z, Time){}

  void computeGain(const Event& event);
  void setGain(const double amp, const double ph=0);
  void setWeight(const Event& event);
  void accumulateOffset(const Event& data, double factor);
  const Complex& getGain() const;
  const Event& getWeight() const;

private:
  Event weight;
  Complex gain;
};


class SpatioTemporalLayer
{
public:
  string exportNeuronsYamlString() const;
  void exportNeuronsYamlFile(const string& fileName) const;

  void evaluate(const Event& event);
  void train(TraceData& td);
  void randomizeNeurons();

  double epsilon(const double time) const;
  double lambda(const double time) const;

  SpatioTemporalNeuron& operator[](size_t pos);
  const SpatioTemporalNeuron& operator[](size_t pos) const;

private:
  vector<SpatioTemporalNeuron> neurons;
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

  /* Vector of N weights, one for each ST neuron */
  vector<Complex> weights;
  Complex gain;
  //int kBMU;
};


class ClassLayer
{
public:
  void evaluate(const SpatioTemporalLayer& stl);

  ClassNeuron& operator[](size_t pos);

private:
  vector<ClassNeuron> neurons;
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

  //TODO: complete training function
  void train(const string& traceFileList);
  void exportYamlFile(const string& traceFile) const;
  void exportCsvFile(const string& mTracesFile) const;
  void loadTraceFileList(const string& traceFileList);
  void exportMTraceYamlFile(const string& mTracesFile) const;

private:
  double computeErrorIncrement(const Complex& gain);

  /* Input */
  TraceData trace;

  vector<TraceData> mTraces;

  /* Neural Network Layers */
  SpatioTemporalLayer stLayer;
  ClassLayer cLayer;

  /* Output */
  vector<double> cumulativeErrors;
};


/************************************************************
 *  Inline Methods
 ***********************************************************/

inline void SpatioTemporalNeuron::setGain(const double amp, const double ph)
{
  gain.amplitude = amp;
  gain.phase = ph;
}


inline void SpatioTemporalNeuron::setWeight(const Event& event)
{
  weight = event;
}


inline void SpatioTemporalNeuron::accumulateOffset(const Event& data, double factor)
{
    weight += ((data - weight) * factor);
}


inline const Complex& SpatioTemporalNeuron::getGain() const
{
  return gain;
}


inline const Event& SpatioTemporalNeuron::getWeight() const
{
  return weight;
}


inline SpatioTemporalNeuron& SpatioTemporalLayer::operator[](size_t pos)
{
  return neurons[pos];
}


inline const SpatioTemporalNeuron& SpatioTemporalLayer::operator[](size_t pos) const
{
  return neurons[pos];
}


inline double SpatioTemporalLayer::epsilon(const double time) const
{
  return 1.0f/(time + 1.0f);
}


inline double SpatioTemporalLayer::lambda(const double time) const
{
  return 1.0f/(time + 1.0f);
}


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
