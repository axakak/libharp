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
      :weight(X, Y, Z, Time){};

  void computeGain(const Event& event);
  void setGain(const double amp, const double ph=0);
  void setWeight(const Event& event);
  Complex getGain() const;
  Event& getWeight();

private:
  Event weight;
  Complex gain;
};


class SpatioTemporalLayer
{
public:
  void evaluate(const Event& event);
  void train(const TraceData& td);
  void randomizeNeurons();
  void estimateDistances(const Event& inputData);
  void adaptWeights(const Event& inputData, int time);

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
  Complex getGain() const;

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
  CRBFNeuralNetwork();
  CRBFNeuralNetwork(const CRBFNeuralNetwork& orig);
  
  void evaluateTrace();

  void train();//TODO: complete training function

private:
  double computeErrorIncrement(const Complex& gain);

  /* Input */
  TraceData trace;

  /* Neural Network Layers */
  SpatioTemporalLayer stLayer;
  ClassLayer cLayer;

  /* Output */
  vector<double> cumulativeErrors;
};


/************************************************************
 *  Inline Methods
 ***********************************************************/

inline Complex SpatioTemporalNeuron::getGain() const
{
  return gain;
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
  return 1/(time + 1);
}

inline double SpatioTemporalLayer::lambda(const double time) const
{
  return 1/(time + 1);
}

inline Complex ClassNeuron::getGain() const
{
  return gain;
}


inline ClassNeuron& ClassLayer::operator[](size_t pos)
{
  return neurons[pos];
}

#endif //C_RBF_H
