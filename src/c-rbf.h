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
  SpatioTemporalNeuron();
  //TODO: init with weights

  void computeGain(const Event& event);
  Complex getGain() const;

private:
  Event weight;
  Complex gain;
};


class SpatioTemporalLayer
{
public:
  void evaluate(const Event& event);
  void train();

  SpatioTemporalNeuron& operator[](size_t pos);

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
  int kBMU;
};


class ClassLayer
{
public:
  void evaluate(const SpatioTemporalLayer& stl, vector<double>cErrors);

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
  CRBFNeuralNetwork(CRBFNeuralNetwork orig);
  
  void evaluateTrace();

  //TODO: write training function
  void train();

private:
  double computeErrorIncrement(const Complex& gain);

  /* Input */
  TraceData trace;

  /* Neural Network Layers */
  SpatioTemporalLayer stLayer;
  ClassLayer cLayer

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


inline Complex ClassNeuron::getGain() const
{
  return gain;
}


inline ClassNeuron& ClassLayer::operator[](size_t pos)
{
  return neurons[pos];
}

#endif //C_RBF_H
