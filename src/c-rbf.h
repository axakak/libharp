//c-rbf.h
#ifndef C_RBF_H
#define C_RBF_H

#include "traceData.h"

#include <cmath>

using namespace std;

struct Complex
{
  double amplitude, phase;
};


class SpatioTemporalNeuron
{
public:
  SpatioTemporalNeuron();
  //TODO: init with weights

  double computeGainAmplitude(Event *event);
  double computeGainPhase(Event *event);

private:
  Event weight;
  Complex gain;
};


class SpatioTemporalLayer
{
public:
  void evaluateEvent(Event& event);

private:
  SpatioTemporalNeuron *neurons;
  int size;
};


class ClassNeuron
{
public:
  int findBestMatchingUnit(SpatioTemporalNeuron *stLayer);

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

private:
  vector<ClassNeuron> neurons;
};


class CRBFNeuralNetwork
{
public:
  CRBFNeuralNetwork();

  CRBFNeuralNetwork(CRBFNeuralNetwork orig);
  
  ~CRBFNeuralNetwork();

  void evaluateTrace();

private:
  /* Input */
  TraceData trace;

  /* Neural Network Layers */
  SpatioTemporalLayer stLayer;
  ClassLayer cLayer

  /* Output */
  vector<double> cumulativeErrors;
};

#endif //C_RBF_H
