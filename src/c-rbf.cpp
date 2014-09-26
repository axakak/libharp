#include "c-rbf.h"

#include <vector>
#include <thread>
#include <random>
#include <algorithm>

/************************************************************
 * Spatio-temporal Objects
 ***********************************************************/

void SpatioTemporalNeuron::computeGain(const Event& event)
{
  double qx, qy, qz;

  qx = event.x - weight.x;
  qy = event.y - weight.y;
  qz = event.z - weight.z;

  gain.amplitude = sqrt( ((qx*qx)+(qy*qy)+(qz*qz)) / 3);

  //TODO: J. Zahradnik paper shows output range < -pi;pi >, double-check
  gain.phase =  event.time - weight.time;
}


void SpatioTemporalNeuron::setGain(const double amp, const double ph)
{
  gain.amplitude = amp;
  gain.phase = ph;
}


void SpatioTemporalNeuron::setWeight(const Event& event)
{
  weight = event;
}


Event& SpatioTemporalNeuron::getWeight()
{
  return weight;
}


// calculate the gain of each neuron for a given event.
void SpatioTemporalLayer::evaluate(const Event& event)
{
  //for each neuron evaluate event gain
  for(auto &neuron : neurons)
    neuron.computeGain(event);
}


/*
void SpatioTemporalLayer::evaluateEventThreaded(const Event& event)
{
  vector<thread> threads;

  for(int i = 0; i < 5; ++i)
    threads.push_back(std::thread(evaluateEvent, event, ));

  for(auto& thread : threads)
    thread.join();
}*/

void SpatioTemporalLayer::train(const TraceData& td)
{
  randomizeNeurons();

  for(int time = 0; time < 100; time++)
  {
    //For each input vector (i.e. training data)
    for(int dataIndex = 0; dataIndex < td.size(); dataIndex++)
    {
      // estimate distances and sort in increasing order
      estimateDistances(td[dataIndex]);

      // perform adaptation step for neural weights, using neural-gas algorithm
      adaptWeights(td[dataIndex], time);
    }
  }
}


void SpatioTemporalLayer::randomizeNeurons()
{
  // assign initial values to the weights with |w| <= 1 & 〈w <= 2pi
  std::random_device rd;// if too slow use as seed to pseudo-random generator
  std::uniform_real_distribution<> sDist(0,1);// range not inclusive [0,1)
  std::uniform_real_distribution<> tDist(0,2*g_pi);// range not inclusive [0,2pi)

  neurons.clear();

  //TODO:LATER: change to dynamic neuron count
  for(int i = 0; i < 200; i++)
    neurons.emplace_back(sDist(rd),sDist(rd),sDist(rd),tDist(rd));
}


void SpatioTemporalLayer::estimateDistances(const Event& inputData)
{
  // estimate distance between input data and all neurons
  evaluate(inputData);

  for(auto &neuron : neurons)
    neuron.setGain(hypot(neuron.getGain().amplitude, neuron.getGain().phase/g_pi));

  // sort neurons in increasing distance from input data
  std::sort(neurons.begin(), neurons.end(),
    [](SpatioTemporalNeuron lhs, SpatioTemporalNeuron rhs)
    {return lhs.getGain().amplitude < rhs.getGain().amplitude;});
}


void SpatioTemporalLayer::adaptWeights(const Event& inputData, int time)
{
  double eps = epsilon(time);
  double lam = lambda(time);
  double factor;
  Event weight;

  // perform adaptation step for the weights
  for(int neuronIndex = 0; neuronIndex < neurons.size(); neuronIndex++)
  {
    factor = eps*exp(-neuronIndex/lam);
    weight = neurons[neuronIndex].getWeight();

    weight.x = weight.x + (factor * (inputData.x - weight.x));
    weight.y = weight.y + (factor * (inputData.y - weight.y));
    weight.z = weight.z + (factor * (inputData.z - weight.z));
    weight.time = weight.time + (factor * (inputData.time - weight.time));

    neurons[neuronIndex].setWeight(weight);
  }
}


/************************************************************
 * Class Objects
 ***********************************************************/

void ClassNeuron::computeGain(const SpatioTemporalLayer& stl)
{
  int kBMU = 0;
  double minMagnitude = 2,
         magnitude = 0,
         re = 0,
         im = 0;

  // find the nearest neuron (best matching unit) from the ST layer
  // and compute the gain for that neuron
  for(int k = 0; k < weights.size(); k++)
  {
    re = omega(stl[k].getGain().amplitude, weights[k].amplitude);
    im = omega(stl[k].getGain().phase, weights[k].phase/g_pi);
    magnitude = hypot(re,im);

    if(magnitude < minMagnitude)
    {
      minMagnitude = magnitude; 
      kBMU = k;
      gain.amplitude = re;
      gain.phase = im;
    }
  }

  gain.phase = copysign(g_pi*gain.phase, stl[kBMU].getGain().phase);
}


double ClassNeuron::omega(double x, double w)
{
  double f = x/log(1-abs(w));

  return 1 - exp(-f*f);
}


void ClassLayer::evaluate(const SpatioTemporalLayer& stl)
{
  // compute the gain for each class neuron
  for(auto& neuron : neurons)
    neuron.computeGain(stl);
}


/************************************************************
 *  Network Objects
 ***********************************************************/

void CRBFNeuralNetwork::evaluateTrace()
{
  if(trace.empty())
  {
    cerr << "error: trace is null, cannot evaluate" << endl;
    exit(1);
  }

  //clear cumulativeErrors
  cumulativeErrors.clear();

  //for each event in trace evaluate event
  for(int i = 0; i < trace.size(); i++)
  {
    //evaluate spatio-temporal layer for given event
    stLayer.evaluate(trace[i]);

    //evaluate class layer
    cLayer.evaluate(stLayer);

    //for each class add error to cumulativeErrors
    for(int k = 0; k < cumulativeErrors.size(); k++)
      cumulativeErrors[k] += computeErrorIncrement(cLayer[k].getGain());
  }
}


double CRBFNeuralNetwork::computeErrorIncrement(const Complex& gain)
{
  double Ys = 1, // spacial priority
         Yt = 1; // temporal priority

  return (1/(trace.size()*hypot(Ys,Yt))) * hypot(Ys*gain.amplitude, Yt*gain.phase/g_pi); 
}


void CRBFNeuralNetwork::train()
{
  //TODO: load all traces into trace for training

  //train Spatio Temporal Layer
  stLayer.train(trace);

  //XXX: CURRENT
  //train Class Layer
  //cLayer.train();
}
