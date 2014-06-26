#include "c-rbf.h"

#include <vector>
#include <thread>

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

  //TODO: J. Zahradnik paper shows output as < -pi;pi >, double-check
  gain.phase =  event.time - weight.time;
}


// calculate the gain of each neuron for a given event.
void SpatioTemporalLayer::evaluate(const Event& event)
{
  int begin = 0,
      end = neurons.size();

  //for each neuron evaluate event gain
  for(int i = begin; i < end; i++)
    neurons[i].computeGain(event);
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

void SpatioTemporalLayer::train()
{
  //randomly init complex valued weights with |w| <= 1.

}

/************************************************************
 * Class Objects
 ***********************************************************/

int ClassNeuron::computeGain(const SpatioTemporalLayer& stl)
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
    im = omega(stl[k].getGain().phase, weights[k].phase/PI);
    magnitude = hypot(re,im);

    if(magnitude < minMagnitude)
    {
      minMagnitude = magnitude; 
      kBMU = k;
      gain.amplitude = re;
      gain.phase = im;
    }
  }

  gain.phase = copysign(PI*gain.phase, stl[kBMU].getGain().phase);
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
    neuron.computeGain();
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

    //evalute class layer
    cLayer.evaluate(stLayer);

    //for each class add error to cumulativeErrors
    for(int k = 0; k < cumulativeErrors.size(); k++)
      cumulativeErrors[k] =+ computeErrorIncrement(cLayer[k].getGain());
  }
}


double CRBFNeuralNetwork::computeErrorIncrement(const Complex& gain)
{
  double Ys = 1, // spacial priority
         Yt = 1; // temporal priority

  return (1/(trace.size()*hypot(Ys,Yt)) * hypot(Ys*gain.amplitude, Yt*gain.phase/PI); 
}


void CRBFNeuralNetwork::train()
{
  //prep traces for training

  //train Spatio Temporal Layer
  stLayer.train();

}
