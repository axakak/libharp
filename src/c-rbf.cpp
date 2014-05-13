#include "c-rbf.h"

#include <vector>

double SpatioTemporalNeuron::computeGainAmplitude(Event *event)
{
  double qx, qy, qz;

  qx = event->x - weight.x;
  qy = event->y - weight.y;
  qz = event->z - weight.z;

  return sqrt( ((qx*qx)+(qy*qy)+(qz*qz)) / 3);
}


double SpatioTemporalNeuron::computeGainPhase(Event *event)
{
  //TODO: J. Zahradnik paper shows output as < -pi;pi >, double-check
  return event->time - weight.time;
}


void SpatioTemporalLayer::evaluateEvent(Event& event)
{
  //for each neuron evaluate event

}


int ClassNeuron::findBestMatchingUnitIndex(/* ST layer gain array */)
{
  /*TODO:consider saving BMU calcuations for later use*/

  int kBMU = 0;
  double minMagnitude = 2,
         magnitude = 0,
         re = 0,
         im = 0;

  for(int k = 0; k < N; k++)
  {
    re = omega(stGainAmp[k],cWeightAmp[k]);
    im = omega(stGainPhase[k],cWeightPhase[k]/pi);
    magnitude = hypot(re,im);

    if(magnitude < minMagnitude)
    {
      minMagnitude = magnitude;
      kBMU = k;
    }
  }

  return kBMU;
}


double ClassNeuron::omega(double x, double w)
{
  double f = x/log(1-abs(w));

  return 1 - exp(-f*f);
}


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
  for(int i = 0, i < trace.size(); i++)
    evaluateEvent(trace.getEvent(i));
}


void CRBFNeuralNetwork::evaluateEvent(Event event)
{
  //evaluate spatio-temporal layer
  stLayer->evaluateEvent(event);

  //evalute class layer, and add error to cumulativeErrors
  cLayer->evaluateSTLayer(stLayer, cumulativeErrors);
}
