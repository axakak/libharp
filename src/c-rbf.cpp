#include "c-rbf.h"

#include <vector>
#include <thread>
#include <random>
#include <algorithm>
#include <fstream>
#include <iostream> 
#include <iomanip>
#include <chrono>

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

  gain.phase =  event.time - weight.time;
 
  // constrain to <-pi;pi>
  /*
  if(gain.phase > g_pi)
    gain.phase -= 2*g_pi;
  else if(gain.phase < -g_pi) 
    gain.phase += 2*g_pi;
  */
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


void SpatioTemporalNeuron::addWeight(const Event& event)
{
  weight += event;
}


void SpatioTemporalNeuron::accumulateOffset(const Event& data, double factor)
{
    weight += ((data - weight) * factor);
}


const Event& SpatioTemporalNeuron::getWeight() const
{
  return weight;
}


YAML::Emitter& operator<< (YAML::Emitter& out, const SpatioTemporalNeuron& v)
{
  out << v.getWeight();

  return out;
}


string SpatioTemporalLayer::exportNeuronsYamlString() const
{
  YAML::Emitter out;

  out << YAML::BeginDoc
      << YAML::BeginMap
      << YAML::Key << "spatio-temporal-neurons"
      << YAML::Value << YAML::BeginSeq
      << YAML::Flow
      << YAML::BeginSeq << "x" << "y" << "z" << "time" << YAML::EndSeq;

  for(auto &neuron : neurons)
    out << neuron;

  out << YAML::EndSeq
      << YAML::EndMap;

  return out.c_str();
}


void SpatioTemporalLayer::exportNeuronsYamlFile(const string& fileName) const
{
  cout << "Exporting spatio-temporal neurons to " << fileName << endl;

  ofstream file(fileName);

  file << "%YAML 1.2" << endl << exportNeuronsYamlString();

  file.close(); 
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


void SpatioTemporalLayer::train(TraceData& td)
{
  cout << "\x1B[36m==>\x1B[0m "
       << "Training spatio-temporal layer using neural-gas algorithm" << endl; 
  
  ofstream file("spatioTemporalTrain.yaml");
  file << "%YAML 1.2";

  chrono::time_point<chrono::system_clock> start, end;
  start = chrono::system_clock::now();

  randomizeNeurons();

  cout << "Adapting spatio-temporal neurons" << endl
       << "  0.0%";

  for(int time = 0, maxTime = 20; time < maxTime; time++)
  {
    double eps = epsilon(time);
    double lam = lambda(time);
    
    file << endl << exportNeuronsYamlString();
    
    // randomly shuffle trace data
    td.shuffleEvents();

    for(int dataIndex = 0; dataIndex < td.size(); dataIndex++)
    {
      // estimate distance between input data and all neurons
      evaluate(td[dataIndex]);

      //convert distance to magnitude
      for(auto &neuron : neurons)
        neuron.setGain(hypot(neuron.getGain().amplitude, neuron.getGain().phase/g_pi));

      // sort neurons in increasing distance from input data
      // TODO: consider using partial_sort
      std::sort(neurons.begin(), neurons.end(),
        [](SpatioTemporalNeuron lhs, SpatioTemporalNeuron rhs)
        {return lhs.getGain().amplitude < rhs.getGain().amplitude;});

      // calculate adaptation step for the weights
      for(int k = 0; k < neurons.size(); k++)
        neurons[k].accumulateOffset(td[dataIndex], eps*exp(-k/lam));
    }

    cout << "\x1B[1K\x1b[10D" << setw(3) << (time/maxTime)*100.0f << "%";
    cout.flush();
    
  }
  
  end = chrono::system_clock::now();
  chrono::duration<double> elapsed_seconds = end-start;

  file << endl << exportNeuronsYamlString();
  file.close(); 

  cout << "\nSpatio-temporal training completed in "
       << elapsed_seconds.count() << "s" << endl;
}


void SpatioTemporalLayer::randomizeNeurons()
{
  cout << "Randomizing spatio-temporal neurons" << endl; 

  // assign initial values to the weights with |w| <= 1 & 〈w <= 2pi
  std::random_device rd;// if too slow, use as seed to pseudo-random generator
  std::uniform_real_distribution<> sDist(0,1);// range not inclusive [0,1)
  std::uniform_real_distribution<> tDist(0,2*g_pi);// range not inclusive [0,2pi)

  neurons.clear();

  //TODO:LATER: use dynamic neuron count (thesis topic)
  for(int i = 0; i < 250; i++)
    neurons.emplace_back(sDist(rd),sDist(rd),sDist(rd),tDist(rd));
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


void CRBFNeuralNetwork::train(const string& traceFileList)
{
  cout << "\x1B[35m==>\x1B[0m " << "Training neural network" << endl; 

  // load all traces into trace for network training
  loadTraceFileList(traceFileList);

  // normalize trace data
  for(auto &mTrace : mTraces)
    mTrace.normalizeEvents();
  
  // export normalized trace data
  exportMTraceYamlFile("normalizedTrainingData.yaml");

  // consolidate traces
  for(auto &mTrace : mTraces)
    trace.insertEvents(mTrace);

  // train Spatio Temporal Layer
  stLayer.train(trace);

  //TODO:train Class Layer
  //cLayer.train();
}


void CRBFNeuralNetwork::exportYamlFile(const string& nnFile) const
{
  cout << "Exporting neural network to " << nnFile << endl;

  ofstream file(nnFile);

  file << "%YAML 1.2" << endl
       << stLayer.exportNeuronsYamlString();// export spatio temporal neurons
  
  //TODO: export class neurons

  file.close();
}


void CRBFNeuralNetwork::exportMTraceYamlFile(const string& fileName) const
{
  cout << "Exporting training trace data to " << fileName << endl;

  ofstream file(fileName);
  file << "%YAML 1.2";

  for(auto &mTrace : mTraces)
    file << endl << mTrace.exportYamlString();

  file.close();
}


void CRBFNeuralNetwork::exportCsvFile(const string& mTracesFile) const
{
  ofstream csvfile(mTracesFile);

  csvfile << "x,y,z,time" << endl;
  
  for(auto &mTrace : mTraces)
    csvfile << mTrace.exportCsvString();

  csvfile.close();
}


void CRBFNeuralNetwork::loadTraceFileList(const string& traceFileList)
{
  cout << "\x1B[36m==>\x1B[0m "
       <<  "Loading traces from " << traceFileList << endl; 

  ifstream fileList(traceFileList);
  string traceFile;

  // for each file in list
  while(getline(fileList, traceFile))
  {
    cout << "Loading trace " << traceFile << endl; 
    mTraces.emplace_back(traceFile);
  }

  //TODO: check all trace data is from the same pattern

  fileList.close();
}
