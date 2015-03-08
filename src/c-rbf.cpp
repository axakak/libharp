#include "c-rbf.h"

#include <thread>
#include <random>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>

//constexpr int g_num_threads = 4;

/************************************************************
 * Spatio-temporal Objects
 ***********************************************************/

void SpatioTemporalNeuron::computeGain(const Event& event)
{
  Event q = event - weight;

  gain.amplitude = sqrt( ((q.x*q.x)+(q.y*q.y)+(q.z*q.z))/3 );
  gain.phase = q.time;

  // constrain to <-pi;pi>
  /*
  if(gain.phase > g_pi)
    gain.phase -= 2*g_pi;
  else if(gain.phase < -g_pi)
    gain.phase += 2*g_pi;
  */
}


void SpatioTemporalNeuron::computeDistance(const Event& event)
{
  //TODO: revaluate the correctness of this method, consolidate with compute gain
  Event q = event - weight;

  q.time /= g_pi;

  distance = (q.x*q.x) + (q.y*q.y) + (q.z*q.z) + (q.time*q.time);
}


void SpatioTemporalNeuron::incramentEdgeAges()
{
  for(auto &edge : edges)
  {
    ++edge.second;
    ++edge.first->edges[this];
  }
}


void SpatioTemporalNeuron::adapt(const Event& event)
{
  //move self towards event
  weight += (event - weight) * 0.1f;

  //move neighbors towards event
  for(auto &edge : edges)
    edge.first->weight += (event - edge.first->weight) * 0.001f;
}


void SpatioTemporalNeuron::connect(SpatioTemporalNeuron* stNeuron)
{
  //set age to zero, if edge does not exist it will be created
  if(stNeuron != this)
  {
    edges[stNeuron] = 0;
    stNeuron->edges[this] = 0;
  }
}


void SpatioTemporalNeuron::disconnect(SpatioTemporalNeuron* stNeuron)
{
  stNeuron->edges.erase(this);
  edges.erase(stNeuron);
}


void SpatioTemporalNeuron::disconnectOld(int maxAge)
{
  for(auto &edge : edges)
  {
    if(edge.second > maxAge)
    {
      edge.first->edges.erase(this);
      edges.erase(edge.first);
    }
  }
}


void SpatioTemporalNeuron::neighbourWithLargestError(const SpatioTemporalNeuron* stNeuron)
{
  stNeuron = edges.begin()->first;

  for(auto &edge : edges)
  {
    if(edge.first->error > stNeuron->error)
    {
      stNeuron = edge.first;
    }
  }
}


void SpatioTemporalNeuron::assignEvent(int group, const Event& event)
{
  eventCluster.insert( {{group, &event}} );
}


void SpatioTemporalNeuron::exportConnectionsYaml(YAML::Emitter& e)
{
  for(auto &edge : edges)
    if(!(index > edge.first->index))
      e  << YAML::Flow
         << YAML::BeginSeq << index << edge.first->index << YAML::EndSeq;
}


string SpatioTemporalLayer::exportNeuronsYamlString()
{
  YAML::Emitter out;

  out << YAML::BeginDoc
      << YAML::BeginMap
      << YAML::Key << "spatio-temporal-neuron-weights"
      << YAML::Value << YAML::BeginSeq
      << YAML::Flow
      << YAML::BeginSeq << "x" << "y" << "z" << "time" << YAML::EndSeq;

  int i = 0;
  for(auto &neuron : neurons)
  {
    out << neuron.getWeight();
    neuron.setIndex(i++);
  }

  out << YAML::EndSeq
      << YAML::Key << "spatio-temporal-neuron-edges"
      << YAML::Value << YAML::BeginSeq
      << YAML::Flow
      << YAML::BeginSeq << "v1" << "v2" << YAML::EndSeq;

  for(auto &neuron : neurons)
    neuron.exportConnectionsYaml(out);

  out << YAML::EndSeq
      << YAML::EndMap;

  return out.c_str();
}


void SpatioTemporalLayer::exportNeuronsYamlFile(const string& fileName)
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


void SpatioTemporalLayer::train(TraceData& td)
{
  cout << "\x1B[36m==>\x1B[0m " << "Training spatio-temporal"
       << " layer using growing neural-gas algorithm" << endl;

  unsigned int time = 0;
  int tdIndex = 0;
  SpatioTemporalNeuron *neuronS1 = 0, *neuronS2 = 0;

  ofstream file("spatioTemporalTrain.yaml");
  file << "%YAML 1.2";

  chrono::time_point<chrono::system_clock> start, end;
  start = chrono::system_clock::now();

  // randomly shuffle trace data
  td.shuffleEvents();

  // setup random number generation
  std::random_device rd;
  std::uniform_int_distribution<int> tdDist(0,td.size());

  //STEP 0: Generate 2 neurons at random positions
  generateRandomNeurons(2);

  do
  {
    //STEP 1: Select an input signal z from training data
    tdIndex = tdDist(rd);

    //STEP 2: Find the 2 nearest neurons, S1 and S2, to the input signal
    computeDistances(td[tdIndex]);

    neuronS1 = &(*(neurons.begin()));
    neuronS2 = &(*(neurons.rbegin()));
    for(auto &neuron : neurons)
    {
      //keep references to the 2 closest neurons
      if(neuron.getDistance() < neuronS1->getDistance())
      {
        neuronS2 = neuronS1;
        neuronS1 = &neuron;
      }
      else if(neuron.getDistance() < neuronS2->getDistance())
      {
        neuronS1 = &neuron;
      }
    }

    //step 3: Increment the age of all the edges emanating from S1
    neuronS1->incramentEdgeAges();

    //step 4: Add the squared distance between the input signal and S1 to S1 error
    neuronS1->accumulateError();

    //step 5: Move S1 and its connected neighbors towards z
    neuronS1->adapt(td[tdIndex]);

    //step 6: If S1 and S2 are connected by an edge, set the age of the
    //edge to 0. If such an age does not exist create it
    neuronS1->connect(neuronS2);

    //step 7: Remove edges with an age larger than a_max. If this results in
    //neurons having no emanating edges, remove them as well.
    neuronS1->disconnectOld(50);

    //TODO: try: neurons.remove_if([](SpatioTemporalNeuron n){ return n.noEdges(); });
    for(auto it = neurons.begin(); it != neurons.end(); ++it)
    {
      if(it->noEdges())
        neurons.erase(it);
    }

    //step 8: If the number of input signals selected is a multiple of lam,
    //insert new neuron
    if(!(time%60000))
    {
      cout << "\x1B[1K\x1B[40D" << "nc:" << neurons.size() <<  "  t:" << time;

      if(!(time%350000))
        file << endl << exportNeuronsYamlString();

      //Determine the neuron q with max accumulated error
      for(auto &neuron : neurons)
      {
        if(neuron.getError() > neuronS1->getError())
          neuronS1 = &neuron;
      }

      cout << "  e:" << setprecision(4) << neuronS1->getError();
      cout.flush();

      //find the neighbour f of q with the largest error
      neuronS1->neighbourWithLargestError(neuronS2);

      //insert new neuron r halfway between q and it's neighbour f with the largest error
      neurons.emplace_back((neuronS1->getWeight() + neuronS2->getWeight())*0.5f);

      //insert edges connecting r with q and f
      neurons.back().connect(neuronS1);
      neurons.back().connect(neuronS2);

      //remove the edge between q and f
      neuronS1->disconnect(neuronS2);

      //decrease the error of q and f by multiplying them with a constant _alpha
      neuronS1->scaleError(0.5f);
      neuronS2->scaleError(0.5f);

      //initialize the error of r with the new error of q
      neurons.back().setError(neuronS1->getError());
    }
    ++time;

    //step 9: Decrease all the error variables by multiplying them with a constant d
    for(auto &neuron : neurons)
      neuron.scaleError(0.995f);

  }while(time < 16020000);//TODO: select better stopping criterion

  end = chrono::system_clock::now();
  chrono::duration<double> elapsed_seconds = end-start;

  file << endl << exportNeuronsYamlString();
  file.close();

  cout << endl << "Spatio-temporal training completed in "
       << elapsed_seconds.count() << "s" << endl;
}


void SpatioTemporalLayer::computeDistances(const Event& event)
{
  // estimate distance between input data and all neurons
  for(auto &neuron : neurons)
    neuron.computeDistance(event);

  //TODO: thread compute distance when neuron count is large
/*
  thread t[g_num_threads];
  //Launch a group of threads
  for(int i = 0; i < g_num_threads; ++i)
    t[i] = thread(call_from_thread, i);

  //Join the threads with the main thread
  for(int i = 0; i < g_num_threads; ++i)
    t[i].join();
*/
}


void SpatioTemporalLayer::generateRandomNeurons(int count)
{
  // assign initial values to the weights with |w| <= 1 & 〈w <= 2pi
  std::random_device rd;// if too slow, use as seed to pseudo-random generator
  std::uniform_real_distribution<> sDist(0,1);// range not inclusive [0,1)
  std::uniform_real_distribution<> tDist(0,2*g_pi);// range not inclusive [0,2pi)

  neurons.clear();

  for(int i = 0; i < count; i++)
    neurons.emplace_back(sDist(rd),sDist(rd),sDist(rd),tDist(rd));
}


void SpatioTemporalLayer::assignEventsToNeurons(const vector<TraceData>& traces)
{
  SpatioTemporalNeuron *nNeuron = &(*(neurons.begin()));
  unsigned int eventCount = 0;

  //for each trace in training set(XXX:TODO: this can be threaded)
  for(auto &trace : traces)
  { //for each event in trace
    for(unsigned int tdIndex = 0; tdIndex < trace.size(); tdIndex++)
    {
      //compute the distance to all st-layer neurons
      computeDistances(trace[tdIndex]);

      //find the st-layer neuron nearest to event
      for(auto &neuron : neurons)
      {
        if(neuron.getDistance() < nNeuron->getDistance())
          nNeuron = &neuron;
      }

      //assign event to st-layer neuron cluster
      nNeuron->assignEvent(trace.getClassificationGroup(), trace[tdIndex]);

      eventCount++;
      cout << "\x1B[1K\x1B[40D" << "events assigned: " << eventCount;
    }
  }
  cout << endl;
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
    //TODO: replace stl [] with somthing that works with lists
    //re = omega(stl[k].getGain().amplitude, weights[k].amplitude);
    //im = omega(stl[k].getGain().phase, weights[k].phase/g_pi);
    magnitude = hypot(re,im);

    if(magnitude < minMagnitude)
    {
      minMagnitude = magnitude;
      kBMU = k;
      gain.amplitude = re;
      gain.phase = im;
    }
  }

  //gain.phase = copysign(g_pi*gain.phase, stl[kBMU].getGain().phase);
}


double ClassNeuron::omega(double x, double w)
{
  double f = x/log(1-abs(w));

  return 1 - exp(-f*f);
}


//XXX: Current top level method, not complete...
void ClassLayer::train(SpatioTemporalLayer& stl, const vector<TraceData>& tdv)
{
  cout << "\x1B[36m==>\x1B[0m " << "Training class layer" << endl;

  chrono::time_point<chrono::system_clock> start, end;
  start = chrono::system_clock::now();

  //TODO: Count the number of classes, if only one class weights are zero
  //countClasses(tdv);

  //seperate training data into N sets, one for each neuron in the st-layer
  stl.assignEventsToNeurons(tdv);

  //XXX: current

  //TODO: for each class calculate the weight of each st-layer neuron

  end = chrono::system_clock::now();
  chrono::duration<double> elapsed_seconds = end-start;

  //TODO: collect garbage

  cout << "Class layer training completed in "
       << elapsed_seconds.count() << "s" << endl;
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
  for(auto &t : traces)
    t.normalizeEvents();

  // export normalized trace data
  exportTracesYamlFile("normalizedTrainingData.yaml");

  // consolidate traces
  for(auto &t : traces)
    trace.insertEvents(t);

  // train Spatio Temporal Layer
  stLayer.train(trace);

  //TODO:train Class Layer
  cLayer.train(stLayer,traces);
}


void CRBFNeuralNetwork::exportYamlFile(const string& nnFile)
{
  cout << "Exporting neural network to " << nnFile << endl;

  ofstream file(nnFile);

  file << "%YAML 1.2" << endl
       << stLayer.exportNeuronsYamlString();// export spatio temporal neurons

  //TODO: export class neurons

  file.close();
}


void CRBFNeuralNetwork::exportTracesYamlFile(const string& fileName) const
{
  cout << "Exporting trace data to " << fileName << endl;

  ofstream file(fileName);
  file << "%YAML 1.2";

  for(auto &trace : traces)
    file << endl << trace.exportYamlString();

  file.close();
}


void CRBFNeuralNetwork::exportCsvFile(const string& tracesFile) const
{
  ofstream csvfile(tracesFile);

  csvfile << "x,y,z,time" << endl;

  for(auto &trace : traces)
    csvfile << trace.exportCsvString();

  csvfile.close();
}


void CRBFNeuralNetwork::loadTraceFileList(const string& traceFileList)
{
  cout << "\x1B[36m==>\x1B[0m "
       <<  "Loading traces from " << traceFileList << endl;

  ifstream fileList(traceFileList);
  string traceFile;
  int eventCount = 0;

  // for each file in list
  while(getline(fileList, traceFile))
  {
    cout << "Loading trace " << traceFile << endl;
    traces.emplace_back(traceFile);
    eventCount += traces.back().size();
  }

  cout << "Loading complete... " << traces.size() << " traces containing "
       << eventCount << " events" << endl;
  //TODO: check all trace data is from the same pattern
  //TODO: create classification group list

  fileList.close();
}


void CRBFNeuralNetwork::loadCRBFNeuralNetworkFile(const string& crbfFile)
{
  cout << "\x1B[36m==>\x1B[0m "
       <<  "Loading C-RBF neural network from " << crbfFile << endl;

  ifstream cFile(crbfFile);

  //TODO: method not complete...
}

/************************************************************
 * YAML parser/emitter
 ***********************************************************/

YAML::Emitter& operator<< (YAML::Emitter& out, const SpatioTemporalNeuron& v)
{
  out << v.getWeight();

  return out;
}
