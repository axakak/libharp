#include "c-rbf.h"

#include <thread>
#include <random>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <unordered_set>

//constexpr int g_num_threads = 4;


/*******************************************************************************
 * Spatio-temporal Objects
 ******************************************************************************/

/* SpatioTemporalNeuron *******************************************************/

Complex SpatioTemporalNeuron::computeGain(const Event& event) const
{
  Event q = event - weight;

  double a = sqrt( ((q.x*q.x)+(q.y*q.y)+(q.z*q.z))/3 );
  double p = q.time;

  // constrain to <-pi;pi>
  /*
  if(gain.phase > g_pi)
    gain.phase -= 2*g_pi;
  else if(gain.phase < -g_pi)
    gain.phase += 2*g_pi;
  */

  return Complex(a, p);
}


double SpatioTemporalNeuron::computeDistance(const Event& event) const
{
  //TODO: revaluate the correctness of this method, consolidate with compute gain
  double x = event.x-weight.x,
         y = event.y-weight.y,
         z = event.z-weight.z,
         time = (event.time-weight.time)/g_pi;

  return ((x*x) + (y*y) + (z*z) + (time*time));
}


void SpatioTemporalNeuron::incramentEdgeAges()
{
  for(auto &edge : edges)
  {
    ++(edge.second);
    ++((edge.first)->edges[this]);
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


bool SpatioTemporalNeuron::disconnectOld(int maxAge)
{
  bool triggered = false;

  for(auto &edge : edges)
  {
    if(edge.second > maxAge)
    {
      edge.first->edges.erase(this);
      edges.erase(edge.first);

      triggered = true;
    }
  }

  return triggered;
}


void SpatioTemporalNeuron::neighbourWithLargestError(const SpatioTemporalNeuron* stNeuron)
{
  stNeuron = edges.begin()->first;

  for(auto &edge : edges)
  {
    if((edge.first)->error > stNeuron->error)
    {
      stNeuron = edge.first;
    }
  }
}


void SpatioTemporalNeuron::exportConnectionsYaml(YAML::Emitter& e)
{
  for(auto &edge : edges)
    if(!(index > edge.first->index))
      e  << YAML::Flow
         << YAML::BeginSeq << index << edge.first->index << YAML::EndSeq;
}

/* SpatioTemporalLayer ********************************************************/

string SpatioTemporalLayer::exportYamlString()
{
  YAML::Emitter out;

  out << YAML::BeginMap
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

  file << "%YAML 1.2" << endl << "---" << exportYamlString();

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
  SpatioTemporalNeuron *neuronS1 = NULL, *neuronS2 = NULL;

  ofstream file("spatioTemporalTrain.yaml");
  file << "%YAML 1.2";

  chrono::time_point<chrono::system_clock> start, end;
  start = chrono::system_clock::now();

  // randomly shuffle trace data
  td.shuffleEvents();

  // setup random number generation
  std::random_device rd;
  std::mt19937 randGen(rd());
  std::uniform_int_distribution<int> tdDist(0,td.size());

  //STEP 0: Initilize 2 neurons at random positions
  initRandomNeurons(2);

  do
  {
    //STEP 1: Select an input signal z from training data
    tdIndex = tdDist(randGen);

    //STEP 2: Find the 2 nearest neurons, S1 and S2, to the input signal
    std::tie(neuronS1, neuronS2) = findNeuronsNearestToEvent(td[tdIndex]);

    //step 3: Increment the age of all the edges emanating from S1
    neuronS1->incramentEdgeAges();

    //step 4: Add the squared distance between the input signal and S1 to S1 error
    neuronS1->accumulateError(neuronS1->computeDistance(td[tdIndex]));

    //step 5: Move S1 and its connected neighbors towards z
    neuronS1->adapt(td[tdIndex]);

    //step 6: If S1 and S2 are connected by an edge, set the age of the
    //edge to 0. If such an age does not exist create it
    neuronS1->connect(neuronS2);

    //step 7: Remove edges with an age larger than a_max. If this results in
    //neurons having no emanating edges, remove them as well.
    //TODO: find a way to check only the neurons that where disconected
    //temp solution: only check all neurons if an edge was disconnected
    if(neuronS1->disconnectOld(50))//remove neurons with no emanating edges
      neurons.remove_if([](SpatioTemporalNeuron& n){ return n.noEdges(); });

    //step 8: If the number of input signals selected is a multiple of lam,
    //insert new neuron
    if(!(time%60000))
    {
      cout << "\x1B[1K\x1B[40D" << "nc:" << neurons.size() <<  "  t:" << time;

      if(!(time%350000))
        file << endl << "---" << endl << exportYamlString();

      //find the neuron q with max accumulated error
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

  file << endl << "---" << endl << exportYamlString();
  file.close();

  cout << endl << "Spatio-temporal training completed in "
       << elapsed_seconds.count() << "s" << endl;
}


void SpatioTemporalLayer::initRandomNeurons(int count)
{
  // assign initial values to the weights with |w| <= 1 & 〈w <= 2pi
  std::random_device rd;// if too slow, use as seed to pseudo-random generator
  std::uniform_real_distribution<> sDist(0,1);// range not inclusive [0,1)
  std::uniform_real_distribution<> tDist(0,2*g_pi);// range not inclusive [0,2pi)

  neurons.clear();

  for(int i = 0; i < count; i++)
    neurons.emplace_back(sDist(rd),sDist(rd),sDist(rd),tDist(rd));
}


SpatioTemporalNeuron*
 SpatioTemporalLayer::findNeuronNearestToEvent(const Event& event)
{
  SpatioTemporalNeuron *n1 = &(*(neurons.begin()));
  double distance1 = n1->computeDistance(event),
         newDistance = 0;

  for(auto &neuron : neurons)
  {
    newDistance = neuron.computeDistance(event);

    if(newDistance < distance1)
    {
      n1 = &neuron;
      distance1 = newDistance;
    }
  }

  return n1;
}


std::pair<SpatioTemporalNeuron*, SpatioTemporalNeuron*>
 SpatioTemporalLayer::findNeuronsNearestToEvent(const Event& event)
{
  SpatioTemporalNeuron *n1 = &(*(neurons.begin())),
                       *n2 = &(*(neurons.rbegin()));
  double distance1 = n1->computeDistance(event),
         distance2 = n2->computeDistance(event),
         newDistance = 0;

  for(auto &neuron : neurons)
  {
    newDistance = neuron.computeDistance(event);

    if(newDistance < distance1)
    {
      n2 = n1;
      n1 = &neuron;
      distance2 = distance1;
      distance1 = newDistance;
    }
    else if(newDistance < distance2 && n1 != &neuron)
    {
      n2 = &neuron;
      distance2 = newDistance;
    }
  }

  return make_pair(n1, n2);
}


/*******************************************************************************
 * Class Objects
 ******************************************************************************/

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


void ClassNeuron::computeWeight(const SpatioTemporalNeuron* stn, unordered_multimap<int, Complex> cgm)
{
  Complex maxGain(0,0);
  //int eventClassGroup;

  //weight is determined by using the event in the current st-neuron cluster
  //which produces the largest st-neuron gain.
  //Weight amplitude and phase are determined independently

  for(auto &eventGain : cgm)
  {
    //if gain is not from event in this class
    if(eventGain.first != classGroup)
    {
      if(eventGain.second.amplitude > maxGain.amplitude)
        maxGain.amplitude = eventGain.second.amplitude;

      //TODO: find max for gain phase
    }
    //else maxGain is 0,0
  }

  //compute weight amplitude from max gain. 1/sqrt(-ln(0.95)) = 4.4153964427018
  maxGain.amplitude = 1-exp(-4.4153964427018 * maxGain.amplitude);

  //assign weight to neuron class link
  weights[stn] = maxGain;
}


void ClassNeuron::exportWeightsYaml(YAML::Emitter& e)
{
  e << YAML::BeginMap
    << YAML::Key << "class-group" << YAML::Value << classGroup
    << YAML::Key << "weights"
    << YAML::Value << YAML::BeginSeq
    << YAML::Flow
    << YAML::BeginSeq << "amplitude" << "phase" << YAML::EndSeq;

    for(auto &weight : weights)
      e << YAML::Flow << YAML::BeginSeq
        << weight.second.amplitude << weight.second.phase << YAML::EndSeq;

  e << YAML::EndSeq
    << YAML::EndMap;
}

/* ClassLayer *****************************************************************/

//XXX: Current top level method, not complete...
void ClassLayer::train(SpatioTemporalLayer& stl, const vector<TraceData>& tdv)
{
  unordered_map<const SpatioTemporalNeuron* ,unordered_multimap<int, const Event*> > eventClusters;

  cout << "\x1B[36m==>\x1B[0m " << "Training class layer" << endl;

  chrono::time_point<chrono::system_clock> start, end;
  start = chrono::system_clock::now();

  //add class neuron for each class in the trainning data
  unordered_set<int> classSet;

  //scan all traces to create class set
  for(auto &trace : tdv)
    classSet.insert(trace.getClassificationGroup());

  //create neuron for each class in set
  for(auto &classGroup : classSet)
    neurons.emplace_back(classGroup);

  cout << "Classifications detected: " << classSet.size() << endl;

  cout << "Seperating events into neuron clusters" << endl;

  //XXX:TODO: thread, give each thread a rang of traces
  //seperate trainning data into N clusters, one for each st-layer neuron
  for(auto &trace : tdv)
  {
    SpatioTemporalNeuron *nNeuron = NULL;
    int group = trace.getClassificationGroup();

    //for each event in trace
    for(auto &event : trace)
    {
      //find the st-layer neuron nearest to event
      nNeuron = stl.findNeuronNearestToEvent(event);

      //assign event to st-layer neuron cluster
      eventClusters[nNeuron].emplace(group, &event);
    }
  }

  cout << "Computing class layer weights" << endl;

  //XXX:TODO: thread, give each thread a range of clusters
  //for each st-neuron cluster
  for(auto &cluster : eventClusters)
  {
    const SpatioTemporalNeuron* stNeuron = cluster.first;
    unordered_multimap<int, Complex> clusterGainMap;

    //compute the gain for each event in cluster
    for(auto &event : cluster.second)
      clusterGainMap.emplace(event.first, stNeuron->computeGain(*(event.second)));

    //for each class neuron determine weight for this st-neuron from gains
    for(auto &neuron : neurons)
      neuron.computeWeight(stNeuron, clusterGainMap);
  }

  end = chrono::system_clock::now();
  chrono::duration<double> elapsed_seconds = end-start;

  cout << "Class layer training completed in "
       << elapsed_seconds.count() << "s" << endl;

  //XXX: current
}


void ClassLayer::evaluate(const SpatioTemporalLayer& stl)
{
  // compute the gain for each class neuron
  for(auto& neuron : neurons)
    neuron.computeGain(stl);
}


string ClassLayer::exportYamlString()
{
  YAML::Emitter out;

  out << YAML::BeginMap
      << YAML::Key << "class-layer"
      << YAML::Value << YAML::BeginMap
      << YAML::Key << "class-neurons"
      << YAML::Value << YAML::BeginSeq;

      for(auto &neuron : neurons)
        neuron.exportWeightsYaml(out);

  out << YAML::EndSeq
      << YAML::EndMap
      << YAML::EndMap;

  return out.c_str();
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

  // consolidate traces//TODO: move into stLayer::train()?
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
       << "---" << endl
       << stLayer.exportYamlString() << endl// export st-neurons
       << cLayer.exportYamlString();// export class neurons

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
