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
  double x = event.x-weight.x,
         y = event.y-weight.y,
         z = event.z-weight.z;

  double amp = sqrt( ( (x*x)+(y*y)+(z*z) )/3 );
  double ph = event.time-weight.time;

  // constrain phase to <-pi;pi>
  if(ph > g_pi)
    ph -= 2*g_pi;
  else if(ph < -g_pi)
    ph += 2*g_pi;

  return Complex(amp, ph);
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


bool SpatioTemporalNeuron::disconnectOld(int aMax)
{
  bool triggered = false;

  for(auto &edge : edges)
  {
    if(edge.second > aMax)
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
      << YAML::BeginSeq << "neuron-index" << "neuron-index" << YAML::EndSeq;

  for(auto &neuron : neurons)
    neuron.exportConnectionsYaml(out);

  out << YAML::EndSeq
      << YAML::EndMap;

  return out.c_str();
}


void SpatioTemporalLayer::loadFile(const string& filename)
{
  YAML::Node doc = YAML::LoadFile(filename);
  const YAML::Node& stnwNode = doc["spatio-temporal-neuron-weights"];

  for(size_t i = 1; i < stnwNode.size(); i++)
    neurons.emplace_back(stnwNode[i].as<Event>());
}


void SpatioTemporalLayer::exportNeuronsYamlFile(const string& filename)
{
  cout << "Exporting spatio-temporal neurons to " << filename << endl;

  ofstream file(filename);

  file << "%YAML 1.2" << endl << "---" << exportYamlString();
  file.close();
}


vector<Complex> SpatioTemporalLayer::evaluate(const Event& event) const
{
  vector<Complex> gains;

  //for each neuron evaluate event gain
  for(auto &neuron : neurons)
    gains.emplace_back(neuron.computeGain(event));

  return gains;
}


void SpatioTemporalLayer::train(vector<TraceData>& tdv, int ageMax, int insertionInterval, int reportCount)
{
  cout << "\x1B[36m==>\x1B[0m " << "Training spatio-temporal"
       << " layer using growing neural-gas algorithm" << endl;

  ofstream file("spatioTemporalTrain.yaml");
  file << "%YAML 1.2";

  chrono::time_point<chrono::system_clock> start, end;
  start = chrono::system_clock::now();

  TraceData td;

  // consolidate traces
  for(auto &t : tdv)
    td.insertEvents(t);

  // randomly shuffle trace data
  td.shuffleEvents();

  // setup random number generation
  std::random_device rd;
  std::mt19937 randGen(rd());
  std::uniform_int_distribution<int> tdDistrobution(0,td.size());

  unsigned int time = 0, maxTime, reportInterval;
  int tdIndex = 0;
  SpatioTemporalNeuron *neuronS1 = NULL, *neuronS2 = NULL;

  // align maxTime and reportInterval with insertionInterval
  maxTime = int(100*td.size()/insertionInterval);
  reportInterval = int(maxTime/(reportCount-2)) * insertionInterval;
  maxTime *= insertionInterval;

  //steps correspond to B. Fritzke, A Growing Neural Gas Network... (1995)
  //STEP 0: Initialize 2 neurons at random positions
  initRandomNeurons(2);

  do
  {
    //STEP 1: Select an input signal z from training data
    tdIndex = tdDistrobution(randGen);

    //STEP 2: Find the 2 nearest neurons, S1 and S2, to the input signal
    std::tie(neuronS1, neuronS2) = findNeuronsNearestToEvent(td[tdIndex]);

    //STEP 3: Increment the age of all the edges emanating from S1
    neuronS1->incramentEdgeAges();

    //STEP 4: Add squared distance between the input signal and S1, to S1 error
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
    if(neuronS1->disconnectOld(ageMax))//remove neurons with no emanating edges
      neurons.remove_if([](SpatioTemporalNeuron& n){ return n.noEdges(); });

    //step 8: If the number of input signals selected is a multiple of lam,
    //insert new neuron
    if(!(time%insertionInterval))
    {
      if(!(time%reportInterval))
        file << endl << "---" << endl << exportYamlString();

      //find the neuron q with max accumulated error
      for(auto &neuron : neurons)
      {
        if(neuron.getError() > neuronS1->getError())
          neuronS1 = &neuron;
      }

      //find adjacent neuron f of q with the largest error
      neuronS1->neighbourWithLargestError(neuronS2);

      //insert new neuron r halfway between q and it's neighbour f with the largest error
      neurons.emplace_back((neuronS1->getWeight() + neuronS2->getWeight())*0.5f);

      //insert edges connecting r with q and f
      neurons.back().connect(neuronS1);
      neurons.back().connect(neuronS2);

      //remove the edge between q and f
      neuronS1->disconnect(neuronS2);

      //decrease the error of q and f by multiplying them with a constant _alpha
      neuronS1->scaleError(0.5f);//TODO: make _alpha an argument
      neuronS2->scaleError(0.5f);

      //initialize the error of r with the new error of q
      neurons.back().setError(neuronS1->getError());

      cout << "\x1B[1K\x1B[40D" << "nc:" << neurons.size() << "  t:" << time
           << "  e:" << setprecision(4) << neuronS1->getError();
      cout.flush();
    }
    ++time;

    //step 9: Decrease all the error variables by multiplying them with a constant d
    for(auto &neuron : neurons)
      neuron.scaleError(0.995f);//TODO: make d and argument

  }while(time < maxTime);//TODO: select better stopping criterion

  end = chrono::system_clock::now();
  chrono::duration<double> elapsed_seconds = end-start;

  file << endl
       << "---" << endl
       << exportYamlString();

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
 SpatioTemporalLayer::findNeuronWithLeastGain(const Event& event)
 {
   SpatioTemporalNeuron *n1 = &neurons.front();
   Complex gain1 = n1->computeGain(event),
           newGain;
   double gainMag1 = hypot(gain1.amplitude, gain1.phase/g_pi),
          newGainMag;


   for(auto &neuron : neurons)
   {
     newGain = neuron.computeGain(event);
     newGainMag = hypot(newGain.amplitude, newGain.phase/g_pi);

     if(newGainMag < gainMag1)
     {
       n1 = &neuron;
       gainMag1 = newGainMag;
     }
   }

   return n1;
 }


SpatioTemporalNeuron*
 SpatioTemporalLayer::findNeuronNearestToEvent(const Event& event)
{
  SpatioTemporalNeuron *n1 = &neurons.front();
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
  SpatioTemporalNeuron *n1 = &neurons.front(),
                       *n2 = &neurons.back();
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

Complex ClassNeuron::computeGain(const vector<Complex>& stlGains) const
{
  size_t k_BMU;
  double minMagnitude = 2,
         magnitude = 0,
         re = 0,
         im = 0;
  Complex gain;

  // find the nearest neuron (best matching unit) from the ST layer
  // and compute the gain for that neuron
  for(size_t k = 0; k < stlGains.size(); k++)
  {
    //re = omega(stlGains[k].amplitude, weights[k].amplitude);
    //im = omega(stlGains[k].phase, weights[k].phase/g_pi);
    magnitude = hypot(re,im);

    if(magnitude < minMagnitude)
    {
      minMagnitude = magnitude;
      k_BMU = k;
      gain.amplitude = re;
      gain.phase = im;
    }
  }

  //gain.phase = copysign(g_pi*gain.phase, stlGains[k_BMU].phase);

  return gain;
}


double ClassNeuron::omega(double x, double w) const
{
  double f = x/log(1-abs(w));

  return 1 - exp(-f*f);
}


void ClassNeuron::computeWeight(const SpatioTemporalNeuron* stn, unordered_multimap<int, Complex> gainMap)
{
  Complex maxGain(0,0);
  //weight is determined by using the event in the current st-neuron cluster
  //which produces the largest st-neuron gain.
  //Weight amplitude and phase are determined independently

  for(auto &eventGain : gainMap)
  {
    //if gain is not from event in this class
    if(eventGain.first != classGroup)
    {
      //find max for gain amplitude in cluster
      if(eventGain.second.amplitude > maxGain.amplitude)
        maxGain.amplitude = eventGain.second.amplitude;

      //find max for gain phase in cluster
      if(abs(eventGain.second.phase) > abs(maxGain.phase))
        maxGain.phase = eventGain.second.phase;
    }//else maxGain is 0,0
  }

  //compute weight amplitude from max gain. 1/sqrt(-ln(0.95)) = 4.4153964427018
  if(maxGain.amplitude)
    maxGain.amplitude = 1-exp(-4.4153964427018 * maxGain.amplitude);

  //compute weight phase from max gain. (1/pi)/sqrt(-ln(0.95)) = 1.40546433913273
  if(maxGain.phase)
    maxGain.phase = copysign((1-exp(-1.40546433913273 * abs(maxGain.phase)))*g_pi, maxGain.phase);

  //assign class neuron weight linkend to stn
  weights[stn] = maxGain;
}


void ClassNeuron::exportWeightsYaml(YAML::Emitter& e) const
{
  e << YAML::BeginMap
    << YAML::Key << "class-group" << YAML::Value << classGroup
    << YAML::Key << "weights"
    << YAML::Value << YAML::BeginSeq
    << YAML::Flow
    << YAML::BeginSeq
    << "st-neuron-index" << "amplitude" << "phase"
    << YAML::EndSeq;

    for(auto &weight : weights)
    {
      //only export nonzero weights
      if(weight.second.amplitude || weight.second.phase)
      {
        e << YAML::Flow
          << YAML::BeginSeq
          << (weight.first)->getIndex()
          << weight.second.amplitude << weight.second.phase
          << YAML::EndSeq;
      }
    }

  e << YAML::EndSeq
    << YAML::EndMap;
}

/* ClassLayer *****************************************************************/

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

  cout << "Classifications detected: " << classSet.size() << endl
       << "Seperating events into neuron clusters" << endl;

  //XXX:TODO: thread, give each thread a rang of traces
  //seperate trainning data into N clusters, one for each st-layer neuron
  for(auto &trace : tdv)
  {
    SpatioTemporalNeuron *nNeuron = NULL;
    int group = trace.getClassificationGroup();

    //for each event in trace
    for(auto &event : trace)
    {
      //find the st-layer neuron with the least amout of gain from event
      nNeuron = stl.findNeuronWithLeastGain(event);

      //assign event to st-layer neuron cluster
      eventClusters[nNeuron].emplace(group, &event);
    }
  }
  //FIXME: eventClusters should always contain all st-neurons?
  // delete st-neurons not clustered to?

  cout << "Cluster count " << eventClusters.size() << endl;
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
}


vector<Complex> ClassLayer::evaluate(const vector<Complex>& stLayerGains)
{
  vector<Complex> gains;

  // compute the gain for each class neuron
  for(auto& neuron : neurons)
    gains.emplace_back(neuron.computeGain(stLayerGains));

  return gains;
}


void ClassLayer::loadFile(const string& filename)
{
  YAML::Node doc = YAML::LoadFile(filename);
  const YAML::Node& clnNode = doc["class-layer"]["class-neurons"];

  for(size_t i = 0; i < clnNode.size(); i++)
  {
    //add class layer neuron for group
    neurons.emplace_back(clnNode[i]["class-group"].as<int>());

    //TODO: load weights into added neuron
    //const YAML::Node& cwNode = clnNode[i]["weights"];
    //for(size_t k = 1; k < cwNode.size(); k++)
      //neurons.push_back(stnwNode[i].as<Event>());
  }
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

vector<double> CRBFNeuralNetwork::evaluateTrace(const string& traceFile)
{
  cout << "\x1B[35m==>\x1B[0m "
       << "Evaluating trace from " << traceFile << endl;

  //TODO: make thread safe
  TraceData trace(traceFile);
  vector<double> cumulativeErrors(cLayer.size(), 0);

  trace.normalizeEvents();

  //for each event in trace evaluate event
  for(auto &event : trace)
  {
    //evaluate spatio-temporal layer for given event, pass gains to class layer
    vector<Complex> clGains = cLayer.evaluate(stLayer.evaluate(event));

    //for each class add error to cumulativeErrors
    for(int k = 0; k < cLayer.size(); k++)
      cumulativeErrors[k] += computeErrorIncrement(clGains[k], trace.size());
  }

  return cumulativeErrors;
}


double CRBFNeuralNetwork::computeErrorIncrement(const Complex& gain, int tdSize) const
{
  double Ys = 1, // spacial priority
         Yt = 1; // temporal priority

  return (1/(tdSize*hypot(Ys,Yt))) * hypot(Ys*gain.amplitude, Yt*gain.phase/g_pi);
}


void CRBFNeuralNetwork::train(const string& traceFileList)
{
  cout << "\x1B[35m==>\x1B[0m " << "Training neural network" << endl;

  vector<TraceData> traces;

  // load all traces into trace for network training
  TraceData::loadTraceDataList(traceFileList, traces);

  // normalize trace data
  for(auto &t : traces)
    t.normalizeEvents();

  // export normalized trace data
  TraceData::exportTracesYamlFile("normalizedTrainingData.yaml", traces);

  // train Spatio Temporal Layer
  stLayer.train(traces, 50, 60000, 10);

  // train Class Layer
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


void CRBFNeuralNetwork::loadCRBFNeuralNetworkFile(const string& crbfFile)
{
  cout << "\x1B[36m==>\x1B[0m "
       <<  "Loading C-RBF neural network from " << crbfFile << endl;

  //TODO: load st-layer
  stLayer.loadFile(crbfFile);

  //TODO: load class-layer
  cLayer.loadFile(crbfFile);
}

/************************************************************
 * YAML parser/emitter
 ***********************************************************/

YAML::Emitter& operator<< (YAML::Emitter& out, const SpatioTemporalNeuron& v)
{
  out << v.getWeight();

  return out;
}


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
