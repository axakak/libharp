#include "c-rbf.h"
//#include <string>

using namespace std;


int main(int argc, char** argv)
{
  CRBFNeuralNetwork neuralNet;

  neuralNet.loadCRBFNeuralNetworkFile(argv[1]);

  neuralNet.evaluateTrace(arv[2]);

  return 0;
}
