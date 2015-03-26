#include "c-rbf.h"

using namespace std;


int main(int argc, char** argv)
{
  CRBFNeuralNetwork neuralNet;

  neuralNet.train(argv[1]);

  neuralNet.exportYamlFile(argv[2]);

  return 0;
}
