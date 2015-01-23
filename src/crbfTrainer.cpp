#include "c-rbf.h"
#include <string>

using namespace std;


int main(int argc, char** argv)
{
  CRBFNeuralNetwork neuralNet;

  neuralNet.train(argv[1]);

  neuralNet.exportYamlFile("crbfNeuralNet.yaml");

  return 0;
}
