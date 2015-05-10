#include <c-rbf.h>

using namespace std;


int main(int argc, char** argv)
{
  CRBFNeuralNetwork neuralNet;

  neuralNet.loadCRBFNeuralNetworkFile(argv[1]);

  neuralNet.evaluateTrace(argv[2]);

  return 0;
}
