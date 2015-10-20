#include <c-rbf.h>

#include <fstream>

using namespace std;


int main(int argc, char** argv)
{
  CRBFNeuralNetwork neuralNet;

  neuralNet.loadCRBFNeuralNetworkFile(argv[1]);

  ifstream fileList(argv[2]);
  string traceFile;

  while(getline(fileList, traceFile))
    neuralNet.evaluateTrace(traceFile);

  return 0;
}
