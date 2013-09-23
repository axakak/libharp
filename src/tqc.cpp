

#include "traceData.h"
#include <iostream>

using namespace std;


int main(int argc, char** argv)
{

  TraceData trace;

  trace.loadYamlFile(argv[1]);

  cout << "patient-id: " << trace.getPaitentID() << endl
       << "total-time: " << trace.getTotalTime() << endl
       << endl;

  return 0;
}

