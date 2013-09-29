

#include "traceData.h"
#include <iostream>

using namespace std;


int main(int argc, char** argv)
{

  TraceData trace;

  trace.loadYamlFile(argv[1]);

  cout << "patient-id: " << trace.getPaitentID() << endl
       << "total-time: " << trace.getTotalTime() << endl
       << "events:" << endl;

  for(int i = 0; i < 20; i++)
    cout << trace.getEvent(i) << endl;

  return 0;
}

