
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

  //trace.normalizeEvents();

  for(int i = 0; i < trace.size(); i+=1000)
    cout << trace.at(i) << endl;

  char exportFile[] = "exportTrace.txt";
    
  trace.exportYamlFile(exportFile);

  return 0;
}

