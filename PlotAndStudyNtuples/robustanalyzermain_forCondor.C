#include "robustanalyzer.hh"
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char* argv[]) {
  
  stringstream ss;
  ss << argv[1];
  int cnt;
  ss >> cnt;
  
  stringstream ssargv2;
  ssargv2 << argv[2];
  int numCores;
  ssargv2 >> numCores;

  try {    
    robustanalyzer rana_data(argv[3], argv[4], numCores, false, false);
    rana_data.analyzersinglefile(cnt);
  }
  catch (char const* exc) {
    cerr<<exc<<endl;
  }
  
  return -1;
}
