#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include "SimpleInputParser.hh"
#include "Control.hh"
using namespace std;

#include <mpi.h>

////////////////////////////////////////////////////////////////////////////////
SimpleInputParser::SimpleInputParser()
{
  inputread_ = false;
}
////////////////////////////////////////////////////////////////////////////////
int SimpleInputParser::readInput(const char* inputfile, Control *ctrl)
{

  int npes, mype;
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);  


  // open input file on pe 0, notify all pes if open was successful
  ifstream input;
  int inputopenerr = -1;
  if (mype == 0)
  {
    input.open(inputfile,ifstream::in);
    inputopenerr = 0;
    if (!input.is_open())
      inputopenerr = 1;
  }  
  MPI_Bcast(&inputopenerr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (inputopenerr == 1)
    return 1;

  int readlines = 1;
  string line;
  while (readlines)
  {
    // read line, broadcast to all pes
    int len;
    if ( mype == 0 )
    {
      if (input)
      {
        getline(input,line);
        len = line.length();
      }
      else
        len = -1;
    }
    MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
    if (len == -1)
    {
      readlines = 0;
    }
    else {
      char* buf = new char[len+1];
      if ( mype == 0 )
      {
        line.copy(buf,string::npos);
        buf[len]=0;
        assert(buf[len]=='\0');
      }
      MPI_Bcast(buf,len+1,MPI_CHAR,0,MPI_COMM_WORLD);
      line = buf;
    }

    // parse line
    transform(line.begin(),line.end(),line.begin(), ::tolower);
    string sub;
    vector<string> words;
    istringstream iss(line);
    while (iss)
    {
      iss >> sub;
      // convert 
      words.push_back(sub);
    }
    words.pop_back();
    
    //// store parsed values in Control struct ////
    if (words.size() > 0)       // skip blank lines
    {
      if (words[0] == "#")
      {
        // ignore lines starting with #
      }
      else if (words[0] == "electrophys")
      {
        ctrl->ephysmodel = words[1];
        if (mype == 0)
          cout << "Input parser:  electrophys = " << words[1] << endl;
      }
      else if (words[0] == "mechanical")
      {
        ctrl->mechmodel = words[1];
        if (mype == 0)
          cout << "Input parser:  mechanical = " << words[1] << endl;
      }
      else if (words[0] == "stimamplitude")
      {
        ctrl->stimAmplitude = atof(words[1].c_str());
        if (mype == 0)
          cout << "Input parser:  stimAmplitude = " << ctrl->stimAmplitude << endl;
      }
      else if (words[0] == "stimlength")
      {
        ctrl->stimLength = atoi(words[1].c_str());
        if (mype == 0)
          cout << "Input parser:  stimLength = " << ctrl->stimLength << endl;
      }
      else if (words[0] == "stiminterval")
      {
        ctrl->stimInterval = atoi(words[1].c_str());
        if (mype == 0)
          cout << "Input parser:  stimInterval = " << ctrl->stimInterval << endl;
      }
      else if (words[0] == "dt")
      {
        ctrl->dt = atof(words[1].c_str());
        if (mype == 0)
          cout << "Input parser:  dt = " << ctrl->dt << endl;
      }
      else if (words[0] == "savestep")
      {
        ctrl->savestep = atoi(words[1].c_str());
        if (mype == 0)
          cout << "Input parser:  savestep = " << ctrl->savestep << endl;
      }
      else if (words[0] == "tend")
      {
        ctrl->laststep = atoi(words[1].c_str());
        if (mype == 0)
          cout << "Input parser:  laststep = " << ctrl->laststep << endl;
      }
      else if (words[0] == "grid")
      {
        if (words.size() < 5) {
          if (mype == 0)
            cout << "Input parser:  syntax error in grid keyword." << endl;
          return 1;
        }
        if (words[1] == "ep")
        {
          ctrl->ngrid_ep.resize(3);
          ctrl->ngrid_ep[0] = atoi(words[2].c_str());
          ctrl->ngrid_ep[1] = atoi(words[3].c_str());
          ctrl->ngrid_ep[2] = atoi(words[4].c_str());
          if (mype == 0)
            cout << "Input parser:  electrophysiology grid = " << ctrl->ngrid_ep[0]
                 << " x " << ctrl->ngrid_ep[1] << " x " << ctrl->ngrid_ep[2] << endl;
        }
        else if (words[1] == "mech")
        {
          ctrl->ngrid_mech.resize(3);
          ctrl->ngrid_mech[0] = atoi(words[2].c_str());
          ctrl->ngrid_mech[1] = atoi(words[3].c_str());
          ctrl->ngrid_mech[2] = atoi(words[4].c_str());
          if (mype == 0)
            cout << "Input parser:  mechanical grid = " << ctrl->ngrid_mech[0]
                 << " x " << ctrl->ngrid_mech[1] << " x " << ctrl->ngrid_mech[2] << endl;
        }
        else
        {
          if (mype == 0)
            cout << "Input parser:  syntax error in grid keyword." << endl;
          return 1;
        }
      }
      else if (words[0] == "proc_grid")
      {
        if (words.size() < 5) {
          if (mype == 0)
            cout << "Input parser:  syntax error in proc_grid keyword." << endl;
          return 1;
        }
        if (words[1] == "ep")
        {
          ctrl->npegrid_ep.resize(3);
          ctrl->npegrid_ep[0] = atoi(words[2].c_str());
          ctrl->npegrid_ep[1] = atoi(words[3].c_str());
          ctrl->npegrid_ep[2] = atoi(words[4].c_str());
          if (mype == 0)
            cout << "Input parser:  electrophysiology process grid = " << ctrl->npegrid_ep[0]
                 << " x " << ctrl->npegrid_ep[1] << " x " << ctrl->npegrid_ep[2] << endl;
        }
        else if (words[1] == "mech")
        {
          ctrl->npegrid_mech.resize(3);
          ctrl->npegrid_mech[0] = atoi(words[2].c_str());
          ctrl->npegrid_mech[1] = atoi(words[3].c_str());
          ctrl->npegrid_mech[2] = atoi(words[4].c_str());
          if (mype == 0)
            cout << "Input parser:  mechanical process grid = " << ctrl->npegrid_mech[0]
                 << " x " << ctrl->npegrid_mech[1] << " x " << ctrl->npegrid_mech[2] << endl;
        }
        else
        {
          if (mype == 0)
            cout << "Input parser:  syntax error in proc_grid keyword." << endl;
          return 1;
        }
      }
      else
      {
        if (mype == 0)
          cout << "WARNING:  unknown keyword " << words[0] << endl;
      }
    }
    words.clear();    
  }
  input.close();
  inputread_ = true;

  return 0;
}
