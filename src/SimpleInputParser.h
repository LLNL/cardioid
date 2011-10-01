#ifndef SIMPLE_INPUT_PARSER_H
#define SIMPLE_INPUT_PARSER_H

#include <string>
#include <vector>
#include "Control.h"
using namespace std;

class SimpleInputParser
{
  private:

  bool inputread_;

  public:

  SimpleInputParser();
  int readInput(char* inputfile, Control* ctrl);
};
#endif
