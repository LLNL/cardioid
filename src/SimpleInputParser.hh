#ifndef SIMPLE_INPUT_PARSER_H
#define SIMPLE_INPUT_PARSER_H

class Control;

class SimpleInputParser
{
  private:

  bool inputread_;

  public:

  SimpleInputParser();
  int readInput(const char* inputfile, Control* ctrl);
};
#endif
