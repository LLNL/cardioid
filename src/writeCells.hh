#ifndef WRITE_CELLS_HH
#define WRITE_CELLS_HH

#include <string>
class Simulate;
class Anatomy;


void writeCells(const Simulate& sim,    const std::string& filename);
void writeCells(const Anatomy& anatomy, const std::string& filename);

#endif
