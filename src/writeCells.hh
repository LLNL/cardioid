#ifndef WRITE_CELLS_HH
#define WRITE_CELLS_HH

#include <string>
class AnatomyReader;


void writeCells(const AnatomyReader& anatomy,
		const std::string& filename);

#endif
