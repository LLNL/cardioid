#include "heartIO.hh"

#include <sstream>
#include <cassert>

using namespace std;


void writeHeader(const string& filename,
		 const HeaderInfo& hdr,
		 int nx, int ny, int nz,
		 const set<int>& cellSet)
{
   FILE* file = fopen(filename.c_str(), "r+");

   stringstream buf;

   unsigned keyLocal;
   memcpy(&keyLocal, "1234", 4);

   buf << 
      "anatomy FILEHEADER {\n"
      "  datatype = "    << hdr.dataType << ";\n"
      "  nfiles = "      << hdr.nFiles << ";"
      "  nrecord = "     << hdr.nRec << ";"
      "  lrec = "        << hdr.lRec << ";"
      "  endian_key = "  << keyLocal << ";\n"
      "  nfields = "     << hdr.nFields << ";\n"
      "  field_names = " << hdr.fieldNames << ";\n"
      "  field_types = " << hdr.fieldTypes << ";\n"
      "  nx = " << nx << "; ny = "<< ny << "; nz = " << nz <<";\n"
      "}"; // do *not* put a new line here!

   assert(buf.str().size() < 2046);
   fprintf(file, "%s", buf.str().c_str());
   fclose(file);
}

string headSpace(int size)
{
   string s(size-2, ' ');
   s.push_back('\n');
   s.push_back('\n');
   return s;
}
