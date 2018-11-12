#include "unionOfStrings.hh"

#include <iterator>
#include <algorithm>

using namespace std;

namespace
{
   string pack(const vector<string>& words);
   vector<string> unpack(const string& buf);
}


void unionOfStrings(vector<string>& words, MPI_Comm comm, int tag)
{
   int nTasks;
   int myRank;
   MPI_Comm_size(comm, &nTasks);
   MPI_Comm_rank(comm, &myRank);

   unsigned nTree = nTasks;

   while (nTree > 1)
   {
      unsigned midpoint = (nTree+1)/2;
      if (myRank < nTree && myRank >= midpoint)
      {
         string buf = pack(words);
	 unsigned dest = myRank - midpoint;
         int bufSize = buf.size();
	 MPI_Send(&bufSize, 1, MPI_INT, dest, tag, comm);
	 MPI_Send(&buf[0], buf.size(), MPI_CHAR, dest, tag, comm);
      }
      else
      {
	 unsigned source = myRank + midpoint;
	 if (source < nTree)
	 {
            string buf;
            int bufSize;
	    MPI_Recv(&bufSize, 1, MPI_INT, source, tag, comm, MPI_STATUS_IGNORE);
	    buf.resize(bufSize);
	    MPI_Recv(&buf[0], buf.size(), MPI_CHAR, source, tag, comm, MPI_STATUS_IGNORE);
            vector<string> remoteWords = unpack(buf);
            vector<string> tmp;
            set_union(words.begin(), words.end(),
                      remoteWords.begin(), remoteWords.end(),
                      back_inserter(tmp));
            words = tmp;
	 }
      }
      nTree = midpoint;
   }
   
   string buf;
   if (myRank == 0)
      buf = pack(words);

   int bufSize = buf.size();
   MPI_Bcast(&bufSize, 1, MPI_INT, 0, comm);
   buf.resize(bufSize);
   MPI_Bcast(&buf[0], buf.size(), MPI_CHAR, 0, comm);

   words = unpack(buf);
}


namespace
{
   /** For now I'm going to be lazy and assume that no word contains a
    * "\n" and use that character as the delimiter of the packed
    * strings.  To really do things right there should be some sort of
    * function that looks through all of the words and picks a
    * delimiter character that does not appear in any word.  For that
    * matter, we could even provide for multi-character delimiters to
    * make sure that we can always find something.  In any case, the
    * delimiter (and its length) gets encoded at the start of the buffer
    * so that the unpack function can find it. */
   string pack(const vector<string>& words)
   {
      string delimiter("\n");
      string result = delimiter;
      for (unsigned ii=0; ii<words.size(); ++ii)
         result += words[ii] + delimiter;
      return result;
   }
}

namespace
{
   /** For now, asssume a single character delimiter.  We need to be a
    * little more clever with the find logic to handle mulit-character
    * delimiters. */
   vector<string> unpack(const string& buf)
   {
      vector<string> result;
      string delimiter = buf.substr(0, 1);
      string::size_type offset = 1;
      while (offset < buf.size())
      {
         string::size_type end = buf.find_first_of(delimiter, offset);
         result.push_back(buf.substr(offset, end-offset));
         offset = end + 1;
      }
      return result;
   }
}
