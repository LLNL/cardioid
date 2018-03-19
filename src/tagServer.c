#include "tagServer.h"
#include <assert.h>

static int getTagUpperBound(void);

/** MPI messages match senders to receivers by three items: The
 *  send/recv rank, the tag, and the MPI_Comm.  This routine provides a
 *  tag that is not in use (at least not by any other client of the
 *  tagServer.
 *
 *  This routine does *not* guarantee that all returned tags are unique,
 *  only that a tag is never reused on the same communicator.  This is
 *  accomplished by processing all tag assignment on rank 0 of the comm.
 *  Since rank 0 of a comm can't change this assures that only unique
 *  tags are returned.  It does not promise that tags will be sequential
 *  since other comms that have the same task as rank 0 will pull from
 *  the same set of tags.
 *
 *  We assume that the application will never ask for more than about
 *  2^31 unique tags.  If this occurs, this routine will crash.
 */
int getUniqueTag(MPI_Comm comm)
{
   static int nextTag = 1000;
   static int tagUpperBound = 0;
   if (tagUpperBound == 0)
      tagUpperBound = getTagUpperBound();

   int myRank;
   MPI_Comm_rank(comm, &myRank);
   int tag;
   if (myRank == 0)
      tag = nextTag++;
   assert(nextTag < tagUpperBound);
   MPI_Bcast(&tag, 1, MPI_INT, 0, comm);
   return tag;
}


int getTagUpperBound(void)
{
   void* v;
   int flag;
   // MPI-2  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &v, &flag);
   MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, &v, &flag);
   int tagUpperBound = *(int*) v;
   return tagUpperBound;
}

