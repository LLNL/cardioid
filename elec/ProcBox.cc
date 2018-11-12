#include "ProcBox.hh"
#include <cassert>
#include <iostream>
//ewd DEBUG
#include "GridPoint.hh"
#include <mpi.h>
//ewd DEBUG
using namespace std;

ProcBox::ProcBox(int peid) : peid_(peid)
{
   addCnt_ = 0;
   planes_.resize(3);
   axisPoints_.resize(3);
   min_.resize(3);
   max_.resize(3);
   for (int ii=0; ii<3; ii++)
   {
      min_[ii] = 999999999;
      max_[ii] = -999999999;
   }
}

ProcBox::~ProcBox()
{
   // delete all Plane objects
   for (int ii=0; ii<3; ++ii)
      for (map<int,Plane*>::iterator it = planes_[ii].begin(); it != planes_[ii].end(); it++)
         if ((*it).second != 0)
            delete (*it).second;
}

int ProcBox::volume()
{
   /*
   int vol = 1;
   for (int ii=0; ii<3; ii++)
      vol *= (maxAxisPoint(ii) - minAxisPoint(ii) + 1);
   if (vol < 0) vol = 0;
   return vol;
   */

   int vol = 1;
   for (int ii=0; ii<3; ii++)
      vol *= (max_[ii] - min_[ii] + 1);
   if (vol < 0) vol = 0;
   return vol;

}

void ProcBox::add3DCoord(int x, int y, int z)
{
   addCnt_++;
   axisPoints_[0].insert(x);
   axisPoints_[1].insert(y);
   axisPoints_[2].insert(z);

   if (x < min_[0]) min_[0] = x;
   if (y < min_[1]) min_[1] = y;
   if (z < min_[2]) min_[2] = z;
   if (x > max_[0]) max_[0] = x;
   if (y > max_[1]) max_[1] = y;
   if (z > max_[2]) max_[2] = z;

   if (planes_[0][x] == 0)
      planes_[0][x] = new Plane(x,x,y,y,z,z);
   if (planes_[1][y] == 0)
      planes_[1][y] = new Plane(x,x,y,y,z,z);
   if (planes_[2][z] == 0)
      planes_[2][z] = new Plane(x,x,y,y,z,z);

   Plane* xp = planes_[0][x];
   if (y < xp->min[1]) xp->min[1] = y;
   if (y > xp->max[1]) xp->max[1] = y;
   if (z < xp->min[2]) xp->min[2] = z;
   if (z > xp->max[2]) xp->max[2] = z;
   Plane* yp = planes_[1][y];
   if (x < xp->min[0]) xp->min[0] = x;
   if (x > xp->max[0]) xp->max[0] = x;
   if (z < xp->min[2]) xp->min[2] = z;
   if (z > xp->max[2]) xp->max[2] = z;
   Plane* zp = planes_[2][z];
   if (x < xp->min[0]) xp->min[0] = x;
   if (x > xp->max[0]) xp->max[0] = x;
   if (y < xp->min[1]) xp->min[1] = y;
   if (y > xp->max[1]) xp->max[1] = y;
}

void ProcBox::planesToCoords(vector<int>& coords)
{
   int x,y,z;
   for (map<int,Plane*>::iterator it = planes_[0].begin(); it != planes_[0].end(); it++)
   {
      x = (*it).first;
      Plane* pp = (*it).second;
      if (pp != 0)
      {
         y = pp->min[1]; z = pp->min[2];
         coords.push_back(x);  coords.push_back(y);  coords.push_back(z);
         y = pp->max[1]; z = pp->min[2];
         coords.push_back(x);  coords.push_back(y);  coords.push_back(z);
         y = pp->min[1]; z = pp->max[2];
         coords.push_back(x);  coords.push_back(y);  coords.push_back(z);
         y = pp->max[1]; z = pp->max[2];
         coords.push_back(x);  coords.push_back(y);  coords.push_back(z);
      }
   }
   for (map<int,Plane*>::iterator it = planes_[1].begin(); it != planes_[1].end(); it++)
   {
      y = (*it).first;
      Plane* pp = (*it).second;
      if (pp != 0)
      {
         x = pp->min[0]; z = pp->min[2];
         coords.push_back(x);  coords.push_back(y);  coords.push_back(z);
         x = pp->max[0]; z = pp->min[2];
         coords.push_back(x);  coords.push_back(y);  coords.push_back(z);
         x = pp->min[0]; z = pp->max[2];
         coords.push_back(x);  coords.push_back(y);  coords.push_back(z);
         x = pp->max[0]; z = pp->max[2];
         coords.push_back(x);  coords.push_back(y);  coords.push_back(z);
      }
   }
   for (map<int,Plane*>::iterator it = planes_[2].begin(); it != planes_[2].end(); it++)
   {
      z = (*it).first;
      Plane* pp = (*it).second;
      if (pp != 0)
      {
         y = pp->min[1]; x = pp->min[0];
         coords.push_back(x);  coords.push_back(y);  coords.push_back(z);
         y = pp->max[1]; x = pp->min[0];
         coords.push_back(x);  coords.push_back(y);  coords.push_back(z);
         y = pp->min[1]; x = pp->max[0];
         coords.push_back(x);  coords.push_back(y);  coords.push_back(z);
         y = pp->max[1]; x = pp->max[0];
         coords.push_back(x);  coords.push_back(y);  coords.push_back(z);
      }
   }
}

void ProcBox::addPlane(int dim, Plane* plptr)
{
   int val = plptr->min[dim];     
   axisPoints_[dim].insert(val);
   
   if (planes_[dim][val] == 0)  // add pointer to planes_
      planes_[dim][val] = plptr;
   else  // update boundaries of current plane
   {
      for (int ii=0; ii<3; ii++)
         if (ii != dim)
         {
            if (plptr->min[ii] < planes_[dim][val]->min[ii])
               planes_[dim][val]->min[ii] = plptr->min[ii];
            if (plptr->max[ii] > planes_[dim][val]->max[ii])
               planes_[dim][val]->max[ii] = plptr->max[ii];
         }
   }

   // update plane boundaries in other dimensions
   for (int ii=0; ii<3; ii++)
   {
      if (ii != dim)
      {
         int pminii = plptr->min[ii];
         int pmaxii = plptr->max[ii];
         for (map<int,Plane*>::iterator it = planes_[ii].begin(); it != planes_[ii].end(); it++)
         {
            int qval = (*it).first;
            Plane* qptr = (*it).second;
            if (qptr != 0)
            {
               if (qval <= pmaxii && qval >= pminii)
               {
                  if (qptr->max[ii] < val)
                     qptr->max[ii] = val;
                  if (qptr->min[ii] > val)
                     qptr->min[ii] = val;
                  /*
                  // only update plane if it intersects in dim jj
                  for (int jj=0; jj<3; jj++)
                  {
                     if (jj != dim && jj != ii)
                     {
                        if (!(plptr->min[jj] > qptr->max[jj]) && !(plptr->max[jj] < qptr->min[jj]))
                           if (qptr->max[ii] < val)
                              qptr->max[ii] = val;
                           if (qptr->min[ii] > val)
                              qptr->min[ii] = val;
                     }
                  }
                  */
               }
            }
         }
      }
   }
}

void ProcBox::deletePlane(int dim, int val)
{
   //ewd:  this function isn't used currently, remove it?

   if (planes_[dim][val] != 0)
      delete planes_[dim][val];
   axisPoints_[dim].erase(val);
}

void ProcBox::removePlane(int dim, int val)
{
   //ewd: need to change min/max of other planes from this value
   //ewd: to next one down the axisPoint list
   for (int ii=0; ii<3; ii++)
   {
      if (ii != dim)
      {
         for (map<int,Plane*>::iterator it = planes_[ii].begin(); it != planes_[ii].end(); it++)
         {
            int plval = (*it).first;
            Plane* plptr = (*it).second;
            if (plptr != 0)
            {
               if (plptr->min[dim] == val && plptr->max[dim] == val)
               {
                  axisPoints_[ii].erase(plval);
                  (*it).second = 0;
               }
               else
               {
                  if (plptr->min[dim] == val)
                  {
                     int nextMin = secondMinAxisPoint(dim);
                     if (nextMin > -1)
                        plptr->min[dim] = nextMin;
                  }
                  if (plptr->max[dim] == val)
                  {
                     int nextMax = secondMaxAxisPoint(dim);
                     if (nextMax > -1)
                        plptr->max[dim] = nextMax;
                  }
               }
            }
         }
      }
   }

   if (planes_[dim][val] != 0)
      planes_[dim][val] = 0;
   axisPoints_[dim].erase(val);
}
    
void ProcBox::updateBoundaries()
{
   /*
   for (int ii=0; ii<3; ii++)
   {
      int tmin = 999999999;
      int tmax = -999999999;
      if (axisPoints_[ii].size() == 0)
      {
         for (int jj=0; jj<3; jj++)
            if (jj != ii)
               for (map<int,Plane*>::iterator it = planes_[jj].begin(); it != planes_[jj].end(); it++)
               {
                  Plane* pp = (*it).second;
                  if (pp != 0)
                  {
                     {
                        if (pp->min[ii] < tmin) tmin = pp->min[ii];
                        if (pp->max[ii] > tmax) tmax = pp->max[ii];
                     }
                  }
               }
      }
      else
      {
         tmin = minAxisPoint(ii);
         tmax = maxAxisPoint(ii);
      }
      min_[ii] = ( tmin > -1 ? tmin : 999999999);
      max_[ii] = ( tmax > -1 ? tmax : -999999999);
   }
   */
   
   // loop through all planes, compute min, max in all directions
   for (int ii=0; ii<3; ii++)
   {
      min_[ii] = 999999999;
      max_[ii] = -999999999;
   }
   for (int ii=0; ii<3; ii++)
      for (map<int,Plane*>::iterator it = planes_[ii].begin(); it != planes_[ii].end(); it++)
      {
         int val = (*it).first;
         Plane* pp = (*it).second;
         if (pp != 0)
         {
            if (val < min_[ii]) min_[ii] = val;
            if (val > max_[ii]) max_[ii] = val;
            for (int jj=0; jj<3; jj++)
               if (jj != ii)
               {
                  if (pp->min[jj] < min_[jj]) min_[jj] = pp->min[jj];
                  if (pp->max[jj] > max_[jj]) max_[jj] = pp->max[jj];
            }
         }
      }

}

int ProcBox::trialMove(int dim, int val, ProcBox& nbrbox)
{
   if (val < 0)
      return 1;
   
   assert(dim < 3 && dim >= 0);
   int myvol = volume();
   int nbrvol = nbrbox.volume();
   int maxvol = ( myvol > nbrvol ? myvol : nbrvol);
   
   Plane* planeptr  = planes_[dim][val];  // save plane pointer 
   Plane nbrcopy;  // save copy of current plane, if necessary (shouldn't be)
   int nbrcopied = 0;

   // move plane from this ProcBox to nbrbox, see if maxvol is reduced
   if (planeptr != 0 && myvol > 0)
   {

      //ewd DEBUG
      {
         int nbr = nbrbox.myRank();
         int npex_ = 32;
         int npey_ = 32;
         int npez_ = 32;
         GridPoint gpt(peid_,npex_,npey_,npez_);
         GridPoint nbrpt(nbr,npex_,npey_,npez_);
         int myRank_;
         MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
         if (myRank_ == 39)
            cout << "TRIAL_START:  dim = " << dim << ", val = " << val 
                 << ", pe " << peid_ << " (" << gpt.x << "," << gpt.y << "," << gpt.z << "), box = " <<
                minPoint(0) << " " << maxPoint(0) << " " <<
                minPoint(1) << " " << maxPoint(1) << " " <<
                minPoint(2) << " " << maxPoint(2) 
                 << ", nbr " << nbr << " (" << nbrpt.x << "," << nbrpt.y << "," << nbrpt.z << "), box = " <<
                nbrbox.minPoint(0) << " " << nbrbox.maxPoint(0) << " " <<
                nbrbox.minPoint(1) << " " << nbrbox.maxPoint(1) << " " <<
                nbrbox.minPoint(2) << " " << nbrbox.maxPoint(2)
                 << ", myvol = " << myvol << ", nbrvol = " << nbrvol << endl;
      }
      //ewd DEBUG


      int tvol = trialVolume(true,dim,planeptr);
      int tnbr = nbrbox.trialVolume(false,dim,planeptr);

      if (tvol > 0 && tnbr > 0 && tvol < maxvol && tnbr < maxvol) // accept move
      {
         removePlane(dim,val);
         nbrbox.addPlane(dim,planeptr);
         updateBoundaries();
         nbrbox.updateBoundaries();

         //ewd DEBUG
         {
            int nbr = nbrbox.myRank();
            int npex_ = 32;
            int npey_ = 32;
            int npez_ = 32;
            GridPoint gpt(peid_,npex_,npey_,npez_);
            GridPoint nbrpt(nbr,npex_,npey_,npez_);
            int myRank_;
            MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
            if (myRank_ == 39)
               cout << "TRIAL_ACCEPTED:  dim = " << dim << ", val = " << val 
                    << ", pe " << peid_ << " (" << gpt.x << "," << gpt.y << "," << gpt.z << "), box = " <<
                   minPoint(0) << " " << maxPoint(0) << " " <<
                   minPoint(1) << " " << maxPoint(1) << " " <<
                   minPoint(2) << " " << maxPoint(2) 
                    << ", nbr " << nbr << " (" << nbrpt.x << "," << nbrpt.y << "," << nbrpt.z << "), box = " <<
                   nbrbox.minPoint(0) << " " << nbrbox.maxPoint(0) << " " <<
                   nbrbox.minPoint(1) << " " << nbrbox.maxPoint(1) << " " <<
                   nbrbox.minPoint(2) << " " << nbrbox.maxPoint(2)
                    << ", myvol = " << tvol << ", nbrvol = " << tnbr << endl;
         }
         //ewd DEBUG

         return 0;
      }
      else // reject move
      {
         //ewd DEBUG
         {
            int nbr = nbrbox.myRank();
            int npex_ = 32;
            int npey_ = 32;
            int npez_ = 32;
            GridPoint gpt(peid_,npex_,npey_,npez_);
            GridPoint nbrpt(nbr,npex_,npey_,npez_);
            int myRank_;
            myvol = volume();
            nbrvol = nbrbox.volume();
            MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
            if (myRank_ == 39)
               cout << "TRIAL_REJECTED1:  dim = " << dim << ", val = " << val 
                    << ", pe " << peid_ << " (" << gpt.x << "," << gpt.y << "," << gpt.z << "), box = " <<
                   minPoint(0) << " " << maxPoint(0) << " " <<
                   minPoint(1) << " " << maxPoint(1) << " " <<
                   minPoint(2) << " " << maxPoint(2) 
                    << ", nbr " << nbr << " (" << nbrpt.x << "," << nbrpt.y << "," << nbrpt.z << "), box = " <<
                   nbrbox.minPoint(0) << " " << nbrbox.maxPoint(0) << " " <<
                   nbrbox.minPoint(1) << " " << nbrbox.maxPoint(1) << " " <<
                   nbrbox.minPoint(2) << " " << nbrbox.maxPoint(2)
                    << ", myvol = " << myvol << ", nbrvol = " << nbrvol << endl;
         }
         //ewd DEBUG
         return 1;
      }
   }
   else // planeptr == 0 || vol == 0
      return 1; 
}
int ProcBox::trialVolume(bool removePlane, int dim, Plane* pptr)
{
   vector<int> tmin(3), tmax(3), pmin(3), pmax(3);
   for (int ii=0; ii<3; ii++)
   {
      tmin[ii] = min_[ii];
      tmax[ii] = max_[ii];
      pmin[ii] = pptr->min[ii];
      pmax[ii] = pptr->max[ii];
   }
   int pval = pmin[dim];

   if (removePlane)
   {
      if (pval == tmin[dim])
         tmin[dim] = secondMinAxisPoint(dim);
      if (pval == tmax[dim])
         tmax[dim] = secondMaxAxisPoint(dim);
      if (tmax[dim] == -1 || tmin[dim] == -1)
         return -1;
   
      for (int ii=0; ii<3; ii++)
      {
         if (ii != dim)
         {
            if (pmin[ii] == tmin[ii])
               tmin[ii] = secondMinAxisPoint(ii);
            if (pmax[ii] == tmax[ii])
               tmax[ii] = secondMaxAxisPoint(ii);
            if (tmax[ii] == -1 || tmin[ii] == -1)
               return -1;
         }
      }
      int tvol = (tmax[2]-tmin[2])*(tmax[1]-tmin[1])*(tmax[0]-tmin[0]);
      return tvol;
   }
   else
   {
      for (int ii=0; ii<3; ii++)
      {
         if (pmin[ii] < tmin[ii])
            tmin[ii] = pmin[ii];
         if (pmax[ii] > tmax[ii])
            tmax[ii] = pmax[ii];
      }
      int tvol = (tmax[2]-tmin[2])*(tmax[1]-tmin[1])*(tmax[0]-tmin[0]);
      return tvol;
   }
}

int ProcBox::maxAxisPoint(int dim)
{
   int val = -1;
   if (axisPoints_[dim].size() > 0)
   {
      set<int>::iterator it;
      it = axisPoints_[dim].end();
      it--;
      val = *it;
   }
   return val;
}

int ProcBox::secondMaxAxisPoint(int dim)
{
   int val = -1;
   if (axisPoints_[dim].size() > 1)
   {
      set<int>::iterator it;
      it = axisPoints_[dim].end();
      it--;
      it--;
      val = *it;
   }
   return val;
}

int ProcBox::minAxisPoint(int dim)
{
   int val = -1;
   if (axisPoints_[dim].size() > 0)
   {
      set<int>::iterator it;
      it = axisPoints_[dim].begin();
      val = *it;
   }
   return val;
}

int ProcBox::secondMinAxisPoint(int dim)
{
   int val = -1;
   if (axisPoints_[dim].size() > 1)
   {
      set<int>::iterator it;
      it = axisPoints_[dim].begin();
      it++;
      val = *it;
   }
   return val;
}

bool ProcBox::containsPoint(int x, int y, int z)
{
   if (x > max_[0]) return false;
   if (x < min_[0]) return false;
   if (y > max_[1]) return false;
   if (y < min_[1]) return false;
   if (z > max_[2]) return false;
   if (z < min_[2]) return false;
   return true;
}
