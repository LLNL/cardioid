#include <iostream>
#include <vector>
#include "GridPoint.hh"
#include "Grid3DStencil.hh"
using std::vector;

Grid3DStencil::Grid3DStencil(int gid, int nx, int ny, int nz) : gid_(gid),
    nx_(nx), ny_(ny), nz_(nz)
{
  const bool faces = true;
  const bool edges = true;
  const bool corners = false;

  GridPoint gpt(gid_,nx_,ny_,nz_);
  if (faces)
  {
    if (gpt.x < nx_-1)
    {
      int newgid = (gpt.x+1) + gpt.y*nx_ + gpt.z*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x > 0 )
    {
      int newgid = (gpt.x-1) + gpt.y*nx_ + gpt.z*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.y < ny_-1)
    {
      int newgid = gpt.x + (gpt.y+1)*nx_ + gpt.z*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.y > 0 )
    {
      int newgid = gpt.x + (gpt.y-1)*nx_ + gpt.z*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.z < nz_-1)
    {
      int newgid = gpt.x + gpt.y*nx_ + (gpt.z+1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.z > 0 )
    {
      int newgid = gpt.x + gpt.y*nx_ + (gpt.z-1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }    
  }
  if (edges)
  {
    if (gpt.x < nx_-1 && gpt.y < ny_-1)
    {
      int newgid = (gpt.x+1) + (gpt.y+1)*nx_ + gpt.z*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x > 0 && gpt.y < ny_-1)
    {
      int newgid = (gpt.x-1) + (gpt.y+1)*nx_ + gpt.z*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x < nx_-1 && gpt.y > 0)
    {
      int newgid = (gpt.x+1) + (gpt.y-1)*nx_ + gpt.z*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x > 0 && gpt.y > 0)
    {
      int newgid = (gpt.x-1) + (gpt.y-1)*nx_ + gpt.z*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x < nx_-1 && gpt.z < nz_-1)
    {
      int newgid = (gpt.x+1) + gpt.y*nx_ + (gpt.z+1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x > 0 && gpt.z < nz_-1)
    {
      int newgid = (gpt.x-1) + gpt.y*nx_ + (gpt.z+1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x < nx_-1 && gpt.z > 0)
    {
      int newgid = (gpt.x+1) + gpt.y*nx_ + (gpt.z-1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x > 0 && gpt.z > 0)
    {
      int newgid = (gpt.x-1) + gpt.y*nx_ + (gpt.z-1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.y < ny_-1 && gpt.z < nz_-1)
    {
      int newgid = gpt.x + (gpt.y+1)*nx_ + (gpt.z+1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.y > 0 && gpt.z < nz_-1)
    {
      int newgid = gpt.x + (gpt.y-1)*nx_ + (gpt.z+1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.y < ny_-1 && gpt.z > 0)
    {
      int newgid = gpt.x + (gpt.y+1)*nx_ + (gpt.z-1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.y > 0 && gpt.z > 0)
    {
      int newgid = gpt.x + (gpt.y-1)*nx_ + (gpt.z-1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
  }
  if (corners)
  {
    if (gpt.x < nx_-1 && gpt.y < ny_-1 && gpt.z < nz_-1)
    {
      int newgid = (gpt.x+1) + (gpt.y+1)*nx_ + (gpt.z+1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x > 0 && gpt.y < ny_-1 && gpt.z < nz_-1)
    {
      int newgid = (gpt.x-1) + (gpt.y+1)*nx_ + (gpt.z+1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x > 0 && gpt.y > 0 && gpt.z < nz_-1)
    {
      int newgid = (gpt.x-1) + (gpt.y-1)*nx_ + (gpt.z+1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x > 0 && gpt.y < ny_-1 && gpt.z > 0)
    {
      int newgid = (gpt.x-1) + (gpt.y+1)*nx_ + (gpt.z-1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x < nx_-1 && gpt.y > 0 && gpt.z < nz_-1)
    {
      int newgid = (gpt.x+1) + (gpt.y-1)*nx_ + (gpt.z+1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x < nx_-1 && gpt.y > 0 && gpt.z > 0)
    {
      int newgid = (gpt.x+1) + (gpt.y-1)*nx_ + (gpt.z-1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x < nx_-1 && gpt.y < ny_-1 && gpt.z > 0)
    {
      int newgid = (gpt.x+1) + (gpt.y+1)*nx_ + (gpt.z-1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
    if (gpt.x > 0 && gpt.y > 0 && gpt.z > 0)
    {
      int newgid = (gpt.x-1) + (gpt.y-1)*nx_ + (gpt.z-1)*nx_*ny_;
      nbr_gids_.push_back(newgid);
    }
  }
  nstencil_ = nbr_gids_.size();
}

Grid3DStencil::~Grid3DStencil(void)
{

}

