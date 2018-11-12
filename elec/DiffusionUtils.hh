#ifndef DIFFUSION_UTILS_HH
#define DIFFUSION_UTILS_HH

class LocalGrid;
class Anatomy;

namespace DiffusionUtils
{
   LocalGrid findBoundingBox(const Anatomy& anatomy, bool verbose);
   LocalGrid findBoundingBox_simd(const Anatomy& anatomy, bool verbose);
}

#endif
