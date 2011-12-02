#ifndef TT06_CELLML_REACTION_HH
#define TT06_CELLML_REACTION_HH

#include "Reaction.hh"
class Anatomy;
class TT06_CellML;


class TT06_CellML_Reaction : public Reaction
{
 public:
   TT06_CellML_Reaction(const Anatomy& anatomy);
   // copy constructor and assignment operator intentionally
   // left unimplemented.
   TT06_CellML_Reaction(const TT06_CellML_Reaction&);
   TT06_CellML_Reaction& operator=(const TT06_CellML_Reaction&);
   ~TT06_CellML_Reaction();
   
   void calc(double dt,
             const std::vector<double>& Vm,
             const std::vector<double>& iStim,
             std::vector<double>& dVm);

 private:
   unsigned nCells_;
   std::vector<int>          ttType_; // maps cellType to ttType
   std::vector<TT06_CellML*> cellModel_;

};

#endif
