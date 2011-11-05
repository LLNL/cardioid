#ifndef TT04_CELLML_REACTION_HH
#define TT04_CELLML_REACTION_HH

#include "Reaction.hh"
class Anatomy;
class TT04_CellML;


class TT04_CellML_Reaction : public Reaction
{
 public:
   TT04_CellML_Reaction(const Anatomy& anatomy);
   // copy constructor and assignment operator intentionally
   // left unimplemented.
   TT04_CellML_Reaction(const TT04_CellML_Reaction&);
   TT04_CellML_Reaction& operator=(const TT04_CellML_Reaction&);
   ~TT04_CellML_Reaction();
   
   void calc(double dt,
	     const std::vector<double>& Vm,
	     const std::vector<double>& iStim,
	     std::vector<double>& dVm);

 private:
   unsigned nCells_;
   std::vector<int>          ttType_; // maps cellType to ttType
   std::vector<TT04_CellML*> cellModel_;

};

#endif
