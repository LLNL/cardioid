/**

   How to convert this code to work for any other model:

   - Search/Replace the model name with your own specific string in the header and source files
   - Add your own code to EDIT_FLAGS and EDIT_PARAMETERS
   - Add your own code to EDIT_PERCELL_FLAGS and EDIT_PERCELL_PARAMETERS
   - Add your own states to EDIT_STATE
   - Add your computation code to the main calc routine, copy pasting frmo matlab.
   
 */


#include "sundnes_et_al_2016_FHN.hh"
#include "object_cc.hh"
#include "mpiUtils.h"
#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>

using namespace std;

#define setDefault(name, value) objectGet(obj, #name, name, TO_STRING(value))

template<typename TTT>
static inline string TO_STRING(const TTT x)
{
   stringstream ss;
   ss << x;
   return ss.str();
}

static const char* interpName[] = {
    NULL
};


   REACTION_FACTORY(sundnes_et_al_2016_FHN)(OBJECT* obj, const double _dt, const int numPoints, const ThreadTeam&)
   {
      sundnes_et_al_2016_FHN::ThisReaction* reaction = new sundnes_et_al_2016_FHN::ThisReaction(numPoints, _dt);

      //override the defaults
      //EDIT_PARAMETERS
      double Vpeak;
      double Vrest;
      double Vthresh;
      double a = 0.130000000000000;
      setDefault(Vrest, -85.0000000000000);
      setDefault(Vpeak, 40.0000000000000);
      double vamp = Vpeak - Vrest;
      setDefault(Vthresh, Vrest + a*vamp);
      reaction->Vpeak = Vpeak;
      reaction->Vrest = Vrest;
      reaction->Vthresh = Vthresh;
      bool reusingInterpolants = false;
      string fitName;
      objectGet(obj, "fit", fitName, "");
      int funcCount = sizeof(reaction->_interpolant)/sizeof(reaction->_interpolant[0])-1; //BGQ_HACKFIX, compiler bug with zero length arrays
      if (fitName != "")
      {
         OBJECT* fitObj = objectFind(fitName, "FIT");
         if (1
         )
         {
            vector<string> functions;
            objectGet(fitObj, "functions", functions);
            OBJECT* funcObj;
            if (functions.size() == funcCount)
            {
               for (int _ii=0; _ii<functions.size(); _ii++)
               {
                  OBJECT* funcObj = objectFind(functions[_ii], "FUNCTION");
                  objectGet(funcObj, "numer", reaction->_interpolant[_ii].numNumer_, "-1");
                  objectGet(funcObj, "denom", reaction->_interpolant[_ii].numDenom_, "-1");
                  objectGet(funcObj, "coeff", reaction->_interpolant[_ii].coeff_);
               }
               reusingInterpolants = true;
            }
         }
      }

      if (!reusingInterpolants)
      {
      reaction->createInterpolants(_dt);

      //save the interpolants
      if (funcCount > 0 && getRank(0) == 0)
      {
         ofstream outfile((string(obj->name) +".fit.data").c_str());
         outfile.precision(16);
         fitName = string(obj->name) + "_fit";
         outfile << obj->name << " REACTION { fit=" << fitName << "; }\n";
         outfile << fitName << " FIT {\n";
         outfile << "   functions = ";
         for (int _ii=0; _ii<funcCount; _ii++) {
            outfile << obj->name << "_interpFunc" << _ii << "_" << interpName[_ii] << " ";
         }
         outfile << ";\n";
         outfile << "}\n";

         for (int _ii=0; _ii<funcCount; _ii++)
         {
            outfile << obj->name << "_interpFunc" << _ii << "_" << interpName[_ii] << " FUNCTION { "
                    << "numer=" << reaction->_interpolant[_ii].numNumer_ << "; "
                    << "denom=" << reaction->_interpolant[_ii].numDenom_ << "; "
                    << "coeff=";
            for (int _jj=0; _jj<reaction->_interpolant[_ii].coeff_.size(); _jj++)
            {
               outfile << reaction->_interpolant[_ii].coeff_[_jj] << " ";
            }
            outfile << "; }\n";
         }
         outfile.close();
      }
      }
#ifdef USE_CUDA
      reaction->constructKernel();
#endif
      return reaction;
   }
#undef setDefault

namespace sundnes_et_al_2016_FHN 
{

void ThisReaction::createInterpolants(const double _dt) {

}

#ifdef USE_CUDA

void generateInterpString(stringstream& ss, const Interpolation& interp, const char* interpVar)
{
   ss <<
   "{\n"
   "   const double _numerCoeff[]={";
    for (int _ii=interp.numNumer_-1; _ii>=0; _ii--)
   {
      if (_ii != interp.numNumer_-1) { ss << ", "; }
      ss << interp.coeff_[_ii];
   }
   ss<< "};\n"
   "   const double _denomCoeff[]={";
   for (int _ii=interp.numDenom_+interp.numNumer_-2; _ii>=interp.numNumer_; _ii--)
   {
      ss << interp.coeff_[_ii] << ", ";
   }
   ss<< "1};\n"
   "   double _inVal = " << interpVar << ";\n"
   "   double _numerator=_numerCoeff[0];\n"
   "   for (int _jj=1; _jj<sizeof(_numerCoeff)/sizeof(_numerCoeff[0]); _jj++)\n"
   "   {\n"
   "      _numerator = _numerCoeff[_jj] + _inVal*_numerator;\n"
   "   }\n"
   "   if (sizeof(_denomCoeff)/sizeof(_denomCoeff[0]) == 1)\n"
   "   {\n"
   "      _ratPoly = _numerator;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      double _denominator=_denomCoeff[0];\n"
   "      for (int _jj=1; _jj<sizeof(_denomCoeff)/sizeof(_denomCoeff[0]); _jj++)\n"
   "      {\n"
   "         _denominator = _denomCoeff[_jj] + _inVal*_denominator;\n"
   "      }\n"
   "      _ratPoly = _numerator/_denominator;\n"
   "   }\n"
   "}"
      ;
}

void ThisReaction::constructKernel()
{

   stringstream ss;
   ss.precision(16);
   ss <<
   "enum StateOffset {\n"

   "   W_off,\n"
   "   NUMSTATES\n"
   "};\n"
   "extern \"C\"\n"
   "__global__ void sundnes_et_al_2016_FHN_kernel(const double* _Vm, const double* _iStim, double* _dVm, double* _state) {\n"
   "const double _dt = " << __cachedDt << ";\n"
   "const int _nCells = " << nCells_ << ";\n"

   "const double Vpeak = " << Vpeak << ";\n"
   "const double Vrest = " << Vrest << ";\n"
   "const double Vthresh = " << Vthresh << ";\n"
   "const int _ii = threadIdx.x + blockIdx.x*blockDim.x;\n"
   "if (_ii >= _nCells) { return; }\n"
   "const double V = _Vm[_ii];\n"
   "double _ratPoly;\n"

   "const double W = _state[_ii+W_off*_nCells];\n"
   "//get the gate updates (diagonalized exponential integrator)\n"
   "//get the other differential updates\n"
   "double b = 0.0130000000000000;\n"
   "double c3 = 1.00000000000000;\n"
   "double W_diff = b*(V - Vrest - W*c3);\n"
   "//get Iion\n"
   "double c1 = 0.260000000000000;\n"
   "double c2 = 0.100000000000000;\n"
   "double vamp = Vpeak - Vrest;\n"
   "double fhn1 = c1/(vamp*vamp);\n"
   "double fhn2 = c2/vamp;\n"
   "double Iion = W*fhn2*(V - Vrest) - fhn1*(-V + Vpeak)*(V - Vrest)*(V - Vthresh);\n"
   "//Do the markov update (1 step rosenbrock with gauss siedel)\n"
   "int _count=0;\n"
   "do\n"
   "{\n"
   "   _count++;\n"
   "} while (_count<50);\n"
   "//EDIT_STATE\n"
   "_state[_ii+W_off*_nCells] += _dt*W_diff;\n"
   "_dVm[_ii] = -Iion;\n"
   "}\n";

   _program_code = ss.str();
   //cout << ss.str();
   nvrtcCreateProgram(&_program,
                      _program_code.c_str(),
                      "sundnes_et_al_2016_FHN_program",
                      0,
                      NULL,
                      NULL);
   nvrtcCompileProgram(_program,
                       0,
                       NULL);
   std::size_t size;
   nvrtcGetPTXSize(_program, &size);
   _ptx.resize(size);
   nvrtcGetPTX(_program, &_ptx[0]);

   cuModuleLoadDataEx(&_module, &_ptx[0], 0, 0, 0);
   cuModuleGetFunction(&_kernel, _module, "sundnes_et_al_2016_FHN_kernel");
}

void ThisReaction::calc(double dt,
                const Managed<ArrayView<double>> Vm_m,
                const Managed<ArrayView<double>> iStim_m,
                Managed<ArrayView<double>> dVm_m)
{
   ArrayView<double> state = stateTransport_.modifyOnDevice();

   {
      int errorCode=-1;
      if (blockSize_ == -1) { blockSize_ = 1024; }
      while(1)
      {
         ConstArrayView<double> Vm = Vm_m.readOnDevice();
         ConstArrayView<double> iStim = iStim_m.readOnDevice();
         ArrayView<double> dVm = dVm_m.modifyOnDevice();
         double* VmRaw = const_cast<double*>(&Vm[0]);
         double* iStimRaw = const_cast<double*>(&iStim[0]);
         double* dVmRaw = &dVm[0];
         double* stateRaw= &state[0];
         void* args[] = { &VmRaw,
                          &iStimRaw,
                          &dVmRaw,
                          &stateRaw};
         int errorCode = cuLaunchKernel(_kernel,
                                        (nCells_+blockSize_-1)/blockSize_, 1, 1,
                                        blockSize_,1,1,
                                        0, NULL,
                                        args, 0);
         if (errorCode == CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES && blockSize_ > 0)
         {
            blockSize_ /= 2;
            continue;
         }
         else if (errorCode != CUDA_SUCCESS)
         {
            printf("Cuda return %d;\n", errorCode);
            assert(0 && "Launch failed");
         }
         else
         {
            break;
         }
         //cuCtxSynchronize();
      }
   }
}

enum StateOffset {
   _W_off,
   NUMSTATES
};

ThisReaction::ThisReaction(const int numPoints, const double __dt)
: nCells_(numPoints)
{
   stateTransport_.setup(PinnedVector<double>(nCells_*NUMSTATES));
   __cachedDt = __dt;
   blockSize_ = -1;
   _program = NULL;
   _module = NULL;
}

ThisReaction::~ThisReaction() {
   if (_program)
   {
      nvrtcDestroyProgram(&_program);
      _program = NULL;
   }
   if (_module)
   {
      cuModuleUnload(_module);
      _module = NULL;
   }
}

#else //USE_CUDA

#define width SIMDOPS_FLOAT64V_WIDTH
#define real simdops::float64v
#define load simdops::load

ThisReaction::ThisReaction(const int numPoints, const double __dt)
: nCells_(numPoints)
{
   state_.resize((nCells_+width-1)/width);
   __cachedDt = __dt;
}

void ThisReaction::calc(double _dt, const VectorDouble32& __Vm,
                       const vector<double>& __iStim , VectorDouble32& __dVm)
{
   //define the constants
   double b = 0.0130000000000000;
   double c1 = 0.260000000000000;
   double c2 = 0.100000000000000;
   double c3 = 1.00000000000000;
   double vamp = Vpeak - Vrest;
   double fhn1 = c1/(vamp*vamp);
   double fhn2 = c2/vamp;
   for (unsigned __jj=0; __jj<(nCells_+width-1)/width; __jj++)
   {
      const int __ii = __jj*width;
      //set Vm
      const real V = load(&__Vm[__ii]);
      const real iStim = load(&__iStim[__ii]);

      //set all state variables
      real W=load(state_[__jj].W);
      //get the gate updates (diagonalized exponential integrator)
      //get the other differential updates
      real W_diff = b*(V - Vrest - W*c3);
      //get Iion
      real Iion = W*fhn2*(V - Vrest) - fhn1*(-V + Vpeak)*(V - Vrest)*(V - Vthresh);
      //Do the markov update (1 step rosenbrock with gauss siedel)
      //EDIT_STATE
      W += _dt*W_diff;
      store(state_[__jj].W, W);
      simdops::store(&__dVm[__ii],-Iion);
   }
}
#endif //USE_CUDA
   
string ThisReaction::methodName() const
{
   return "sundnes_et_al_2016_FHN";
}

#ifdef USE_CUDA
void ThisReaction::initializeMembraneVoltage(ArrayView<double> __Vm)
#else //USE_CUDA
void ThisReaction::initializeMembraneVoltage(VectorDouble32& __Vm)
#endif //USE_CUDA
{
   assert(__Vm.size() >= nCells_);

#ifdef USE_CUDA
#define READ_STATE(state,index) (stateData[_##state##_off*nCells_+index])
   ArrayView<double> stateData = stateTransport_;
#else //USE_CUDA
#define READ_STATE(state,index) (state_[index/width].state[index % width])
   state_.resize((nCells_+width-1)/width);
#endif //USE_CUDA


   double V_init = Vrest;
   double V = V_init;
   double W_init = 1;
   double W = W_init;
   for (int iCell=0; iCell<nCells_; iCell++)
   {
      READ_STATE(W,iCell) = W;
   }

   __Vm.assign(__Vm.size(), V_init);
}

enum varHandles
{
   W_handle,
   NUMHANDLES
};

const string ThisReaction::getUnit(const std::string& varName) const
{
   if(0) {}
   else if (varName == "W") { return "1"; }
   return "INVALID";
}

int ThisReaction::getVarHandle(const std::string& varName) const
{
   if (0) {}
   else if (varName == "W") { return W_handle; }
   return -1;
}

void ThisReaction::setValue(int iCell, int varHandle, double value) 
{
#ifdef USE_CUDA
   ArrayView<double> stateData = stateTransport_;
#endif //USE_CUDA



   if (0) {}
   else if (varHandle == W_handle) { READ_STATE(W,iCell) = value; }
}


double ThisReaction::getValue(int iCell, int varHandle) const
{
#ifdef USE_CUDA
   ConstArrayView<double> stateData = stateTransport_;
#endif //USE_CUDA


   if (0) {}
   else if (varHandle == W_handle) { return READ_STATE(W,iCell); }
   return NAN;
}

double ThisReaction::getValue(int iCell, int varHandle, double V) const
{
#ifdef USE_CUDA
   ConstArrayView<double> stateData = stateTransport_;
#endif //USE_CUDA


   const double W=READ_STATE(W,iCell);
   if (0) {}
   else if (varHandle == W_handle)
   {
      return W;
   }
   return NAN;
}

void ThisReaction::getCheckpointInfo(vector<string>& fieldNames,
                                     vector<string>& fieldUnits) const
{
   fieldNames.clear();
   fieldUnits.clear();
   fieldNames.push_back("W");
   fieldUnits.push_back(getUnit("W"));
}

}
