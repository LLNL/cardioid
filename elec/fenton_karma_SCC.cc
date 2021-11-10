/**

   How to convert this code to work for any other model:

   - Search/Replace the model name with your own specific string in the header and source files
   - Add your own code to EDIT_FLAGS and EDIT_PARAMETERS
   - Add your own code to EDIT_PERCELL_FLAGS and EDIT_PERCELL_PARAMETERS
   - Add your own states to EDIT_STATE
   - Add your computation code to the main calc routine, copy pasting frmo matlab.
   
 */


#include "fenton_karma_SCC.hh"
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


   REACTION_FACTORY(fenton_karma_SCC)(OBJECT* obj, const double _dt, const int numPoints, const ThreadTeam&)
   {
      fenton_karma_SCC::ThisReaction* reaction = new fenton_karma_SCC::ThisReaction(numPoints, _dt);

      //override the defaults
      //EDIT_PARAMETERS
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

namespace fenton_karma_SCC 
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

   "   v_off,\n"
   "   w_off,\n"
   "   NUMSTATES\n"
   "};\n"
   "extern \"C\"\n"
   "__global__ void fenton_karma_SCC_kernel(const int* _indexArray, const double* _Vm, const double* _iStim, double* _dVm, double* _state) {\n"
   "const double _dt = " << __cachedDt << ";\n"
   "const int _nCells = " << nCells_ << ";\n"

   "const int _ii = threadIdx.x + blockIdx.x*blockDim.x;\n"
   "if (_ii >= _nCells) { return; }\n"
   "const double V = _Vm[_indexArray[_ii]];\n"
   "double _ratPoly;\n"

   "const double v = _state[_ii+v_off*_nCells];\n"
   "const double w = _state[_ii+w_off*_nCells];\n"
   "//get the gate updates (diagonalized exponential integrator)\n"
   "//get the other differential updates\n"
   "double V_init = -85;\n"
   "double V_fi = 15;\n"
   "double u = (V - V_init)/(V_fi - V_init);\n"
   "double u_c = 0.13;\n"
   "double __melodee_temp_000 = u < u_c;\n"
   "double __melodee_temp_001;\n"
   "if (__melodee_temp_000)\n"
   "{\n"
   "   __melodee_temp_001 = 0;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_001 = 1;\n"
   "}\n"
   "double p = __melodee_temp_001;\n"
   "double u_v = 0;\n"
   "double __melodee_temp_002 = u < u_v;\n"
   "double __melodee_temp_003;\n"
   "if (__melodee_temp_002)\n"
   "{\n"
   "   __melodee_temp_003 = 0;\n"
   "}\n"
   "else\n"
   "{\n"
   "   __melodee_temp_003 = 1;\n"
   "}\n"
   "double q = __melodee_temp_003;\n"
   "double tau_v1_minus = 18.199999999999999;\n"
   "double tau_v2_minus = 18.199999999999999;\n"
   "double tau_v_plus = 10;\n"
   "double tau_v_minus = q*tau_v1_minus + tau_v2_minus*(1 - q);\n"
   "double v_diff = -p*v/tau_v_plus + (1 - p)*(1 - v)/tau_v_minus;\n"
   "double tau_w_minus = 80;\n"
   "double tau_w_plus = 1020;\n"
   "double w_diff = -p*w/tau_w_plus + (1 - p)*(1 - w)/tau_w_minus;\n"
   "//get Iion\n"
   "double Cm = 1;\n"
   "double g_fi_max = 5.7999999999999998;\n"
   "double tau_d = Cm/g_fi_max;\n"
   "double J_fi = -p*v*(1 - u)*(u - u_c)/tau_d;\n"
   "double k = 10;\n"
   "double tau_si = 127;\n"
   "double u_csi = 0.84999999999999998;\n"
   "double _expensive_functions = tanh(k*(u - u_csi));\n"
   "double J_si = -1.0/2.0*w*(_expensive_functions + 1)/tau_si;\n"
   "double tau_0 = 12.5;\n"
   "double tau_r = 130;\n"
   "double J_so = p/tau_r + u*(1 - p)/tau_0;\n"
   "double Iion = (V_fi - V_init)*(J_fi + J_si + J_so);\n"
   "//Do the markov update (1 step rosenbrock with gauss siedel)\n"
   "int _count=0;\n"
   "do\n"
   "{\n"
   "   _count++;\n"
   "} while (_count<50);\n"
   "//EDIT_STATE\n"
   "_state[_ii+v_off*_nCells] += _dt*v_diff;\n"
   "_state[_ii+w_off*_nCells] += _dt*w_diff;\n"
   "_dVm[_indexArray[_ii]] = -Iion;\n"
   "}\n";

   _program_code = ss.str();
   //cout << ss.str();
   nvrtcCreateProgram(&_program,
                      _program_code.c_str(),
                      "fenton_karma_SCC_program",
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
   cuModuleGetFunction(&_kernel, _module, "fenton_karma_SCC_kernel");
}

void ThisReaction::calc(double dt,
                ro_mgarray_ptr<int> indexArray_m,
                ro_mgarray_ptr<double> Vm_m,
                ro_mgarray_ptr<double> iStim_m,
                wo_mgarray_ptr<double> dVm_m)
{
   if (nCells_ == 0) { return; }

   {
      int errorCode=-1;
      if (blockSize_ == -1) { blockSize_ = 1024; }
      while(1)
      {
         const int* indexArrayRaw = indexArray_m.useOn(GPU).raw();
         const double* VmRaw = Vm_m.useOn(GPU).raw();
         const double* iStimRaw = iStim_m.useOn(GPU).raw();
         double* dVmRaw = dVm_m.useOn(GPU).raw();
         double* stateRaw= stateTransport_.readwrite(GPU).raw();
         void* args[] = { &indexArrayRaw,
                          &VmRaw,
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
   _v_off,
   _w_off,
   NUMSTATES
};

ThisReaction::ThisReaction(const int numPoints, const double __dt)
: nCells_(numPoints)
{
   stateTransport_.resize(nCells_*NUMSTATES);
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

ThisReaction::~ThisReaction() {}

void ThisReaction::calc(double _dt,
                ro_mgarray_ptr<int> ___indexArray,
                ro_mgarray_ptr<double> ___Vm,
                ro_mgarray_ptr<double> ,
                wo_mgarray_ptr<double> ___dVm)
{
   ro_array_ptr<int>    __indexArray = ___indexArray.useOn(CPU);
   ro_array_ptr<double> __Vm = ___Vm.useOn(CPU);
   wo_array_ptr<double> __dVm = ___dVm.useOn(CPU);

   //define the constants
   double V_init = -85;
   double V_fi = 15;
   double Cm = 1;
   double u_c = 0.13;
   double u_v = 0;
   double tau_v1_minus = 18.199999999999999;
   double tau_v2_minus = 18.199999999999999;
   double tau_v_plus = 10;
   double g_fi_max = 5.7999999999999998;
   double tau_d = Cm/g_fi_max;
   double tau_w_minus = 80;
   double tau_w_plus = 1020;
   double k = 10;
   double tau_si = 127;
   double u_csi = 0.84999999999999998;
   double tau_0 = 12.5;
   double tau_r = 130;
   for (unsigned __jj=0; __jj<(nCells_+width-1)/width; __jj++)
   {
      const int __ii = __jj*width;
      //set Vm
      double __Vm_local[width];
      {
         int cursor = 0;
         for (int __kk=0; __kk<width; __kk++)
         {
            __Vm_local[__kk] = __Vm[__indexArray[__ii+cursor]];
            if (__ii+__kk < nCells_) { cursor++; }
         }
      }
      const real V = load(&__Vm_local[0]);
      //set all state variables
      real v=load(state_[__jj].v);
      real w=load(state_[__jj].w);
      //get the gate updates (diagonalized exponential integrator)
      //get the other differential updates
      real u = (V - V_init)/(V_fi - V_init);
      real __melodee_temp_000 = u < u_c;
      real __melodee_temp_001;
      if (__melodee_temp_000)
      {
         __melodee_temp_001 = 0;
      }
      else
      {
         __melodee_temp_001 = 1;
      }
      real p = __melodee_temp_001;
      real __melodee_temp_002 = u < u_v;
      real __melodee_temp_003;
      if (__melodee_temp_002)
      {
         __melodee_temp_003 = 0;
      }
      else
      {
         __melodee_temp_003 = 1;
      }
      real q = __melodee_temp_003;
      real tau_v_minus = q*tau_v1_minus + tau_v2_minus*(1 - q);
      real v_diff = -p*v/tau_v_plus + (1 - p)*(1 - v)/tau_v_minus;
      real w_diff = -p*w/tau_w_plus + (1 - p)*(1 - w)/tau_w_minus;
      //get Iion
      real J_fi = -p*v*(1 - u)*(u - u_c)/tau_d;
      real _expensive_functions = tanh(k*(u - u_csi));
      real J_si = -1.0/2.0*w*(_expensive_functions + 1)/tau_si;
      real J_so = p/tau_r + u*(1 - p)/tau_0;
      real Iion = (V_fi - V_init)*(J_fi + J_si + J_so);
      //Do the markov update (1 step rosenbrock with gauss siedel)
      //EDIT_STATE
      v += _dt*v_diff;
      w += _dt*w_diff;
      store(state_[__jj].v, v);
      store(state_[__jj].w, w);
      double __dVm_local[width];
      simdops::store(&__dVm_local[0],-Iion);
      {
         int cursor = 0;
         for (int __kk=0; __kk<width && __ii+__kk<nCells_; __kk++)
         {
            __dVm[__indexArray[__ii+__kk]] = __dVm_local[__kk];
         }
      }
   }
}
#endif //USE_CUDA
   
string ThisReaction::methodName() const
{
   return "fenton_karma_SCC";
}

void ThisReaction::initializeMembraneVoltage(ro_mgarray_ptr<int> __indexArray_m, wo_mgarray_ptr<double> __Vm_m)
{
   assert(__Vm_m.size() >= nCells_);

   wo_array_ptr<double> __Vm = __Vm_m.useOn(CPU);
   ro_array_ptr<int> __indexArray = __indexArray_m.useOn(CPU);
#ifdef USE_CUDA
#define READ_STATE(state,index) (stateData[_##state##_off*nCells_+index])
   wo_array_ptr<double> stateData = stateTransport_.useOn(CPU);
#else //USE_CUDA
#define READ_STATE(state,index) (state_[index/width].state[index % width])
   state_.resize((nCells_+width-1)/width);
#endif //USE_CUDA


   double V_init = -85;
   double V = V_init;
   double v_init = 1;
   double v = v_init;
   double w_init = 1;
   double w = w_init;
   for (int iCell=0; iCell<nCells_; iCell++)
   {
      READ_STATE(v,iCell) = v;
      READ_STATE(w,iCell) = w;
      __Vm[__indexArray[iCell]] = V_init;
   }
}

enum varHandles
{
   v_handle,
   w_handle,
   NUMHANDLES
};

const string ThisReaction::getUnit(const std::string& varName) const
{
   if(0) {}
   else if (varName == "v") { return "1"; }
   else if (varName == "w") { return "1"; }
   return "INVALID";
}

int ThisReaction::getVarHandle(const std::string& varName) const
{
   if (0) {}
   else if (varName == "v") { return v_handle; }
   else if (varName == "w") { return w_handle; }
   return -1;
}

void ThisReaction::setValue(int iCell, int varHandle, double value) 
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readwrite(CPU);
#endif //USE_CUDA



   if (0) {}
   else if (varHandle == v_handle) { READ_STATE(v,iCell) = value; }
   else if (varHandle == w_handle) { READ_STATE(w,iCell) = value; }
}


double ThisReaction::getValue(int iCell, int varHandle) const
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readonly(CPU);
#endif //USE_CUDA


   if (0) {}
   else if (varHandle == v_handle) { return READ_STATE(v,iCell); }
   else if (varHandle == w_handle) { return READ_STATE(w,iCell); }
   return NAN;
}

double ThisReaction::getValue(int iCell, int varHandle, double V) const
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readonly(CPU);
#endif //USE_CUDA


   const double v=READ_STATE(v,iCell);
   const double w=READ_STATE(w,iCell);
   if (0) {}
   else if (varHandle == v_handle)
   {
      return v;
   }
   else if (varHandle == w_handle)
   {
      return w;
   }
   return NAN;
}

void ThisReaction::getCheckpointInfo(vector<string>& fieldNames,
                                     vector<string>& fieldUnits) const
{
   fieldNames.clear();
   fieldUnits.clear();
   fieldNames.push_back("v");
   fieldUnits.push_back(getUnit("v"));
   fieldNames.push_back("w");
   fieldUnits.push_back(getUnit("w"));
}

}
