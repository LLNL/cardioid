/**

   How to convert this code to work for any other model:

   - Search/Replace the model name with your own specific string in the header and source files
   - Add your own code to EDIT_FLAGS and EDIT_PARAMETERS
   - Add your own code to EDIT_PERCELL_FLAGS and EDIT_PERCELL_PARAMETERS
   - Add your own states to EDIT_STATE
   - Add your computation code to the main calc routine, copy pasting frmo matlab.
   
 */


#include "Grandi.hh"
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
   "_d_RLA",
   "_d_RLB",
   "_expensive_functions_045",
   "_expensive_functions_046",
   "_expensive_functions_047",
   "_expensive_functions_048",
   "_expensive_functions_049",
   "_expensive_functions_050",
   "_expensive_functions_051",
   "_expensive_functions_052",
   "_expensive_functions_053",
   "_expensive_functions_054",
   "_expensive_functions_060",
   "_expensive_functions_061",
   "_expensive_functions_062",
   "_expensive_functions_063",
   "_expensive_functions_065",
   "_expensive_functions_067",
   "_f_RLA",
   "_f_RLB",
   "_hL_RLB",
   "_h_RLA",
   "_h_RLB",
   "_j_RLA",
   "_j_RLB",
   "_mL_RLA",
   "_mL_RLB",
   "_m_RLA",
   "_m_RLB",
   "_xkr_RLA",
   "_xkr_RLB",
   "_xks_RLA",
   "_xks_RLB",
   "_xkur_RLA",
   "_xkur_RLB",
   "_xtf_RLA",
   "_xtf_RLB",
   "_ykur_RLA",
   "_ykur_RLB",
   "_ytf_RLA",
   "_ytf_RLB",
   "fnak",
   "kp_kp",
   "rkr",
   "IK1",
    NULL
};


   REACTION_FACTORY(Grandi)(OBJECT* obj, const double _dt, const int numPoints, const ThreadTeam&)
   {
      Grandi::ThisReaction* reaction = new Grandi::ThisReaction(numPoints, _dt);

      //override the defaults
      //EDIT_PARAMETERS
      double AF;
      double ISO;
      double RA;
      setDefault(AF, 1);
      setDefault(RA, 1);
      setDefault(ISO, 1);
      reaction->AF = AF;
      reaction->ISO = ISO;
      reaction->RA = RA;
      bool reusingInterpolants = false;
      string fitName;
      objectGet(obj, "fit", fitName, "");
      int funcCount = sizeof(reaction->_interpolant)/sizeof(reaction->_interpolant[0])-1; //BGQ_HACKFIX, compiler bug with zero length arrays
      if (fitName != "")
      {
         OBJECT* fitObj = objectFind(fitName, "FIT");
         double _fit_dt; objectGet(fitObj, "dt", _fit_dt, "nan");
         double _fit_AF; objectGet(fitObj, "AF", _fit_AF, "nan");
         double _fit_ISO; objectGet(fitObj, "ISO", _fit_ISO, "nan");
         if (1
            && _fit_dt == _dt
            && _fit_AF == reaction->AF
            && _fit_ISO == reaction->ISO
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
            outfile << "   dt = " << _dt << ";\n";
            outfile << "   AF = " << reaction->AF << ";\n";
            outfile << "   ISO = " << reaction->ISO << ";\n";
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

namespace Grandi 
{

void ThisReaction::createInterpolants(const double _dt) {

   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_040 = exp(-0.5*ISO - 0.16666666666666666*v);
         double dss = 1.0/(0.22313016014842982*_expensive_functions_040 + 1.0);
         double _expensive_functions_041 = exp(-0.5*ISO - 0.16666666666666666*v);
         double taud = 1.0*dss*(-0.22313016014842982*_expensive_functions_041 + 1.0)/(0.10500000000000001*ISO + 0.035000000000000003*v + 0.31500000000000006);
         double _expensive_functions_078 = exp(-_dt/taud);
         double _d_RLA = _expensive_functions_078 - 1;
         _outputs[_ii] = _d_RLA;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[0].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _d_RLA: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_040 = exp(-0.5*ISO - 0.16666666666666666*v);
         double dss = 1.0/(0.22313016014842982*_expensive_functions_040 + 1.0);
         double _d_RLB = -dss;
         _outputs[_ii] = _d_RLB;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[1].create(_inputs,_outputs, relError,1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _d_RLB: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double _expensive_functions_045 = exp(2.0*FoRT*v);
         _outputs[_ii] = _expensive_functions_045;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[2].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_045: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double _expensive_functions_046 = exp(2.0*FoRT*v);
         _outputs[_ii] = _expensive_functions_046;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[3].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_046: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double _expensive_functions_047 = exp(2.0*FoRT*v);
         _outputs[_ii] = _expensive_functions_047;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[4].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_047: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double _expensive_functions_048 = exp(2.0*FoRT*v);
         _outputs[_ii] = _expensive_functions_048;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[5].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_048: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double _expensive_functions_049 = exp(FoRT*v);
         _outputs[_ii] = _expensive_functions_049;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[6].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_049: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double _expensive_functions_050 = exp(FoRT*v);
         _outputs[_ii] = _expensive_functions_050;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[7].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_050: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double _expensive_functions_051 = exp(FoRT*v);
         _outputs[_ii] = _expensive_functions_051;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[8].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_051: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double _expensive_functions_052 = exp(FoRT*v);
         _outputs[_ii] = _expensive_functions_052;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[9].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_052: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double _expensive_functions_053 = exp(FoRT*v);
         _outputs[_ii] = _expensive_functions_053;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[10].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_053: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double _expensive_functions_054 = exp(FoRT*v);
         _outputs[_ii] = _expensive_functions_054;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[11].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_054: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double nu = 0.34999999999999998;
         double _expensive_functions_060 = exp(FoRT*nu*v);
         _outputs[_ii] = _expensive_functions_060;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[12].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_060: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double nu = 0.34999999999999998;
         double _expensive_functions_061 = exp(FoRT*nu*v);
         _outputs[_ii] = _expensive_functions_061;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[13].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_061: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double nu = 0.34999999999999998;
         double _expensive_functions_062 = exp(FoRT*v*(nu - 1.0));
         _outputs[_ii] = _expensive_functions_062;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[14].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_062: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double nu = 0.34999999999999998;
         double _expensive_functions_063 = exp(FoRT*v*(nu - 1.0));
         _outputs[_ii] = _expensive_functions_063;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[15].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_063: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double nu = 0.34999999999999998;
         double _expensive_functions_065 = exp(FoRT*v*(nu - 1.0));
         _outputs[_ii] = _expensive_functions_065;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[16].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_065: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double nu = 0.34999999999999998;
         double _expensive_functions_067 = exp(FoRT*v*(nu - 1.0));
         _outputs[_ii] = _expensive_functions_067;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[17].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _expensive_functions_067: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_044 = exp(-((0.1011*ISO + 0.033700000000000001*v + 0.84250000000000003)*(0.1011*ISO + 0.033700000000000001*v + 0.84250000000000003)));
         double tauf = 1.0/(0.019699999999999999*_expensive_functions_044 + 0.02);
         double _expensive_functions_079 = exp(-_dt/tauf);
         double _f_RLA = _expensive_functions_079 - 1;
         _outputs[_ii] = _f_RLA;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[18].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _f_RLA: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_042 = exp(0.42857142857142855*ISO + 0.14285714285714285*v);
         double _expensive_functions_043 = exp(-0.15000000000000002*ISO - 0.050000000000000003*v);
         double fss = 0.20000000000000001/(12.182493960703473*_expensive_functions_043 + 1) + 1.0/(72.654424207165462*_expensive_functions_042 + 1.0);
         double _f_RLB = -fss;
         _outputs[_ii] = _f_RLB;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[19].create(_inputs,_outputs, relError,1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _f_RLB: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_014 = exp(0.16393442622950821*v);
         double hlinf = 1.0/(3011752.7821238632*_expensive_functions_014 + 1.0);
         double _hL_RLB = -hlinf;
         _outputs[_ii] = _hL_RLB;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[20].create(_inputs,_outputs, relError,1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _hL_RLB: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double __melodee_temp_000 = v >= -40;
         double ah;
         double aj;
         double bh;
         double bj;
         if (__melodee_temp_000)
         {
            ah = 0.0;
            double _expensive_functions_010 = exp(-0.0900900900900901*v);
            bh = 0.77000000000000002/(0.049758141083938695*_expensive_functions_010 + 0.13);
         }
         else
         {
            double _expensive_functions_010 = exp(-0.14705882352941177*v);
            ah = 4.4312679295805147e-7*_expensive_functions_010;
            double _expensive_functions_011 = exp(0.34849999999999998*v);
            double _expensive_functions_012 = exp(0.079000000000000001*v);
            bh = 310000.0*_expensive_functions_011 + 2.7000000000000002*_expensive_functions_012;
         }
         double tauh = 1.0/(ah + bh);
         double _expensive_functions_080 = exp(-_dt/tauh);
         double _h_RLA = _expensive_functions_080 - 1;
         _outputs[_ii] = _h_RLA;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[21].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _h_RLA: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_010 = exp(0.13458950201884254*v);
         double hss = 1.0*(1.0/(15212.593285654404*_expensive_functions_010 + 1.0)/(15212.593285654404*_expensive_functions_010 + 1.0));
         double _h_RLB = -hss;
         _outputs[_ii] = _h_RLB;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[22].create(_inputs,_outputs, relError,1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _h_RLB: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double __melodee_temp_000 = v >= -40;
         double ah;
         double aj;
         double bh;
         double bj;
         if (__melodee_temp_000)
         {
            aj = 0.0;
            double _expensive_functions_011 = exp(-0.10000000000000001*v);
            double _expensive_functions_012 = exp(0.057000000000000002*v);
            bj = 0.59999999999999998*_expensive_functions_012/(0.040762203978366204*_expensive_functions_011 + 1.0);
         }
         else
         {
            double _expensive_functions_013 = exp(0.311*v);
            double _expensive_functions_014 = exp(0.24440000000000001*v);
            double _expensive_functions_015 = exp(-0.043909999999999998*v);
            aj = (-25428.0*_expensive_functions_014 - 6.9480000000000002e-6*_expensive_functions_015)*(v + 37.780000000000001)/(50262745825.953987*_expensive_functions_013 + 1.0);
            double _expensive_functions_016 = exp(-0.13780000000000001*v);
            double _expensive_functions_017 = exp(-0.01052*v);
            bj = 0.024240000000000001*_expensive_functions_017/(0.003960868339904256*_expensive_functions_016 + 1.0);
         }
         double tauj = 1.0/(aj + bj);
         double _expensive_functions_082 = exp(-_dt/tauj);
         double _j_RLA = _expensive_functions_082 - 1;
         _outputs[_ii] = _j_RLA;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[23].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _j_RLA: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_011 = exp(0.13458950201884254*v);
         double jss = 1.0*(1.0/(15212.593285654404*_expensive_functions_011 + 1.0)/(15212.593285654404*_expensive_functions_011 + 1.0));
         double _j_RLB = -jss;
         _outputs[_ii] = _j_RLB;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[24].create(_inputs,_outputs, relError,1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _j_RLB: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_012 = exp(-0.10000000000000001*v);
         double aml = (0.32000000000000001*v + 15.081600000000002)/(-0.0089778037306972435*_expensive_functions_012 + 1.0);
         double _expensive_functions_013 = exp(-0.090909090909090912*v);
         double bml = 0.080000000000000002*_expensive_functions_013;
         double _expensive_functions_084 = exp(_dt*(-aml - bml));
         double _mL_RLA = _expensive_functions_084 - 1;
         _outputs[_ii] = _mL_RLA;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[25].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _mL_RLA: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_012 = exp(-0.10000000000000001*v);
         double aml = (0.32000000000000001*v + 15.081600000000002)/(-0.0089778037306972435*_expensive_functions_012 + 1.0);
         double _expensive_functions_013 = exp(-0.090909090909090912*v);
         double bml = 0.080000000000000002*_expensive_functions_013;
         double _mL_RLB = 1.0*aml/(-aml - bml);
         _outputs[_ii] = _mL_RLB;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[26].create(_inputs,_outputs, relError,1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _mL_RLB: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_008 = exp(-((0.064350064350064351*v + 2.9465894465894467)*(0.064350064350064351*v + 2.9465894465894467)));
         double _expensive_functions_009 = exp(-((0.019561815336463225*v - 0.094346635367762138)*(0.019561815336463225*v - 0.094346635367762138)));
         double taum = 0.12920000000000001*_expensive_functions_008 + 0.064869999999999997*_expensive_functions_009;
         double _expensive_functions_083 = exp(-_dt/taum);
         double _m_RLA = _expensive_functions_083 - 1;
         _outputs[_ii] = _m_RLA;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[27].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _m_RLA: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_007 = exp(-0.11074197120708749*v);
         double mss = 1.0*(1.0/(0.0018422115811651339*_expensive_functions_007 + 1.0)/(0.0018422115811651339*_expensive_functions_007 + 1.0));
         double _m_RLB = -mss;
         _outputs[_ii] = _m_RLB;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[28].create(_inputs,_outputs, relError,1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _m_RLB: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_020 = exp(0.050000000000000003*v);
         double _expensive_functions_021 = exp(-0.1111111111111111*v);
         double _expensive_functions_022 = exp(0.1111111111111111*v);
         double tauxr = 3300.0/((0.086774329473929268*_expensive_functions_021 + 1.0)*(3.394723187098903*_expensive_functions_022 + 1.0)) + 230.0/(7.3890560989306504*_expensive_functions_020 + 1.0);
         double _expensive_functions_085 = exp(-_dt/tauxr);
         double _xkr_RLA = _expensive_functions_085 - 1;
         _outputs[_ii] = _xkr_RLA;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[29].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _xkr_RLA: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_019 = exp(-0.20000000000000001*v);
         double xrss = 1.0/(0.1353352832366127*_expensive_functions_019 + 1.0);
         double _xkr_RLB = -xrss;
         _outputs[_ii] = _xkr_RLB;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[30].create(_inputs,_outputs, relError,1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _xkr_RLB: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_025 = exp(-2.8328611898017*ISO - 0.070821529745042494*v);
         double tauxs = 990.10000000000002/(0.84154040886810166*_expensive_functions_025 + 1.0);
         double _expensive_functions_086 = exp(-_dt/tauxs);
         double _xks_RLA = _expensive_functions_086 - 1;
         _outputs[_ii] = _xks_RLA;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[31].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _xks_RLA: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_024 = exp(-2.807017543859649*ISO - 0.070175438596491224*v);
         double xsss = 1.0/(0.76592833836464869*_expensive_functions_024 + 1.0);
         double _xks_RLB = -xsss;
         _outputs[_ii] = _xks_RLB;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[32].create(_inputs,_outputs, relError,1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _xks_RLB: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_032 = exp(0.083333333333333329*v);
         double tauxkur = 0.5 + 9.0/(1.5168967963882134*_expensive_functions_032 + 1.0);
         double _expensive_functions_087 = exp(-_dt/tauxkur);
         double _xkur_RLA = _expensive_functions_087 - 1;
         _outputs[_ii] = _xkur_RLA;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[33].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _xkur_RLA: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_031 = exp(-0.11627906976744186*v);
         double xkurss = 1.0/(0.49774149722499028*_expensive_functions_031 + 1.0);
         double _xkur_RLB = -xkurss;
         _outputs[_ii] = _xkur_RLB;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[34].create(_inputs,_outputs, relError,1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _xkur_RLB: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_028 = exp(-0.0011111111111111111*(v*v));
         double tauxtf = 3.5*_expensive_functions_028 + 1.5;
         double _expensive_functions_088 = exp(-_dt/tauxtf);
         double _xtf_RLA = _expensive_functions_088 - 1;
         _outputs[_ii] = _xtf_RLA;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[35].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _xtf_RLA: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_027 = exp(-0.090909090909090912*v);
         double xtss = 1.0/(0.9131007162822623*_expensive_functions_027 + 1.0);
         double _xtf_RLB = -xtss;
         _outputs[_ii] = _xtf_RLB;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[36].create(_inputs,_outputs, relError,1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _xtf_RLB: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_034 = exp(0.10000000000000001*v);
         double tauykur = 3050.0 + 590.0/(403.42879349273511*_expensive_functions_034 + 1.0);
         double _expensive_functions_089 = exp(-_dt/tauykur);
         double _ykur_RLA = _expensive_functions_089 - 1;
         _outputs[_ii] = _ykur_RLA;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[37].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _ykur_RLA: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_033 = exp(0.10000000000000001*v);
         double ykurss = 1.0/(2.1170000166126748*_expensive_functions_033 + 1.0);
         double _ykur_RLB = -ykurss;
         _outputs[_ii] = _ykur_RLB;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[38].create(_inputs,_outputs, relError,1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _ykur_RLB: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_030 = exp(-((0.062961587135688515*v + 3.3023352452668626)*(0.062961587135688515*v + 3.3023352452668626)));
         double tauytf = 25.635000000000002*_expensive_functions_030 + 24.140000000000001;
         double _expensive_functions_090 = exp(-_dt/tauytf);
         double _ytf_RLA = _expensive_functions_090 - 1;
         _outputs[_ii] = _ytf_RLA;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[39].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _ytf_RLA: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_029 = exp(0.086956521739130432*v);
         double ytss = 1.0/(33.843235113007339*_expensive_functions_029 + 1.0);
         double _ytf_RLB = -ytss;
         _outputs[_ii] = _ytf_RLB;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[40].create(_inputs,_outputs, relError,1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for _ytf_RLB: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double R = 8314.0;
         double Frdy = 96485.0;
         double Temp = 310.0;
         double FoRT = Frdy/(R*Temp);
         double Nao = 140.0;
         double _expensive_functions_015 = exp(0.01485884101040119*Nao);
         double sigma = 0.14285714285714285*_expensive_functions_015 - 0.14285714285714285;
         double _expensive_functions_016 = exp(-0.10000000000000001*FoRT*v);
         double _expensive_functions_017 = exp(-FoRT*v);
         double fnak = 1.0/(0.1245*_expensive_functions_016 + 0.036499999999999998*_expensive_functions_017*sigma + 1.0);
         _outputs[_ii] = fnak;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[41].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for fnak: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_026 = exp(-0.16722408026755853*v);
         double kp_kp = 1.0/(1786.4755653786237*_expensive_functions_026 + 1.0);
         _outputs[_ii] = kp_kp;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[42].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for kp_kp: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double v = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = v;
         double _expensive_functions_023 = exp(0.041666666666666664*v);
         double rkr = 1.0/(21.831051418620834*_expensive_functions_023 + 1.0);
         _outputs[_ii] = rkr;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[43].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for rkr: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
   {
      int _numPoints = (100 - -100)/1e-2;
      vector<double> _inputs(_numPoints);
      vector<double> _outputs(_numPoints);
      for (int _ii=0; _ii<_numPoints; _ii++)
      {
         double vek = -100 + (100 - -100)*(_ii+0.5)/_numPoints;
         _inputs[_ii] = vek;
         double Ko = 5.4000000000000004;
         double _expensive_functions_035 = exp(0.23849999999999999*vek);
         double aki = 1.02/(7.3545425104644605e-7*_expensive_functions_035 + 1);
         double _expensive_functions_036 = exp(-0.51429999999999998*vek);
         double _expensive_functions_037 = exp(0.080320000000000003*vek);
         double _expensive_functions_038 = exp(0.061749999999999999*vek);
         double bki = (0.76262400650630813*_expensive_functions_037 + 1.1534056351865558e-16*_expensive_functions_038)/(0.086772294157693317*_expensive_functions_036 + 1.0);
         double kiss = aki/(aki + bki);
         double _expensive_functions_039 = sqrt(Ko);
         double IK1 = 0.43033148291193518*_expensive_functions_039*kiss*vek*(0.052499999999999998*AF + 0.052499999999999998);
         _outputs[_ii] = IK1;
      }
      double relError = 0.0001;
      double actualTolerance = _interpolant[44].create(_inputs,_outputs, relError,0.1);
      if (actualTolerance > relError  && getRank(0) == 0)
      {
         cerr << "Warning: Could not meet tolerance for IK1: " 
              << actualTolerance << " > " << relError
              << " target" << endl;
      }
   }
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

   "   CaM_off,\n"
   "   Cai_off,\n"
   "   Caj_off,\n"
   "   Casl_off,\n"
   "   Casr_off,\n"
   "   Ki_off,\n"
   "   Myc_off,\n"
   "   Mym_off,\n"
   "   NaBj_off,\n"
   "   NaBsl_off,\n"
   "   Nai_off,\n"
   "   Naj_off,\n"
   "   Nasl_off,\n"
   "   RyRi_off,\n"
   "   RyRo_off,\n"
   "   RyRr_off,\n"
   "   SLHj_off,\n"
   "   SLHsl_off,\n"
   "   SLLj_off,\n"
   "   SLLsl_off,\n"
   "   SRB_off,\n"
   "   TnCHc_off,\n"
   "   TnCHm_off,\n"
   "   TnCL_off,\n"
   "   d_off,\n"
   "   f_off,\n"
   "   fcaBj_off,\n"
   "   fcaBsl_off,\n"
   "   h_off,\n"
   "   hL_off,\n"
   "   j_off,\n"
   "   m_off,\n"
   "   mL_off,\n"
   "   xkr_off,\n"
   "   xks_off,\n"
   "   xkur_off,\n"
   "   xtf_off,\n"
   "   ykur_off,\n"
   "   ytf_off,\n"
   "   NUMSTATES\n"
   "};\n"
   "extern \"C\"\n"
   "__global__ void Grandi_kernel(const double* _Vm, const double* _iStim, double* _dVm, double* _state) {\n"
   "const double _dt = " << __cachedDt << ";\n"
   "const int _nCells = " << nCells_ << ";\n"

   "const double AF = " << AF << ";\n"
   "const double ISO = " << ISO << ";\n"
   "const double RA = " << RA << ";\n"
   "const int _ii = threadIdx.x + blockIdx.x*blockDim.x;\n"
   "if (_ii >= _nCells) { return; }\n"
   "const double V = _Vm[_ii];\n"
   "double _ratPoly;\n"

   "const double CaM = _state[_ii+CaM_off*_nCells];\n"
   "const double Cai = _state[_ii+Cai_off*_nCells];\n"
   "const double Caj = _state[_ii+Caj_off*_nCells];\n"
   "const double Casl = _state[_ii+Casl_off*_nCells];\n"
   "const double Casr = _state[_ii+Casr_off*_nCells];\n"
   "const double Ki = _state[_ii+Ki_off*_nCells];\n"
   "const double Myc = _state[_ii+Myc_off*_nCells];\n"
   "const double Mym = _state[_ii+Mym_off*_nCells];\n"
   "const double NaBj = _state[_ii+NaBj_off*_nCells];\n"
   "const double NaBsl = _state[_ii+NaBsl_off*_nCells];\n"
   "const double Nai = _state[_ii+Nai_off*_nCells];\n"
   "const double Naj = _state[_ii+Naj_off*_nCells];\n"
   "const double Nasl = _state[_ii+Nasl_off*_nCells];\n"
   "const double RyRi = _state[_ii+RyRi_off*_nCells];\n"
   "const double RyRo = _state[_ii+RyRo_off*_nCells];\n"
   "const double RyRr = _state[_ii+RyRr_off*_nCells];\n"
   "const double SLHj = _state[_ii+SLHj_off*_nCells];\n"
   "const double SLHsl = _state[_ii+SLHsl_off*_nCells];\n"
   "const double SLLj = _state[_ii+SLLj_off*_nCells];\n"
   "const double SLLsl = _state[_ii+SLLsl_off*_nCells];\n"
   "const double SRB = _state[_ii+SRB_off*_nCells];\n"
   "const double TnCHc = _state[_ii+TnCHc_off*_nCells];\n"
   "const double TnCHm = _state[_ii+TnCHm_off*_nCells];\n"
   "const double TnCL = _state[_ii+TnCL_off*_nCells];\n"
   "const double d = _state[_ii+d_off*_nCells];\n"
   "const double f = _state[_ii+f_off*_nCells];\n"
   "const double fcaBj = _state[_ii+fcaBj_off*_nCells];\n"
   "const double fcaBsl = _state[_ii+fcaBsl_off*_nCells];\n"
   "const double h = _state[_ii+h_off*_nCells];\n"
   "const double hL = _state[_ii+hL_off*_nCells];\n"
   "const double j = _state[_ii+j_off*_nCells];\n"
   "const double m = _state[_ii+m_off*_nCells];\n"
   "const double mL = _state[_ii+mL_off*_nCells];\n"
   "const double xkr = _state[_ii+xkr_off*_nCells];\n"
   "const double xks = _state[_ii+xks_off*_nCells];\n"
   "const double xkur = _state[_ii+xkur_off*_nCells];\n"
   "const double xtf = _state[_ii+xtf_off*_nCells];\n"
   "const double ykur = _state[_ii+ykur_off*_nCells];\n"
   "const double ytf = _state[_ii+ytf_off*_nCells];\n"
   "//get the gate updates (diagonalized exponential integrator)\n"
   "double v = V;\n"
   "double tauhl = 600.0;\n"
   ""; generateInterpString(ss,_interpolant[0], "v"); ss << "\n"
   "double _d_RLA = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[1], "v"); ss << "\n"
   "double _d_RLB = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[18], "v"); ss << "\n"
   "double _f_RLA = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[19], "v"); ss << "\n"
   "double _f_RLB = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[21], "v"); ss << "\n"
   "double _h_RLA = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[22], "v"); ss << "\n"
   "double _h_RLB = _ratPoly;\n"
   "double _expensive_functions_081 = exp(-_dt/tauhl);\n"
   "double _hL_RLA = _expensive_functions_081 - 1;\n"
   ""; generateInterpString(ss,_interpolant[20], "v"); ss << "\n"
   "double _hL_RLB = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[23], "v"); ss << "\n"
   "double _j_RLA = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[24], "v"); ss << "\n"
   "double _j_RLB = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[27], "v"); ss << "\n"
   "double _m_RLA = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[28], "v"); ss << "\n"
   "double _m_RLB = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[25], "v"); ss << "\n"
   "double _mL_RLA = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[26], "v"); ss << "\n"
   "double _mL_RLB = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[29], "v"); ss << "\n"
   "double _xkr_RLA = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[30], "v"); ss << "\n"
   "double _xkr_RLB = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[31], "v"); ss << "\n"
   "double _xks_RLA = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[32], "v"); ss << "\n"
   "double _xks_RLB = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[33], "v"); ss << "\n"
   "double _xkur_RLA = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[34], "v"); ss << "\n"
   "double _xkur_RLB = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[35], "v"); ss << "\n"
   "double _xtf_RLA = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[36], "v"); ss << "\n"
   "double _xtf_RLB = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[37], "v"); ss << "\n"
   "double _ykur_RLA = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[38], "v"); ss << "\n"
   "double _ykur_RLB = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[39], "v"); ss << "\n"
   "double _ytf_RLA = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[40], "v"); ss << "\n"
   "double _ytf_RLB = _ratPoly;\n"
   "//get the other differential updates\n"
   "double R = 8314.0;\n"
   "double Frdy = 96485.0;\n"
   "double Temp = 310.0;\n"
   "double FoRT = Frdy/(R*Temp);\n"
   "double Cmem = 1.0999999999999999e-10;\n"
   "double Qpow = 0.10000000000000001*Temp - 31.0;\n"
   "double pi = 3.1415926535897932;\n"
   "double cellLength = 100.0;\n"
   "double cellRadius = 10.25;\n"
   "double Vcell = 1.0000000000000001e-15*(cellRadius*cellRadius)*cellLength*pi;\n"
   "double Vmyo = 0.65000000000000002*Vcell;\n"
   "double Vsr = 0.035000000000000003*Vcell;\n"
   "double Vsl = 0.02*Vcell;\n"
   "double Vjunc = 0.00053900000000000009*Vcell;\n"
   "double Jca_juncsl = 8.2413054227789685e-13;\n"
   "double Jca_slmyo = 3.7242560798480505e-12;\n"
   "double Jna_juncsl = 1.8312782322060799e-14;\n"
   "double Jna_slmyo = 1.6386279222197947e-12;\n"
   "double Fjunc = 0.11;\n"
   "double Fsl = -Fjunc + 1;\n"
   "double Fjunc_CaL = 0.90000000000000002;\n"
   "double Fsl_CaL = -Fjunc_CaL + 1;\n"
   "double Ko = 5.4000000000000004;\n"
   "double Nao = 140.0;\n"
   "double Cao = 1.8;\n"
   "double Mgi = 1.0;\n"
   "double _expensive_functions = log(Nao/Naj);\n"
   "double ENa_junc = 1.0*_expensive_functions/FoRT;\n"
   "double _expensive_functions_001 = log(Nao/Nasl);\n"
   "double ENa_sl = 1.0*_expensive_functions_001/FoRT;\n"
   "double _expensive_functions_002 = log(Ko/Ki);\n"
   "double EK = 1.0*_expensive_functions_002/FoRT;\n"
   "double pNaK = 0.018329999999999999;\n"
   "double _expensive_functions_003 = log((Ko + Nao*pNaK)/(Ki + Nai*pNaK));\n"
   "double EKs = 1.0*_expensive_functions_003/FoRT;\n"
   "double _expensive_functions_004 = log(Cao/Caj);\n"
   "double ECa_junc = 0.5*_expensive_functions_004/FoRT;\n"
   "double _expensive_functions_005 = log(Cao/Casl);\n"
   "double ECa_sl = 0.5*_expensive_functions_005/FoRT;\n"
   "double vek = -EK + v;\n"
   "double veks = -EKs + v;\n"
   "double GNa = -2.3000000000000003*AF + 23.0;\n"
   "double INa_junc = (m*m*m)*Fjunc*GNa*h*j*(-ENa_junc + v);\n"
   "double INa_sl = (m*m*m)*Fsl*GNa*h*j*(-ENa_sl + v);\n"
   "double GNaL = -0.0025000000000000001*AF + 0.0025000000000000001;\n"
   "double INaL_junc = (mL*mL*mL)*Fjunc*GNaL*hL*(-ENa_junc + v);\n"
   "double INaL_sl = (mL*mL*mL)*Fsl*GNaL*hL*(-ENa_sl + v);\n"
   "double GNaB = 0.00059699999999999998;\n"
   "double INaBk_junc = Fjunc*GNaB*(-ENa_junc + v);\n"
   "double INaBk_sl = Fsl*GNaB*(-ENa_sl + v);\n"
   "double KmNaip = -2.75*ISO + 11.0;\n"
   "double KmKo = 1.5;\n"
   ""; generateInterpString(ss,_interpolant[41], "v"); ss << "\n"
   "double fnak = _ratPoly;\n"
   "double IbarNaK = 1.26;\n"
   "double INAK_junc = Fjunc*IbarNaK*Ko*fnak/((KmKo + Ko)*(((KmNaip/Naj)*(KmNaip/Naj)*(KmNaip/Naj)*(KmNaip/Naj)) + 1.0));\n"
   "double INAK_sl = Fsl*IbarNaK*Ko*fnak/((KmKo + Ko)*(((KmNaip/Nasl)*(KmNaip/Nasl)*(KmNaip/Nasl)*(KmNaip/Nasl)) + 1.0));\n"
   "double INAK = INAK_junc + INAK_sl;\n"
   "double _expensive_functions_018 = sqrt(Ko);\n"
   "double gkr = 0.015061601901917732*_expensive_functions_018;\n"
   ""; generateInterpString(ss,_interpolant[43], "v"); ss << "\n"
   "double rkr = _ratPoly;\n"
   "double IKr = gkr*rkr*vek*xkr;\n"
   "double gks_junc = 0.0035000000000000001*AF + 0.0070000000000000001*ISO + 0.0035000000000000001;\n"
   "double gks_sl = 0.0035000000000000001*AF + 0.0070000000000000001*ISO + 0.0035000000000000001;\n"
   "double IKs_junc = (xks*xks)*Fjunc*gks_junc*veks;\n"
   "double IKs_sl = (xks*xks)*Fsl*gks_sl*veks;\n"
   "double IKs = IKs_junc + IKs_sl;\n"
   ""; generateInterpString(ss,_interpolant[42], "v"); ss << "\n"
   "double kp_kp = _ratPoly;\n"
   "double gkp = 0.002;\n"
   "double IKp_junc = Fjunc*gkp*kp_kp*vek;\n"
   "double IKp_sl = Fsl*gkp*kp_kp*vek;\n"
   "double IKp = IKp_junc + IKp_sl;\n"
   "double GtoFast = -0.11549999999999999*AF + 0.16500000000000001;\n"
   "double Ito = GtoFast*vek*xtf*ytf;\n"
   "double Gkur = 0.044999999999999998*(-0.5*AF + 1.0)*(2.0*ISO + 1.0)*(0.20000000000000001*RA + 1.0);\n"
   "double IKur = Gkur*vek*xkur*ykur;\n"
   ""; generateInterpString(ss,_interpolant[44], "vek"); ss << "\n"
   "double IK1 = _ratPoly;\n"
   "double fcaBj_diff = 1.7*Caj*(-fcaBj + 1.0) - 0.011900000000000001*fcaBj;\n"
   "double fcaBsl_diff = 1.7*Casl*(-fcaBsl + 1) - 0.011900000000000001*fcaBsl;\n"
   "double pNa = 7.4999999999999993e-9*(-0.5*AF + 1.0)*(0.5*ISO + 1.0);\n"
   "double pCa = 0.00027*(-0.5*AF + 1.0)*(0.5*ISO + 1.0);\n"
   "double pK = 1.35e-7*(-0.5*AF + 1.0)*(0.5*ISO + 1.0);\n"
   "double Q10CaL = 1.8;\n"
   ""; generateInterpString(ss,_interpolant[2], "v"); ss << "\n"
   "double _expensive_functions_045 = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[3], "v"); ss << "\n"
   "double _expensive_functions_046 = _ratPoly;\n"
   "double ibarca_j = 4.0*FoRT*Frdy*pCa*v*(0.34100000000000003*Caj*_expensive_functions_046 - 0.34100000000000003*Cao)/(_expensive_functions_045 - 1.0);\n"
   ""; generateInterpString(ss,_interpolant[4], "v"); ss << "\n"
   "double _expensive_functions_047 = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[5], "v"); ss << "\n"
   "double _expensive_functions_048 = _ratPoly;\n"
   "double ibarca_sl = 4.0*FoRT*Frdy*pCa*v*(-0.34100000000000003*Cao + 0.34100000000000003*Casl*_expensive_functions_048)/(_expensive_functions_047 - 1.0);\n"
   ""; generateInterpString(ss,_interpolant[6], "v"); ss << "\n"
   "double _expensive_functions_049 = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[7], "v"); ss << "\n"
   "double _expensive_functions_050 = _ratPoly;\n"
   "double ibark = FoRT*Frdy*pK*v*(0.75*Ki*_expensive_functions_050 - 0.75*Ko)/(_expensive_functions_049 - 1.0);\n"
   ""; generateInterpString(ss,_interpolant[8], "v"); ss << "\n"
   "double _expensive_functions_051 = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[9], "v"); ss << "\n"
   "double _expensive_functions_052 = _ratPoly;\n"
   "double ibarna_j = FoRT*Frdy*pNa*v*(0.75*Naj*_expensive_functions_052 - 0.75*Nao)/(_expensive_functions_051 - 1.0);\n"
   ""; generateInterpString(ss,_interpolant[10], "v"); ss << "\n"
   "double _expensive_functions_053 = _ratPoly;\n"
   ""; generateInterpString(ss,_interpolant[11], "v"); ss << "\n"
   "double _expensive_functions_054 = _ratPoly;\n"
   "double ibarna_sl = FoRT*Frdy*pNa*v*(-0.75*Nao + 0.75*Nasl*_expensive_functions_054)/(_expensive_functions_053 - 1.0);\n"
   "double _expensive_functions_055 = pow(Q10CaL, Qpow);\n"
   "double ICa_junc = 0.45000000000000001*Fjunc_CaL*_expensive_functions_055*d*f*ibarca_j*(-fcaBj + 1);\n"
   "double _expensive_functions_056 = pow(Q10CaL, Qpow);\n"
   "double ICa_sl = 0.45000000000000001*Fsl_CaL*_expensive_functions_056*d*f*ibarca_sl*(-fcaBsl + 1);\n"
   "double _expensive_functions_057 = pow(Q10CaL, Qpow);\n"
   "double ICaNa_junc = 0.45000000000000001*Fjunc_CaL*_expensive_functions_057*d*f*ibarna_j*(-fcaBj + 1);\n"
   "double _expensive_functions_058 = pow(Q10CaL, Qpow);\n"
   "double ICaNa_sl = 0.45000000000000001*Fsl_CaL*_expensive_functions_058*d*f*ibarna_sl*(-fcaBsl + 1);\n"
   "double _expensive_functions_059 = pow(Q10CaL, Qpow);\n"
   "double ICaK = 0.45000000000000001*_expensive_functions_059*d*f*ibark*(Fjunc_CaL*(-fcaBj + 1) + Fsl_CaL*(-fcaBsl + 1));\n"
   "double KmCai = 0.0035899999999999999;\n"
   "double KmCao = 1.3;\n"
   "double KmNai = 12.289999999999999;\n"
   "double KmNao = 87.5;\n"
   "double ksat = 0.27000000000000002;\n"
   "double Kdact = 0.00038400000000000001;\n"
   "double Q10NCX = 1.5700000000000001;\n"
   "double Ka_junc = 1.0/(((Kdact/Caj)*(Kdact/Caj)) + 1.0);\n"
   "double Ka_sl = 1.0/(((Kdact/Casl)*(Kdact/Casl)) + 1.0);\n"
   ""; generateInterpString(ss,_interpolant[12], "v"); ss << "\n"
   "double _expensive_functions_060 = _ratPoly;\n"
   "double s1_junc = (Naj*Naj*Naj)*Cao*_expensive_functions_060;\n"
   ""; generateInterpString(ss,_interpolant[13], "v"); ss << "\n"
   "double _expensive_functions_061 = _ratPoly;\n"
   "double s1_sl = (Nasl*Nasl*Nasl)*Cao*_expensive_functions_061;\n"
   ""; generateInterpString(ss,_interpolant[14], "v"); ss << "\n"
   "double _expensive_functions_062 = _ratPoly;\n"
   "double s2_junc = (Nao*Nao*Nao)*Caj*_expensive_functions_062;\n"
   "double s3_junc = (KmNao*KmNao*KmNao)*Caj*(Caj/KmCai + 1.0) + (Nao*Nao*Nao)*Caj + (Naj*Naj*Naj)*Cao + (Nao*Nao*Nao)*KmCai*(((Naj/KmNai)*(Naj/KmNai)*(Naj/KmNai)) + 1) + (Naj*Naj*Naj)*KmCao;\n"
   ""; generateInterpString(ss,_interpolant[15], "v"); ss << "\n"
   "double _expensive_functions_063 = _ratPoly;\n"
   "double s2_sl = (Nao*Nao*Nao)*Casl*_expensive_functions_063;\n"
   "double s3_sl = (Nasl*Nasl*Nasl)*Cao + (KmNao*KmNao*KmNao)*Casl*(Casl/KmCai + 1.0) + (Nao*Nao*Nao)*Casl + (Nao*Nao*Nao)*KmCai*(((Nasl/KmNai)*(Nasl/KmNai)*(Nasl/KmNai)) + 1.0) + (Nasl*Nasl*Nasl)*KmCao;\n"
   "double IbarNCX = 1.26*AF + 3.1499999999999999;\n"
   "double _expensive_functions_064 = pow(Q10NCX, Qpow);\n"
   ""; generateInterpString(ss,_interpolant[16], "v"); ss << "\n"
   "double _expensive_functions_065 = _ratPoly;\n"
   "double Incx_junc = Fjunc*IbarNCX*Ka_junc*_expensive_functions_064*(s1_junc - s2_junc)/(s3_junc*(_expensive_functions_065*ksat + 1.0));\n"
   "double _expensive_functions_066 = pow(Q10NCX, Qpow);\n"
   ""; generateInterpString(ss,_interpolant[17], "v"); ss << "\n"
   "double _expensive_functions_067 = _ratPoly;\n"
   "double Incx_sl = Fsl*IbarNCX*Ka_sl*_expensive_functions_066*(s1_sl - s2_sl)/(s3_sl*(_expensive_functions_067*ksat + 1.0));\n"
   "double Q10SLCaP = 2.3500000000000001;\n"
   "double KmPCa = 0.00050000000000000001;\n"
   "double IbarSLCaP = 0.047100000000000003;\n"
   "double Caj_pow = pow(Caj, 1.6000000000000001);\n"
   "double Casl_pow = pow(Casl, 1.6000000000000001);\n"
   "double _expensive_functions_068 = pow(Q10SLCaP, Qpow);\n"
   "double _expensive_functions_069 = pow(KmPCa, 1.6000000000000001);\n"
   "double Ipca_junc = Caj_pow*Fjunc*IbarSLCaP*_expensive_functions_068/(Caj_pow + _expensive_functions_069);\n"
   "double _expensive_functions_070 = pow(Q10SLCaP, Qpow);\n"
   "double _expensive_functions_071 = pow(KmPCa, 1.6000000000000001);\n"
   "double Ipca_sl = Casl_pow*Fsl*IbarSLCaP*_expensive_functions_070/(Casl_pow + _expensive_functions_071);\n"
   "double GCaB = 0.00060643000000000003;\n"
   "double Icabk_junc = Fjunc*GCaB*(-ECa_junc + v);\n"
   "double Icabk_sl = Fsl*GCaB*(-ECa_sl + v);\n"
   "double Q10SRCaP = 2.6000000000000001;\n"
   "double Vmax_SRCaP = -0.0026557*AF + 0.0053114;\n"
   "double Kmf = -0.00030750000000000005*ISO + 0.0006150000000000001;\n"
   "double Kmr = 1.7;\n"
   "double hillSRCaP = 1.7869999999999999;\n"
   "double ks = 25.0;\n"
   "double koCa = 20.0*AF + 10.0*ISO*(-AF + 1.0) + 10.0;\n"
   "double kom = 0.059999999999999998;\n"
   "double kiCa = 0.5;\n"
   "double kim = 0.0050000000000000001;\n"
   "double ec50SR = 0.45000000000000001;\n"
   "double MaxSR = 15.0;\n"
   "double MinSR = 1.0;\n"
   "double _expensive_functions_072 = pow(ec50SR/Casr, 2.5);\n"
   "double kCaSR = MaxSR - (MaxSR - MinSR)/(_expensive_functions_072 + 1.0);\n"
   "double koSRCa = koCa/kCaSR;\n"
   "double kiSRCa = kCaSR*kiCa;\n"
   "double RI = -RyRi - RyRo - RyRr + 1.0;\n"
   "double RyRr_diff = -(Caj*Caj)*RyRr*koSRCa - Caj*RyRr*kiSRCa + RI*kim + RyRo*kom;\n"
   "double RyRo_diff = (Caj*Caj)*RyRr*koSRCa - Caj*RyRo*kiSRCa + RyRi*kim - RyRo*kom;\n"
   "double RyRi_diff = (Caj*Caj)*RI*koSRCa + Caj*RyRo*kiSRCa - RyRi*kim - RyRi*kom;\n"
   "double JSRCarel = RyRo*ks*(-Caj + Casr);\n"
   "double _expensive_functions_073 = pow(Q10SRCaP, Qpow);\n"
   "double _expensive_functions_074 = pow(Cai/Kmf, hillSRCaP);\n"
   "double _expensive_functions_075 = pow(Casr/Kmr, hillSRCaP);\n"
   "double _expensive_functions_076 = pow(Cai/Kmf, hillSRCaP);\n"
   "double _expensive_functions_077 = pow(Casr/Kmr, hillSRCaP);\n"
   "double Jserca = 1.0*Vmax_SRCaP*_expensive_functions_073*(_expensive_functions_076 - _expensive_functions_077)/(_expensive_functions_074 + _expensive_functions_075 + 1);\n"
   "double JSRleak = (1.3370000000000001e-6*AF + 5.3480000000000003e-6)*(-Caj + Casr);\n"
   "double Bmax_Naj = 7.5609999999999999;\n"
   "double Bmax_Nasl = 1.6499999999999999;\n"
   "double koff_na = 0.001;\n"
   "double kon_na = 0.0001;\n"
   "double NaBj_diff = -NaBj*koff_na + Naj*kon_na*(Bmax_Naj - NaBj);\n"
   "double NaBsl_diff = -NaBsl*koff_na + Nasl*kon_na*(Bmax_Nasl - NaBsl);\n"
   "double Bmax_TnClow = 0.070000000000000007;\n"
   "double koff_tncl = 0.0097999999999999997*ISO + 0.019599999999999999;\n"
   "double kon_tncl = 32.700000000000003;\n"
   "double Bmax_TnChigh = 0.14000000000000001;\n"
   "double koff_tnchca = 3.1999999999999999e-5;\n"
   "double kon_tnchca = 2.3700000000000001;\n"
   "double koff_tnchmg = 0.0033300000000000001;\n"
   "double kon_tnchmg = 0.0030000000000000001;\n"
   "double Bmax_CaM = 0.024;\n"
   "double koff_cam = 0.23799999999999999;\n"
   "double kon_cam = 34.0;\n"
   "double Bmax_myosin = 0.14000000000000001;\n"
   "double koff_myoca = 0.00046000000000000001;\n"
   "double kon_myoca = 13.800000000000001;\n"
   "double koff_myomg = 5.7000000000000003e-5;\n"
   "double kon_myomg = 0.015699999999999999;\n"
   "double Bmax_SR = 0.017100000000000001;\n"
   "double koff_sr = 0.059999999999999998;\n"
   "double kon_sr = 100.0;\n"
   "double TnCL_diff = Cai*kon_tncl*(Bmax_TnClow - TnCL) - TnCL*koff_tncl;\n"
   "double TnCHc_diff = Cai*kon_tnchca*(Bmax_TnChigh - TnCHc - TnCHm) - TnCHc*koff_tnchca;\n"
   "double TnCHm_diff = Mgi*kon_tnchmg*(Bmax_TnChigh - TnCHc - TnCHm) - TnCHm*koff_tnchmg;\n"
   "double CaM_diff = -CaM*koff_cam + Cai*kon_cam*(Bmax_CaM - CaM);\n"
   "double Myc_diff = Cai*kon_myoca*(Bmax_myosin - Myc - Mym) - Myc*koff_myoca;\n"
   "double Mym_diff = Mgi*kon_myomg*(Bmax_myosin - Myc - Mym) - Mym*koff_myomg;\n"
   "double SRB_diff = Cai*kon_sr*(Bmax_SR - SRB) - SRB*koff_sr;\n"
   "double JCaB_cytsol = CaM_diff + Myc_diff + Mym_diff + SRB_diff + TnCHc_diff + TnCHm_diff + TnCL_diff;\n"
   "double Bmax_SLlowsl = 0.037400000000000003*Vmyo/Vsl;\n"
   "double Bmax_SLlowj = 0.00046000000000000001*Vmyo/Vjunc;\n"
   "double koff_sll = 1.3;\n"
   "double kon_sll = 100.0;\n"
   "double Bmax_SLhighsl = 0.0134*Vmyo/Vsl;\n"
   "double Bmax_SLhighj = 0.000165*Vmyo/Vjunc;\n"
   "double koff_slh = 0.029999999999999999;\n"
   "double kon_slh = 100.0;\n"
   "double SLLj_diff = Caj*kon_sll*(Bmax_SLlowj - SLLj) - SLLj*koff_sll;\n"
   "double SLLsl_diff = Casl*kon_sll*(Bmax_SLlowsl - SLLsl) - SLLsl*koff_sll;\n"
   "double SLHj_diff = Caj*kon_slh*(Bmax_SLhighj - SLHj) - SLHj*koff_slh;\n"
   "double SLHsl_diff = Casl*kon_slh*(Bmax_SLhighsl - SLHsl) - SLHsl*koff_slh;\n"
   "double JCaB_junc = SLHj_diff + SLLj_diff;\n"
   "double JCaB_sl = SLHsl_diff + SLLsl_diff;\n"
   "double Bmax_Csqn = 0.14000000000000001*Vmyo/Vsr;\n"
   "double koff_csqn = 65.0;\n"
   "double kon_csqn = 100.0;\n"
   "double KmCsqnb = koff_csqn/kon_csqn;\n"
   "double Buff_Csqnb = 1.0/(Bmax_Csqn*KmCsqnb*(1.0/(Casr + KmCsqnb)/(Casr + KmCsqnb)) + 1.0);\n"
   "double Casr_diff = Buff_Csqnb*(-JSRCarel - JSRleak*Vmyo/Vsr + Jserca);\n"
   "double INa_tot_junc = ICaNa_junc + 3.0*INAK_junc + INaBk_junc + INaL_junc + INa_junc + 3.0*Incx_junc;\n"
   "double INa_tot_sl = ICaNa_sl + 3.0*INAK_sl + INaBk_sl + INaL_sl + INa_sl + 3.0*Incx_sl;\n"
   "double Naj_diff = -Cmem*INa_tot_junc/(Frdy*Vjunc) + Jna_juncsl*(-Naj + Nasl)/Vjunc - NaBj_diff;\n"
   "double Nasl_diff = -Cmem*INa_tot_sl/(Frdy*Vsl) + Jna_juncsl*(Naj - Nasl)/Vsl + Jna_slmyo*(Nai - Nasl)/Vsl - NaBsl_diff;\n"
   "double Nai_diff = Jna_slmyo*(-Nai + Nasl)/Vmyo;\n"
   "double IK_tot = ICaK + IK1 + IKp + IKr + IKs + IKur - 2.0*INAK + Ito;\n"
   "double Ki_diff = -Cmem*IK_tot/(Frdy*Vmyo);\n"
   "double ICa_tot_junc = ICa_junc + Icabk_junc - 2.0*Incx_junc + Ipca_junc;\n"
   "double ICa_tot_sl = ICa_sl + Icabk_sl - 2.0*Incx_sl + Ipca_sl;\n"
   "double Caj_diff = -0.5*Cmem*ICa_tot_junc/(Frdy*Vjunc) - JCaB_junc + JSRCarel*Vsr/Vjunc + JSRleak*Vmyo/Vjunc + Jca_juncsl*(-Caj + Casl)/Vjunc;\n"
   "double Casl_diff = -0.5*Cmem*ICa_tot_sl/(Frdy*Vsl) - JCaB_sl + Jca_juncsl*(Caj - Casl)/Vsl + Jca_slmyo*(Cai - Casl)/Vsl;\n"
   "double Cai_diff = -JCaB_cytsol + Jca_slmyo*(-Cai + Casl)/Vmyo - Jserca*Vsr/Vmyo;\n"
   "//get Iion\n"
   "double Cli = 15.0;\n"
   "double Clo = 150.0;\n"
   "double _expensive_functions_006 = log(Cli/Clo);\n"
   "double ECl = 1.0*_expensive_functions_006/FoRT;\n"
   "double KdClCa = 0.10000000000000001;\n"
   "double GClCa = 0.054800000000000001;\n"
   "double IClCa_junc = Fjunc*GClCa*(-ECl + v)/(1.0 + KdClCa/Caj);\n"
   "double IClCa_sl = Fsl*GClCa*(-ECl + v)/(1.0 + KdClCa/Casl);\n"
   "double IClCa = IClCa_junc + IClCa_sl;\n"
   "double GClB = 0.0089999999999999993;\n"
   "double IClbk = GClB*(-ECl + v);\n"
   "double INa_tot = INa_tot_junc + INa_tot_sl;\n"
   "double ICl_tot = IClCa + IClbk;\n"
   "double ICa_tot = ICa_tot_junc + ICa_tot_sl;\n"
   "double Iion = ICa_tot + ICl_tot + IK_tot + INa_tot;\n"
   "//Do the markov update (1 step rosenbrock with gauss siedel)\n"
   "double _mi_new_RyRi = RyRi_diff;\n"
   "double _mi_new_RyRo = RyRo_diff;\n"
   "double _mi_new_RyRr = RyRr_diff;\n"
   "int _count=0;\n"
   "do\n"
   "{\n"
   "   double _mi_old_RyRi = _mi_new_RyRi;\n"
   "   double _mi_old_RyRo = _mi_new_RyRo;\n"
   "   double _mi_old_RyRr = _mi_new_RyRr;\n"
   "   double _d_RI_wrt_RyRi = -1;\n"
   "   double _d_RI_wrt_RyRo = -1;\n"
   "   double _d_RI_wrt_RyRr = -1;\n"
   "   _mi_new_RyRi = (RyRi_diff + _dt*((Caj*Caj)*_d_RI_wrt_RyRr*_mi_old_RyRr*koSRCa + _mi_old_RyRo*((Caj*Caj)*_d_RI_wrt_RyRo*koSRCa + Caj*kiSRCa)))/(-_dt*((Caj*Caj)*_d_RI_wrt_RyRi*koSRCa - kim - kom) + 1);\n"
   "   _mi_new_RyRo = (RyRo_diff + _dt*((Caj*Caj)*_mi_old_RyRr*koSRCa + _mi_new_RyRi*kim))/(-_dt*(-Caj*kiSRCa - kom) + 1);\n"
   "   _mi_new_RyRr = (RyRr_diff + _dt*(_d_RI_wrt_RyRi*_mi_new_RyRi*kim + _mi_new_RyRo*(_d_RI_wrt_RyRo*kim + kom)))/(-_dt*(-(Caj*Caj)*koSRCa - Caj*kiSRCa + _d_RI_wrt_RyRr*kim) + 1);\n"
   "   _count++;\n"
   "} while (_count<50);\n"
   "//EDIT_STATE\n"
   "_state[_ii+CaM_off*_nCells] += _dt*CaM_diff;\n"
   "_state[_ii+Cai_off*_nCells] += _dt*Cai_diff;\n"
   "_state[_ii+Caj_off*_nCells] += _dt*Caj_diff;\n"
   "_state[_ii+Casl_off*_nCells] += _dt*Casl_diff;\n"
   "_state[_ii+Casr_off*_nCells] += _dt*Casr_diff;\n"
   "_state[_ii+Ki_off*_nCells] += _dt*Ki_diff;\n"
   "_state[_ii+Myc_off*_nCells] += _dt*Myc_diff;\n"
   "_state[_ii+Mym_off*_nCells] += _dt*Mym_diff;\n"
   "_state[_ii+NaBj_off*_nCells] += _dt*NaBj_diff;\n"
   "_state[_ii+NaBsl_off*_nCells] += _dt*NaBsl_diff;\n"
   "_state[_ii+Nai_off*_nCells] += _dt*Nai_diff;\n"
   "_state[_ii+Naj_off*_nCells] += _dt*Naj_diff;\n"
   "_state[_ii+Nasl_off*_nCells] += _dt*Nasl_diff;\n"
   "_state[_ii+SLHj_off*_nCells] += _dt*SLHj_diff;\n"
   "_state[_ii+SLHsl_off*_nCells] += _dt*SLHsl_diff;\n"
   "_state[_ii+SLLj_off*_nCells] += _dt*SLLj_diff;\n"
   "_state[_ii+SLLsl_off*_nCells] += _dt*SLLsl_diff;\n"
   "_state[_ii+SRB_off*_nCells] += _dt*SRB_diff;\n"
   "_state[_ii+TnCHc_off*_nCells] += _dt*TnCHc_diff;\n"
   "_state[_ii+TnCHm_off*_nCells] += _dt*TnCHm_diff;\n"
   "_state[_ii+TnCL_off*_nCells] += _dt*TnCL_diff;\n"
   "_state[_ii+fcaBj_off*_nCells] += _dt*fcaBj_diff;\n"
   "_state[_ii+fcaBsl_off*_nCells] += _dt*fcaBsl_diff;\n"
   "_state[_ii+d_off*_nCells] += _d_RLA*(d+_d_RLB);\n"
   "_state[_ii+f_off*_nCells] += _f_RLA*(f+_f_RLB);\n"
   "_state[_ii+h_off*_nCells] += _h_RLA*(h+_h_RLB);\n"
   "_state[_ii+hL_off*_nCells] += _hL_RLA*(hL+_hL_RLB);\n"
   "_state[_ii+j_off*_nCells] += _j_RLA*(j+_j_RLB);\n"
   "_state[_ii+m_off*_nCells] += _m_RLA*(m+_m_RLB);\n"
   "_state[_ii+mL_off*_nCells] += _mL_RLA*(mL+_mL_RLB);\n"
   "_state[_ii+xkr_off*_nCells] += _xkr_RLA*(xkr+_xkr_RLB);\n"
   "_state[_ii+xks_off*_nCells] += _xks_RLA*(xks+_xks_RLB);\n"
   "_state[_ii+xkur_off*_nCells] += _xkur_RLA*(xkur+_xkur_RLB);\n"
   "_state[_ii+xtf_off*_nCells] += _xtf_RLA*(xtf+_xtf_RLB);\n"
   "_state[_ii+ykur_off*_nCells] += _ykur_RLA*(ykur+_ykur_RLB);\n"
   "_state[_ii+ytf_off*_nCells] += _ytf_RLA*(ytf+_ytf_RLB);\n"
   "_state[_ii+RyRi_off*_nCells] += _dt*_mi_new_RyRi;\n"
   "_state[_ii+RyRo_off*_nCells] += _dt*_mi_new_RyRo;\n"
   "_state[_ii+RyRr_off*_nCells] += _dt*_mi_new_RyRr;\n"
   "_dVm[_ii] = -Iion;\n"
   "}\n";

   _program_code = ss.str();
   //cout << ss.str();
   nvrtcCreateProgram(&_program,
                      _program_code.c_str(),
                      "Grandi_program",
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
   cuModuleGetFunction(&_kernel, _module, "Grandi_kernel");
}

void ThisReaction::calc(double dt,
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
         const double* VmRaw = Vm_m.useOn(GPU).raw();
         const double* iStimRaw = iStim_m.useOn(GPU).raw();
         double* dVmRaw = dVm_m.useOn(GPU).raw();
         double* stateRaw= stateTransport_.readwrite(GPU).raw();
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
   _CaM_off,
   _Cai_off,
   _Caj_off,
   _Casl_off,
   _Casr_off,
   _Ki_off,
   _Myc_off,
   _Mym_off,
   _NaBj_off,
   _NaBsl_off,
   _Nai_off,
   _Naj_off,
   _Nasl_off,
   _RyRi_off,
   _RyRo_off,
   _RyRr_off,
   _SLHj_off,
   _SLHsl_off,
   _SLLj_off,
   _SLLsl_off,
   _SRB_off,
   _TnCHc_off,
   _TnCHm_off,
   _TnCL_off,
   _d_off,
   _f_off,
   _fcaBj_off,
   _fcaBsl_off,
   _h_off,
   _hL_off,
   _j_off,
   _m_off,
   _mL_off,
   _xkr_off,
   _xks_off,
   _xkur_off,
   _xtf_off,
   _ykur_off,
   _ytf_off,
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
                ro_mgarray_ptr<double> ___Vm,
                ro_mgarray_ptr<double> ___iStim,
                wo_mgarray_ptr<double> ___dVm)
{
   ro_array_ptr<double> __Vm = ___Vm.useOn(CPU);
   ro_array_ptr<double> __iStim = ___iStim.useOn(CPU);
   wo_array_ptr<double> __dVm = ___dVm.useOn(CPU);

   //define the constants
   double R = 8314.0;
   double Frdy = 96485.0;
   double Temp = 310.0;
   double FoRT = Frdy/(R*Temp);
   double Cmem = 1.0999999999999999e-10;
   double Qpow = 0.10000000000000001*Temp - 31.0;
   double pi = 3.1415926535897932;
   double cellLength = 100.0;
   double cellRadius = 10.25;
   double Vcell = 1.0000000000000001e-15*(cellRadius*cellRadius)*cellLength*pi;
   double Vmyo = 0.65000000000000002*Vcell;
   double Vsr = 0.035000000000000003*Vcell;
   double Vsl = 0.02*Vcell;
   double Vjunc = 0.00053900000000000009*Vcell;
   double Jca_juncsl = 8.2413054227789685e-13;
   double Jca_slmyo = 3.7242560798480505e-12;
   double Jna_juncsl = 1.8312782322060799e-14;
   double Jna_slmyo = 1.6386279222197947e-12;
   double Fjunc = 0.11;
   double Fsl = -Fjunc + 1;
   double Fjunc_CaL = 0.90000000000000002;
   double Fsl_CaL = -Fjunc_CaL + 1;
   double Cli = 15.0;
   double Clo = 150.0;
   double Ko = 5.4000000000000004;
   double Nao = 140.0;
   double Cao = 1.8;
   double Mgi = 1.0;
   double pNaK = 0.018329999999999999;
   double _expensive_functions_006 = log(Cli/Clo);
   double ECl = 1.0*_expensive_functions_006/FoRT;
   double GNa = -2.3000000000000003*AF + 23.0;
   double tauhl = 600.0;
   double GNaL = -0.0025000000000000001*AF + 0.0025000000000000001;
   double GNaB = 0.00059699999999999998;
   double KmNaip = -2.75*ISO + 11.0;
   double KmKo = 1.5;
   double IbarNaK = 1.26;
   double _expensive_functions_018 = sqrt(Ko);
   double gkr = 0.015061601901917732*_expensive_functions_018;
   double gks_junc = 0.0035000000000000001*AF + 0.0070000000000000001*ISO + 0.0035000000000000001;
   double gks_sl = 0.0035000000000000001*AF + 0.0070000000000000001*ISO + 0.0035000000000000001;
   double gkp = 0.002;
   double GtoFast = -0.11549999999999999*AF + 0.16500000000000001;
   double Gkur = 0.044999999999999998*(-0.5*AF + 1.0)*(2.0*ISO + 1.0)*(0.20000000000000001*RA + 1.0);
   double KdClCa = 0.10000000000000001;
   double GClCa = 0.054800000000000001;
   double GClB = 0.0089999999999999993;
   double pNa = 7.4999999999999993e-9*(-0.5*AF + 1.0)*(0.5*ISO + 1.0);
   double pCa = 0.00027*(-0.5*AF + 1.0)*(0.5*ISO + 1.0);
   double pK = 1.35e-7*(-0.5*AF + 1.0)*(0.5*ISO + 1.0);
   double Q10CaL = 1.8;
   double _expensive_functions_055 = pow(Q10CaL, Qpow);
   double _expensive_functions_056 = pow(Q10CaL, Qpow);
   double _expensive_functions_057 = pow(Q10CaL, Qpow);
   double _expensive_functions_058 = pow(Q10CaL, Qpow);
   double _expensive_functions_059 = pow(Q10CaL, Qpow);
   double KmCai = 0.0035899999999999999;
   double KmCao = 1.3;
   double KmNai = 12.289999999999999;
   double KmNao = 87.5;
   double ksat = 0.27000000000000002;
   double Kdact = 0.00038400000000000001;
   double Q10NCX = 1.5700000000000001;
   double IbarNCX = 1.26*AF + 3.1499999999999999;
   double _expensive_functions_064 = pow(Q10NCX, Qpow);
   double _expensive_functions_066 = pow(Q10NCX, Qpow);
   double Q10SLCaP = 2.3500000000000001;
   double KmPCa = 0.00050000000000000001;
   double IbarSLCaP = 0.047100000000000003;
   double _expensive_functions_068 = pow(Q10SLCaP, Qpow);
   double _expensive_functions_069 = pow(KmPCa, 1.6000000000000001);
   double _expensive_functions_070 = pow(Q10SLCaP, Qpow);
   double _expensive_functions_071 = pow(KmPCa, 1.6000000000000001);
   double GCaB = 0.00060643000000000003;
   double Q10SRCaP = 2.6000000000000001;
   double Vmax_SRCaP = -0.0026557*AF + 0.0053114;
   double Kmf = -0.00030750000000000005*ISO + 0.0006150000000000001;
   double Kmr = 1.7;
   double hillSRCaP = 1.7869999999999999;
   double ks = 25.0;
   double koCa = 20.0*AF + 10.0*ISO*(-AF + 1.0) + 10.0;
   double kom = 0.059999999999999998;
   double kiCa = 0.5;
   double kim = 0.0050000000000000001;
   double ec50SR = 0.45000000000000001;
   double MaxSR = 15.0;
   double MinSR = 1.0;
   double _d_RI_wrt_RyRi = -1;
   double _d_RI_wrt_RyRo = -1;
   double _d_RI_wrt_RyRr = -1;
   double _expensive_functions_073 = pow(Q10SRCaP, Qpow);
   double Bmax_Naj = 7.5609999999999999;
   double Bmax_Nasl = 1.6499999999999999;
   double koff_na = 0.001;
   double kon_na = 0.0001;
   double Bmax_TnClow = 0.070000000000000007;
   double koff_tncl = 0.0097999999999999997*ISO + 0.019599999999999999;
   double kon_tncl = 32.700000000000003;
   double Bmax_TnChigh = 0.14000000000000001;
   double koff_tnchca = 3.1999999999999999e-5;
   double kon_tnchca = 2.3700000000000001;
   double koff_tnchmg = 0.0033300000000000001;
   double kon_tnchmg = 0.0030000000000000001;
   double Bmax_CaM = 0.024;
   double koff_cam = 0.23799999999999999;
   double kon_cam = 34.0;
   double Bmax_myosin = 0.14000000000000001;
   double koff_myoca = 0.00046000000000000001;
   double kon_myoca = 13.800000000000001;
   double koff_myomg = 5.7000000000000003e-5;
   double kon_myomg = 0.015699999999999999;
   double Bmax_SR = 0.017100000000000001;
   double koff_sr = 0.059999999999999998;
   double kon_sr = 100.0;
   double Bmax_SLlowsl = 0.037400000000000003*Vmyo/Vsl;
   double Bmax_SLlowj = 0.00046000000000000001*Vmyo/Vjunc;
   double koff_sll = 1.3;
   double kon_sll = 100.0;
   double Bmax_SLhighsl = 0.0134*Vmyo/Vsl;
   double Bmax_SLhighj = 0.000165*Vmyo/Vjunc;
   double koff_slh = 0.029999999999999999;
   double kon_slh = 100.0;
   double Bmax_Csqn = 0.14000000000000001*Vmyo/Vsr;
   double koff_csqn = 65.0;
   double kon_csqn = 100.0;
   double KmCsqnb = koff_csqn/kon_csqn;
   double _expensive_functions_081 = exp(-_dt/tauhl);
   double _hL_RLA = _expensive_functions_081 - 1;
   for (unsigned __jj=0; __jj<(nCells_+width-1)/width; __jj++)
   {
      const int __ii = __jj*width;
      //set Vm
      const real V = load(&__Vm[__ii]);
      const real iStim = load(&__iStim[__ii]);

      //set all state variables
      real CaM=load(state_[__jj].CaM);
      real Cai=load(state_[__jj].Cai);
      real Caj=load(state_[__jj].Caj);
      real Casl=load(state_[__jj].Casl);
      real Casr=load(state_[__jj].Casr);
      real Ki=load(state_[__jj].Ki);
      real Myc=load(state_[__jj].Myc);
      real Mym=load(state_[__jj].Mym);
      real NaBj=load(state_[__jj].NaBj);
      real NaBsl=load(state_[__jj].NaBsl);
      real Nai=load(state_[__jj].Nai);
      real Naj=load(state_[__jj].Naj);
      real Nasl=load(state_[__jj].Nasl);
      real RyRi=load(state_[__jj].RyRi);
      real RyRo=load(state_[__jj].RyRo);
      real RyRr=load(state_[__jj].RyRr);
      real SLHj=load(state_[__jj].SLHj);
      real SLHsl=load(state_[__jj].SLHsl);
      real SLLj=load(state_[__jj].SLLj);
      real SLLsl=load(state_[__jj].SLLsl);
      real SRB=load(state_[__jj].SRB);
      real TnCHc=load(state_[__jj].TnCHc);
      real TnCHm=load(state_[__jj].TnCHm);
      real TnCL=load(state_[__jj].TnCL);
      real d=load(state_[__jj].d);
      real f=load(state_[__jj].f);
      real fcaBj=load(state_[__jj].fcaBj);
      real fcaBsl=load(state_[__jj].fcaBsl);
      real h=load(state_[__jj].h);
      real hL=load(state_[__jj].hL);
      real j=load(state_[__jj].j);
      real m=load(state_[__jj].m);
      real mL=load(state_[__jj].mL);
      real xkr=load(state_[__jj].xkr);
      real xks=load(state_[__jj].xks);
      real xkur=load(state_[__jj].xkur);
      real xtf=load(state_[__jj].xtf);
      real ykur=load(state_[__jj].ykur);
      real ytf=load(state_[__jj].ytf);
      //get the gate updates (diagonalized exponential integrator)
      real v = V;
      real _d_RLA = _interpolant[0].eval(v);
      real _d_RLB = _interpolant[1].eval(v);
      real _f_RLA = _interpolant[18].eval(v);
      real _f_RLB = _interpolant[19].eval(v);
      real _h_RLA = _interpolant[21].eval(v);
      real _h_RLB = _interpolant[22].eval(v);
      real _hL_RLB = _interpolant[20].eval(v);
      real _j_RLA = _interpolant[23].eval(v);
      real _j_RLB = _interpolant[24].eval(v);
      real _m_RLA = _interpolant[27].eval(v);
      real _m_RLB = _interpolant[28].eval(v);
      real _mL_RLA = _interpolant[25].eval(v);
      real _mL_RLB = _interpolant[26].eval(v);
      real _xkr_RLA = _interpolant[29].eval(v);
      real _xkr_RLB = _interpolant[30].eval(v);
      real _xks_RLA = _interpolant[31].eval(v);
      real _xks_RLB = _interpolant[32].eval(v);
      real _xkur_RLA = _interpolant[33].eval(v);
      real _xkur_RLB = _interpolant[34].eval(v);
      real _xtf_RLA = _interpolant[35].eval(v);
      real _xtf_RLB = _interpolant[36].eval(v);
      real _ykur_RLA = _interpolant[37].eval(v);
      real _ykur_RLB = _interpolant[38].eval(v);
      real _ytf_RLA = _interpolant[39].eval(v);
      real _ytf_RLB = _interpolant[40].eval(v);
      //get the other differential updates
      real _expensive_functions = log(Nao/Naj);
      real ENa_junc = 1.0*_expensive_functions/FoRT;
      real _expensive_functions_001 = log(Nao/Nasl);
      real ENa_sl = 1.0*_expensive_functions_001/FoRT;
      real _expensive_functions_002 = log(Ko/Ki);
      real EK = 1.0*_expensive_functions_002/FoRT;
      real _expensive_functions_003 = log((Ko + Nao*pNaK)/(Ki + Nai*pNaK));
      real EKs = 1.0*_expensive_functions_003/FoRT;
      real _expensive_functions_004 = log(Cao/Caj);
      real ECa_junc = 0.5*_expensive_functions_004/FoRT;
      real _expensive_functions_005 = log(Cao/Casl);
      real ECa_sl = 0.5*_expensive_functions_005/FoRT;
      real vek = -EK + v;
      real veks = -EKs + v;
      real INa_junc = (m*m*m)*Fjunc*GNa*h*j*(-ENa_junc + v);
      real INa_sl = (m*m*m)*Fsl*GNa*h*j*(-ENa_sl + v);
      real INaL_junc = (mL*mL*mL)*Fjunc*GNaL*hL*(-ENa_junc + v);
      real INaL_sl = (mL*mL*mL)*Fsl*GNaL*hL*(-ENa_sl + v);
      real INaBk_junc = Fjunc*GNaB*(-ENa_junc + v);
      real INaBk_sl = Fsl*GNaB*(-ENa_sl + v);
      real fnak = _interpolant[41].eval(v);
      real INAK_junc = Fjunc*IbarNaK*Ko*fnak/((KmKo + Ko)*(((KmNaip/Naj)*(KmNaip/Naj)*(KmNaip/Naj)*(KmNaip/Naj)) + 1.0));
      real INAK_sl = Fsl*IbarNaK*Ko*fnak/((KmKo + Ko)*(((KmNaip/Nasl)*(KmNaip/Nasl)*(KmNaip/Nasl)*(KmNaip/Nasl)) + 1.0));
      real INAK = INAK_junc + INAK_sl;
      real rkr = _interpolant[43].eval(v);
      real IKr = gkr*rkr*vek*xkr;
      real IKs_junc = (xks*xks)*Fjunc*gks_junc*veks;
      real IKs_sl = (xks*xks)*Fsl*gks_sl*veks;
      real IKs = IKs_junc + IKs_sl;
      real kp_kp = _interpolant[42].eval(v);
      real IKp_junc = Fjunc*gkp*kp_kp*vek;
      real IKp_sl = Fsl*gkp*kp_kp*vek;
      real IKp = IKp_junc + IKp_sl;
      real Ito = GtoFast*vek*xtf*ytf;
      real IKur = Gkur*vek*xkur*ykur;
      real IK1 = _interpolant[44].eval(vek);
      real fcaBj_diff = 1.7*Caj*(-fcaBj + 1.0) - 0.011900000000000001*fcaBj;
      real fcaBsl_diff = 1.7*Casl*(-fcaBsl + 1) - 0.011900000000000001*fcaBsl;
      real _expensive_functions_045 = _interpolant[2].eval(v);
      real _expensive_functions_046 = _interpolant[3].eval(v);
      real ibarca_j = 4.0*FoRT*Frdy*pCa*v*(0.34100000000000003*Caj*_expensive_functions_046 - 0.34100000000000003*Cao)/(_expensive_functions_045 - 1.0);
      real _expensive_functions_047 = _interpolant[4].eval(v);
      real _expensive_functions_048 = _interpolant[5].eval(v);
      real ibarca_sl = 4.0*FoRT*Frdy*pCa*v*(-0.34100000000000003*Cao + 0.34100000000000003*Casl*_expensive_functions_048)/(_expensive_functions_047 - 1.0);
      real _expensive_functions_049 = _interpolant[6].eval(v);
      real _expensive_functions_050 = _interpolant[7].eval(v);
      real ibark = FoRT*Frdy*pK*v*(0.75*Ki*_expensive_functions_050 - 0.75*Ko)/(_expensive_functions_049 - 1.0);
      real _expensive_functions_051 = _interpolant[8].eval(v);
      real _expensive_functions_052 = _interpolant[9].eval(v);
      real ibarna_j = FoRT*Frdy*pNa*v*(0.75*Naj*_expensive_functions_052 - 0.75*Nao)/(_expensive_functions_051 - 1.0);
      real _expensive_functions_053 = _interpolant[10].eval(v);
      real _expensive_functions_054 = _interpolant[11].eval(v);
      real ibarna_sl = FoRT*Frdy*pNa*v*(-0.75*Nao + 0.75*Nasl*_expensive_functions_054)/(_expensive_functions_053 - 1.0);
      real ICa_junc = 0.45000000000000001*Fjunc_CaL*_expensive_functions_055*d*f*ibarca_j*(-fcaBj + 1);
      real ICa_sl = 0.45000000000000001*Fsl_CaL*_expensive_functions_056*d*f*ibarca_sl*(-fcaBsl + 1);
      real ICaNa_junc = 0.45000000000000001*Fjunc_CaL*_expensive_functions_057*d*f*ibarna_j*(-fcaBj + 1);
      real ICaNa_sl = 0.45000000000000001*Fsl_CaL*_expensive_functions_058*d*f*ibarna_sl*(-fcaBsl + 1);
      real ICaK = 0.45000000000000001*_expensive_functions_059*d*f*ibark*(Fjunc_CaL*(-fcaBj + 1) + Fsl_CaL*(-fcaBsl + 1));
      real Ka_junc = 1.0/(((Kdact/Caj)*(Kdact/Caj)) + 1.0);
      real Ka_sl = 1.0/(((Kdact/Casl)*(Kdact/Casl)) + 1.0);
      real _expensive_functions_060 = _interpolant[12].eval(v);
      real s1_junc = (Naj*Naj*Naj)*Cao*_expensive_functions_060;
      real _expensive_functions_061 = _interpolant[13].eval(v);
      real s1_sl = (Nasl*Nasl*Nasl)*Cao*_expensive_functions_061;
      real _expensive_functions_062 = _interpolant[14].eval(v);
      real s2_junc = (Nao*Nao*Nao)*Caj*_expensive_functions_062;
      real s3_junc = (KmNao*KmNao*KmNao)*Caj*(Caj/KmCai + 1.0) + (Nao*Nao*Nao)*Caj + (Naj*Naj*Naj)*Cao + (Nao*Nao*Nao)*KmCai*(((Naj/KmNai)*(Naj/KmNai)*(Naj/KmNai)) + 1) + (Naj*Naj*Naj)*KmCao;
      real _expensive_functions_063 = _interpolant[15].eval(v);
      real s2_sl = (Nao*Nao*Nao)*Casl*_expensive_functions_063;
      real s3_sl = (Nasl*Nasl*Nasl)*Cao + (KmNao*KmNao*KmNao)*Casl*(Casl/KmCai + 1.0) + (Nao*Nao*Nao)*Casl + (Nao*Nao*Nao)*KmCai*(((Nasl/KmNai)*(Nasl/KmNai)*(Nasl/KmNai)) + 1.0) + (Nasl*Nasl*Nasl)*KmCao;
      real _expensive_functions_065 = _interpolant[16].eval(v);
      real Incx_junc = Fjunc*IbarNCX*Ka_junc*_expensive_functions_064*(s1_junc - s2_junc)/(s3_junc*(_expensive_functions_065*ksat + 1.0));
      real _expensive_functions_067 = _interpolant[17].eval(v);
      real Incx_sl = Fsl*IbarNCX*Ka_sl*_expensive_functions_066*(s1_sl - s2_sl)/(s3_sl*(_expensive_functions_067*ksat + 1.0));
      real Caj_pow = pow(Caj, 1.6000000000000001);
      real Casl_pow = pow(Casl, 1.6000000000000001);
      real Ipca_junc = Caj_pow*Fjunc*IbarSLCaP*_expensive_functions_068/(Caj_pow + _expensive_functions_069);
      real Ipca_sl = Casl_pow*Fsl*IbarSLCaP*_expensive_functions_070/(Casl_pow + _expensive_functions_071);
      real Icabk_junc = Fjunc*GCaB*(-ECa_junc + v);
      real Icabk_sl = Fsl*GCaB*(-ECa_sl + v);
      real _expensive_functions_072 = pow(ec50SR/Casr, 2.5);
      real kCaSR = MaxSR - (MaxSR - MinSR)/(_expensive_functions_072 + 1.0);
      real koSRCa = koCa/kCaSR;
      real kiSRCa = kCaSR*kiCa;
      real RI = -RyRi - RyRo - RyRr + 1.0;
      real RyRr_diff = -(Caj*Caj)*RyRr*koSRCa - Caj*RyRr*kiSRCa + RI*kim + RyRo*kom;
      real RyRo_diff = (Caj*Caj)*RyRr*koSRCa - Caj*RyRo*kiSRCa + RyRi*kim - RyRo*kom;
      real RyRi_diff = (Caj*Caj)*RI*koSRCa + Caj*RyRo*kiSRCa - RyRi*kim - RyRi*kom;
      real JSRCarel = RyRo*ks*(-Caj + Casr);
      real _expensive_functions_074 = pow(Cai/Kmf, hillSRCaP);
      real _expensive_functions_075 = pow(Casr/Kmr, hillSRCaP);
      real _expensive_functions_076 = pow(Cai/Kmf, hillSRCaP);
      real _expensive_functions_077 = pow(Casr/Kmr, hillSRCaP);
      real Jserca = 1.0*Vmax_SRCaP*_expensive_functions_073*(_expensive_functions_076 - _expensive_functions_077)/(_expensive_functions_074 + _expensive_functions_075 + 1);
      real JSRleak = (1.3370000000000001e-6*AF + 5.3480000000000003e-6)*(-Caj + Casr);
      real NaBj_diff = -NaBj*koff_na + Naj*kon_na*(Bmax_Naj - NaBj);
      real NaBsl_diff = -NaBsl*koff_na + Nasl*kon_na*(Bmax_Nasl - NaBsl);
      real TnCL_diff = Cai*kon_tncl*(Bmax_TnClow - TnCL) - TnCL*koff_tncl;
      real TnCHc_diff = Cai*kon_tnchca*(Bmax_TnChigh - TnCHc - TnCHm) - TnCHc*koff_tnchca;
      real TnCHm_diff = Mgi*kon_tnchmg*(Bmax_TnChigh - TnCHc - TnCHm) - TnCHm*koff_tnchmg;
      real CaM_diff = -CaM*koff_cam + Cai*kon_cam*(Bmax_CaM - CaM);
      real Myc_diff = Cai*kon_myoca*(Bmax_myosin - Myc - Mym) - Myc*koff_myoca;
      real Mym_diff = Mgi*kon_myomg*(Bmax_myosin - Myc - Mym) - Mym*koff_myomg;
      real SRB_diff = Cai*kon_sr*(Bmax_SR - SRB) - SRB*koff_sr;
      real JCaB_cytsol = CaM_diff + Myc_diff + Mym_diff + SRB_diff + TnCHc_diff + TnCHm_diff + TnCL_diff;
      real SLLj_diff = Caj*kon_sll*(Bmax_SLlowj - SLLj) - SLLj*koff_sll;
      real SLLsl_diff = Casl*kon_sll*(Bmax_SLlowsl - SLLsl) - SLLsl*koff_sll;
      real SLHj_diff = Caj*kon_slh*(Bmax_SLhighj - SLHj) - SLHj*koff_slh;
      real SLHsl_diff = Casl*kon_slh*(Bmax_SLhighsl - SLHsl) - SLHsl*koff_slh;
      real JCaB_junc = SLHj_diff + SLLj_diff;
      real JCaB_sl = SLHsl_diff + SLLsl_diff;
      real Buff_Csqnb = 1.0/(Bmax_Csqn*KmCsqnb*(1.0/(Casr + KmCsqnb)/(Casr + KmCsqnb)) + 1.0);
      real Casr_diff = Buff_Csqnb*(-JSRCarel - JSRleak*Vmyo/Vsr + Jserca);
      real INa_tot_junc = ICaNa_junc + 3.0*INAK_junc + INaBk_junc + INaL_junc + INa_junc + 3.0*Incx_junc;
      real INa_tot_sl = ICaNa_sl + 3.0*INAK_sl + INaBk_sl + INaL_sl + INa_sl + 3.0*Incx_sl;
      real Naj_diff = -Cmem*INa_tot_junc/(Frdy*Vjunc) + Jna_juncsl*(-Naj + Nasl)/Vjunc - NaBj_diff;
      real Nasl_diff = -Cmem*INa_tot_sl/(Frdy*Vsl) + Jna_juncsl*(Naj - Nasl)/Vsl + Jna_slmyo*(Nai - Nasl)/Vsl - NaBsl_diff;
      real Nai_diff = Jna_slmyo*(-Nai + Nasl)/Vmyo;
      real IK_tot = ICaK + IK1 + IKp + IKr + IKs + IKur - 2.0*INAK + Ito;
      real Ki_diff = -Cmem*IK_tot/(Frdy*Vmyo);
      real ICa_tot_junc = ICa_junc + Icabk_junc - 2.0*Incx_junc + Ipca_junc;
      real ICa_tot_sl = ICa_sl + Icabk_sl - 2.0*Incx_sl + Ipca_sl;
      real Caj_diff = -0.5*Cmem*ICa_tot_junc/(Frdy*Vjunc) - JCaB_junc + JSRCarel*Vsr/Vjunc + JSRleak*Vmyo/Vjunc + Jca_juncsl*(-Caj + Casl)/Vjunc;
      real Casl_diff = -0.5*Cmem*ICa_tot_sl/(Frdy*Vsl) - JCaB_sl + Jca_juncsl*(Caj - Casl)/Vsl + Jca_slmyo*(Cai - Casl)/Vsl;
      real Cai_diff = -JCaB_cytsol + Jca_slmyo*(-Cai + Casl)/Vmyo - Jserca*Vsr/Vmyo;
      //get Iion
      real IClCa_junc = Fjunc*GClCa*(-ECl + v)/(1.0 + KdClCa/Caj);
      real IClCa_sl = Fsl*GClCa*(-ECl + v)/(1.0 + KdClCa/Casl);
      real IClCa = IClCa_junc + IClCa_sl;
      real IClbk = GClB*(-ECl + v);
      real INa_tot = INa_tot_junc + INa_tot_sl;
      real ICl_tot = IClCa + IClbk;
      real ICa_tot = ICa_tot_junc + ICa_tot_sl;
      real Iion = ICa_tot + ICl_tot + IK_tot + INa_tot;
      //Do the markov update (1 step rosenbrock with gauss siedel)
      real _mi_new_RyRi = RyRi_diff;
      real _mi_new_RyRo = RyRo_diff;
      real _mi_new_RyRr = RyRr_diff;
      int _count=0;
      real _error;
      do
      {
         real _mi_old_RyRi = _mi_new_RyRi;
         real _mi_old_RyRo = _mi_new_RyRo;
         real _mi_old_RyRr = _mi_new_RyRr;
         _mi_new_RyRi = (RyRi_diff + _dt*((Caj*Caj)*_d_RI_wrt_RyRr*_mi_old_RyRr*koSRCa + _mi_old_RyRo*((Caj*Caj)*_d_RI_wrt_RyRo*koSRCa + Caj*kiSRCa)))/(-_dt*((Caj*Caj)*_d_RI_wrt_RyRi*koSRCa - kim - kom) + 1);
         _mi_new_RyRo = (RyRo_diff + _dt*((Caj*Caj)*_mi_old_RyRr*koSRCa + _mi_new_RyRi*kim))/(-_dt*(-Caj*kiSRCa - kom) + 1);
         _mi_new_RyRr = (RyRr_diff + _dt*(_d_RI_wrt_RyRi*_mi_new_RyRi*kim + _mi_new_RyRo*(_d_RI_wrt_RyRo*kim + kom)))/(-_dt*(-(Caj*Caj)*koSRCa - Caj*kiSRCa + _d_RI_wrt_RyRr*kim) + 1);
         _error = 0;
         _error += (_mi_old_RyRi-_mi_new_RyRi)*(_mi_old_RyRi-_mi_new_RyRi);
         _error += (_mi_old_RyRo-_mi_new_RyRo)*(_mi_old_RyRo-_mi_new_RyRo);
         _error += (_mi_old_RyRr-_mi_new_RyRr)*(_mi_old_RyRr-_mi_new_RyRr);
         _count++;
      } while (simdops::any(_error > 1e-100) && _count<50);
      //EDIT_STATE
      CaM += _dt*CaM_diff;
      Cai += _dt*Cai_diff;
      Caj += _dt*Caj_diff;
      Casl += _dt*Casl_diff;
      Casr += _dt*Casr_diff;
      Ki += _dt*Ki_diff;
      Myc += _dt*Myc_diff;
      Mym += _dt*Mym_diff;
      NaBj += _dt*NaBj_diff;
      NaBsl += _dt*NaBsl_diff;
      Nai += _dt*Nai_diff;
      Naj += _dt*Naj_diff;
      Nasl += _dt*Nasl_diff;
      SLHj += _dt*SLHj_diff;
      SLHsl += _dt*SLHsl_diff;
      SLLj += _dt*SLLj_diff;
      SLLsl += _dt*SLLsl_diff;
      SRB += _dt*SRB_diff;
      TnCHc += _dt*TnCHc_diff;
      TnCHm += _dt*TnCHm_diff;
      TnCL += _dt*TnCL_diff;
      fcaBj += _dt*fcaBj_diff;
      fcaBsl += _dt*fcaBsl_diff;
      d += _d_RLA*(d+_d_RLB);
      f += _f_RLA*(f+_f_RLB);
      h += _h_RLA*(h+_h_RLB);
      hL += _hL_RLA*(hL+_hL_RLB);
      j += _j_RLA*(j+_j_RLB);
      m += _m_RLA*(m+_m_RLB);
      mL += _mL_RLA*(mL+_mL_RLB);
      xkr += _xkr_RLA*(xkr+_xkr_RLB);
      xks += _xks_RLA*(xks+_xks_RLB);
      xkur += _xkur_RLA*(xkur+_xkur_RLB);
      xtf += _xtf_RLA*(xtf+_xtf_RLB);
      ykur += _ykur_RLA*(ykur+_ykur_RLB);
      ytf += _ytf_RLA*(ytf+_ytf_RLB);
      RyRi += _dt*_mi_new_RyRi;
      RyRo += _dt*_mi_new_RyRo;
      RyRr += _dt*_mi_new_RyRr;
      store(state_[__jj].CaM, CaM);
      store(state_[__jj].Cai, Cai);
      store(state_[__jj].Caj, Caj);
      store(state_[__jj].Casl, Casl);
      store(state_[__jj].Casr, Casr);
      store(state_[__jj].Ki, Ki);
      store(state_[__jj].Myc, Myc);
      store(state_[__jj].Mym, Mym);
      store(state_[__jj].NaBj, NaBj);
      store(state_[__jj].NaBsl, NaBsl);
      store(state_[__jj].Nai, Nai);
      store(state_[__jj].Naj, Naj);
      store(state_[__jj].Nasl, Nasl);
      store(state_[__jj].RyRi, RyRi);
      store(state_[__jj].RyRo, RyRo);
      store(state_[__jj].RyRr, RyRr);
      store(state_[__jj].SLHj, SLHj);
      store(state_[__jj].SLHsl, SLHsl);
      store(state_[__jj].SLLj, SLLj);
      store(state_[__jj].SLLsl, SLLsl);
      store(state_[__jj].SRB, SRB);
      store(state_[__jj].TnCHc, TnCHc);
      store(state_[__jj].TnCHm, TnCHm);
      store(state_[__jj].TnCL, TnCL);
      store(state_[__jj].d, d);
      store(state_[__jj].f, f);
      store(state_[__jj].fcaBj, fcaBj);
      store(state_[__jj].fcaBsl, fcaBsl);
      store(state_[__jj].h, h);
      store(state_[__jj].hL, hL);
      store(state_[__jj].j, j);
      store(state_[__jj].m, m);
      store(state_[__jj].mL, mL);
      store(state_[__jj].xkr, xkr);
      store(state_[__jj].xks, xks);
      store(state_[__jj].xkur, xkur);
      store(state_[__jj].xtf, xtf);
      store(state_[__jj].ykur, ykur);
      store(state_[__jj].ytf, ytf);
      simdops::store(&__dVm.raw()[__ii],-Iion);
   }
}
#endif //USE_CUDA
   
string ThisReaction::methodName() const
{
   return "Grandi";
}

void ThisReaction::initializeMembraneVoltage(wo_mgarray_ptr<double> __Vm_m)
{
   assert(__Vm_m.size() >= nCells_);

   wo_array_ptr<double> __Vm = __Vm_m.useOn(CPU);
#ifdef USE_CUDA
#define READ_STATE(state,index) (stateData[_##state##_off*nCells_+index])
   wo_array_ptr<double> stateData = stateTransport_.useOn(CPU);
#else //USE_CUDA
#define READ_STATE(state,index) (state_[index/width].state[index % width])
   state_.resize((nCells_+width-1)/width);
#endif //USE_CUDA


   double V_init = -87.840000000000003;
   double V = V_init;
   double NaBj_init = 3.5;
   double NaBj = NaBj_init;
   double NaBsl_init = 0.80000000000000004;
   double NaBsl = NaBsl_init;
   double Naj_init = 9.1359999999999992;
   double Naj = Naj_init;
   double Nasl_init = 9.1359999999999992;
   double Nasl = Nasl_init;
   double Nai_init = 9.1359999999999992;
   double Nai = Nai_init;
   double Ki_init = 120.0;
   double Ki = Ki_init;
   double Casr_init = 0.01;
   double Casr = Casr_init;
   double Caj_init = 0.00017000000000000001;
   double Caj = Caj_init;
   double Casl_init = 0.0001;
   double Casl = Casl_init;
   double Cai_init = 0.0001;
   double Cai = Cai_init;
   double m_init = 0.0;
   double m = m_init;
   double h_init = 1.0;
   double h = h_init;
   double j_init = 1.0;
   double j = j_init;
   double mL_init = 0.0;
   double mL = mL_init;
   double hL_init = 1.0;
   double hL = hL_init;
   double xkr_init = 0.0;
   double xkr = xkr_init;
   double xks_init = 0.0;
   double xks = xks_init;
   double xtf_init = 0.0;
   double xtf = xtf_init;
   double ytf_init = 1.0;
   double ytf = ytf_init;
   double xkur_init = 0.0;
   double xkur = xkur_init;
   double ykur_init = 1.0;
   double ykur = ykur_init;
   double d_init = 0.0;
   double d = d_init;
   double f_init = 1.0;
   double f = f_init;
   double fcaBj_init = 0.025000000000000001;
   double fcaBj = fcaBj_init;
   double fcaBsl_init = 0.014999999999999999;
   double fcaBsl = fcaBsl_init;
   double RyRr_init = 1.0;
   double RyRr = RyRr_init;
   double RyRo_init = 0.0;
   double RyRo = RyRo_init;
   double RyRi_init = 0.0;
   double RyRi = RyRi_init;
   double TnCL_init = 0.01;
   double TnCL = TnCL_init;
   double TnCHc_init = 0.10000000000000001;
   double TnCHc = TnCHc_init;
   double TnCHm_init = 0.01;
   double TnCHm = TnCHm_init;
   double CaM_init = 0.00029999999999999997;
   double CaM = CaM_init;
   double Myc_init = 0.0012999999999999999;
   double Myc = Myc_init;
   double Mym_init = 0.14000000000000001;
   double Mym = Mym_init;
   double SRB_init = 0.002;
   double SRB = SRB_init;
   double SLLj_init = 0.01;
   double SLLj = SLLj_init;
   double SLLsl_init = 0.10000000000000001;
   double SLLsl = SLLsl_init;
   double SLHj_init = 0.0073000000000000001;
   double SLHj = SLHj_init;
   double SLHsl_init = 0.072999999999999995;
   double SLHsl = SLHsl_init;
   for (int iCell=0; iCell<nCells_; iCell++)
   {
      READ_STATE(CaM,iCell) = CaM;
      READ_STATE(Cai,iCell) = Cai;
      READ_STATE(Caj,iCell) = Caj;
      READ_STATE(Casl,iCell) = Casl;
      READ_STATE(Casr,iCell) = Casr;
      READ_STATE(Ki,iCell) = Ki;
      READ_STATE(Myc,iCell) = Myc;
      READ_STATE(Mym,iCell) = Mym;
      READ_STATE(NaBj,iCell) = NaBj;
      READ_STATE(NaBsl,iCell) = NaBsl;
      READ_STATE(Nai,iCell) = Nai;
      READ_STATE(Naj,iCell) = Naj;
      READ_STATE(Nasl,iCell) = Nasl;
      READ_STATE(RyRi,iCell) = RyRi;
      READ_STATE(RyRo,iCell) = RyRo;
      READ_STATE(RyRr,iCell) = RyRr;
      READ_STATE(SLHj,iCell) = SLHj;
      READ_STATE(SLHsl,iCell) = SLHsl;
      READ_STATE(SLLj,iCell) = SLLj;
      READ_STATE(SLLsl,iCell) = SLLsl;
      READ_STATE(SRB,iCell) = SRB;
      READ_STATE(TnCHc,iCell) = TnCHc;
      READ_STATE(TnCHm,iCell) = TnCHm;
      READ_STATE(TnCL,iCell) = TnCL;
      READ_STATE(d,iCell) = d;
      READ_STATE(f,iCell) = f;
      READ_STATE(fcaBj,iCell) = fcaBj;
      READ_STATE(fcaBsl,iCell) = fcaBsl;
      READ_STATE(h,iCell) = h;
      READ_STATE(hL,iCell) = hL;
      READ_STATE(j,iCell) = j;
      READ_STATE(m,iCell) = m;
      READ_STATE(mL,iCell) = mL;
      READ_STATE(xkr,iCell) = xkr;
      READ_STATE(xks,iCell) = xks;
      READ_STATE(xkur,iCell) = xkur;
      READ_STATE(xtf,iCell) = xtf;
      READ_STATE(ykur,iCell) = ykur;
      READ_STATE(ytf,iCell) = ytf;
   }

   __Vm.assign(__Vm.size(), V_init);
}

enum varHandles
{
   CaM_handle,
   Cai_handle,
   Caj_handle,
   Casl_handle,
   Casr_handle,
   Ki_handle,
   Myc_handle,
   Mym_handle,
   NaBj_handle,
   NaBsl_handle,
   Nai_handle,
   Naj_handle,
   Nasl_handle,
   RyRi_handle,
   RyRo_handle,
   RyRr_handle,
   SLHj_handle,
   SLHsl_handle,
   SLLj_handle,
   SLLsl_handle,
   SRB_handle,
   TnCHc_handle,
   TnCHm_handle,
   TnCL_handle,
   d_handle,
   f_handle,
   fcaBj_handle,
   fcaBsl_handle,
   h_handle,
   hL_handle,
   j_handle,
   m_handle,
   mL_handle,
   xkr_handle,
   xks_handle,
   xkur_handle,
   xtf_handle,
   ykur_handle,
   ytf_handle,
   NUMHANDLES
};

const string ThisReaction::getUnit(const std::string& varName) const
{
   if(0) {}
   else if (varName == "CaM") { return "1"; }
   else if (varName == "Cai") { return "mM"; }
   else if (varName == "Caj") { return "mM"; }
   else if (varName == "Casl") { return "mM"; }
   else if (varName == "Casr") { return "mM"; }
   else if (varName == "Ki") { return "mM"; }
   else if (varName == "Myc") { return "1"; }
   else if (varName == "Mym") { return "1"; }
   else if (varName == "NaBj") { return "mM"; }
   else if (varName == "NaBsl") { return "mM"; }
   else if (varName == "Nai") { return "mM"; }
   else if (varName == "Naj") { return "mM"; }
   else if (varName == "Nasl") { return "mM"; }
   else if (varName == "RyRi") { return "1"; }
   else if (varName == "RyRo") { return "1"; }
   else if (varName == "RyRr") { return "1"; }
   else if (varName == "SLHj") { return "1"; }
   else if (varName == "SLHsl") { return "1"; }
   else if (varName == "SLLj") { return "1"; }
   else if (varName == "SLLsl") { return "1"; }
   else if (varName == "SRB") { return "1"; }
   else if (varName == "TnCHc") { return "1"; }
   else if (varName == "TnCHm") { return "1"; }
   else if (varName == "TnCL") { return "1"; }
   else if (varName == "d") { return "1"; }
   else if (varName == "f") { return "1"; }
   else if (varName == "fcaBj") { return "1"; }
   else if (varName == "fcaBsl") { return "1"; }
   else if (varName == "h") { return "1"; }
   else if (varName == "hL") { return "1"; }
   else if (varName == "j") { return "1"; }
   else if (varName == "m") { return "1"; }
   else if (varName == "mL") { return "1"; }
   else if (varName == "xkr") { return "1"; }
   else if (varName == "xks") { return "1"; }
   else if (varName == "xkur") { return "1"; }
   else if (varName == "xtf") { return "1"; }
   else if (varName == "ykur") { return "1"; }
   else if (varName == "ytf") { return "1"; }
   return "INVALID";
}

int ThisReaction::getVarHandle(const std::string& varName) const
{
   if (0) {}
   else if (varName == "CaM") { return CaM_handle; }
   else if (varName == "Cai") { return Cai_handle; }
   else if (varName == "Caj") { return Caj_handle; }
   else if (varName == "Casl") { return Casl_handle; }
   else if (varName == "Casr") { return Casr_handle; }
   else if (varName == "Ki") { return Ki_handle; }
   else if (varName == "Myc") { return Myc_handle; }
   else if (varName == "Mym") { return Mym_handle; }
   else if (varName == "NaBj") { return NaBj_handle; }
   else if (varName == "NaBsl") { return NaBsl_handle; }
   else if (varName == "Nai") { return Nai_handle; }
   else if (varName == "Naj") { return Naj_handle; }
   else if (varName == "Nasl") { return Nasl_handle; }
   else if (varName == "RyRi") { return RyRi_handle; }
   else if (varName == "RyRo") { return RyRo_handle; }
   else if (varName == "RyRr") { return RyRr_handle; }
   else if (varName == "SLHj") { return SLHj_handle; }
   else if (varName == "SLHsl") { return SLHsl_handle; }
   else if (varName == "SLLj") { return SLLj_handle; }
   else if (varName == "SLLsl") { return SLLsl_handle; }
   else if (varName == "SRB") { return SRB_handle; }
   else if (varName == "TnCHc") { return TnCHc_handle; }
   else if (varName == "TnCHm") { return TnCHm_handle; }
   else if (varName == "TnCL") { return TnCL_handle; }
   else if (varName == "d") { return d_handle; }
   else if (varName == "f") { return f_handle; }
   else if (varName == "fcaBj") { return fcaBj_handle; }
   else if (varName == "fcaBsl") { return fcaBsl_handle; }
   else if (varName == "h") { return h_handle; }
   else if (varName == "hL") { return hL_handle; }
   else if (varName == "j") { return j_handle; }
   else if (varName == "m") { return m_handle; }
   else if (varName == "mL") { return mL_handle; }
   else if (varName == "xkr") { return xkr_handle; }
   else if (varName == "xks") { return xks_handle; }
   else if (varName == "xkur") { return xkur_handle; }
   else if (varName == "xtf") { return xtf_handle; }
   else if (varName == "ykur") { return ykur_handle; }
   else if (varName == "ytf") { return ytf_handle; }
   return -1;
}

void ThisReaction::setValue(int iCell, int varHandle, double value) 
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readwrite(CPU);
#endif //USE_CUDA



   if (0) {}
   else if (varHandle == CaM_handle) { READ_STATE(CaM,iCell) = value; }
   else if (varHandle == Cai_handle) { READ_STATE(Cai,iCell) = value; }
   else if (varHandle == Caj_handle) { READ_STATE(Caj,iCell) = value; }
   else if (varHandle == Casl_handle) { READ_STATE(Casl,iCell) = value; }
   else if (varHandle == Casr_handle) { READ_STATE(Casr,iCell) = value; }
   else if (varHandle == Ki_handle) { READ_STATE(Ki,iCell) = value; }
   else if (varHandle == Myc_handle) { READ_STATE(Myc,iCell) = value; }
   else if (varHandle == Mym_handle) { READ_STATE(Mym,iCell) = value; }
   else if (varHandle == NaBj_handle) { READ_STATE(NaBj,iCell) = value; }
   else if (varHandle == NaBsl_handle) { READ_STATE(NaBsl,iCell) = value; }
   else if (varHandle == Nai_handle) { READ_STATE(Nai,iCell) = value; }
   else if (varHandle == Naj_handle) { READ_STATE(Naj,iCell) = value; }
   else if (varHandle == Nasl_handle) { READ_STATE(Nasl,iCell) = value; }
   else if (varHandle == RyRi_handle) { READ_STATE(RyRi,iCell) = value; }
   else if (varHandle == RyRo_handle) { READ_STATE(RyRo,iCell) = value; }
   else if (varHandle == RyRr_handle) { READ_STATE(RyRr,iCell) = value; }
   else if (varHandle == SLHj_handle) { READ_STATE(SLHj,iCell) = value; }
   else if (varHandle == SLHsl_handle) { READ_STATE(SLHsl,iCell) = value; }
   else if (varHandle == SLLj_handle) { READ_STATE(SLLj,iCell) = value; }
   else if (varHandle == SLLsl_handle) { READ_STATE(SLLsl,iCell) = value; }
   else if (varHandle == SRB_handle) { READ_STATE(SRB,iCell) = value; }
   else if (varHandle == TnCHc_handle) { READ_STATE(TnCHc,iCell) = value; }
   else if (varHandle == TnCHm_handle) { READ_STATE(TnCHm,iCell) = value; }
   else if (varHandle == TnCL_handle) { READ_STATE(TnCL,iCell) = value; }
   else if (varHandle == d_handle) { READ_STATE(d,iCell) = value; }
   else if (varHandle == f_handle) { READ_STATE(f,iCell) = value; }
   else if (varHandle == fcaBj_handle) { READ_STATE(fcaBj,iCell) = value; }
   else if (varHandle == fcaBsl_handle) { READ_STATE(fcaBsl,iCell) = value; }
   else if (varHandle == h_handle) { READ_STATE(h,iCell) = value; }
   else if (varHandle == hL_handle) { READ_STATE(hL,iCell) = value; }
   else if (varHandle == j_handle) { READ_STATE(j,iCell) = value; }
   else if (varHandle == m_handle) { READ_STATE(m,iCell) = value; }
   else if (varHandle == mL_handle) { READ_STATE(mL,iCell) = value; }
   else if (varHandle == xkr_handle) { READ_STATE(xkr,iCell) = value; }
   else if (varHandle == xks_handle) { READ_STATE(xks,iCell) = value; }
   else if (varHandle == xkur_handle) { READ_STATE(xkur,iCell) = value; }
   else if (varHandle == xtf_handle) { READ_STATE(xtf,iCell) = value; }
   else if (varHandle == ykur_handle) { READ_STATE(ykur,iCell) = value; }
   else if (varHandle == ytf_handle) { READ_STATE(ytf,iCell) = value; }
}


double ThisReaction::getValue(int iCell, int varHandle) const
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readonly(CPU);
#endif //USE_CUDA


   if (0) {}
   else if (varHandle == CaM_handle) { return READ_STATE(CaM,iCell); }
   else if (varHandle == Cai_handle) { return READ_STATE(Cai,iCell); }
   else if (varHandle == Caj_handle) { return READ_STATE(Caj,iCell); }
   else if (varHandle == Casl_handle) { return READ_STATE(Casl,iCell); }
   else if (varHandle == Casr_handle) { return READ_STATE(Casr,iCell); }
   else if (varHandle == Ki_handle) { return READ_STATE(Ki,iCell); }
   else if (varHandle == Myc_handle) { return READ_STATE(Myc,iCell); }
   else if (varHandle == Mym_handle) { return READ_STATE(Mym,iCell); }
   else if (varHandle == NaBj_handle) { return READ_STATE(NaBj,iCell); }
   else if (varHandle == NaBsl_handle) { return READ_STATE(NaBsl,iCell); }
   else if (varHandle == Nai_handle) { return READ_STATE(Nai,iCell); }
   else if (varHandle == Naj_handle) { return READ_STATE(Naj,iCell); }
   else if (varHandle == Nasl_handle) { return READ_STATE(Nasl,iCell); }
   else if (varHandle == RyRi_handle) { return READ_STATE(RyRi,iCell); }
   else if (varHandle == RyRo_handle) { return READ_STATE(RyRo,iCell); }
   else if (varHandle == RyRr_handle) { return READ_STATE(RyRr,iCell); }
   else if (varHandle == SLHj_handle) { return READ_STATE(SLHj,iCell); }
   else if (varHandle == SLHsl_handle) { return READ_STATE(SLHsl,iCell); }
   else if (varHandle == SLLj_handle) { return READ_STATE(SLLj,iCell); }
   else if (varHandle == SLLsl_handle) { return READ_STATE(SLLsl,iCell); }
   else if (varHandle == SRB_handle) { return READ_STATE(SRB,iCell); }
   else if (varHandle == TnCHc_handle) { return READ_STATE(TnCHc,iCell); }
   else if (varHandle == TnCHm_handle) { return READ_STATE(TnCHm,iCell); }
   else if (varHandle == TnCL_handle) { return READ_STATE(TnCL,iCell); }
   else if (varHandle == d_handle) { return READ_STATE(d,iCell); }
   else if (varHandle == f_handle) { return READ_STATE(f,iCell); }
   else if (varHandle == fcaBj_handle) { return READ_STATE(fcaBj,iCell); }
   else if (varHandle == fcaBsl_handle) { return READ_STATE(fcaBsl,iCell); }
   else if (varHandle == h_handle) { return READ_STATE(h,iCell); }
   else if (varHandle == hL_handle) { return READ_STATE(hL,iCell); }
   else if (varHandle == j_handle) { return READ_STATE(j,iCell); }
   else if (varHandle == m_handle) { return READ_STATE(m,iCell); }
   else if (varHandle == mL_handle) { return READ_STATE(mL,iCell); }
   else if (varHandle == xkr_handle) { return READ_STATE(xkr,iCell); }
   else if (varHandle == xks_handle) { return READ_STATE(xks,iCell); }
   else if (varHandle == xkur_handle) { return READ_STATE(xkur,iCell); }
   else if (varHandle == xtf_handle) { return READ_STATE(xtf,iCell); }
   else if (varHandle == ykur_handle) { return READ_STATE(ykur,iCell); }
   else if (varHandle == ytf_handle) { return READ_STATE(ytf,iCell); }
   return NAN;
}

double ThisReaction::getValue(int iCell, int varHandle, double V) const
{
#ifdef USE_CUDA
   auto stateData = stateTransport_.readonly(CPU);
#endif //USE_CUDA


   const double CaM=READ_STATE(CaM,iCell);
   const double Cai=READ_STATE(Cai,iCell);
   const double Caj=READ_STATE(Caj,iCell);
   const double Casl=READ_STATE(Casl,iCell);
   const double Casr=READ_STATE(Casr,iCell);
   const double Ki=READ_STATE(Ki,iCell);
   const double Myc=READ_STATE(Myc,iCell);
   const double Mym=READ_STATE(Mym,iCell);
   const double NaBj=READ_STATE(NaBj,iCell);
   const double NaBsl=READ_STATE(NaBsl,iCell);
   const double Nai=READ_STATE(Nai,iCell);
   const double Naj=READ_STATE(Naj,iCell);
   const double Nasl=READ_STATE(Nasl,iCell);
   const double RyRi=READ_STATE(RyRi,iCell);
   const double RyRo=READ_STATE(RyRo,iCell);
   const double RyRr=READ_STATE(RyRr,iCell);
   const double SLHj=READ_STATE(SLHj,iCell);
   const double SLHsl=READ_STATE(SLHsl,iCell);
   const double SLLj=READ_STATE(SLLj,iCell);
   const double SLLsl=READ_STATE(SLLsl,iCell);
   const double SRB=READ_STATE(SRB,iCell);
   const double TnCHc=READ_STATE(TnCHc,iCell);
   const double TnCHm=READ_STATE(TnCHm,iCell);
   const double TnCL=READ_STATE(TnCL,iCell);
   const double d=READ_STATE(d,iCell);
   const double f=READ_STATE(f,iCell);
   const double fcaBj=READ_STATE(fcaBj,iCell);
   const double fcaBsl=READ_STATE(fcaBsl,iCell);
   const double h=READ_STATE(h,iCell);
   const double hL=READ_STATE(hL,iCell);
   const double j=READ_STATE(j,iCell);
   const double m=READ_STATE(m,iCell);
   const double mL=READ_STATE(mL,iCell);
   const double xkr=READ_STATE(xkr,iCell);
   const double xks=READ_STATE(xks,iCell);
   const double xkur=READ_STATE(xkur,iCell);
   const double xtf=READ_STATE(xtf,iCell);
   const double ykur=READ_STATE(ykur,iCell);
   const double ytf=READ_STATE(ytf,iCell);
   if (0) {}
   else if (varHandle == CaM_handle)
   {
      return CaM;
   }
   else if (varHandle == Cai_handle)
   {
      return Cai;
   }
   else if (varHandle == Caj_handle)
   {
      return Caj;
   }
   else if (varHandle == Casl_handle)
   {
      return Casl;
   }
   else if (varHandle == Casr_handle)
   {
      return Casr;
   }
   else if (varHandle == Ki_handle)
   {
      return Ki;
   }
   else if (varHandle == Myc_handle)
   {
      return Myc;
   }
   else if (varHandle == Mym_handle)
   {
      return Mym;
   }
   else if (varHandle == NaBj_handle)
   {
      return NaBj;
   }
   else if (varHandle == NaBsl_handle)
   {
      return NaBsl;
   }
   else if (varHandle == Nai_handle)
   {
      return Nai;
   }
   else if (varHandle == Naj_handle)
   {
      return Naj;
   }
   else if (varHandle == Nasl_handle)
   {
      return Nasl;
   }
   else if (varHandle == RyRi_handle)
   {
      return RyRi;
   }
   else if (varHandle == RyRo_handle)
   {
      return RyRo;
   }
   else if (varHandle == RyRr_handle)
   {
      return RyRr;
   }
   else if (varHandle == SLHj_handle)
   {
      return SLHj;
   }
   else if (varHandle == SLHsl_handle)
   {
      return SLHsl;
   }
   else if (varHandle == SLLj_handle)
   {
      return SLLj;
   }
   else if (varHandle == SLLsl_handle)
   {
      return SLLsl;
   }
   else if (varHandle == SRB_handle)
   {
      return SRB;
   }
   else if (varHandle == TnCHc_handle)
   {
      return TnCHc;
   }
   else if (varHandle == TnCHm_handle)
   {
      return TnCHm;
   }
   else if (varHandle == TnCL_handle)
   {
      return TnCL;
   }
   else if (varHandle == d_handle)
   {
      return d;
   }
   else if (varHandle == f_handle)
   {
      return f;
   }
   else if (varHandle == fcaBj_handle)
   {
      return fcaBj;
   }
   else if (varHandle == fcaBsl_handle)
   {
      return fcaBsl;
   }
   else if (varHandle == h_handle)
   {
      return h;
   }
   else if (varHandle == hL_handle)
   {
      return hL;
   }
   else if (varHandle == j_handle)
   {
      return j;
   }
   else if (varHandle == m_handle)
   {
      return m;
   }
   else if (varHandle == mL_handle)
   {
      return mL;
   }
   else if (varHandle == xkr_handle)
   {
      return xkr;
   }
   else if (varHandle == xks_handle)
   {
      return xks;
   }
   else if (varHandle == xkur_handle)
   {
      return xkur;
   }
   else if (varHandle == xtf_handle)
   {
      return xtf;
   }
   else if (varHandle == ykur_handle)
   {
      return ykur;
   }
   else if (varHandle == ytf_handle)
   {
      return ytf;
   }
   return NAN;
}

void ThisReaction::getCheckpointInfo(vector<string>& fieldNames,
                                     vector<string>& fieldUnits) const
{
   fieldNames.clear();
   fieldUnits.clear();
   fieldNames.push_back("CaM");
   fieldUnits.push_back(getUnit("CaM"));
   fieldNames.push_back("Cai");
   fieldUnits.push_back(getUnit("Cai"));
   fieldNames.push_back("Caj");
   fieldUnits.push_back(getUnit("Caj"));
   fieldNames.push_back("Casl");
   fieldUnits.push_back(getUnit("Casl"));
   fieldNames.push_back("Casr");
   fieldUnits.push_back(getUnit("Casr"));
   fieldNames.push_back("Ki");
   fieldUnits.push_back(getUnit("Ki"));
   fieldNames.push_back("Myc");
   fieldUnits.push_back(getUnit("Myc"));
   fieldNames.push_back("Mym");
   fieldUnits.push_back(getUnit("Mym"));
   fieldNames.push_back("NaBj");
   fieldUnits.push_back(getUnit("NaBj"));
   fieldNames.push_back("NaBsl");
   fieldUnits.push_back(getUnit("NaBsl"));
   fieldNames.push_back("Nai");
   fieldUnits.push_back(getUnit("Nai"));
   fieldNames.push_back("Naj");
   fieldUnits.push_back(getUnit("Naj"));
   fieldNames.push_back("Nasl");
   fieldUnits.push_back(getUnit("Nasl"));
   fieldNames.push_back("RyRi");
   fieldUnits.push_back(getUnit("RyRi"));
   fieldNames.push_back("RyRo");
   fieldUnits.push_back(getUnit("RyRo"));
   fieldNames.push_back("RyRr");
   fieldUnits.push_back(getUnit("RyRr"));
   fieldNames.push_back("SLHj");
   fieldUnits.push_back(getUnit("SLHj"));
   fieldNames.push_back("SLHsl");
   fieldUnits.push_back(getUnit("SLHsl"));
   fieldNames.push_back("SLLj");
   fieldUnits.push_back(getUnit("SLLj"));
   fieldNames.push_back("SLLsl");
   fieldUnits.push_back(getUnit("SLLsl"));
   fieldNames.push_back("SRB");
   fieldUnits.push_back(getUnit("SRB"));
   fieldNames.push_back("TnCHc");
   fieldUnits.push_back(getUnit("TnCHc"));
   fieldNames.push_back("TnCHm");
   fieldUnits.push_back(getUnit("TnCHm"));
   fieldNames.push_back("TnCL");
   fieldUnits.push_back(getUnit("TnCL"));
   fieldNames.push_back("d");
   fieldUnits.push_back(getUnit("d"));
   fieldNames.push_back("f");
   fieldUnits.push_back(getUnit("f"));
   fieldNames.push_back("fcaBj");
   fieldUnits.push_back(getUnit("fcaBj"));
   fieldNames.push_back("fcaBsl");
   fieldUnits.push_back(getUnit("fcaBsl"));
   fieldNames.push_back("h");
   fieldUnits.push_back(getUnit("h"));
   fieldNames.push_back("hL");
   fieldUnits.push_back(getUnit("hL"));
   fieldNames.push_back("j");
   fieldUnits.push_back(getUnit("j"));
   fieldNames.push_back("m");
   fieldUnits.push_back(getUnit("m"));
   fieldNames.push_back("mL");
   fieldUnits.push_back(getUnit("mL"));
   fieldNames.push_back("xkr");
   fieldUnits.push_back(getUnit("xkr"));
   fieldNames.push_back("xks");
   fieldUnits.push_back(getUnit("xks"));
   fieldNames.push_back("xkur");
   fieldUnits.push_back(getUnit("xkur"));
   fieldNames.push_back("xtf");
   fieldUnits.push_back(getUnit("xtf"));
   fieldNames.push_back("ykur");
   fieldUnits.push_back(getUnit("ykur"));
   fieldNames.push_back("ytf");
   fieldUnits.push_back(getUnit("ytf"));
}

}
