// Classname: TFitConstraintMGaus
// Author: Jan E. Sundermann, Verena Klose (TU Dresden)


//________________________________________________________________
//
// TFitConstraintMGaus::
// --------------------
//
// Fit constraint: mass conservation ( m_i - m_j - alpha * MassConstraint == 0 )
//

using namespace std;

#include <iostream>
#include <iomanip>
#include <sstream>
#include "TFitConstraintMGaus.h"

#include "TClass.h"

//ClassImp(TFitConstraintMGaus)

//----------------
// Constructor --
//----------------
TFitConstraintMGaus::TFitConstraintMGaus() :
   TFitConstraintM() {

   init();

}

TFitConstraintMGaus::TFitConstraintMGaus(vector<TAbsFitParticle*>* ParList1, vector<
      TAbsFitParticle*>* ParList2, Double_t Mass, Double_t Width, Double_t normalization) :
   TFitConstraintM(ParList1, ParList2, Mass, normalization) {

   init();
   setMassConstraint(Mass, Width);

}

TFitConstraintMGaus::TFitConstraintMGaus(const TString &name, const TString &title, vector<
      TAbsFitParticle*>* ParList1, vector<TAbsFitParticle*>* ParList2, Double_t Mass,
      Double_t Width, Double_t normalization) :
   TFitConstraintM(name, title, ParList1, ParList2, Mass, normalization) {

   init();
   setMassConstraint(Mass, Width);

}

void TFitConstraintMGaus::init() {

   _nPar = 1;
   _iniparameters.ResizeTo(1, 1);
   _iniparameters(0, 0) = 1.;
   _parameters.ResizeTo(1, 1);
   _parameters = _iniparameters;

}

//--------------
// Destructor --
//--------------
TFitConstraintMGaus::~TFitConstraintMGaus() {

}

//--------------
// Operations --
//--------------

void TFitConstraintMGaus::setMassConstraint(Double_t Mass, Double_t Width) {

   _TheMassConstraint = Mass;
   _width = Width;
   setCovMatrix(0);
   _covMatrix(0, 0) = (Width * Width) / (Mass * Mass);

}

Double_t TFitConstraintMGaus::getInitValue() {
   // Get initial value of constraint (before the fit)

   Double_t InitValue = CalcMass(&_ParList1, true) - CalcMass(&_ParList2, true) - _iniparameters(0,
         0) * _TheMassConstraint;

   return InitValue / _normalization;

}

Double_t TFitConstraintMGaus::getCurrentValue() {
   // Get value of constraint after the fit

   Double_t CurrentValue = CalcMass(&_ParList1, false) - CalcMass(&_ParList2, false) - _parameters(
         0, 0) * _TheMassConstraint;
   return CurrentValue / _normalization;

}

TMatrixD* TFitConstraintMGaus::getDerivativeAlpha() {
   // Calculate dF/dAlpha = -1 * M

   TMatrixD* DerivativeMatrix = new TMatrixD(1, 1);
   DerivativeMatrix->Zero();

   (*DerivativeMatrix)(0, 0) = -1. * _TheMassConstraint;

   return DerivativeMatrix;

}

TString TFitConstraintMGaus::getInfoString() {
   // Collect information to be used for printout

   stringstream info;
   info << scientific << setprecision(6);

   info << "__________________________" << endl << endl;
   info << "OBJ: " << IsA()->GetName() << "\t" << GetName() << "\t" << GetTitle() << endl;

   info << "initial value: " << getInitValue() << endl;
   info << "current value: " << getCurrentValue() << endl;
   info << "mean mass: " << _TheMassConstraint << endl;
   info << "width: " << _width << endl;
   info << "initial mass: " << _iniparameters(0, 0) * _TheMassConstraint << endl;
   info << "current mass: " << _parameters(0, 0) * _TheMassConstraint << endl;

   return info.str();

}

void TFitConstraintMGaus::print() {
   // Print constraint contents

   std::cout << "\n" << this->getInfoString();

}

