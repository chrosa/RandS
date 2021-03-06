// Classname: TAbsFitConstraint
// Author: Jan E. Sundermann, Verena Klose (TU Dresden)


//________________________________________________________________
//
// TAbsFitConstraint::
// --------------------
//
// Abstract base class for fit constraints
//

using namespace std;

#include "TAbsFitConstraint.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include "TClass.h"

//ClassImp(TAbsFitConstraint)


TAbsFitConstraint::TAbsFitConstraint(Double_t normalization) :
   TNamed("NoName", "NoTitle"), _covMatrix(), _covMatrixFit(), _covMatrixDeltaAlpha(),
         _iniparameters(), _parameters()

{
   _nPar = 0;
   _normalization = normalization;
}

TAbsFitConstraint::TAbsFitConstraint(const TString &name, const TString &title,
      Double_t normalization) :
   TNamed(name, title), _covMatrix(), _covMatrixFit(), _covMatrixDeltaAlpha(), _iniparameters(),
         _parameters()

{
   _nPar = 0;
   _normalization = normalization;
}

TAbsFitConstraint::~TAbsFitConstraint() {

}

void TAbsFitConstraint::reset() {
   // Reset parameters to initial values

   _parameters = _iniparameters;
   setCovMatrixFit(0);

}

void TAbsFitConstraint::setCovMatrix(const TMatrixD* theCovMatrix) {
   // Set measured alpha covariance matrix

   _covMatrix.ResizeTo(_nPar, _nPar);
   if (theCovMatrix == 0) {
      _covMatrix.Zero();
   } else if (theCovMatrix->GetNcols() == _nPar && theCovMatrix->GetNrows() == _nPar) {
      _covMatrix = (*theCovMatrix);
   } else {
      std::cout << "\n" << GetName()
            << "::setCovMatrix - Measured alpha covariance matrix needs to be a " << _nPar << "x"
            << _nPar << " matrix.";
   }

}

void TAbsFitConstraint::setCovMatrixFit(const TMatrixD* theCovMatrixFit) {
   // Set the fitted covariance matrix

   _covMatrixFit.ResizeTo(_nPar, _nPar);
   if (theCovMatrixFit == 0) {
      _covMatrixFit.Zero();
   } else if (theCovMatrixFit->GetNcols() == _nPar && theCovMatrixFit->GetNrows() == _nPar) {
      _covMatrixFit = (*theCovMatrixFit);
   } else {
      std::cout << "\n" << GetName()
            << "::setCovMatrixFit - Fitted covariance matrix needs to be a " << _nPar << "x"
            << _nPar << " matrix.";
   }

}

void TAbsFitConstraint::calcCovMatrixDeltaAlpha() {
   // Calculates V(deltaAlpha) ==  V(alpha_meas) - V(alpha_fit)

   _covMatrixDeltaAlpha.ResizeTo(_nPar, _nPar);
   _covMatrixDeltaAlpha = _covMatrix;
   if (_covMatrixFit.GetNrows() == _nPar && _covMatrixFit.GetNcols() == _nPar)
      _covMatrixDeltaAlpha -= _covMatrixFit;
   else
      std::cout << "\n" << GetName()
            << "::calcCovMatrixDeltaAlpha - _covMatrixFit probably not set.";
}

void TAbsFitConstraint::applyDeltaAlpha(TMatrixD* corrMatrix) {
   // Apply corrections to the parameters wrt. to the
   // initial parameters alpha* = alpha + delta(alpha)

   _parameters = _iniparameters;
   _parameters += (*corrMatrix);

}

void TAbsFitConstraint::setParIni(const TMatrixD* parini) {
   // Set initial parameter values (before the fit)

   if (parini == 0)
      return;
   else if (parini->GetNrows() == _iniparameters.GetNrows() && parini->GetNcols()
         == _iniparameters.GetNcols())
      _iniparameters = (*parini);
   else {
      std::cout << "\n" << GetName() << "::setParIni - Matrices don't fit.";
      return;
   }

}

const TMatrixD* TAbsFitConstraint::getCovMatrixDeltaAlpha() {
   // Returns covariance matrix delta(alpha)

   calcCovMatrixDeltaAlpha();
   return &_covMatrixDeltaAlpha;

}

TString TAbsFitConstraint::getInfoString() {
   // Collect information to be used for printout

   std::stringstream info;
   info << scientific << setprecision(6);

   info << "__________________________" << endl << endl;
   info << "OBJ: " << IsA()->GetName() << "\t" << GetName() << "\t" << GetTitle() << endl;

   info << "initial value: " << getInitValue() << endl;
   info << "current value: " << getCurrentValue() << endl;

   return info.str();

}

void TAbsFitConstraint::print() {
   // Print constraint contents

   std::cout << "\n" << this->getInfoString();

}

