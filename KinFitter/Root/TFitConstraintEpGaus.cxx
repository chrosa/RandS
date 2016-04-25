// Classname: TFitConstraintEpGaus
// Author: Holger Enderle (Uni Hamburg)


//________________________________________________________________
// 
// TFitConstraintEpGaus::
// --------------------
//
// Fit constraint: energy and momentum conservation
//

#include "KinFitter/TFitConstraintEpGaus.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include "TClass.h"

//----------------
// Constructor --
//----------------
TFitConstraintEpGaus::TFitConstraintEpGaus()
  : TFitConstraintEp()
{

  init();

}

TFitConstraintEpGaus::TFitConstraintEpGaus(std::vector<TAbsFitParticle*>* particles, 
					   TFitConstraintEp::component thecomponent, 
					   Double_t constraint,
					   Double_t Width,
					   Double_t normalization)
  :TFitConstraintEp(particles, thecomponent, constraint, normalization)
{

  init();
  setConstraint( constraint, Width );

}

TFitConstraintEpGaus::TFitConstraintEpGaus(const TString &name, const TString &title,
					   std::vector<TAbsFitParticle*>* particles, 
					   TFitConstraintEp::component thecomponent, 
					   Double_t constraint,
					   Double_t Width,
					   Double_t normalization)
  : TFitConstraintEp( name, title, particles, thecomponent, constraint, normalization)
{

  init();
  setConstraint( constraint, Width );

}

void
TFitConstraintEpGaus::init() {

  _nPar = 1;
  _iniparameters.ResizeTo(1,1);
  _iniparameters(0,0) = 0.;
  _parameters.ResizeTo(1,1);
  _parameters = _iniparameters;

}


//--------------
// Destructor --
//--------------
TFitConstraintEpGaus::~TFitConstraintEpGaus() {

}

//--------------
// Operations --
//--------------

void TFitConstraintEpGaus::setConstraint(Double_t constraint, Double_t Width) { 
  
  _constraint = constraint;
  _width = Width;
  setCovMatrix( 0 );
  _covMatrix(0,0) = (Width*Width);

}

TMatrixD* TFitConstraintEpGaus::getDerivativeAlpha() {
  // Calculate dF/dAlpha

  TMatrixD* DerivativeMatrix = new TMatrixD(1,1);
  DerivativeMatrix->Zero();
  (*DerivativeMatrix)(0,0) = -1.;
  return DerivativeMatrix;
}


Double_t TFitConstraintEpGaus::getInitValue() {
  // Get initial value of constraint (before the fit)

  Double_t InitValue(0) ; 
  UInt_t Npart = _particles.size();
  for (unsigned int i=0;i<Npart;i++) {
    const TLorentzVector* FourVec = _particles[i]->getIni4Vec();
    InitValue += (*FourVec)[(int) _component];
  }
  InitValue -= _constraint;
  InitValue -= _iniparameters(0,0);
  return InitValue / _normalization;
}

Double_t TFitConstraintEpGaus::getCurrentValue() {
  // Get value of constraint after the fit

  Double_t CurrentValue(0);
  UInt_t Npart = _particles.size();
  for (unsigned int i=0;i<Npart;i++) {
    const TLorentzVector* FourVec = _particles[i]->getCurr4Vec();
    CurrentValue += (*FourVec)[(int) _component];
  }
  CurrentValue -= _constraint;
  CurrentValue -= _parameters(0,0);
  return CurrentValue / _normalization;
}

TString TFitConstraintEpGaus::getInfoString() {
  // Collect information to be used for printout

  std::stringstream info;
  info << std::scientific << std::setprecision(6);

  info << "__________________________" << std::endl
       << std::endl;
  info << "OBJ: " << IsA()->GetName() << "\t" << GetName() << "\t" << GetTitle() << std::endl;

  info << "initial value: " << getInitValue() << std::endl;
  info << "current value: " << getCurrentValue() << std::endl;
  info << "component: " << _component << std::endl;
  info << "constraint: " << _constraint << std::endl;
  info << "width: " << _width << std::endl;
  info << "initial constraint: " << _iniparameters(0,0)  << std::endl;
  info << "current constraint: " << _parameters(0,0)  << std::endl;

  return info.str();

}

void TFitConstraintEpGaus::print() {
  // Print constraint contents

  cout << "\n" << this->getInfoString();

}
