using namespace std;

#ifndef TFitConstraintEpGaus_hh
#define TFitConstraintEpGaus_hh

#include "TFitConstraintEp.h"




class TFitConstraintEpGaus: public TFitConstraintEp {

public :

  enum component {
    pX, pY, pZ, E
  };

  TFitConstraintEpGaus( );

  TFitConstraintEpGaus(  std::vector<TAbsFitParticle*>* particles, 
			 TFitConstraintEp::component thecomponent, Double_t constraint = 0.,
			 Double_t Width = 0, Double_t normalization = 1);

  TFitConstraintEpGaus(  const TString &name, const TString &title,
		     std::vector<TAbsFitParticle*>* particles, TFitConstraintEp::component thecomponent, 
		     Double_t constraint = 0., Double_t Width = 0, Double_t normalization = 1);
  virtual ~TFitConstraintEpGaus();
    






  // returns derivative df/dP with P=(p,E) and f the constraint f=0.
  // The matrix contains one row (df/dp, df/dE).
  virtual TMatrixD* getDerivativeAlpha();
  virtual Double_t getInitValue();
  virtual Double_t getCurrentValue();

  void setConstraint(Double_t constraint, Double_t Width);

  virtual TString getInfoString();
  virtual void print(); 
 
protected :
  
  Double_t _width;

  void init();


};

#endif
