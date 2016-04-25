#ifndef TAbsFitConstraint_hh
#define TAbsFitConstraint_hh

#include <vector>
#include "TAbsFitParticle.h"
#include "TMatrixD.h"
#include "TNamed.h"
#include "TString.h"

class TAbsFitConstraint: public TNamed {

   public:

      TAbsFitConstraint(Double_t normalization = 1.);
      TAbsFitConstraint(const TString &name, const TString &title, Double_t normalization = 1.);
      virtual ~TAbsFitConstraint();

      // returns derivative df/dP with P=(p,E) and f the constraint f=0.
      // The matrix contains one row (df/dp, df/dE).
      virtual TMatrixD* getDerivative(TAbsFitParticle* particle) = 0;
      virtual Double_t getInitValue() = 0;
      virtual Double_t getCurrentValue() = 0;

      // new ---  additional parameters
      Int_t getNPar() {
         return _nPar;
      }

      virtual TMatrixD* getDerivativeAlpha() {
         return 0;
      }

      virtual const TMatrixD* getCovMatrix() const {
         return &_covMatrix;
      }
      virtual void setCovMatrix(const TMatrixD* theCovMatrix);

      virtual const TMatrixD* getCovMatrixFit() const {
         return &_covMatrixFit;
      }
      virtual void setCovMatrixFit(const TMatrixD* theCovMatrixFit);

      virtual const TMatrixD* getCovMatrixDeltaAlpha();

      const TMatrixD* getParIni() {
         return &_iniparameters;
      }
      void setParIni(const TMatrixD* parini);
      virtual void applyDeltaAlpha(TMatrixD* corrMatrix);
      const TMatrixD* getParCurr() {
         return &_parameters;
      }
      Double_t getNormalization() {
         return _normalization;
      }

      virtual TString getInfoString();
      virtual void print();
      virtual void reset();

   protected:

      void calcCovMatrixDeltaAlpha();

      Int_t _nPar;
      Double_t _normalization;

      TMatrixD _covMatrix; // covariance matrix
      TMatrixD _covMatrixFit; // fitted covariance matrix
      TMatrixD _covMatrixDeltaAlpha; // V(deltaAlpha) == V(alpha_meas) - V(alpha_fit)
      TMatrixD _iniparameters; // initialized parameters (parameters values before the fit)
      TMatrixD _parameters; // fitted parameters

      //ClassDef(TAbsFitConstraint, 1)   // Abstract base class for fit constraints

};

#endif
