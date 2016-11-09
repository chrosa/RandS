#include <iostream>
#include "TKinFitter.h"
#include "TFitParticleEtEtaPhi.h"
#include "TFitConstraintM.h"

#include "TLorentzVector.h"
#include "TMatrixD.h"

int main()
{
   std::cout << "Do nothing\n" << endl;

   // The fit
   TKinFitter* myFit = new TKinFitter();

   // The particles before fitting
   TLorentzVector* lvec1 = new TLorentzVector(0,0,0,38);
   TLorentzVector* lvec2 = new TLorentzVector(0,0,0,41);

   // The covariance Matrix
   TMatrixD* covMat1 = new TMatrixD(3,3);
   TMatrixD* covMat2 = new TMatrixD(3,3);
   
   (*covMat1)(0,0)=5;
   (*covMat1)(1,1)=5;
   (*covMat1)(2,2)=5;

   (*covMat2)(0,0)=5;
   (*covMat2)(1,1)=5;
   (*covMat2)(2,2)=5;

   TFitParticleEtEtaPhi* jet1 = new TFitParticleEtEtaPhi("","",lvec1, covMat1);
   TFitParticleEtEtaPhi* jet2 = new TFitParticleEtEtaPhi("","",lvec2, covMat2);

   TFitConstraintM* massConstr = new TFitConstraintM();
   massConstr->addParticle1( jet1 );
   massConstr->addParticle1( jet2 );
   massConstr->setMassConstraint( 80. );
   
   myFit->addConstraint(massConstr);
   myFit->fit();

}
