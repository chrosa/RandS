# Makefile for TKinFitter code
# christian.sander@desy.de

C=g++
LD=g++

ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)

CFLAGS = -fPIC -Wno-deprecated -O2
LFLAGS = -lz

RCXX=$(CFLAGS) $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS)

SRC=TAbsFitConstraint.cc TFitConstraintMGaus.cc TFitParticleEScaledMomDev.cc TFitParticleMCCart.cc TFitParticleMomDev.cc TAbsFitParticle.cc TFitParticleCart.cc TFitParticleESpher.cc TFitParticleMCMomDev.cc TFitParticleSpher.cc TFitConstraintEpGaus.cc TFitConstraintEp.cc TFitParticleECart.cc TFitParticleEtEtaPhi.cc TFitParticleMCPInvSpher.cc TKinFitter.cc TFitConstraintM.cc TFitParticleEMomDev.cc TFitParticleEtThetaPhi.cc TFitParticleMCSpher.cc

%.o: %.cc
		$(C) $(RCXX) -c $<

all: libkinfitter.so

libkinfitter.so: $(SRC:.cc=.o)
	${LD} $(SRC:.cc=.o) ${RLXX} -shared -o libkinfitter.so 
	@echo '-> shared library created.'
		

clean:
	@rm -f *~
	@rm -f *.o 

cleanlib:
	@rm -f libkinfitter.so*


