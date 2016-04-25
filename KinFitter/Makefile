# Makefile for Global Fit code
# benedikt.mura@desy.de

C=g++
LD=g++
#F77=g77
SPECIALFLAGS= #--thread_safe --exceptions
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs) #-lMinuit

#add -DNDEBUG to suppress debugging output

#CFLAGS = $(SPECIALFLAGS) -I. -I./include -I$(SRT_PUBLIC_CONTEXT)/include -I$(ROOTSYS)/include -DNDEBUG -Wno-deprecated -g -m32 -Wa,--32
#LFLAGS = $(SPECIALFLAGS) -L../../lib/$(SRT_SUBDIR)/ -lz -DNDEBUG -g -m32 -lg2c
CFLAGS = $(SPECIALFLAGS) -fPIC -I. -I./include -I$(SRT_PUBLIC_CONTEXT)/include -I$(ROOTSYS)/include -DNDEBUG -Wno-deprecated -g -Wa,
LFLAGS = $(SPECIALFLAGS) -L../../lib/$(SRT_SUBDIR)/ -lz -DNDEBUG -g #-lg2c

RCXX=$(CFLAGS) $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS)

SRC=TAbsFitConstraint.cc TFitConstraintMGaus.cc TFitParticleEScaledMomDev.cc TFitParticleMCCart.cc TFitParticleMomDev.cc TAbsFitParticle.cc TFitParticleCart.cc TFitParticleESpher.cc TFitParticleMCMomDev.cc TFitParticleSpher.cc TFitConstraintEpGaus.cc TFitConstraintEp.cc TFitParticleECart.cc TFitParticleEtEtaPhi.cc TFitParticleMCPInvSpher.cc TKinFitter.cc TFitConstraintM.cc TFitParticleEMomDev.cc TFitParticleEtThetaPhi.cc TFitParticleMCSpher.cc
#SRC= test.cc TAbsFitConstraint.cc TAbsFitParticle.cc TKinFitter.cc

%.o: %.cc
		$(C) $(RCXX) -c $<

#all: FitExec

libkinfitter.so: $(SRC:.cc=.o)
	        ${LD} $(SRC:.cc=.o) ${RLXX} -shared -Wl,-soname,libkinfitter.so -o libkinfitter.so 
		@echo '-> shared library created.'
		
#test.o: test.cc
#		$(C) $(RCXX) -c test.cc
		
#FitExec: test.o libkinfitter.so
#		${LD} $< $(RLXX) -L. -lkinfitter -o FitExec	
#		@echo '-> test executable created.'

clean:
		@rm -f *~
		@rm -f *.o 
#		@rm -f FitExec

cleanlib:
		@rm -f libkinfitter.so*


