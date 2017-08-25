C=g++
LD=g++

ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)

CFLAGS =
LFLAGS =

RCXX=$(CFLAGS) $(ROOTCFLAGS) -O2
RLXX=$(LFLAGS) $(ROOTLIBS)

SRC0= RandS.o MyJet.o MyElectron.o MyPhoton.o MyTau.o MyMuon.o SmearFunction.o
SRC1= MyResolution.o MyJet.o
SRC2= Prediction.o
SRC3= MyMETStudies.o MyJet.o MyElectron.o MyPhoton.o MyTau.o MyMuon.o
SRC4= MyTriggerStudies.o MyJet.o
SRC5= MyABCDStudies.o MyJet.o MyElectron.o MyPhoton.o MyTau.o MyMuon.o

%.o: %.C %.h
	$(C) $(RCXX) -c $<

all: GetPrediction RunMyResolution RunMyMETStudies RunMyTriggerStudies RunMyABCDStudies RunRandS

RunRandS: $(SRC0) RunRandS.o
	$(LD) $(SRC0) RunRandS.o $(RLXX) -l TreePlayer -L./TKinFitter -lkinfitter -o RunRandS.x
	@echo '-> RunRandS.x created!'

RunMyResolution: $(SRC1) RunMyResolution.o
	$(LD) $(SRC1) RunMyResolution.o $(RLXX) -l TreePlayer -o RunMyResolution.x
	@echo '-> RunMyResolution.x created!'

RunMyMETStudies: $(SRC3) RunMyMETStudies.o
	$(LD) $(SRC3) RunMyMETStudies.o $(RLXX) -l TreePlayer -o RunMyMETStudies.x
	@echo '-> RunMyMETStudies.x created!'

RunMyABCDStudies: $(SRC5) RunMyABCDStudies.o
	$(LD) $(SRC5) RunMyABCDStudies.o $(RLXX) -l TreePlayer -o RunMyABCDStudies.x
	@echo '-> RunMyABCDStudies.x created!'

RunMyTriggerStudies: $(SRC4) RunMyTriggerStudies.o
	$(LD) $(SRC4) RunMyTriggerStudies.o $(RLXX) -l TreePlayer -o RunMyTriggerStudies.x
	@echo '-> RunMyTriggerStudies.x created!'

GetPrediction: $(SRC2) GetPrediction.o
	$(LD) $(SRC2) GetPrediction.o $(RLXX) -o GetPrediction.x
	@echo '-> GetPrediction.x created!'

RunRandS.o: RunRandS.C
	$(C) $(RCXX) -c $<

RunMyResolution.o: RunMyResolution.C
	$(C) $(RCXX) -c $<

RunMyMETStudies.o: RunMyMETStudies.C
	$(C) $(RCXX) -c $<

RunMyABCDStudies.o: RunMyABCDStudies.C
	$(C) $(RCXX) -c $<

RunMyTriggerStudies.o: RunMyTriggerStudies.C
	$(C) $(RCXX) -c $<

GetPrediction.o: GetPrediction.C
	$(C) $(RCXX) -c $<

clean:
	@rm -f *~
	@rm -f *.o 
	@rm -f *.x
