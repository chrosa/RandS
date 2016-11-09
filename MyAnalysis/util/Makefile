C=g++
LD=g++

ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)

CFLAGS =
LFLAGS =

RCXX=$(CFLAGS) $(ROOTCFLAGS) -O2
RLXX=$(LFLAGS) $(ROOTLIBS)

SRC0= RandS.o MyJet.o SmearFunction.o
SRC1= MyResolution.o MyJet.o
SRC2= Prediction.o

%.o: %.C %.h
	$(C) $(RCXX) -c $<

all: GetPrediction RunMyResolution RunRandS

RunRandS: $(SRC0) RunRandS.o
	$(LD) $(SRC0) RunRandS.o $(RLXX) -l TreePlayer -L./TKinFitter -lkinfitter -o RunRandS.x
	@echo '-> RunRandS.x created!'

RunMyResolution: $(SRC1) RunMyResolution.o
	$(LD) $(SRC1) RunMyResolution.o $(RLXX) -l TreePlayer -o RunMyResolution.x
	@echo '-> RunMyResolution.x created!'

GetPrediction: $(SRC2) GetPrediction.o
	$(LD) $(SRC2) GetPrediction.o $(RLXX) -o GetPrediction.x
	@echo '-> GetPrediction.x created!'

RunRandS.o: RunRandS.C
	$(C) $(RCXX) -c $<

RunMyResolution.o: RunMyResolution.C
	$(C) $(RCXX) -c $<

GetPrediction.o: GetPrediction.C
	$(C) $(RCXX) -c $<

clean:
	@rm -f *~
	@rm -f *.o 
	@rm -f *.x
