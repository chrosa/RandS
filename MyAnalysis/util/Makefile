C=g++
LD=g++

ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)

CFLAGS =
LFLAGS =

RCXX=$(CFLAGS) $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS)

SRC1= Prediction.o

SRC2= BinPrediction.o RA2bBin.o

%.o: %.C %.h
		$(C) $(RCXX) -c $<

all: GetPrediction

GetPrediction.o: GetPrediction.C
		$(C) $(RCXX) -c $<

GetPrediction: $(SRC1) GetPrediction.o
		$(LD) $(SRC1) GetPrediction.o $(RLXX) -o GetPrediction.x
		@echo '-> GetPrediction.x created!'

clean:
		@rm -f *~
		@rm -f *.o 
		@rm -f *.x
