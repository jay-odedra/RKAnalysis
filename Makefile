ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lTMVA

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)

CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR) -I.
OUTLIB	         = $(INCLUDEDIR)/lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/


$(OUTLIB)BParkBase.o: $(INCLUDEDIR)/src/BParkBase.C \
	$(INCLUDEDIR)/src/BPark.cc \
	$(INCLUDEDIR)/src/TestMva.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)BParkBase.o $<
$(OUTLIB)BPark.o: $(INCLUDEDIR)/src/BPark.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)BPark.o $<
$(OUTLIB)TestMva.o: $(INCLUDEDIR)/src/TestMva.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -L../FastForest/build/ -lfastforest -o $(OUTLIB)TestMva.o $<
$(OUTLIB)BParkBaseNew.o: $(INCLUDEDIR)/src/BParkBaseNew.C \
	$(INCLUDEDIR)/src/BParkBaseNew.C \
	$(INCLUDEDIR)/src/SkimmerWithKStar.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)BParkBaseNew.o $<
#$(OUTLIB)BParkBaseNewWithKStar.o: $(INCLUDEDIR)/src/BParkBaseNewWithKStar.C \
#	$(INCLUDEDIR)/src/BParkBaseNewWithKStar.C \
#	$(INCLUDEDIR)/src/SkimmerWithKStar.cc 
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)BParkBaseNewWithKStar.o $<
$(OUTLIB)SkimmerWithKStar.o: $(INCLUDEDIR)/src/SkimmerWithKStar.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -L../FastForest/build/ -lfastforest -o $(OUTLIB)SkimmerWithKStar.o $<

# ==================== BParkApp =============================================
SkimmerWithKStar:  $(INCLUDEDIR)/src/BParkApp.C \
	$(OUTLIB)BParkBase.o \
	$(OUTLIB)BPark.o \
	$(OUTLIB)TestMva.o \
	$(OUTLIB)BParkBaseNew.o \
	$(OUTLIB)SkimmerWithKStar.o 
	$(CXX) $(CXXFLAGS) -ldl -L../FastForest/build/ -lfastforest -o SkimmerWithKStar.exe $(OUTLIB)/*.o $(GLIBS) $(LDFLAGS) $ $<
SkimmerWithKStar.clean:
	rm -f SkimmerWithKStar.exe

# ==================== reduced trees =============================================

clean:
	rm -f $(OUTLIB)*.o
	rm -f SkimmerWithKStar.exe

all:  SkimmerWithKStar
