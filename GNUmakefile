CXX = g++
LIBS = 

INCDIR = include
SRCDIR = src
OBJDIR = obj
EXECDIR = .
#TARGET = mainFitProf_FreeSys
#TARGET = mainFitReweightPMPrint
#TARGET = mainFitSysPrint
#TARGET = mainFit
#TARGET = mainFitWithScanKs2
TARGET = mainFitTest

DEPS = vectorUtils.h utils.h organizer.h experiment.h fitKdePdf.h mockData.h fitGlobal.h
SRC = vectorUtils.cc utils.cc organizer.cc experiment.cc fitKdePdf.cc mockData.cc fitGlobal.cc \
			mainFit.cc

LOCAL_INC = -I$(INCDIR) -I$(ROOTSYS)/include

CXXFLAGS = $(LOCAL_INC) -fpermissive

OBJS = $(SRC:%.cc=$(OBJDIR)/%.o) 

all: $(TARGET)

# Compile c++ objects
$(OBJDIR)/%.o : $(SRCDIR)/%.cc
	$(CXX) $^ -c -o $@ $(shell root-config --libs --cflags) $(CXXFLAGS)

$(TARGET): $(OBJS)
	$(CXX) $^ -o $@ $(shell root-config --libs --cflags) $(CXXFLAGS)
	

clean: 
	$(RM) $(OBJDIR)/*.o core.* $(TARGET) 

emptyrule:: all
