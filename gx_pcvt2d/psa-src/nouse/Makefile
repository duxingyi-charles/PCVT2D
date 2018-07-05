###################### user-configuration section #######################
INC := -I/opt/local/include
LIB := -L/opt/local/lib

# 32 bit
#MARCH := -m32

# 64 bit
MARCH := -m64
#########################################################################

CXX := g++
CXXFLAGS := -Wall -fopenmp
OPTFLAGS := -O3 $(MARCH)

LIB += -lfftw3 -lm -lcairo

CXXFILES := main.cpp fileio.cpp measurements.cpp spectralanalysis.cpp
OBJS := $(patsubst %.cpp, %.o, $(CXXFILES))
TARGET := psa

.PHONY: all clean cleanall rebuild

all: $(TARGET)

$(TARGET): $(OBJS)
	@$(CXX) $(CXXFLAGS) $(OPTFLAGS) -o $(TARGET) $(LIB) $(OBJS)

%.o : %.cpp
	@$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(INC) -o $@ -c $<

clean:
	rm -rf *.o

cleanall:	clean
	rm -rf $(TARGET)

rebuild:
	$(MAKE) cleanall
	$(MAKE) all
