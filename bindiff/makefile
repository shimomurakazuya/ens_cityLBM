CPPVER := 14

CXX := g++
CXXFLAGS := -O2 -std=c++$(CPPVER)
CXXFLAGS += -fopenmp

LD := $(CXX)
LDFLAGS := -std=c++$(CPPVER)
LDFLAGS += -lz
LDFLAGS += -fopenmp

RUN:= bindiff.out
.PHONY: all clean

all: $(RUN)

$(RUN): main.o
	$(LD) -o $@ $^ $(LDFLAGS)

.cpp.o:
	$(CXX) $< -c $(CXXFLAGS)

clean:
	rm -f $(RUN)
	rm -f *.o

