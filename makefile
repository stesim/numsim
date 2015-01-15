
CXX := mpic++ #g++
SRC := $(wildcard *.cpp)
OBJ := $(SRC:%.cpp=%.o)
CXXDEPENDFILE := .depend

#ComSpec only exists on Windows
ifdef ComSpec 
   RM = del /Q
   TARGET := blatt1.exe
else
   TARGET := blatt1
   ifeq ($(shell uname), Linux)
      RM = rm -f
   endif
endif


# Wall: Warn all warnings
# march=native optimize for the current processor
# O3 : optimize max. (only -Ofast is more but it enables Assoziativität bei floats..
# g : nearly no optimization but make debugging possible
# -D_DEBUG : definiert das _DEBUG flag 
# -flto link time optimization
#-flto -fuse-linker-plugin# linking errors if you have old gcc i.e. ld older 2.21
# fwhole-program:Assume that the current compilation unit represents the whole program being compiled. All public functions and variables with the exception of main and those merged by attribute externally_visible become static functions and in effect are optimized more aggressively by interprocedural optimizers. This option should ***not be used in combination with -flto***. Instead relying on a linker plugin should provide safer and more precise information. 
# -Ofast: allow e.g. interchange of 1/c*a*(c+d) = a + 1/c*a*d etc. (=> Assoziativität=AN) => ggf. andere Ergebnisse
# siehe auch http://www.phoronix.com/scan.php?page=article&item=gcc_48_og&num=1
# Ofast ist ähnlich zu -fno-signed-zeros -freciprocal-math -fno-trapping-math -fassociative-math 
# -fopenmp: enable #pragma omp parallel for
#CXXFLAGS := -Wall -D_DEBUG -g -std=c++11 -D_CPP11 #-pedantic
CXXFLAGS := -DCOMM_MPI -Wall -march=native -O3 -flto -fuse-linker-plugin -DNDEBUG -std=c++11 -fno-signed-zeros -freciprocal-math -fno-trapping-math -fassociative-math #-flto -fuse-linker-plugin -fopenmp
LIBS := -pthread
INCLUDE := -I.

$(TARGET): $(CXXDEPENDFILE) $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) $(LIBS) -o $(TARGET)

$(OBJ): %.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(CXXDEPENDFILE): $(SRC) $(SRCPROTO)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -MM $(SRC) > $(CXXDEPENDFILE)

ifneq ($(MAKECMDGOALS),clean)
-include $(CXXDEPENDFILE)
endif

run: $(TARGET)
	./$(TARGET)

debug: $(TARGET)
	gdb $(TARGET)

clean:
	$(RM) $(OBJ) $(TARGET) $(CXXDEPENDFILE)

.PHONY: run debug clean
