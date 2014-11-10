TARGET := blatt1
CXX := g++
SRC := $(wildcard *.cpp)
OBJ := $(SRC:%.cpp=%.o)
CXXDEPENDFILE := .depend

#ComSpec only exists on Windows
ifdef ComSpec 
   RM = del /Q
   TARGET +=.exe
else
   ifeq ($(shell uname), Linux)
      RM = rm -f
   endif
endif


# Wall: Warn all warnings
# march=native optimize for the current processor
# O3 : optimize max. (only -Ofast is more but it enables Assoziativität bei floats..
# g : nearly no optimization but make debugging possible
# -D_DEBUG : definiert das _DEBUG flag 
# flto link time optimization
#-flto -fwhole-program -fuse-linker-plugin# may cause linking errors(!)
# -Ofast: allow e.g. interchange of 1/c*a*(c+d) = a + 1/c*a*d etc. (=> Assoziativität=AN) => ggf. andere Ergebnisse
# siehe auch http://www.phoronix.com/scan.php?page=article&item=gcc_48_og&num=1
CXXFLAGS := -Wall -march=native -pedantic -D_DEBUG -O3 #-g
LIBS :=

$(TARGET): $(CXXDEPENDFILE) $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) $(LIBS) -o $(TARGET)

$(OBJ): %.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(CXXDEPENDFILE): $(SRC) $(SRCPROTO)
	$(CXX) $(CXXFLAGS) -MM $(SRC) > $(CXXDEPENDFILE)

-include $(CXXDEPENDFILE)

run: $(TARGET)
	./$(TARGET)

debug: $(TARGET)
	gdb $(TARGET)

clean:
	$(RM) $(OBJ) $(TARGET) $(CXXDEPENDFILE)

.PHONY: run debug clean
