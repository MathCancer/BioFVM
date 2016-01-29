PROGRAM_NAME := project

CC := g++

# change this to your own CPU archicture. 
# Here is a list for gcc 4.9.0
# https://gcc.gnu.org/onlinedocs/gcc-4.9.0/gcc/i386-and-x86-64-Options.html#i386-and-x86-64-Options 

ARCH := core2 # a reasonably safe default for most CPUs since 2007
# ARCH := corei7
# ARCH := corei7-avx # earlier i7 
# ARCH := core-avx-i # i7 ivy bridge or newer 
# ARCH := core-avx2 # i7 with Haswell or newer
# ARCH := nehalem
# ARCH := westmere
# ARCH := sandybridge
# ARCH := ivybridge
# ARCH := haswell
# ARCH := broadwell
# ARCH := bonnell
# ARCH := silvermont
# ARCH := nocona #64-bit pentium 4 or later 

CFLAGS := -march=$(ARCH) -O3 -s -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11

BioFVM_OBJECTS := BioFVM_vector.o BioFVM_matlab.o BioFVM_utilities.o BioFVM_mesh.o \
BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_basic_agent.o \
BioFVM_agent_container.o BioFVM_MultiCellDS.o 

pugixml_OBJECTS := pugixml.o

COMPILE_COMMAND := $(CC) $(CFLAGS) 

all: $(BioFVM_OBJECTS) $(pugixml_OBJECTS) main.o
	$(COMPILE_COMMAND) -o $(PROGRAM_NAME) $(BioFVM_OBJECTS) $(pugixml_OBJECTS) main.o

template: $(BioFVM_OBJECTS) $(pugixml_OBJECTS) template_project.cpp
	$(COMPILE_COMMAND) -o template $(BioFVM_OBJECTS) $(pugixml_OBJECTS) template_project.cpp
	
main.o: main.cpp
	$(COMPILE_COMMAND) -c main.cpp

BioFVM_vector.o: BioFVM_vector.cpp
	$(COMPILE_COMMAND) -c BioFVM_vector.cpp 

BioFVM_agent_container.o: BioFVM_agent_container.cpp
	$(COMPILE_COMMAND) -c BioFVM_agent_container.cpp 
	
BioFVM_mesh.o: BioFVM_mesh.cpp
	$(COMPILE_COMMAND) -c BioFVM_mesh.cpp 

BioFVM_microenvironment.o: BioFVM_microenvironment.cpp
	$(COMPILE_COMMAND) -c BioFVM_microenvironment.cpp 

BioFVM_solvers.o: BioFVM_solvers.cpp
	$(COMPILE_COMMAND) -c BioFVM_solvers.cpp 

BioFVM_utilities.o: BioFVM_utilities.cpp
	$(COMPILE_COMMAND) -c BioFVM_utilities.cpp 
	
BioFVM_basic_agent.o: BioFVM_basic_agent.cpp
	$(COMPILE_COMMAND) -c BioFVM_basic_agent.cpp 
	
BioFVM_matlab.o: BioFVM_matlab.cpp
	$(COMPILE_COMMAND) -c BioFVM_matlab.cpp

BioFVM_MultiCellDS.o: BioFVM_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c BioFVM_MultiCellDS.cpp

pugixml.o: pugixml.cpp
	$(COMPILE_COMMAND) -c pugixml.cpp
	
clean:
	rm -f *.o
	rm -f $(PROGRAM_NAME)*
	
zip:
	zip $$(date +%b_%d_%Y_%H%M).zip *.cpp *.h *akefile* *.xml *.tex *.bib *hanges*.txt config/*.xml *.txt
	zip VERSION_$(VERSION).zip *.cpp *.h *akefile* *.xml *.tex *.bib *hanges*.txt *.txt
	mv *.zip archives/
