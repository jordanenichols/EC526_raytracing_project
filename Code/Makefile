.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================
TARGET1 =  main
C_OBJS1 = main.o

MY_INCLUDES =


CCX = g++
CXXFLAGS = -std=c++11 -g -O2 -fopenmp $(INC)
#-std=c++11x

#============================================================
all: $(TARGET1) $(TARGET2) 

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<

$(TARGET1) :   $(C_OBJS1)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@

# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites 
#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES)


clean:
	rm -f $(TARGET1) $(C_OBJS1) *.~*~ core *.ppm

tar: $(ALL_SOURCES) $(NOTES)
	tar cvf MG.tar $(ALL_SOURCES)  $(NOTES)


