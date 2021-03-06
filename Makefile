FC = gfortran
FCFLAGS = -g -fcheck=all -Warray-bounds

# source files and objects
SRCS = Tools.f90 Functions.f90 Solver.f90 main.f90

# dependances
#main.f90: Tools.o

# program name
PROGRAM = runUzawa

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod
