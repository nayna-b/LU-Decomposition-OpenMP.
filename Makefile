# Compiler and flags
CXX = icpc
NOWARN = -wd3180
OPT = -O2 -g
DEBUG = -O0 -g
OMP = -fopenmp
NUMA = -lnuma
LDFLAGS = -lrt

# Executable names
EXEC = lu-omp-num
OBJ = $(EXEC) $(EXEC)-debug $(EXEC)-serial

# Matrix size to use when running
MATRIX_SIZE = 8000
W := $(shell grep processor /proc/cpuinfo | wc -l)

# Thread checker and viewer
CHECKER = inspxe-cl -collect=ti3 -knob scope=extreme -result-dir=check_races
VIEWER = inspxe-gui

all: $(OBJ)

# Build optimized parallel version with NUMA support
$(EXEC): $(EXEC).cpp
	$(CXX) $(OPT) $(OMP) $(NUMA) -o $(EXEC) $(EXEC).cpp $(LDFLAGS)

# Build debug version
$(EXEC)-debug: $(EXEC).cpp
	$(CXX) $(DEBUG) $(OMP) $(NUMA) -o $(EXEC)-debug $(EXEC).cpp $(LDFLAGS)

# Build serial version
$(EXEC)-serial: $(EXEC).cpp
	$(CXX) $(OPT) $(NOWARN) -o $(EXEC)-serial $(EXEC).cpp $(LDFLAGS) -liomp5

# Run with all available workers
runp: $(EXEC)
	@echo "Running with $(W) threads..."
	./$(EXEC) $(MATRIX_SIZE) $(W)

# Run serial version with 1 thread
runs: $(EXEC)-serial
	@echo "Running serial version..."
	./$(EXEC)-serial $(MATRIX_SIZE) 1

# Run thread checker (Intel Inspector)
check: $(EXEC)
	$(CHECKER) ./$(EXEC) 300 4

# Open Intel Inspector GUI
view:
	$(VIEWER) check_races*

# Clean all generated files
clean:
	/bin/rm -rf $(OBJ) check_races*

