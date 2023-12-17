TARGET  := hh

BIN_DIR := .
SRC_DIR := src
BUILD_DIR := build

SRC := $(shell find $(SRC_DIR) -type f -name *.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
DEP := $(OBJ:.o=.d)

INC_DIR := $(shell find $(SRC_DIR) -type d)
INC_FLAG := $(addprefix -I,$(INC_DIR))

CXX    := mpic++

ARMADILLO_LIB := -larmadillo
LAPACK_LIB    := -llapack
BLAS_LIB      := -lblas
GMP_LIB		  := -lgmp
GMPXX_LIB     := -lgmpxx
GSL_LIB		  := -lgsl
BOOST_LIB     := -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system
ARPACK_LIB   := -larpack

ALL_LIBS     := $(ARMADILLO_LIB) $(LAPACK_LIB) $(BLAS_LIB)  $(ARPACK_LIB) ##$(GMP_LIB) $(GMPXX_LIB) $(GSL_LIB)

CXX_WARNINGS := -Wall -pedantic-errors
CXX_DEBUG    := #-g -DDEBUG
CXX_OPTIM    := -O3
CXX_DEP      := -MMD

CXXFLAGS    := -std=c++11 $(CXX_DEBUG) $(CXX_WARNINGS) $(CXX_OPTIM) $(CXX_DEP) $(INC_FLAG)

$(BIN_DIR)/$(TARGET): $(OBJ)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ $(ALL_LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	$(RM) -rf $(BUILD_DIR)/* $(BIN_DIR)/$(TARGET) $(BUILD_DIR) *.txt *.dat

-include $(DEP)
