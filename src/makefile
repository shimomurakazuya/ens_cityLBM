PROGRAM_NAME := run.exe


OUT_DIR         := ../bin
SOURCE_DIR      := ./
HEADER_DIR      := ../include
OBJ_DIR         := ../obj
TEST_DIR        := ../test_srcs
VALIDATION_DIR  := ../validation_srcs
FUNCTION_DIR    := ../function_srcs
INDEX_DIR       := ../index_srcs
LIBRARY_SRC_DIR := ../library_srcs
INITIALIZE_DIR  := ../initialize_srcs
DEFINE_DIR      := ../defines
RESULT_DIR      := ../io
RUN_DIR         := ../run

CPP  = icpc -O3  -std=c++14 -fopenmp
NVCC = nvcc -O3 -DUSE_NVCC -std=c++14 -restrict -Xcompiler -fopenmp #--expt-relaxed-constexpr
#NVCC = nvcc -O3  -std=c++11 -restrict -Xcompiler -fopenmp #--expt-relaxed-constexpr


# CFLAGS
CFLAGS   = -use_fast_math # --ptxas-options=-v
#CFLAGS  += -I/home/8/hpci000852/software/hdf5/1.8.11/include # HDF5 tsubame #
#CFLAGS  += -I/apps/t3/sles12sp2/free/hdf5/1.10.1/gnu/include
#CFLAGS  += -I/home/g5/a170085/software/hdf5/1.8.11/include # HDF5 icex

#CFLAGS  += -I/lustre/app/acc/cuda/9.0.176/include # nvcc reedbush
#CFLAGS  += -I/lustre/app/mvapich2-gdr/ofed4.2/gnu/include # MPI reedbush
#CFLAGS  += -I/apps/t3/sles12sp2/free/openmpi/1.10.7/gnu-cuda8.0/include # MPI tsubame gcc
#CFLAGS  += -I/apps/mvapich2-gdr/2.3/gcc4.8.5_cuda9.2/include # MPI abci
#CFLAGS  += -I/apps/openmpi/2.1.5/gcc4.8.5_cuda9.2.88.1/include # MPI abci
CFLAGS  += -I/opt/mvapich2/2.3/gcc-7.3.0.lp/include # MPI DGX2
#CFLAGS  += -I/opt/openmpi/2.1.5/gcc-7.3.0.lp/include # MPI DGX2
#CFLAGS  += -I/apps/t3/sles12sp2/free/openmpi/2.1.2/gnu/include # MPI tsubame gcc
#CFLAGS  += -I/apps/t3/sles12sp2/isv/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/include # MPI tsubame icc
#CFLAGS  += -I/home/app/intel/compilers_and_libraries_2017.2.174/linux/mpi/intel64/include # MPI icex #

CFLAGS += --default-stream per-thread
CFLAGS += -Xptxas -v -lineinfo
#CFLAGS += -maxrregcount 31
#CFLAGS += -maxrregcount 63
CFLAGS += -maxrregcount 127
CFLAGS += --gpu-architecture=sm_70
CFLAGS += --nvlink-options -Werror


# LDFLAGS
LDFLAGS  = -lm -lstdc++
LDFLAGS += -Xcompiler -fopenmp
#LDFLAGS += -lhdf5 -lhdf5_hl -L/home/8/hpci000852/software/hdf5/1.8.11/lib # HDF5 tsubame #
#LDFLAGS += -lhdf5 -lhdf5_hl -L/apps/t3/sles12sp2/free/hdf5/1.10.1/gnu/lib # HDF5 tsubame #
#LDFLAGS += -lhdf5 -lhdf5_hl -L/home/g5/a170085/software/hdf5/1.8.11/lib # HDF5 icex #


#LDFLAGS += -L/lustre/app/mvapich2-gdr/ofed4.2/gnu/lib64 # MPI reedbush
#LDFLAGS += -L/lustre/app/mvapich2-gdr/ofed4.2/cuda9.0/gdrcopy/gnu/lib64 # MPI reedbush
#LDFLAGS += -L/lustre/app/pgi/linux86-64/17.1/lib # MPI reedbush
#LDFLAGS += -L/apps/t3/sles12sp2/free/openmpi/1.10.7/gnu-cuda8.0/lib # MPI tsubame gcc
#LDFLAGS += -L/apps/mvapich2-gdr/2.3/gcc4.8.5_cuda9.2/lib64 # MPI abci
#LDFLAGS += -L/apps/openmpi/2.1.5/gcc4.8.5_cuda9.2.88.1/lib # MPI abci
LDFLAGS += -L/opt/mvapich2/2.3/gcc-7.3.0.lp/lib # MPI DGX2
#LDFLAGS += -L/opt/openmpi/2.1.5/gcc-7.3.0.lp/lib # MPI DGX2
#LDFLAGS += -L/apps/t3/sles12sp2/free/openmpi/2.1.2/gnu/lib # MPI tsubame gcc
#LDFLAGS += -L/apps/t3/sles12sp2/isv/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/lib # MPI tsubame icc
#LDFLAGS += -L/home/app/intel/compilers_and_libraries_2017.2.174/linux/mpi/intel64/lib # MPI icex #
#LDFLAGS += -lmpi -lmpicxx # reedbush
LDFLAGS += -lmpi # abci
#LDFLAGS += -lmpi -lmpi_cxx # tsubame


#LDFLAGS += -L/lustre/app/acc/cuda/9.0.176/lib64 # nvcc reedbush
#LDFLAGS += -L/apps/t3/sles12sp2/cuda/8.0.61/lib64 # nvcc tsubame
#LDFLAGS += -L/apps/cuda/9.2.88.1/lib64 # nvcc abci
LDFLAGS += -L/opt/cuda/9.2/lib64 # nvcc DGX2
LDFLAGS += -lcudart
LDFLAGS += -rdc=true
LDFLAGS += -Xnvlink -w
LDFLAGS += --gpu-architecture=sm_70



LIBS     =
INCLUDE  = -I$(HEADER_DIR)
INCLUDE += -I$(DEFINE_DIR)
INCLUDE += -I$(FUNCTION_DIR)
INCLUDE += -I$(INDEX_DIR)
INCLUDE += -I$(LIBRARY_SRC_DIR)
INCLUDE += -I$(INITIALIZE_DIR)
INCLUDE += -I$(TEST_DIR)
INCLUDE += -I$(VALIDATION_DIR)
TARGET   = $(OUT_DIR)/$(PROGRAM_NAME)


ifeq "$(strip $(OBJ_DIR))" ""
  OBJ_DIR= .
endif
CPP_SOURCES  := $(notdir $(wildcard $(SOURCE_DIR)/*.cpp)) $(notdir $(wildcard $(FUNCTION_DIR)/*.cpp)) $(notdir $(wildcard $(INDEX_DIR)/*.cpp)) $(notdir $(wildcard $(LIBRARY_SRC_DIR)/*.cpp)) $(notdir $(wildcard $(INITIALIZE_DIR)/*.cpp)) $(notdir $(wildcard $(TEST_DIR)/*.cpp)) $(notdir $(wildcard $(VALIDATION_DIR)/*.cpp))
NVCC_SOURCES := $(notdir $(wildcard $(SOURCE_DIR)/*.cu)) $(notdir $(wildcard $(FUNCTION_DIR)/*.cu))  $(notdir $(wildcard $(INDEX_DIR)/*.cu)) $(notdir $(wildcard $(LIBRARY_SRC_DIR)/*.cu)) $(notdir $(wildcard $(INITIALIZE_DIR)/*.cu)) $(notdir $(wildcard $(TEST_DIR)/*.cu)) $(notdir $(wildcard $(VALIDATION_DIR)/*.cu))
OBJECTS  = $(addprefix $(OBJ_DIR)/, $(CPP_SOURCES:.cpp=.o) $(NVCC_SOURCES:.cu=.o))

DEPENDS  = $(OBJECTS:.o=.d)
vpath %.cu  $(SOURCE_DIR)
vpath %.cpp $(SOURCE_DIR)
vpath %.h   $(HEADER_DIR)
vpath %.hpp $(HEADER_DIR)
vpath %.h   $(DEFINE_DIR)
vpath %.cu  $(FUNCTION_DIR)
vpath %.cpp $(FUNCTION_DIR)
vpath %.h   $(FUNCTION_DIR)
vpath %.hpp $(FUNCTION_DIR)
vpath %.cu  $(INDEX_DIR)
vpath %.cpp $(INDEX_DIR)
vpath %.h   $(INDEX_DIR)
vpath %.hpp $(INDEX_DIR)
vpath %.cu  $(LIBRARY_SRC_DIR)
vpath %.cpp $(LIBRARY_SRC_DIR)
vpath %.h   $(LIBRARY_SRC_DIR)
vpath %.hpp $(LIBRARY_SRC_DIR)
vpath %.cu  $(INITIALIZE_DIR)
vpath %.cpp $(INITIALIZE_DIR)
vpath %.h   $(INITIALIZE_DIR)
vpath %.hpp $(INITIALIZE_DIR)
vpath %.cu  $(TEST_DIR)
vpath %.cpp $(TEST_DIR)
vpath %.h   $(TEST_DIR)
vpath %.hpp $(TEST_DIR)
vpath %.cu  $(VALIDATION_DIR)
vpath %.cpp $(VALIDATION_DIR)
vpath %.h   $(VALIDATION_DIR)
vpath %.hpp $(VALIDATION_DIR)


$(TARGET): $(OBJECTS) $(LIBS)
	@mkdir -p $(OUT_DIR)
	$(NVCC) -o $@ $^ $(LDFLAGS)


$(OBJ_DIR)/%.o: %.cpp
	@[ -d $(OBJ_DIR) ]|| mkdir -p $(OBJ_DIR)
	$(CPP) $(CFLAGS) $(INCLUDE) -o $@ -c $<


$(OBJ_DIR)/%.o: %.cu
	@[ -d $(OBJ_DIR) ]|| mkdir -p $(OBJ_DIR)
	$(NVCC) $(CFLAGS) $(INCLUDE) -o $@ -c $<


all: $(TARGET)

-include $(DEPENDS)


clean:
	rm -fr $(OUT_DIR)/*.exe
	rm -fr $(OUT_DIR)/*.optrpt
	rm -fr $(OBJ_DIR)/*.o
	rm -fr $(OBJ_DIR)/*.optrpt
	rm -rf $(DEPEND_DIR)

resultclean:
	rm -fr $(OUT_DIR)/*.exe
	rm -fr $(OUT_DIR)/*.optrpt
	rm -fr $(OBJ_DIR)/*.o
	rm -fr $(OBJ_DIR)/*.optrpt
	rm -rf $(DEPEND_DIR)
	rm -rf $(RESULT_DIR)/output/*
	rm -rf $(RESULT_DIR)/iofiles/*.dat
	rm -rf $(RESULT_DIR)/iofiles/*.h5
	rm -rf $(RESULT_DIR)/iofiles/*.xmf
	rm -rf $(RESULT_DIR)/iofiles/*.vtu
	rm -rf $(RESULT_DIR)/iofiles/*.pvtu
	rm -rf $(RUN_DIR)/core

tagfiles:
	ctags -R --langmap=c:+.hpp --langmap=c:+.cu --exclude=bin --exclude=recodes --exclude=log  --exclude=obj  --exclude=io  --exclude=run ../*

