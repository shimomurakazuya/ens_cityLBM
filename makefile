# Usage 
# export DEVICE={bdw,p100,A64FX}
# In case DEVICE=bdw, it is needed to export OMPI_CXX=icpc
# make
ifeq (,$(DEVICE))
$(error "\$(DEVICE) is not defined. please set, e.g., export DEVICE=v100")
endif
DEVICES = $(DEVICE)
## note: compiler infomation of $(DEVICES) is in `config/devices/*`
include config/devices/$(DEVICES)

# problem selection
TEST ?= oklahoma
#TEST ?= cavity
## note: defined macros of the test are in `config/tests/*`
include config/tests/$(TEST)

# Source directories
OUT_DIR         := bin
SRC_DIRS        := src include validation_srcs function_srcs index_srcs library_srcs initialize_srcs defines util
OBJ_DIR         := obj
RESULT_DIR      := io
RUN_DIR         := run

# Finalize $(CXXFLAGS)
ifeq ("$(COMPILER)","nvcc")
citylbm_addopt += USE_NVCC
else
## revert GPU-featured macros
citylbm_addopt := $(filter-out NO_FIELD_WRITEDAT, $(citylbm_addopt))
endif

# g++ / icpc
MPICXX_CXX ?= $(firstword $(shell mpicxx -show))
ifneq (, $(findstring g++, $(MPICXX_CXX)))
MPICXXFLAGS_EXTRA += -fopenmp
else
ifneq (, $(findstring nvc++, $(MPICXX_CXX)))
MPICXXFLAGS_EXTRA += -fopenmp
else
ifneq (, $(findstring hipcc, $(MPICXX_CXX)))
MPICXXFLAGS_EXTRA += -fopenmp
else
ifneq (, $(findstring icpx,$(MPICXX_CXX)))
MPICXXFLAGS_EXTRA += -qopenmp
MPICXXFLAGS_EXTRA += -ipo
else
ifneq (, $(findstring icpc,$(MPICXX_CXX)))
MPICXXFLAGS_EXTRA += -qopenmp
MPICXXFLAGS_EXTRA += -ipo -qopt-report=5
#MPICXXFLAGS_EXTRA += -xHost
else
$(error "unknown compiler: $(MPICXX_CXX)")
endif
endif
endif
endif
endif

ifeq ("$(COMPILER)","nvcc")
CXXFLAGS += -Xcompiler "$(MPICXXFLAGS_EXTRA)"
LDFLAGS  += -Xcompiler "$(MPICXXFLAGS_EXTRA)"
else
CXXFLAGS += $(MPICXXFLAGS_EXTRA)
LDFLAGS  += $(MPICXXFLAGS_EXTRA)
endif

# Finalize CXXFLAGS
CXXFLAGS += $(addprefix -D,$(citylbm_addopt))
CXXFLAGS += $(addprefix -I,$(SRC_DIRS))
EXTRA += submodule

# Disabled files
NO_MAKE_LIST += $(filter-out $(VVTARGET), $(basename $(notdir $(wildcard validation_srcs/*.cpp))))
NO_MAKE_LIST += CreateTreeMPI
NO_MAKE_LIST += _tmp_FuncMapData
$(info NO_MAKE_LIST=$(NO_MAKE_LIST))

# Source files
ALL_SRCS := $(basename $(notdir $(wildcard $(addsuffix /*.cu, $(SRC_DIRS))) $(wildcard $(addsuffix /*.cpp, $(SRC_DIRS))) ))
SRCS := $(filter-out $(NO_MAKE_LIST), $(ALL_SRCS))
OBJS = $(addprefix $(OBJ_DIR)/, $(addsuffix .o, $(SRCS)))

# Add paths
vpath %.cu  $(SRC_DIRS)
vpath %.cpp $(SRC_DIRS)
vpath %.h   $(SRC_DIRS)
vpath %.hpp $(SRC_DIRS)

.PHONY: all build flaginfo clean tagfiles install_pbvr_files submodule

# Targets
all: $(TARGET) $(EXTRA)

build: $(TARGET)

$(TARGET): $(OBJS)
	@mkdir -p $(OUT_DIR)
	$(LINKER) -o $@ $^ $(LDFLAGS)

# Compiler .cpp files to create object files
$(OBJ_DIR)/%.o : %.cpp
	@[ -d $(OBJ_DIR) ]|| mkdir -p $(OBJ_DIR)
	$(COMPILER) $(CXXFLAGS) -o $@ -c $<

$(OBJ_DIR)/%.o : %.cu
	@[ -d $(OBJ_DIR) ]|| mkdir -p $(OBJ_DIR)
	$(COMPILER) $(CXXFLAGS) -o $@ -c $<

clean:
	find $(OUT_DIR) -name "*.exe"     | xargs -r rm -fv
	find $(OUT_DIR) -name "*.optrpt"  | xargs -r rm -fv
	find $(OUT_DIR) -name "citylbm.*" | xargs -r rm -fv
	find $(OBJ_DIR) -name "*.o"       | xargs -r rm -fv
	find $(OBJ_DIR) -name "*.optrpt"  | xargs -r rm -fv
	find $(OBJ_DIR) -name "*.d"       | xargs -r rm -fv
	find ./  -name "*.lst"            | xargs -r rm -fv

resultclean:
	find $(RESULT_DIR)/output -name "*.csv"                 | xargs -r rm -fv
	find $(RESULT_DIR)/output -name "MLUPS.txt"             | xargs -r rm -fv
	find $(RESULT_DIR)/output -mindepth 2 -type d           | xargs -r rmdir
	find $(RESULT_DIR)/output -mindepth 1 -type d           | xargs -r rmdir
	find $(RESULT_DIR)/iofiles/restart -name "timestep.txt" | xargs -r rm -fv
	find $(RESULT_DIR)/iofiles/restart -name "*.dat"        | xargs -r rm -fv
	find $(RESULT_DIR)/iofiles/restart -name "*.dat.gz"     | xargs -r rm -fv
	find $(RESULT_DIR)/iofiles/restart -mindepth 1 -type  d | xargs -r rmdir
	find $(RESULT_DIR)/iofiles -name "*.dat"                | xargs -r rm -fv
	find $(RESULT_DIR)/iofiles -name "*.h5"                 | xargs -r rm -fv
	find $(RESULT_DIR)/iofiles -name "*.xmf"                | xargs -r rm -fv
	find $(RESULT_DIR)/iofiles -name "*.vtu"                | xargs -r rm -fv
	find $(RESULT_DIR)/iofiles -name "*.pvtu"               | xargs -r rm -fv
	rm -fv $(RUN_DIR)/core

tagfiles:
	ctags -R --langmap=c:+.hpp --langmap=c:+.cu $(SRC_DIRS)

submodule:
	git submodule init
	git submodule update

