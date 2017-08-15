
# Debug build flags
ifeq ($(dbg),1)
      NVCCFLAGS += -g -G
      TARGET := debug
else
      TARGET := release
endif

ALL_CCFLAGS :=
ALL_CCFLAGS += $(NVCCFLAGS)
ALL_CCFLAGS += $(EXTRA_NVCCFLAGS)
ALL_CCFLAGS += $(addprefix -Xcompiler ,$(CCFLAGS))
ALL_CCFLAGS += $(addprefix -Xcompiler ,$(EXTRA_CCFLAGS))

ALL_LDFLAGS :=
ALL_LDFLAGS += $(ALL_CCFLAGS)
ALL_LDFLAGS += $(addprefix -Xlinker ,$(LDFLAGS))
ALL_LDFLAGS += $(addprefix -Xlinker ,$(EXTRA_LDFLAGS))


# Gencode arguments
ifeq ($(OS_ARCH),armv7l)
SMS ?= 20 30 32 35 37 50
else
SMS ?= 11 20 30 35 37 50
endif

ifeq ($(SMS),)
$(info >>> WARNING - no SM architectures have been specified - waiving sample <<<)
SAMPLE_ENABLED := 0
endif

ifeq ($(GENCODE_FLAGS),)
# Generate SASS code for each SM architecture listed in $(SMS)
$(foreach sm,$(SMS),$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))
endif

# Generate PTX code from the highest SM architecture in $(SMS) to guarantee forward-compatibility
HIGHEST_SM := $(lastword $(sort $(SMS)))
ifneq ($(HIGHEST_SM),)
GENCODE_FLAGS += -gencode arch=compute_$(HIGHEST_SM),code=compute_$(HIGHEST_SM)
endif

CC = nvcc

INCLUDES = -Iopenssl_built/include

LFLAGS = -L/openssl_built/lib

LIBS = -lcrypto -lssl

SRCS = main.cu GCD.cu

OBJS = $(SRCS:.cu=.o)

MAIN = GCD_finder


.PHONY: depend clean

all:    $(MAIN)
	@echo  Program has been compiled

main.o: main.cu
	$(CC) $(NVCCFLAGS) $(INCLUDES) $(ALL_LDFLAGS) $(GENCODE_FLAGS) -c $<  -o $@

GCD.o: GCD.cu
	$(CC) $(NVCCFLAGS) $(INCLUDES) $(ALL_LDFLAGS) $(GENCODE_FLAGS) -c $<  -o $@

$(MAIN): main.o GCD.o
	$(CC) $(NVCCFLAGS) $(INCLUDES) $(GENCODE_FLAGS) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

run: build
	./$(MAIN)

clean:
	$(RM) *.o *~ $(MAIN)

