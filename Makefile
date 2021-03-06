################################################################################
#
# Copyright 1993-2013 NVIDIA Corporation.  All rights reserved.
#
# NOTICE TO USER:
#
# This source code is subject to NVIDIA ownership rights under U.S. and
# international Copyright laws.
#
# NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
# CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
# IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH
# REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
# IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL,
# OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
# OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
# OR PERFORMANCE OF THIS SOURCE CODE.
#
# U.S. Government End Users.  This source code is a "commercial item" as
# that term is defined at 48 C.F.R. 2.101 (OCT 1995), consisting  of
# "commercial computer software" and "commercial computer software
# documentation" as such terms are used in 48 C.F.R. 12.212 (SEPT 1995)
# and is provided to the U.S. Government only as a commercial end item.
# Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through
# 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the
# source code with only those rights set forth herein.
#
################################################################################
#
# Makefile project only supported on Mac OS X and Linux Platforms)
#
################################################################################

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

GCC = gcc 

CC = nvcc -ccbin $(GCC)

INCLUDES = -Iopenssl_built/include

LFLAGS = -Lopenssl_built/lib

LIBS = -lcrypto -lssl

MAIN_FILE = main

SRCS =  $(MAIN_FILE).cu cuda_bignum.cu test.cu files_manager.cu device_cuda_bignum.cu

OBJS = $(SRCS:.cu=.o)

MAIN = GCD_RSA


.PHONY: depend clean

all:    $(MAIN)
	@echo  Program has been compiled

$(MAIN_FILE).o: $(MAIN_FILE).cu
	$(CC) $(NVCCFLAGS) $(INCLUDES) $(ALL_LDFLAGS) $(GENCODE_FLAGS) -c $<  -o $@

cuda_bignum.o: cuda_bignum.cu
	$(CC) $(NVCCFLAGS) $(INCLUDES) $(ALL_LDFLAGS) $(GENCODE_FLAGS) -c $<  -o $@

device_cuda_bignum.o: device_cuda_bignum.cu
	$(CC) $(NVCCFLAGS) $(INCLUDES) $(ALL_LDFLAGS) $(GENCODE_FLAGS) -c $<  -o $@

test.o: test.cu
	$(CC) $(NVCCFLAGS) $(INCLUDES) $(ALL_LDFLAGS) $(GENCODE_FLAGS) -c $<  -o $@

files_manager.o: files_manager.cu
	$(CC) $(NVCCFLAGS) $(INCLUDES) $(ALL_LDFLAGS) $(GENCODE_FLAGS) -c $<  -o $@

$(MAIN): $(MAIN_FILE).o test.o cuda_bignum.o files_manager.o device_cuda_bignum.o
	$(CC) $(NVCCFLAGS) $(INCLUDES) $(GENCODE_FLAGS) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

run: build
	./$(MAIN)

clean:
	$(RM) *.o *.a *~ $(MAIN)

