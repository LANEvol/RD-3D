CUDA_PATH    ?= /usr/local/cuda


GCC := g++

NVCC := $(CUDA_PATH)/bin/nvcc -ccbin $(GCC)





# Common includes and paths for CUDA
INCLUDES  := -I $(CUDA_PATH)/samples/common/inc 
LIBRARIES :=  -L /usr/lib -L /usr/local/lib -lcurand

# Target rules
all: build

build: 3D-RD

3D-RD.o:3D-RD.cu
	$(EXEC) $(NVCC) $(INCLUDES)  -std=c++11 -gencode=arch=compute_70,code=compute_70 --prec-div=true --prec-sqrt=true --fmad=false  -o $@ -c $< 

3D-RD: 3D-RD.o
	$(EXEC) $(NVCC) -gencode=arch=compute_70,code=compute_70 --prec-div=true --prec-sqrt=true --fmad=false -o $@ $+ $(LIBRARIES)

run: build
	$(EXEC) ./3D-RD

clean:
	rm -f 3D-RD 3D-RD.o

clobber: clean
