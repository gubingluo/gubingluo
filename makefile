CUDA_INSTALL_PATH = /usr/local/cuda
GCC_INSTALL_PATH = /usr

NVCC = $(CUDA_INSTALL_PATH)/bin/nvcc
GCC = $(GCC_INSTALL_PATH)/bin/gcc

LDFLAGS = -L$(CUDA_INSTALL_PATH)/lib64

LIB1 = -lcudart -lcurand
LIB2 = -lfftw3 -lm
LIB3 = -lstdc++

CFILES1 = RTM_GPU_velstr.cpp
CFILES2 = function.cpp
CFILES3 = Interp_function.cpp
CFILES4 = IO_function.cpp
CFILES5 = sgyhead_function.cpp
CFILES6 = Velocity_function.cpp
CFILES7 = FD_function.cpp
CFILES8 = image_function.cpp
CUFILE1 = rtm_model.cu
CUFILE2 = rtm_real.cu
CUFILE3 = kernel.cu

OBJECTS = RTM_GPU_velstr.o rtm_model.o rtm_real.o kernel.o function.o Interp_function.o IO_function.o FD_function.o image_function.o sgyhead_function.o Velocity_function.o
EXECNAME = rtm_gpu

all:
	$(GCC) -c $(CFILES1) 
	$(GCC) -c $(CFILES2)
	$(GCC) -c $(CFILES3)
	$(GCC) -c $(CFILES4)
	$(GCC) -c $(CFILES5)
	$(GCC) -c $(CFILES6)
	$(GCC) -c $(CFILES7)
	$(GCC) -c $(CFILES8)
	$(NVCC) -c $(CUFILE1)
	$(NVCC) -c $(CUFILE2)
	$(NVCC) -c $(CUFILE3)
	$(GCC) -fopenmp -O3 -o $(EXECNAME) $(OBJECTS) $(LDFLAGS) $(LIB1) $(LIB2) $(LIB3)

clean:
	rm -rf *.o
