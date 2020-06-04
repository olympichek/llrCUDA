CUDA = /usr/local/cuda
CUC = $(CUDA)/bin/nvcc
CULIB = $(CUDA)/lib64
CUINC = $(CUDA)/include
# CUFLAGS = -O$(OptLevel) --compiler-options=-Wall -I$(CUINC)
llrCUDA: lprime.o lmenu.o gwdbldbl.o gwypnum.o giants.o giantext.o
	g++ -fPIC -O3 -g -o llrCUDA lprime.o lmenu.o gwdbldbl.o giants.o giantext.o gwypnum.o -L$(CULIB) -lgmp -lcufft -lcudart -lm -lstdc++
gwypnum.o: gwypnum.cu
	$(CUC) -g -O3 -I$(CUINC) -I. -Wno-deprecated-gpu-targets gwypnum.cu -c
lprime.o: lprime.cu
	$(CUC) -g -O -I$(CUINC) -I. -Wno-deprecated-gpu-targets lprime.cu -c
lmenu.o: lmenu.cu
	$(CUC) -g -O3 -I$(CUINC) -I. -Wno-deprecated-gpu-targets lmenu.cu -c
gwdbldbl.o: gwdbldbl.cpp
	$(CUC) -g -O3 -I$(CUINC) -I. -Wno-deprecated-gpu-targets gwdbldbl.cpp -c
giants.o: giants.cu
	$(CUC) -g -O3 -I. -I/home/jpenne/include -w -Wno-deprecated-gpu-targets giants.cu -c
giantext.o: giantext.cu
	$(CUC) -g -O3 -I$(CUINC) -I. -Wno-deprecated-gpu-targets giantext.cu -c
clean:
	-rm *.o llrCUDA
