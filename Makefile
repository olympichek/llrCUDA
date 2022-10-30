CUDA = /usr/local/cuda
CUC = $(CUDA)/bin/nvcc
CULIB = $(CUDA)/lib64
CUINC = $(CUDA)/include
# CUFLAGS = -O$(OptLevel) --compiler-options=-Wall -I$(CUINC)
llrCUDA: lprime.o lmenu.o gwdbldbl.o gwypnum.o giants.o giantext.o
	$(CUC) -O3 -arch sm_61 -o llrCUDA lprime.o lmenu.o gwdbldbl.o giants.o giantext.o gwypnum.o -L$(CULIB) -lgmp -lcufft -lcudart -lm -lstdc++
sllrCUDA: lprime.o lmenu.o gwdbldbl.o gwypnum.o giants.o giantext.o
	$(CUC) -O3  -arch sm_61 lprime.o lmenu.o gwdbldbl.o giants.o giantext.o gwypnum.o ./libcudart_static.a ./libcufft_static.a ./libculibos.a -L. -lgmp  -lm -lstdc++ -o sllrCUDA
gwypnum.o: gwypnum.cu
	$(CUC) -O3 -I$(CUINC) -I. -Wno-deprecated-gpu-targets gwypnum.cu -c
lprime.o: lprime.cu
	$(CUC) -O3 -I$(CUINC) -I. -Wno-deprecated-gpu-targets lprime.cu -c
lmenu.o: lmenu.cu
	$(CUC) -O3 -I$(CUINC) -I. -Wno-deprecated-gpu-targets lmenu.cu -c
gwdbldbl.o: gwdbldbl.cpp
	$(CUC) -O3 -I$(CUINC) -I. -Wno-deprecated-gpu-targets gwdbldbl.cpp -c
giants.o: giants.cu
	$(CUC) -O3 -I. -I/home/jpenne/include -w -Wno-deprecated-gpu-targets giants.cu -c
giantext.o: giantext.cu
	$(CUC) -O3 -I$(CUINC) -I. -Wno-deprecated-gpu-targets giantext.cu -c
clean:
	-rm *.o llrCUDA sllrCUDA
