echo compiling..
mpic++ -O2 -I/work/pi_gkhanna_uri_edu/gkhanna_uri_edu/include -c main.c proj.c
ar r lib.a main.o proj.o
rm main.o proj.o

nvcc -O2 -c -w body.cu -I./inc -I/work/pi_gkhanna_uri_edu/gkhanna_uri_edu/include \
 -I/$OPENMPI_ROOT/include \
   -gencode arch=compute_70,code=sm_70 --fmad=false

ar r lib.a body.o
rm body.o

echo linking...
mpic++ -o runme.bin\
	lib.a\
	-L/$CUDA_ROOT/lib64 \
        -L/work/pi_gkhanna_uri_edu/gkhanna_uri_edu/lib -lqd  -lcudart

echo cleaning.. 
rm lib.a
