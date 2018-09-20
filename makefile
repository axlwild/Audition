all:
	gcc -std=gnu99 -o jack_filtrado jack_fft_filtrado.c -ljack -lfftw3 -lm
#	gcc -std=gnu99 -o hann_clase hann_clase.c -ljack -lfftw3 -lm

in_to_out:
	gcc -std=gnu99 -o jack_in_to_out jack_in_to_out.c -ljack -lm
hann:
	gcc -std=gnu99 -o hann_clase hann_clase.c -ljack -lfftw3 -lm

#Transformada de Fourier sin hacer procesamiento
fft:
	gcc -std=gnu99 -o fourier jack_fft.c -ljack -lfftw3 -lm
#Eigenvalores
eigen:
	g++ -I/usr/include/eigen3 eigen_examples.cpp -o eigen_examples
#readmicwavs
readmicwavs:
	g++ -g -O2 -I/usr/include/ ReadMicWavs.cpp -o ReadMicWavs -ljack -lsndfile

