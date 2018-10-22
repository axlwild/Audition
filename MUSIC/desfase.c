/**

 gcc -std=gnu99 -o desfase.ej desfase.c -ljack -lfftw3 -lm
 make -C ~/Escritorio/baudline_1.08_linux_x86_64/


SEGUN ESTE EJEMPLO FOURIER SÏ DA EL COMPONENTE CENTRAL
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <jack/jack.h>

// Include FFTW header
#include <complex.h> //needs to be included before fftw3.h for compatibility
#include <fftw3.h>


fftw_plan i_forward, o_inverse;
double complex *i_fft, *i_time, *o_fft, *o_time;




void shiftInTime(double complex *buff, int buff_size, int time){
	//printf("Dentro de la funcion\n");
	//printf("Antes\n");
	for(int i = 0; i < buff_size; ++i){
		//printf("%.1f + %.1f\n", creal(buff[i]),cimag(buff[i]));	
	}
	//printf("Después\n");
	for(int i = 0; i < buff_size; ++i){
		buff[i] *= cexp(-I*2* M_PI*time*i);
	//	printf("%.1f + %.1f\n", creal(buff[i]),cimag(buff[i]));	
	}
	return;
}


int main(int argc, char const *argv[])
{
	//Para pasar a Fourier
	int buffer_size = 10;
	double complex equis = 4 + 3I;
	i_fft = (double complex *) fftw_malloc(sizeof(double complex) * buffer_size);
	i_time = (double complex *) fftw_malloc(sizeof(double complex) * buffer_size);
	o_fft = (double complex *) fftw_malloc(sizeof(double complex) * buffer_size);
	o_time = (double complex *) fftw_malloc(sizeof(double complex) * buffer_size);
	
	i_forward = fftw_plan_dft_1d(buffer_size, i_time, i_fft , FFTW_FORWARD, FFTW_MEASURE);
	o_inverse = fftw_plan_dft_1d(buffer_size, o_fft , o_time, FFTW_BACKWARD, FFTW_MEASURE);


	jack_default_audio_sample_t buff[] = {1,0,0,0,1,0,0,0,0,1};
	
	for (int i = 0; i < 10; ++i)
	{
		i_time[i] = buff[i];
		printf("%.0f",buff[i]);
	}
	printf("\n");
	fftw_execute(i_forward);
	for(int i=0 ; i<10; i++){
		o_fft[i] = i_fft[i];
		printf("%.1f + %.1f\n", creal(o_fft[i]),cimag(o_fft[i]));
	}
	printf("Antes del shift\n");
	for(int i=0 ; i<buffer_size; i++){
		printf("%.1f + %.1f\n", creal(o_fft[i]),cimag(o_fft[i]));	
	}
	shiftInTime(o_fft,buffer_size,1);
	printf("Después del shift\n");
	for(int i=0 ; i<buffer_size; i++){
		printf("%.1f + %.1f\n", creal(o_fft[i]),cimag(o_fft[i]));	
	}
	
	fftw_execute(o_inverse);
	
	for(int i=0 ; i<10; i++){
		buff[i] = abs((creal(o_time[i])/buffer_size));
		printf("%.0f ",buff[i]);
	}
	printf("\n");
	printf("%f + %f\n",creal(cexp(-I*M_PI)),cimag(cexp(-I*M_PI)));
	//printf("%.1f\n", cabs(equis));


}
















