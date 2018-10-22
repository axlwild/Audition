/*
	Cálculo de múltiples direcciones de arribo con MUSIC
	Primero correr este:
	./ReadMicWavs music input ../corpusAIRA/corpus44100/clean-2source 2
	Luego el corpus aira:
	./ReadMicWavs music input_ rutadelwav
*/

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <jack/jack.h>
#include <Eigen/Eigen>
#include <complex.h>
#include <fftw3.h>

#define REAL 0
#define IMAG 1
#define N_MUESTRAS 20
#define FMIN 0.0
#define FMAX 4000.0

//Variables para transformada de una señal.
fftw_complex *x_t, *x_f;	//Entradas de un micrófono.
fftw_complex *y_t, *y_f;	//Entradas en otro micrófono.
fftw_plan x_forward, y_forward;


						//NMUESTRAS
fftw_complex **Xacum;	//Matriz con los acumulados de las FFTs de un micrófono.

						//NMUESTRAS
fftw_complex **Yacum;	//Matriz con los acumulados de las FFTs de otro micrófono.



//Variables para manejo de jack.
jack_port_t *input_1;
jack_port_t *input_2;
jack_port_t *output_port;
jack_client_t *client;

//Índices para las frecuencias a considerar
int i_min, i_max;

//Unidad imaginaria auxiliar.
std::complex<double> img(0.0, 1.0);

//Datos del servidor.
double sample_rate;
int nframes;

double MUSIC(Eigen::MatrixXcd X) {
	/*	PARA LAS SIGUIENTES LÍNEAS ESTOY HACIENDO LA SUPOSICIÓN DE QUE
		EL MENOR EIGENVECTOR (OBTENIDO POR EL ÍNDICE DEL MENOR EIGENVALOR)
		ES EL DE RUIDO Y QUE NO HAY OTRO.
		FALTARÍA IMPLEMENTAR QUE REALMENTE CALCULE CUÁL ES EL VECTOR DE RUIDO.
	*/

	Eigen::MatrixXcd R(2, 2);		//Matriz de covarianzas.
	Eigen::VectorXcd eval(2);		//Vector de eigen valores.
	Eigen::MatrixXcd evec(2, 2);	//Matriz de eigen vectores.
	Eigen::MatrixXcd Q(2, 1);		//Matriz con eigenvector ruidoso.
	Eigen::MatrixXcd Q_aux(2, 2);	//Resultado de multiplicar Q por su hermitiana * Q.adjoint()
	Eigen::MatrixXcd a_pot(2, 180);	//Vector de direcciones probables.
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(R);	//Para sacar eigens
	int i_min_evec;					//Almacena el índice el menor eigenvalor.
	double angulo_max;				//Almacena el ángulo de dirección más probable.
	double valor, valor_max = 0;	//Variables auxiliares para comparar valores de MUSIC.

	R = (X * X.adjoint())/N_MUESTRAS;		//Se obtiene la matriz de covarianza R
	
	eval = es.eigenvalues();				//Se obtienen los eigen valores de R.
	evec = es.eigenvectors();
	if (eval(0).real() > eval(1).real())	//Se obtiene el índice del mínimo eigenvalor
		i_min_evec = 0;
	else
		i_min_evec = 1;

	Q(0, 0) = evec(REAL, i_min_evec);		//Se asigna la parte real del eigenvector de ruido.
	Q(1, 0) = evec(IMAG, i_min_evec);		//Se asigna la parte imaginaria del eigenvector de ruido.

	a_pot = Eigen::MatrixXcd::Ones(2, 180);	//Se crea la matriz de direcciones potenciales.
	//El primer renglón es de puros unos, el segundo es de direcciones potenciales.
	for (int i = -90; i < 90; i++){
		a_pot(1, i+90) = exp(img*2.0*M_PI*sin(i*M_PI/180));
	}

	Q_aux = Q * Q.adjoint();				//Se calcula Q * Q' para acelerar cálculos.

	//MUSIC
	for (int i = 0; i < 180; i++){
		//Se obtiene cada valor del espectrograma de music.
		valor = real(1.0 / (a_pot.col(i).adjoint() * Q_aux * a_pot.col(i)) (0,0) );
		//std::cout << valor << std::endl;
		if (valor > valor_max && valor >= 20.0) {
			//Se obtiene el ángulo correspondiente al mayor valor que se haya encontrado.
			valor_max = valor;
			angulo_max = i-90;
		}
	}
	//std::cout << valor_max << " - " << angulo_max << std::endl;
	return angulo_max;
}

int jack_callback (jack_nframes_t nframes, void *arg){
	jack_default_audio_sample_t *x, *y, *out;	//Variables para entrada y salida de señales.
	Eigen::MatrixXcd X(2, N_MUESTRAS);			//Matriz X con señales transformadas.
	int angulo_max = 0, angulo;

	int i, freq, t;

	x = (jack_default_audio_sample_t *) jack_port_get_buffer (input_1, nframes);
	y = (jack_default_audio_sample_t *) jack_port_get_buffer (input_2, nframes);
	out = (jack_default_audio_sample_t *)jack_port_get_buffer (output_port, nframes);

	//Se recorren los datos de las matrices de acumulados un renglón hacia abajo.
	
	/*
	
	for (t = N_MUESTRAS-1; t > 0; t--) {
		for (freq = 0; freq < i_max; freq++) {
			Xacum[t][freq][REAL] = Xacum[t-1][freq][REAL];
			Xacum[t][freq][IMAG] = Xacum[t-1][freq][IMAG];
			Yacum[t][freq][REAL] = Yacum[t-1][freq][REAL];
			Yacum[t][freq][IMAG] = Yacum[t-1][freq][IMAG];
		}
	}
	
	*/

	
	for (t = N_MUESTRAS-1; t > 0; t--) {	
		Xacum[t] = Xacum[t-1];
		Yacum[t] = Yacum[t-1];
	}
	
	

	
	//Se preparan las señales de entrada para hacerles FFT
	for (i = 0; i < nframes; i++) 
		x_t[i][REAL] = x[i];
	for (i = 0; i < nframes; i++) 
		y_t[i][REAL] = y[i];

	//Se hace FFT sobre las señales de entrada
	fftw_execute(x_forward);
	fftw_execute(y_forward);

	//Se pasan las señales actuales en frecuencia a las matrices de acumulados.
	for (freq = 0; freq < i_max; freq++) {
		Xacum[0][freq][REAL] = x_f[freq][REAL];
		Xacum[0][freq][IMAG] = x_f[freq][IMAG];
		Yacum[0][freq][REAL] = y_f[freq][REAL];
		Yacum[0][freq][IMAG] = y_f[freq][IMAG];
	}

	//Para cada columna (frecuencia) de la matriz de acumulados.
	for (freq = 0; freq < i_max; freq++) {
		//Se pasa la columna a la matriz de datos como un renglón.
		for (t = 0; t < N_MUESTRAS; t++) {
			X(0, t) = std::complex<double>(Xacum[t][freq][REAL], Xacum[t][freq][IMAG]);
			X(1, t) = std::complex<double>(Yacum[t][freq][REAL], Yacum[t][freq][IMAG]);			
		}
		//Se obtiene el mayor ángulo por MUSIC.
		angulo = MUSIC(X);
		if (angulo > angulo_max)
			angulo_max = angulo;
	}

	//Se imprime el ándulo sólo si es diferente de 0. Si es 0, no detectó nada.
	if (angulo_max != 0)
		std::cout << angulo_max << std::endl;

	//Se manda la señal hacia la salida.
	for (i = 0; i < nframes; i++) 
		out[i] = x[i];

	return 0;
}

void jack_shutdown (void *arg){
	exit (1);
}


int main (int argc, char *argv[]) {
	const char *client_name = "music";
	jack_options_t options = JackNoStartServer;
	jack_status_t status;
	
	/* open a client connection to the JACK server */
	client = jack_client_open (client_name, options, &status);
	if (client == NULL){
		/* if connection failed, say why */
		printf ("jack_client_open() failed, status = 0x%2.0x\n", status);
		if (status & JackServerFailed) {
			printf ("Unable to connect to JACK server.\n");
		}
		exit (1);
	}
	
	/* if connection was successful, check if the name we proposed is not in use */
	if (status & JackNameNotUnique){
		client_name = jack_get_client_name(client);
		printf ("Warning: other agent with our name is running, `%s' has been assigned to us.\n", client_name);
	}
	
	jack_set_process_callback (client, jack_callback, 0);
	jack_on_shutdown (client, jack_shutdown, 0);
	
	printf ("Sample rate: %d\n", jack_get_sample_rate (client));
	printf ("Window size: %d\n", jack_get_buffer_size (client));
	sample_rate = (double) jack_get_sample_rate(client);
	nframes = jack_get_buffer_size (client);

	//Se calculan los índices de las frecuencias a analizar:
	i_min = (int) FMIN/(sample_rate/(nframes));
	i_max = (int) FMAX/(sample_rate/(nframes));
	
	//Reservando buffers para FFT
	x_t = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nframes);
	y_t = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nframes);
	x_f = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nframes);
	y_f = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nframes);


	Xacum = (fftw_complex **) fftw_malloc(sizeof(fftw_complex*)*N_MUESTRAS);
	Yacum = (fftw_complex **) fftw_malloc(sizeof(fftw_complex*)*N_MUESTRAS);

	//Reservando matrices de acumulados
	for (int i = 0; i < N_MUESTRAS; i++) {
		Xacum[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * i_max);
		Yacum[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * i_max);
	}
	
	x_forward = fftw_plan_dft_1d(nframes, x_t, x_f, FFTW_FORWARD, FFTW_MEASURE);
	y_forward = fftw_plan_dft_1d(nframes, y_t, y_f, FFTW_FORWARD, FFTW_MEASURE);

	/* create the agent input ports */
	input_1 = jack_port_register (client, "input_1", JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput, 0);
	input_2 = jack_port_register (client, "input_2", JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput, 0);

	
	/* create the agent output ports */
	output_port = jack_port_register (client, "output_port", JACK_DEFAULT_AUDIO_TYPE,JackPortIsOutput, 0);
	
	/* check that both ports were created succesfully */
	if (input_1 == NULL || input_2 == NULL || output_port == NULL) {
		printf("Could not create agent ports. Have we reached the maximum amount of JACK agent ports?\n");
		exit (1);
	}
	
	if (jack_activate (client)) {
		printf ("Cannot activate client.");
		exit (1);
	}
	
	printf ("Agent activated.\n");

	printf ("Connecting ports... ");
	const char **serverports_names;

	//ENTRADAS
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsOutput);
	if (serverports_names == NULL) {
		printf("No available physical capture (server output) ports.\n");
		exit (1);
	}
	if (jack_connect (client, serverports_names[0], jack_port_name (input_1)) || jack_connect (client, serverports_names[1], jack_port_name (input_2))) {
		printf("Cannot connect input port.\n");
		exit (1);
	}
	free (serverports_names);
	
	//SALIDAS
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsInput);
	if (serverports_names == NULL) {
		printf("No available physical playback (server input) ports.\n");
		exit (1);
	}
	if (jack_connect (client, jack_port_name (output_port), serverports_names[0])) {
		printf ("Cannot connect output ports.\n");
		exit (1);
	}
	free (serverports_names);
	
	printf ("done.\n");

	sleep (-1);
	
	jack_client_close (client);
	exit (0);
}
