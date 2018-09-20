/*
 * A simple 1-input to 1-output JACK client.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <jack/jack.h>
#define CANALES 2

jack_port_t **input_port;
jack_port_t **output_port;
jack_client_t *client;

/**
 * The process callback for this JACK application is called in a
 * special realtime thread once for each audio cycle.
 *
 * This client does nothing more than copy data from its input
 * port to its output port. It will exit when stopped by 
 * the user (e.g. using Ctrl-C on a unix-ish operating system)
 */

//Recibe el tamaño de la ventana y cualquier argumentos (si lo requiere)
int jack_callback (jack_nframes_t nframes, void *arg){
	jack_default_audio_sample_t **in, **out;
	int i, j;
	
	in = (jack_default_audio_sample_t **) malloc (sizeof(jack_default_audio_sample_t *)*CANALES);
	out = (jack_default_audio_sample_t **) malloc (sizeof(jack_default_audio_sample_t *)*CANALES);
	/*
	//Asigna la salida a la entrada correspondiente
	for (j = 0; j < CANALES; j++) {
		in[j] = jack_port_get_buffer (input_port[j], nframes);
		out[j] = jack_port_get_buffer (output_port[j], nframes);

		for (int i = 0; i < nframes; i++) {
			out[j][i] = in[j][i];
		}
	}
	*/
	
	
	return 0;
}


/**
 * JACK calls this shutdown_callback if the server ever shuts down or
 * decides to disconnect the client.
 */
void jack_shutdown (void *arg){
	exit (1);
}


int main (int argc, char *argv[]) {
	const char *client_name = "in_to_out";		
	jack_options_t options = JackNoStartServer;
	jack_status_t status;
	int i;
	
	/* open a client connection to the JACK server */
	//Le decimos el nombre del cliente, las opciones de inicio del servidor 
	//(en este caso, le decimos que no inicie un servidor, sino que use configuraciones 
	//del servidor) ya levantado y cuál es el etado de inicio del cliente.
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
	
	/* tell the JACK server to call 'jack_callback()' whenever there is work to be done. */
	//Recibe el cliente, la funcion a hacer y el segundo parámetro de esa función.
	jack_set_process_callback (client, jack_callback, 0);
	
	
	/* tell the JACK server to call 'jack_shutdown()' if it ever shuts down,
	   either entirely, or if it just decides to stop calling us. */
	jack_on_shutdown (client, jack_shutdown, 0);
	
	
	/* display the current sample rate. */
	printf ("Engine sample rate: %d\n", jack_get_sample_rate (client));

	/* display the current window size. */
	printf ("Window size: %d\n", jack_get_buffer_size (client));
	
	input_port = (jack_port_t **) malloc(sizeof(jack_port_t *)*CANALES);
	output_port = (jack_port_t **) malloc(sizeof(jack_port_t *)*CANALES);

	char nombre[50];
	for (i = 0; i < CANALES; i++) {
		sprintf(nombre, "input_%d", i+1);
		input_port[i] = jack_port_register (client, nombre, JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput, 0);
		sprintf(nombre, "output_%d", i+1);
		output_port[i] = jack_port_register (client, nombre, JACK_DEFAULT_AUDIO_TYPE,JackPortIsOutput, 0);
	}
	
	/* check that both ports were created succesfully */
	if ((input_port == NULL) || (output_port == NULL)) {
		printf("Could not create agent ports. Have we reached the maximum amount of JACK agent ports?\n");
		exit (1);
	}
	
	/* Tell the JACK server that we are ready to roll.
	   Our jack_callback() callback will start running now. */
	if (jack_activate (client)) {
		printf ("Cannot activate client.");
		exit (1);
	}
	
	printf ("Agent activated.\n");
	
	/* Connect the ports.  You can't do this before the client is
	 * activated, because we can't make connections to clients
	 * that aren't running.  Note the confusing (but necessary)
	 * orientation of the driver backend ports: playback ports are
	 * "input" to the backend, and capture ports are "output" from
	 * it.
	 */
	printf ("Connecting ports... ");
	 
	/* Assign our input port to a server output port*/
	// Find possible output server port names
	const char **serverports_names;
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsOutput);
	if (serverports_names == NULL) {
		printf("No available physical capture (server output) ports.\n");
		exit (1);
	}

	for (i = 0; i < CANALES; i++) {
		if (jack_connect (client, serverports_names[i], jack_port_name (input_port[i]))) {
			printf("Cannot connect input port.\n");
			exit (1);
		}
	}
	// Connect the first available to our input port
	
	// free serverports_names variable for reuse in next part of the code
	free (serverports_names);
	
	
	/* Assign our output port to a server input port*/
	// Find possible input server port names
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsInput);
	if (serverports_names == NULL) {
		printf("No available physical playback (server input) ports.\n");
		exit (1);
	}
	for (int i = 0; i < CANALES; i++) {
		if (jack_connect (client, jack_port_name (output_port[i]), serverports_names[i])) {
			printf ("Cannot connect output ports.\n");
			exit (1);
		}
	}
	// free serverports_names variable, we're not going to use it again
	free (serverports_names);
	
	
	printf ("done.\n");
	/* keep running until stopped by the user */
	sleep (-1);
	
	
	/* this is never reached but if the program
	   had some other way to exit besides being killed,
	   they would be important to call.
	*/
	jack_client_close (client);
	exit (0);
}
