#pragma once
#include <time.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_set>
#include "Graph.h"
#include "DisjointSet.h"
#include "neighborhood_fns.h"

typedef unsigned long long int ulongint;

const int max_iterations = 1000;
const double tolerance = 0.000001;

int MessagePassing(std::string fname, int r, double z_real_min, double z_real_max,
		double z_imag, double delta_z, std::string outfname){
	// Open output file
	std::ofstream outfile;
	outfile.open(outfname);

	std::srand(time(NULL));

	//std::complex<double> hi(0.5,-10.1);
	std::complex<double> message_initial_value(0.5,-10.1);
	std::complex<double> z(z_real_min,z_imag);

	Graph network(fname);
	ulongint n = network.number_of_nodes();

	std::cout << "Loaded a network with " << n << " nodes and ";
	std::cout << network.number_of_edges() << " edges.\n";

	Neighborhood *Fnode;
	Fnode = new Neighborhood[n];
	init_node_F(Fnode, network, r);

	// Count number of messages:
	ulongint number_of_messages = 0;
	for (ulongint i=0; i<n; ++i) number_of_messages += Fnode[i].num_nodes;

	// initialize H (the message vector), and message lookups
	std::vector<std::complex<double>> H; // vector with message values
	std::vector<ulongint> message_to_edge; // vector with egde id for messages
	std::unordered_map<ulongint,ulongint> edge_to_message; // location of message for edge
	ulongint k = 0;
	for (ulongint i=0; i<n; i++){
		for (unsigned int x=0; x<Fnode[i].num_nodes; ++x){
			ulongint j = Fnode[i].global_name[x]; //node
			ulongint u = i*n + j; //edge id
			message_to_edge.push_back(u);
			edge_to_message[u] = k;
			k += 1;
			H.push_back(message_initial_value);
		}
	}
	for (ulongint i=0; i<n; ++i) Fnode[i].update_message_names(edge_to_message, n);

	// Initialize neighborhoods for messages
	Neighborhood *F;
	F = new Neighborhood[number_of_messages];
	init_message_F(F, network, r, message_to_edge, edge_to_message);

	std::vector<ulongint> message_order;
	for (ulongint i=0; i<number_of_messages; ++i) message_order.push_back(i);

	double current_density = -10;
	double old_density = -50;

	std::cout << "Starting message passing..." << std::endl;

	while (z.real()<z_real_max){

		double percent_complete = (((z-z_real_min)/delta_z)*100.0/((z_real_max-z_real_min)/delta_z)).real();
		printf("%2.1f%%, ", percent_complete);
		fflush(stdout);
	
		for (int itr=0; itr<max_iterations; ++itr){
			for (int x=0; x<4; ++x){
				std::random_shuffle(message_order.begin(), message_order.end());
				update_messages(F, network, H, z, message_order);
			}
			current_density = density(Fnode, network, H, z);
			if(std::abs(current_density-old_density)<tolerance) break;
			old_density = current_density;
		}
		outfile << z.real() << " " << current_density << std::endl;
		z = z+delta_z;
	}

	std::cout << "Complete." << std::endl;

	return 0;
}
