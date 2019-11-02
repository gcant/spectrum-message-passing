#include <string>
#include <iostream>
#include "include/MessagePassing.h"
#include <unistd.h>

int main(int argc, char* argv[]){

	std::string fname = argv[1];
	int r=1;
	std::string outfname = "out.txt";
	std::vector<double> z;
	double z_real_min=-1;
	double z_real_max=1;
	double z_imag=0.01;
	int num_points=100;
	double delta_z = (z_real_max-z_real_min)/num_points;

	int c=0;

	while ((c = getopt (argc, argv, "i:r:o:z:")) != -1){
		switch(c){
			case 'i':
				fname = optarg;
				break;
			case 'r':
				r = std::stoi(optarg);
				break;
			case 'o':
				outfname = optarg;
				break;
			case 'z':
				std::stringstream ss(optarg);
				while( ss.good() ){
					std::string substr;
					getline( ss, substr, ',' );
					z.push_back( std::stof(substr) );
				}
				z_real_min = z[0];
				z_real_max = z[1];
				delta_z = (z_real_max-z_real_min)/z[2]; 
				z_imag = z[3];
				break;
		}
	}


	std::cout << "Computing r=" << r << " message passing for spectrum of ";
	std::cout << fname << ", for z = " << z_real_min << ":" << z_real_max;
	std::cout << " with Im(z)=" << z_imag << std::endl;

	MessagePassing(fname, r, z_real_min, z_real_max, z_imag, delta_z, outfname);

	return 0;
}
