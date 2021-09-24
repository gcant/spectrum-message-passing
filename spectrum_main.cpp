#include <iostream>
#include "neighborhoods.h"

const double TOL = 1.0e-11;

double run_MP(std::vector<Neighborhood> &H,
    std::unordered_map<int,std::unordered_map<int,Neighborhood>> &H_diff,
    WGraph &G,
    COMPLEX &z)
{
  int num_nodes = G.number_of_nodes();
  for (int s=0; s<100; ++s) {
    double Delta = 0.0;
    #pragma omp parallel for reduction(+:Delta)
    for (int i=0; i<num_nodes; ++i) {
      for (int j : H[i].nodes) {
        Delta += H_diff[i][j].update_value(z,H_diff);
      }
    }
    if (Delta < TOL) break;
  }

  #pragma omp parallel for
  for (int i=0; i<num_nodes; ++i) {
    H[i].update_value(z,H_diff);
  }

  double rho = 0.0;
  for (int i=0; i<num_nodes; ++i){
    rho += (1./(z-H[i].c_value)).imag();
  }
  rho = -rho/(num_nodes*3.14159265);
  return rho;
}


int main(int argc, char* argv[]) {
  WGraph G = read_edgelist(argv[1]);
  int num_nodes = G.number_of_nodes();
  int r=std::stoi(argv[2]);
  double eps = std::stof(argv[3]);
  double x_min = std::stof(argv[4]);
  double x_max = std::stof(argv[5]);
  int num_pts = std::stof(argv[6]);
  double dx = (x_max-x_min)/num_pts;

  std::vector<Neighborhood> H(G.number_of_nodes());
  for (int i : G.nodes())
    H[i].init(find_neighborhood_edges(G,i,r),i,G);

  std::unordered_map<int,std::unordered_map<int,Neighborhood>> H_diff;
  for (int i : G.nodes()) {
    for (int j : H[i].nodes) {
      auto edges = difference(G,find_neighborhood_edges(G,j,r), find_neighborhood_edges(G,i,r));
      H_diff[i][j].init(edges,j,G);
    }
  }

  for (double x=x_min; x<=x_max; x+=dx) {
    COMPLEX z = {x,eps};
    std::cout << x << " " << run_MP(H,H_diff,G,z) << std::endl;
  }

  return 0;
}
