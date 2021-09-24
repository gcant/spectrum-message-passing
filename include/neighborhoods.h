#ifndef NEIGHBORHOOD_HEADER 
#define NEIGHBORHOOD_HEADER

#include "WeightedGraph.h"
#include "find_neighborhoods.h"
#include <complex>
#include <Eigen/Dense>


typedef std::pair<int,int> EDGE;
typedef std::vector<EDGE> EDGE_LIST;
class Neighborhood;
typedef std::unordered_map<int,std::unordered_map<int,Neighborhood>> N_map;

typedef std::complex<double> COMPLEX;


class Neighborhood {
  public:
    int f_node;
    EDGE_LIST edges;
    EDGE_LIST local_edges;
    std::vector<double> edge_weights;
    std::vector<int> nodes;
    COMPLEX c_value=0.9;
    double self_loop;
    Eigen::Matrix<COMPLEX, Eigen::Dynamic, 1> u;
    Eigen::Matrix<COMPLEX, Eigen::Dynamic, 1> v;

    COMPLEX compute_value(COMPLEX const&, N_map&);
    double update_value(COMPLEX const&, N_map&);
    void init(EDGE_LIST const&, int, WGraph &);
};


void Neighborhood::init(EDGE_LIST const &neighborhood_edges, int node, WGraph &G) {
  f_node = node;
  edges = neighborhood_edges;
  self_loop = G(node);

  std::unordered_set<int> node_set;
  node_set.insert(node);
  for (auto edge : edges) {
    node_set.insert(edge.first);
    node_set.insert(edge.second);
  }

  node_set.erase(node);
  nodes.resize(node_set.size());
  std::unordered_map<int,int> local_names;
  int temp=0;
  for (auto j : node_set) {
    nodes[temp] = j;
    local_names[j] = temp;
    ++temp;
  }

  u.resize(nodes.size(),1);
  v.resize(nodes.size(),1);
  for (int i=0; i<nodes.size(); ++i) {
    u(i) = 0.0;
    v(i) = 0.0;
  }

  for (int s=0; s<edges.size(); ++s) {
    if ((edges[s].first!=node) && (edges[s].second!=node)) {
      local_edges.push_back({ local_names[edges[s].first], local_names[edges[s].second] });
      edge_weights.push_back(G(edges[s].first,edges[s].second));
    }
    if (edges[s].first==node) {
      int j = edges[s].second;
      v(local_names[j]) = G(node,j);
      u(local_names[j]) = G(j,node);
    }
    if (edges[s].second==node) {
      int j = edges[s].first;
      v(local_names[j]) = G(node,j);
      u(local_names[j]) = G(j,node);
    }
  }
}

COMPLEX Neighborhood::compute_value(COMPLEX const &z, N_map &N_diff) {

  int size = nodes.size();
  if (size==0) return self_loop;

  Eigen::Matrix<COMPLEX, Eigen::Dynamic, Eigen::Dynamic> M;
  M.setIdentity(size,size);
  for (int j=0; j<size; ++j) {
    M(j,j) = (z - N_diff.at(f_node).at(nodes[j]).c_value) ;
  }

  for (int s=0; s<local_edges.size(); ++s) {
    int j = local_edges[s].first;
    int k = local_edges[s].second;
    M(k,j) = M(j,k) = -edge_weights[s];
  }

  return v.dot(M.colPivHouseholderQr().solve(u)) + self_loop;
}

double Neighborhood::update_value(COMPLEX const &z, N_map &N_diff) {
  auto new_value = compute_value(z,N_diff);
  double diff = (double) abs(c_value - new_value);
  c_value = new_value;
  return diff;
}


#endif
