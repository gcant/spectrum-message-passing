#pragma once
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>

class WGraph {
  protected:
    int m; //number of edges
    std::unordered_set<int> node_set;
    std::unordered_map<int, std::unordered_map<int,double>> edge_weights;
    std::unordered_map<int, double> self_loops;
  public:
    WGraph();
    int number_of_nodes() const;
    int number_of_edges() const;
    std::unordered_set<int> nodes() const;
    std::unordered_map<int,double> neighbors(int) const;
    int degree(int) const;
    bool has_node(int) const;
    bool has_edge(int, int) const;
    double operator()(int, int) const;
    double operator()(int) const;
    void add_node(int);
    void remove_node(int);
    void add_edge(int, int, double);
    void remove_edge(int, int);
    bool has_self_loop(int) const;
    void add_self_loop(int, double);
};


WGraph::WGraph(){m = 0;}

int WGraph::number_of_nodes() const {return node_set.size();}

int WGraph::number_of_edges() const {return m;}

std::unordered_set<int> WGraph::nodes() const {return node_set;}

std::unordered_map<int,double> WGraph::neighbors(int i) const {return edge_weights.at(i);}

int WGraph::degree(int i) const {return edge_weights.at(i).size();}

bool WGraph::has_node(int i) const {
  if (node_set.count(i))
    return true;
  else
    return false;
}

bool WGraph::has_edge(int i, int j) const {
  if (has_node(i) && edge_weights.at(i).count(j)) 
    return true;
  else
    return false;
}

double WGraph::operator()(int i, int j) const {
  return has_edge(i,j) ? edge_weights.at(i).at(j) : 0.0;
}

double WGraph::operator()(int i) const {
  return has_self_loop(i) ? self_loops.at(i) : 0.0;
}

void WGraph::add_node(int i){
  if (not(has_node(i))){
    node_set.insert(i);
    edge_weights[i];
  }
}

void WGraph::remove_node(int i){
  m -= degree(i);
  node_set.erase(i);
  for (const auto & [ j, w ] : edge_weights[i]) edge_weights[j].erase(i);
  edge_weights.erase(i);
}

void WGraph::add_edge(int i, int j, double w){
  if (i==j) add_self_loop(i,w);
  else {
    if (not(has_edge(i,j))) m += 1;
    add_node(i);
    add_node(j);
    edge_weights[i][j]=w;
    edge_weights[j][i]=w;
  }
}

void WGraph::remove_edge(int i, int j){
  if (has_edge(i,j)){
    edge_weights[i].erase(j);
    edge_weights[j].erase(i);
    m -= 1;
  }
};

void WGraph::add_self_loop(int i, double w){
  add_node(i);
  self_loops[i]=w;
}

bool WGraph::has_self_loop(int i) const {
  if (self_loops.count(i))
    return true;
  else
    return false;
}


WGraph read_edgelist(std::string f_name){

  std::ifstream in_file(f_name);
  WGraph G;
  int i,j;
  double w;

  while (in_file >> i >> j >> w){
    G.add_edge(i,j,w);
    in_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  in_file.close();

  return G;
};

