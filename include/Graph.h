#pragma once
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

// Simple class for networks. Assumes there are no degree 0 nodes,
// and nodes are labelled by integers, starting at 0.

typedef unsigned long long int ulongint;
typedef std::unordered_set<unsigned long long int> int_set;

class Graph {
	private:
		ulongint n;
		ulongint m;
		int_set node_set;
		std::unordered_map<ulongint, int_set> edge_set;
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> W;
	public:
		Graph(){};
		Graph(std::string fname){ load_edgelist(fname); };
		void load_edgelist(std::string fname);
		ulongint number_of_nodes(){return n;};
		ulongint number_of_edges(){return m;};
		bool has_edge(ulongint const &, ulongint const &);
		int_set nodes(){return node_set;};
		int_set neighbors(ulongint i){return edge_set[i];};
		int_set common_neighbors(ulongint const &, ulongint const &);
		ulongint uid(ulongint const &, ulongint const &);
		double degree(ulongint);
		std::pair<ulongint,ulongint> uid(ulongint const &);
                double weight(ulongint, ulongint);
};


void Graph::load_edgelist(std::string fname){
	
	// Load graph from file fname.
	std::ifstream file;
	file.open(fname);
	ulongint i,j;
        double w;
	while (file >> i >> j >> w){
		if (i!=j){
			node_set.insert(i);
			node_set.insert(j);
			edge_set[i].insert(j);
			edge_set[j].insert(i);
		}
	}
	file.close();

	// Set n = number of nodes; m = number of edges.
	n = node_set.size();
	m = 0;
	for (ulongint i : node_set) m += edge_set[i].size();
	m = m/2;

	// set weights in weight matrix
	W.resize(n,n);
	file.open(fname);
	while (file >> i >> j >> w){
		W(i,j) = W(j,i) = w;
	}
	file.close();

};

bool Graph::has_edge(ulongint const &i, ulongint const &j){
	return (bool) edge_set.at(i).count(j);
};

int_set Graph::common_neighbors(ulongint const &i, ulongint const &j){
	int_set shared_neighbors;
	for (ulongint k : neighbors(i)){
		if (has_edge(j, k)) shared_neighbors.insert(k);
	};
	return shared_neighbors;
};

ulongint Graph::uid(ulongint const &i, ulongint const &j){
	if(i<j) return n*i + j;
	else return n*j + i;
};

std::pair<ulongint,ulongint> Graph::uid(ulongint const &u){
	std::pair<ulongint,ulongint> nodes;
	nodes = { std::min(u/n,u%n), std::max(u/n,u%n) };
	return nodes;
};

double Graph::degree(ulongint i){
	return (double) edge_set.at(i).size();
};

double Graph::weight(ulongint i, ulongint j){
        return W(i,j);
};
