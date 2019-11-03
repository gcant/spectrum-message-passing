#pragma once
#include <complex>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <Eigen/Dense>

typedef unsigned long long int ulongint;
typedef std::unordered_set<unsigned long long int> int_set;
typedef std::complex<double> complexd;

// Find the uid for each edge in the neighborhood of i, at loop length r.
//
int_set find_node_neighborhood_edges(Graph &network, ulongint const &i, int const &r){
    int_set ans;
    for (ulongint j : network.neighbors(i)){
        ans.insert(network.uid(i,j));
        if (r>0){
            for (ulongint k : network.neighbors(i)){
                if (j<k){
                    if (network.has_edge(j,k)) ans.insert(network.uid(j,k));
                    if (r>1){
                        for (ulongint l : network.common_neighbors(j,k)){
                            if (l!=i){
                                ans.insert(network.uid(l,k));
                                ans.insert(network.uid(l,j));
                            }
                        }
                    }
                }
            }
        }
    }
    return ans;
};


// set difference for unordered set of unsigned long long ints
int_set set_difference(int_set S1, int_set S2){
	int_set diff;
	for (ulongint i : S1)
		if(not(S2.count(i))) diff.insert(i);
	return diff;
};


// find message neighborhood, defined to be the set difference of
// edges between the two neighborhoods
int_set find_message_neighborhood_edges(Graph &network, ulongint const &i,
		ulongint const &j, int const &r){
	int_set Ni = find_node_neighborhood_edges(network, i, r);
	int_set Nj = find_node_neighborhood_edges(network, j, r);
	return set_difference(Nj,Ni);
};


// Structure for neighborhood data.
// direct and indirect are vectors of nodes that are either directly
// connected to the central node, or in the neighborhood but not directly connected.
// edges is a vector of edges between nodes in the neighborhood (this vector would
// be empty for locally tree graphs)
struct DirectIndirectEdges{
	std::vector<ulongint> direct;
	std::vector<ulongint> indirect;
	std::vector<std::pair<ulongint,ulongint>> edges;
};


// Find which nodes are directly connected to the center node, which
// are indirectly connected and all edges between neighbors
DirectIndirectEdges construct_neighborhood_from_edges(Graph &network,
		int_set &edges_in, ulongint center_node){

	int_set direct;
	int_set indirect;
	std::vector<std::pair<ulongint,ulongint>> edges_out;

	// find all nodes in edges_in
	int_set all_nodes;
	for (ulongint u : edges_in){
		auto nodes = network.uid(u);
		all_nodes.insert(nodes.first);
		all_nodes.insert(nodes.second);
	}

	// find all nodes that are in the same component as the center node
	DisjointSet component;
	component.init(all_nodes);
	for (ulongint u : edges_in){
		auto nodes = network.uid(u);
		component.join_sets(nodes.first, nodes.second);
	}
	int_set allowed_nodes = component.nodes_in_same_set(center_node);
	allowed_nodes.insert(center_node);
	
	// Create a new set of the edges in the component of center node
	int_set new_edges_in;
	for (ulongint u : edges_in)
		if (allowed_nodes.count(network.uid(u).first))
			new_edges_in.insert(u);
	

	// create sets for direct and indirectly connected edges
	for (ulongint u : new_edges_in){
		auto nodes= network.uid(u);
		ulongint k1 = nodes.first, k2 = nodes.second;
		if(k1==center_node) direct.insert(k2);
		else if(k2==center_node) direct.insert(k1);
		else edges_out.push_back(nodes);
	}
	for (ulongint u : new_edges_in){
		auto nodes = network.uid(u);
		ulongint k1 = nodes.first, k2 = nodes.second;
		if(not(direct.count(k1))) indirect.insert(k1);
		if(not(direct.count(k2))) indirect.insert(k2);
	}
	direct.erase(center_node);
	indirect.erase(center_node);
	
	DirectIndirectEdges ans;
	for (ulongint i : direct) ans.direct.push_back(i);
	for (ulongint i : indirect) ans.indirect.push_back(i);
	ans.edges = edges_out;

	return ans;

};

// class for neighborhoods, stores which edges and nodes are present,
// and calculates its message
class Neighborhood{
	public:
	        void initialize(ulongint const &, std::vector<ulongint> &, 
	        	std::vector<ulongint> &,
			std::vector<std::pair<ulongint,ulongint>> &, double, Graph &);
		ulongint central_node;
		ulongint M; //number of edges in neighborhood
		ulongint num_nodes;
		ulongint num_direct;
		double self_loop;
		Eigen::Matrix<complexd, Eigen::Dynamic, 1> n;
		ulongint *global_name;
		ulongint *message_name;
		std::pair<ulongint,ulongint> *edges;
		complexd *weights;
		complexd compute_message(std::vector<complexd> const &H,
				complexd const &z);
		void update_message_names(std::unordered_map<ulongint,ulongint> &,
				ulongint const &);
};

void Neighborhood::initialize(ulongint const &center_node,
	std::vector<ulongint> &direct_nodes_in, std::vector<ulongint> &indirect_nodes_in,
	std::vector<std::pair<ulongint,ulongint>> &edges_in, double self_loop_in, Graph &network){

        self_loop = self_loop_in;
	central_node = center_node;
	M = edges_in.size();
	num_direct = direct_nodes_in.size();
	num_nodes = num_direct + indirect_nodes_in.size();
	n.resize(num_nodes,1);
	global_name = new ulongint[num_nodes];
	message_name = new ulongint[num_nodes];
	edges = new std::pair<ulongint,ulongint>[M];
	weights = new complexd[M];
	std::unordered_map<ulongint,ulongint> name_map; // global to local names

	for (ulongint i=0; i<num_direct; ++i){
		global_name[i] = direct_nodes_in[i];
		n(i) = network.weight(center_node, direct_nodes_in[i]);
		name_map[global_name[i]] = i;
	}
	for (ulongint i=num_direct; i<num_nodes; ++i){
		global_name[i] = indirect_nodes_in[i-num_direct];
		n[i] = 0;
		name_map[global_name[i]] = i;
	}

	for (ulongint i=0; i<M; ++i){
		std::pair<ulongint,ulongint> edge;
		edge.first = name_map.at(edges_in[i].first);
		edge.second = name_map.at(edges_in[i].second);
		edges[i] = edge;
		weights[i] = network.weight(edges_in[i].first, edges_in[i].second);
	}

};


inline complexd Neighborhood::compute_message(std::vector<complexd> const &H,
		complexd const &z){

	if (num_nodes==0) return self_loop;

	Eigen::Matrix<complexd, Eigen::Dynamic, Eigen::Dynamic> B;
	B.setIdentity(num_nodes,num_nodes);

	for (ulongint i=0; i<num_nodes; ++i)
		B(i,i) = z-H[message_name[i]];

	for (ulongint i=0; i<M; ++i){
		ulongint ii = edges[i].first;
		ulongint jj = edges[i].second;
		B(ii,jj) = B(jj,ii) = -weights[i];
	}
	return n.dot(B.colPivHouseholderQr().solve(n)) + self_loop;
}


void Neighborhood::update_message_names(std::unordered_map<ulongint,ulongint> &m_names,
		ulongint const &N){

    for (ulongint i=0; i<num_nodes; ++i)
        message_name[i] = m_names.at(central_node*N + global_name[i]);

};

void init_node_F(Neighborhood *Fnode, Graph &network, int r){
	for (ulongint i : network.nodes()){
		int_set edges = find_node_neighborhood_edges(network, i, r);
		DirectIndirectEdges temp = construct_neighborhood_from_edges(network, edges, i);
		Fnode[i].initialize(i,temp.direct,temp.indirect,temp.edges,network.weight(i,i),network);
	}
};

void init_message_F(Neighborhood *F, Graph &network, int r, std::vector<ulongint> const &m_id,
		std::unordered_map<ulongint,ulongint> &edge_to_message){
	#pragma omp parallel for	
	for (ulongint k=0; k<m_id.size(); ++k){
		ulongint i = m_id[k]/network.number_of_nodes();
		ulongint j = m_id[k]%network.number_of_nodes();
		int_set edges = find_message_neighborhood_edges(network, i, j, r);
		DirectIndirectEdges temp = construct_neighborhood_from_edges(network, edges, j);
		F[k].initialize(j,temp.direct,temp.indirect,temp.edges,network.weight(j,j),network);
		F[k].update_message_names(edge_to_message, network.number_of_nodes());
	}
};

void update_messages(Neighborhood *F, Graph &network,
		std::vector<complexd> &H, complexd const &z,
		std::vector<ulongint> const &message_order){

	ulongint N = message_order.size();
	#pragma omp parallel for
	for (ulongint k=0; k<N; ++k){
		complexd new_message = F[message_order[k]].compute_message(H,z);
		H[message_order[k]] = new_message;
	}
};

double density(Neighborhood *Fnode, Graph &network,
		std::vector<complexd> &H, complexd const &z){
	
	ulongint n = network.number_of_nodes();
	double rho = 0.0;
	#pragma omp parallel for reduction(+: rho)
	for (ulongint i=0; i<n; ++i){
		rho += (1./(z-Fnode[i].compute_message(H,z))).imag();
	}
	return -rho/(n*3.14159265);
}
