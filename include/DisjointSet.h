#pragma once
#include <unordered_set>
#include <unordered_map>
#include <tuple>

// Class for disjoint set data structure.
// Initialize class with an unordered set of nodes using init method.
typedef unsigned long long int ulongint;

class DisjointSet {
	private:
		std::unordered_map<ulongint, ulongint> parent;
	public:
		void init(std::unordered_set<ulongint> const &);
		std::unordered_map<ulongint, ulongint> cluster_size;
		void join_sets(ulongint const &, ulongint const &);
		ulongint &operator[](ulongint const &);
		std::unordered_set<ulongint> nodes_in_same_set(ulongint const &);
};


void DisjointSet::init(std::unordered_set<ulongint> const &nodes){
	for (ulongint i : nodes){
		parent[i]=i;
		cluster_size[i]=1;
	}
};

ulongint &DisjointSet::operator[](ulongint const &index){
	//Find the root.
	ulongint root;
	root=parent[index];
	while(root!=parent[root]) root = parent[root];

	//Compress the path.
	ulongint x,y;
	x = parent[index];
	parent[index]=root;
	while(x!=parent[x]){
		y=parent[x];
		parent[x]=root;
		x=y;
	};
	return parent[index];
};

void DisjointSet::join_sets(ulongint const &x, ulongint const &y){
	ulongint root_x, root_y;
	ulongint large, small;
	root_x = operator[](x);
	root_y = operator[](y);
	if(root_x!=root_y){
		large = root_x; small = root_y;
		if (cluster_size[root_x]<cluster_size[root_y]){
			large = root_y;
			small = root_x;
		};
		cluster_size[large] += cluster_size[small];
		parent[small] = large;
	}
};


std::unordered_set<ulongint> DisjointSet::nodes_in_same_set(ulongint const &i){
	std::unordered_set<ulongint> ans;
	ulongint root_i, root_j;
	root_i = operator[](i);
	for (std::pair<ulongint,ulongint> j : parent){
		root_j = operator[](j.first);
		if (root_i == root_j) ans.insert(j.first);
	}
	return ans;
};
