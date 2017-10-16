// BoostLibraryDijkstraExample.cpp: define el punto de entrada de la aplicación de consola.
//

#include "stdafx.h"
#include <boost/config.hpp>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <chrono>


#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/indirect_cmp.hpp>
#include <boost/range/irange.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/depth_first_search.hpp>
//#include <boost/pending/integer_range.hpp>
//#include <boost/pending/indirect_cmp.hpp>
#include <boost/graph/visitors.hpp>


using namespace boost;
clock_t clkStart;
clock_t clkFinish;


template <class VisitorList>
struct edge_categorizer : public dfs_visitor<VisitorList> {
	typedef dfs_visitor<VisitorList> Base;

	edge_categorizer(const VisitorList& v = null_visitor()) : Base(v) { }

	template <class Edge, class Graph>
	void tree_edge(Edge e, Graph& G) {
		cout << "Tree edge: " << source(e, G) <<
			" --> " << target(e, G) << endl;
		Base::tree_edge(e, G);
	}
	template <class Edge, class Graph>
	void back_edge(Edge e, Graph& G) {
		cout << "Back edge: " << source(e, G)
			<< " --> " << target(e, G) << endl;
		Base::back_edge(e, G);
	}
	template <class Edge, class Graph>
	void forward_or_cross_edge(Edge e, Graph& G) {
		cout << "Forward or cross edge: " << source(e, G)
			<< " --> " << target(e, G) << endl;
		Base::forward_or_cross_edge(e, G);
	}
	template <class Edge, class Graph>
	void finish_edge(Edge e, Graph& G) {
		cout << "Finish edge: " << source(e, G) <<
			" --> " << target(e, G) << endl;
		Base::finish_edge(e, G);
	}
};
template <class VisitorList>
edge_categorizer<VisitorList>
categorize_edges(const VisitorList& v) {
	return edge_categorizer<VisitorList>(v);
}


template < typename TimeMap > class bfs_time_visitor :public default_bfs_visitor {
	typedef typename property_traits < TimeMap >::value_type T;
public:
	bfs_time_visitor(TimeMap tmap, T & t) :m_timemap(tmap), m_time(t) { }
	template < typename Vertex, typename Graph >
	void discover_vertex(Vertex u, const Graph & g) const
	{
		put(m_timemap, u, m_time++);
	}
	TimeMap m_timemap;
	T & m_time;
};

template < typename TimeMap > class dfs_time_visitor :public default_dfs_visitor {
	typedef typename property_traits < TimeMap >::value_type T;
public:
	dfs_time_visitor(TimeMap dmap, TimeMap fmap, T & t)
		: m_dtimemap(dmap), m_ftimemap(fmap), m_time(t) {
	}
	template < typename Vertex, typename Graph >
	void discover_vertex(Vertex u, const Graph & g) const
	{
		put(m_dtimemap, u, m_time++);
	}
	template < typename Vertex, typename Graph >
	void finish_vertex(Vertex u, const Graph & g) const
	{
		put(m_ftimemap, u, m_time++);
	}
	TimeMap m_dtimemap;
	TimeMap m_ftimemap;
	T & m_time;
};

//void DFS() {
//	typedef adjacency_list<> Graph;
//
//	Graph G(5);
//	add_edge(0, 2, G);
//	add_edge(1, 1, G);
//	add_edge(1, 3, G);
//	add_edge(2, 1, G);
//	add_edge(2, 3, G);
//	add_edge(3, 1, G);
//	add_edge(3, 4, G);
//	add_edge(4, 0, G);
//	add_edge(4, 1, G);
//
//	typedef graph_traits<Graph>::vertices_size_type size_type;
//
//	std::vector<size_type> d(num_vertices(G));
//	std::vector<size_type> f(num_vertices(G));
//	int t = 0;
//	depth_first_search(G, visitor(categorize_edges(
//		make_pair(stamp_times(&d[0], t, on_discover_vertex()),
//			stamp_times(&f[0], t, on_finish_vertex())))));
//
//	std::vector<size_type>::iterator i, j;
//	for (i = d.begin(), j = f.begin(); i != d.end(); ++i, ++j)
//		std::cout << *i << " " << *j << std::endl;
//
//}

void DFS() {
	// Select the graph type we wish to use
	typedef adjacency_list < vecS, vecS, directedS > graph_t;
	typedef graph_traits < graph_t >::vertices_size_type size_type;
	// Set up the vertex names
	enum
	{
		v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, N
	};
	char name[] = { '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e' };
	// Specify the edges in the graph
	typedef std::pair < int, int >E;
	E edge_array[] = { E(1, 3), E(1, 4), E(3, 2), E(3, 5),
		E(3, 'a'), E(4, 5), E(4, 7), E(4, 8),
		E(2, 5), E(5, 6), E('a', 3), E('a', 6),
		E(7, 4), E(8, 7), E(8, 9), E(6, 'd'),
		E(9, 'a'), E(9, 'c'), E('d', 'e'), E('c', 9),
		E('c', 'b'), E('c', 'e'), E('e', 'd'), E('b', 'c')
	};
#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	graph_t g(N);
	for (std::size_t j = 0; j < sizeof(edge_array) / sizeof(E); ++j)
		add_edge(edge_array[j].first, edge_array[j].second, g);
#else
	graph_t g(edge_array, edge_array + sizeof(edge_array) / sizeof(E), N);
#endif

	// discover time and finish time properties
	std::vector < size_type > dtime(num_vertices(g));
	std::vector < size_type > ftime(num_vertices(g));
	typedef
		iterator_property_map<std::vector<size_type>::iterator,
		property_map<graph_t, vertex_index_t>::const_type>
		time_pm_type;
	time_pm_type dtime_pm(dtime.begin(), get(vertex_index, g));
	time_pm_type ftime_pm(ftime.begin(), get(vertex_index, g));
	size_type t = 0;
	dfs_time_visitor < time_pm_type >vis(dtime_pm, ftime_pm, t);

	auto start = std::chrono::system_clock::now();
	depth_first_search(g, visitor(vis));
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "DFS elapsed time: " << elapsed_seconds.count() << "s\n";

	// use std::sort to order the vertices by their discover time
	std::vector < size_type > discover_order(N);
	integer_range < size_type > r(0, N);
	std::copy(r.begin(), r.end(), discover_order.begin());
	std::sort(discover_order.begin(), discover_order.end(),
		indirect_cmp < time_pm_type, std::less < size_type > >(dtime_pm));
	std::cout << "order of discovery: ";
	int i;
	for (i = 0; i < N; ++i)
		std::cout << name[discover_order[i]] << " ";

	std::vector < size_type > finish_order(N);
	std::copy(r.begin(), r.end(), finish_order.begin());
	std::sort(finish_order.begin(), finish_order.end(),
		indirect_cmp < time_pm_type, std::less < size_type > >(ftime_pm));
	std::cout << std::endl << "order of finish: ";
	for (i = 0; i < N; ++i)
		std::cout << name[finish_order[i]] << " ";
	std::cout << std::endl;
}

/*void BFS()
{
	using namespace boost;
	// Select the graph type we wish to use
	typedef adjacency_list < vecS, vecS, undirectedS > graph_t;
	// Set up the vertex IDs and names
	enum { v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, N };
	const char *name = "bfs";
	// Specify the edges in the graph
	typedef std::pair < int, int >E;
	E edge_array[] = { E(v1, v3), E(v1, v4), E(v3, v2), E(v3, v5),
						E(v3, v10), E(v4, v5), E(v4, v7), E(v4, v8),
						E(v2, v5), E(v5, v6), E(v10, v3), E(v10, v6),
						E(v7, v4), E(v8, v7), E(v8, v9), E(v6, v13),
						E(v9, v10), E(v9, v12), E(v13, v14), E(v12, v9),
						E(v12, v11), E(v12, v14), E(v14, v13), E(v11, v12)
	};
	// Create the graph object
	const int n_edges = sizeof(edge_array) / sizeof(E);
#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	// VC++ has trouble with the edge iterator constructor
	graph_t g(N);
	for (std::size_t j = 0; j < n_edges; ++j)
		add_edge(edge_array[j].first, edge_array[j].second, g);
#else
	typedef graph_traits<graph_t>::vertices_size_type v_size_t;
	graph_t g(edge_array, edge_array + n_edges, v_size_t(N));
#endif

	// Typedefs
	typedef graph_traits < graph_t >::vertex_descriptor Vertex;
	typedef graph_traits < graph_t >::vertices_size_type Size;
	typedef Size* Iiter;

	// a vector to hold the discover time property for each vertex
	std::vector < Size > dtime(num_vertices(g));

	Size time = 0;
	bfs_time_visitor < Size * >vis(&dtime[0], time);

	auto start = std::chrono::system_clock::now();
	breadth_first_search(g, vertex(s, g), visitor(vis));
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "BFS elapsed time: " << elapsed_seconds.count() << "s\n";

	// Use std::sort to order the vertices by their discover time
	std::vector<graph_traits<graph_t>::vertices_size_type > discover_order(N);
	integer_range < int >range(0, N);
	std::copy(range.begin(), range.end(), discover_order.begin());
	std::sort(discover_order.begin(), discover_order.end(),
		indirect_cmp < Iiter, std::less < Size > >(&dtime[0]));

	std::cout << "order of discovery: ";
	for (int i = 0; i < N; ++i)
		std::cout << name[discover_order[i]] << " ";
	std::cout << std::endl;
}*/

void Dijkstra()
{
	typedef adjacency_list < listS, vecS, directedS,
		no_property, property < edge_weight_t, int > > graph_t;
	typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
	typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
	typedef std::pair<int, int> Edge;

	const int num_nodes = 5;
	enum nodes { v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14 };
	char name[] = "Dijkstra";
	Edge edge_array[] = { Edge(v1, v3), Edge(v1, v4), Edge(v3, v2), Edge(v3, v5),
		Edge(v3, v10), Edge(v4, v5), Edge(v4, v7), Edge(v4, v8),
		Edge(v2, v5), Edge(v5, v6), Edge(v10, v3), Edge(v10, v6),
		Edge(v7, v4), Edge(v8, v7), Edge(v8, v9), Edge(v6, v13),
		Edge(v9, v10), Edge(v9, v12), Edge(v13, v14), Edge(v12, v9),
		Edge(v12, v11), Edge(v12, v14), Edge(v14, v13), Edge(v11, v12)
	};
	int weights[] = { 8,8,7,8,4,1,3,2,7,9,10,6,6,3,3,4,2,4,6,2,8,9,2,6 };
	int num_arcs = sizeof(edge_array) / sizeof(Edge);
#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	graph_t g(num_nodes);
	property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
	for (std::size_t j = 0; j < num_arcs; ++j) {
		edge_descriptor e; bool inserted;
		boost::tie(e, inserted) = add_edge(edge_array[j].first, edge_array[j].second, g);
		weightmap[e] = weights[j];
	}
#else
	graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);
	property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
#endif
	std::vector<vertex_descriptor> p(num_vertices(g));
	std::vector<int> d(num_vertices(g));
	vertex_descriptor s = vertex(v1, g);

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	// VC++ has trouble with the named parameters mechanism
	property_map<graph_t, vertex_index_t>::type indexmap = get(vertex_index, g);
	dijkstra_shortest_paths(g, s, &p[0], &d[0], weightmap, indexmap,
		std::less<int>(), closed_plus<int>(),
		(std::numeric_limits<int>::max)(), 0,
		default_dijkstra_visitor());
#else
	auto start = std::chrono::system_clock::now();
	dijkstra_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&d[0]));
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "Dijkstra elapsed time: " << elapsed_seconds.count() << "s\n";
#endif

	std::cout << "distances and parents:" << std::endl;
	graph_traits < graph_t >::vertex_iterator vi, vend;
	for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) {
		std::cout << "distance(" << name[*vi] << ") = " << d[*vi] << ", ";
		std::cout << "parent(" << name[*vi] << ") = " << name[p[*vi]] << std::
			endl;
	}
	std::cout << std::endl;

	std::ofstream dot_file("figs/dijkstra-eg.dot");

	dot_file << "digraph D {\n"
		<< "  rankdir=LR\n"
		<< "  size=\"4,3\"\n"
		<< "  ratio=\"fill\"\n"
		<< "  edge[style=\"bold\"]\n" << "  node[shape=\"circle\"]\n";

	graph_traits < graph_t >::edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
		graph_traits < graph_t >::edge_descriptor e = *ei;
		graph_traits < graph_t >::vertex_descriptor
			u = source(e, g), v = target(e, g);
		dot_file << name[u] << " -> " << name[v]
			<< "[label=\"" << get(weightmap, e) << "\"";
		if (p[v] == u)
			dot_file << ", color=\"black\"";
		else
			dot_file << ", color=\"grey\"";
		dot_file << "]";
	}
	dot_file << "}";
}

void Prim() {
	using namespace boost;
	typedef adjacency_list < vecS, vecS, undirectedS,
		property<vertex_distance_t, int>, property < edge_weight_t, int > > Graph;
	typedef std::pair < int, int >E;
	const int num_nodes = 14;
	E edges[] = { E(1, 3), E(1, 4), E(3, 2), E(3, 5),
		E(3, 10), E(4, 5), E(4, 7), E(4, 8),
		E(2, 5), E(5, 6), E(10, 3), E(10, 6),
		E(7, 4), E(8, 7), E(8, 9), E(6, 13),
		E(9, 10), E(9, 12), E(13, 14), E(12, 9),
		E(12, 11), E(12, 14), E(14, 13), E(11, 12)
	};
	int weights[] = { 8,8,7,8,4,1,3,2,7,9,10,6,6,3,3,4,2,4,6,2,8,9,2,6 };
#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	Graph g(num_nodes);
	property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
	for (std::size_t j = 0; j < sizeof(edges) / sizeof(E); ++j) {
		graph_traits<Graph>::edge_descriptor e; bool inserted;
		boost::tie(e, inserted) = add_edge(edges[j].first, edges[j].second, g);
		weightmap[e] = weights[j];
	}
#else
	Graph g(edges, edges + sizeof(edges) / sizeof(E), weights, num_nodes);
	property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
#endif
	std::vector < graph_traits < Graph >::vertex_descriptor >
		p(num_vertices(g));

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	property_map<Graph, vertex_distance_t>::type distance = get(vertex_distance, g);
	property_map<Graph, vertex_index_t>::type indexmap = get(vertex_index, g);
	prim_minimum_spanning_tree
	(g, *vertices(g).first, &p[0], distance, weightmap, indexmap,
		default_dijkstra_visitor());
#else
	auto start = std::chrono::system_clock::now();
	prim_minimum_spanning_tree(g, &p[0]);
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "Prim elapsed time: " << elapsed_seconds.count() << "s\n";
#endif

	for (std::size_t i = 0; i != p.size(); ++i)
		if (p[i] != i)
			std::cout << "parent[" << i << "] = " << p[i] << std::endl;
		else
			std::cout << "parent[" << i << "] = no parent" << std::endl;
}

void Kruskal() {
	using namespace boost;
	typedef adjacency_list < vecS, vecS, undirectedS,
		no_property, property < edge_weight_t, int > > Graph;
	typedef graph_traits < Graph >::edge_descriptor Edge;
	typedef graph_traits < Graph >::vertex_descriptor Vertex;
	typedef std::pair<int, int> E;

	const int num_nodes = 14;
	E edge_array[] = { E(1, 3), E(1, 4), E(3, 2), E(3, 5),
		E(3, 10), E(4, 5), E(4, 7), E(4, 8),
		E(2, 5), E(5, 6), E(10, 3), E(10, 6),
		E(7, 4), E(8, 7), E(8, 9), E(6, 13),
		E(9, 10), E(9, 12), E(13, 14), E(12, 9),
		E(12, 11), E(12, 14), E(14, 13), E(11, 12)
	};
	int weights[] = { 8,8,7,8,4,1,3,2,7,9,10,6,6,3,3,4,2,4,6,2,8,9,2,6 };
	std::size_t num_edges = sizeof(edge_array) / sizeof(E);
#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	Graph g(num_nodes);
	property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
	for (std::size_t j = 0; j < num_edges; ++j) {
		Edge e; bool inserted;
		boost::tie(e, inserted) = add_edge(edge_array[j].first, edge_array[j].second, g);
		weightmap[e] = weights[j];
	}
#else
	Graph g(edge_array, edge_array + num_edges, weights, num_nodes);
#endif
	property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
	std::vector < Edge > spanning_tree;

	auto start = std::chrono::system_clock::now();
	kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "Kruskal elapsed time: " << elapsed_seconds.count() << "s\n";

	std::cout << "Print the edges in the MST:" << std::endl;
	for (std::vector < Edge >::iterator ei = spanning_tree.begin();
		ei != spanning_tree.end(); ++ei) {
		std::cout << source(*ei, g) << " <--> " << target(*ei, g)
			<< " with weight of " << weight[*ei]
			<< std::endl;
	}

	std::ofstream fout("figs/kruskal-eg.dot");
	fout << "graph A {\n"
		<< " rankdir=LR\n"
		<< " size=\"3,3\"\n"
		<< " ratio=\"filled\"\n"
		<< " edge[style=\"bold\"]\n" << " node[shape=\"circle\"]\n";
	graph_traits<Graph>::edge_iterator eiter, eiter_end;
	for (boost::tie(eiter, eiter_end) = edges(g); eiter != eiter_end; ++eiter) {
		fout << source(*eiter, g) << " -- " << target(*eiter, g);
		if (std::find(spanning_tree.begin(), spanning_tree.end(), *eiter)
			!= spanning_tree.end())
			fout << "[color=\"black\", label=\"" << get(edge_weight, g, *eiter)
			<< "\"];\n";
		else
			fout << "[color=\"gray\", label=\"" << get(edge_weight, g, *eiter)
			<< "\"];\n";
	}
	fout << "}\n";
}

int main() {
	//DFS();
	//BFS();
	//Dijkstra();
	Prim();
	Kruskal();
	return 0;
}

