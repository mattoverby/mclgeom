// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_GRAPHCOLOR_HPP
#define MCL_GEOM_GRAPHCOLOR_HPP

#include <Eigen/Sparse>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>
#include <mutex>
#include <unordered_set>
#include <cstdlib>
#include <random>

namespace mcl {

struct GraphNode {
    int index = -1;
    int color = -1;
    bool conflict = false;
    std::vector<int> neighbors;
    std::unordered_set<int> palette;
};

/// @brief Stochastic graph coloring from:
/// "Fast Distributed Algorithms for Brooks-Vizing Colourings" by Grable and Panconesi (2000)
template<typename T>
inline void graph_color(const Eigen::SparseMatrix<T> &adjacency, std::vector<std::vector<int>> &colors);

/// @brief Creates a directed graph from a symmetric adjacency matrix A (only the upper-triangular is used).
/// Any values in A are considered an edge, no matter the sign/value.
template<typename T>
inline void make_nodes(const Eigen::SparseMatrix<T> &adjacency, std::vector<GraphNode> &nodes);

/// @brief Creates a directed graph from a symmetric adjacency matrix A (only the upper-triangular is used).
/// Any values in A are considered an edge, no matter the sign/value.
template<typename T>
inline void color_nodes(std::vector<GraphNode> &nodes);

/// @brief Given a colored graph, bin them into color groups
template <typename T>
inline void map_colors(const std::vector<GraphNode> &nodes, std::vector<std::vector<int>> &colors);

//
// Implementation
//

template <typename T>
void graph_color(const Eigen::SparseMatrix<T> &adjacency, std::vector<std::vector<int>> &colors)
{
    std::vector<GraphNode> nodes;
    make_nodes<T>(adjacency, nodes);
    color_nodes<T>(nodes);
    map_colors<T>(nodes, colors);
}

template <typename T>
void make_nodes(const Eigen::SparseMatrix<T> &adjacency, std::vector<GraphNode> &nodes)
{
    Eigen::SparseMatrix<T> mat = adjacency;
    mat.makeCompressed();

    int num_nodes = mat.rows();
    nodes.clear();
    nodes.resize(num_nodes);

    // Iterate over outer dimension (columns if column-major, rows if row-major)
    for (int k = 0; k < mat.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            if (i < j && std::abs(it.value()) > T(0)) {
                nodes[i].neighbors.push_back(j);
            }
        }
    }
}

template <typename T>
void color_nodes(std::vector<GraphNode> &nodes)
{
    int init_palette = 6;
    std::vector<int> nodeq(nodes.size());
    std::iota(nodeq.begin(), nodeq.end(), 0);

    // Initialize palette
    tbb::parallel_for(size_t(0), nodes.size(), [&](size_t i) {
        auto& node = nodes[i];
        node.index = i;

        // Initialize palette
        for (int j = 0; j < init_palette; ++j) {
            node.palette.insert(j);
        }

        // Validate directed graph
        // for (int nbr : node.neighbors) {
        //     if (nbr == node.idx)
        //         throw std::runtime_error("Error: Node has self as neighbor");
        //     if (nbr < node.idx) {
        //         std::cerr << "Error: Not a directed graph — node "
        //                   << node.idx << " -> " << nbr << "\n";
        //         throw std::runtime_error("Graph must be upper-triangular directed.");
        //     }
        // }
    });

    int max_iter = nodes.size();
    std::mutex palette_mutex; // For atomic palette merges

    for (int iter = 0; !nodeq.empty() && iter < max_iter; ++iter) {

        // Assign Random Colors
        tbb::parallel_for(size_t(0), nodeq.size(), [&](size_t i) {
            auto& node = nodes[nodeq[i]];
            static thread_local std::mt19937 gen(std::random_device{}());
            std::uniform_int_distribution<> dist(0, node.palette.size() - 1);
            int offset = dist(gen);
            node.color = *std::next(node.palette.begin(), offset);
        });

        // Conflict Detection
        tbb::parallel_for(size_t(0), nodeq.size(), [&](size_t i) {
            auto& node = nodes[nodeq[i]];
            node.conflict = false;
            for (int nbr_idx : node.neighbors) {
                const auto& nbr = nodes[nbr_idx];
                if (node.color == nbr.color) {
                    node.conflict = true;
                    break;
                }
            }
        });

        // Each thread collects palette removals locally, then merges
        tbb::concurrent_vector<std::pair<int,int>> removals;
        tbb::parallel_for(size_t(0), nodeq.size(), [&](size_t i) {
            auto& node = nodes[nodeq[i]];
            const int color = node.color;
            for (int nbr_idx : node.neighbors) {
                const auto& nbr = nodes[nbr_idx];
                // Remove neighbor's color from self palette
                if (!nbr.conflict) {
                    node.palette.erase(nbr.color);
                }
                // Buffer removal of own color from neighbors palette
                if (!node.conflict) {
                    removals.emplace_back(nbr_idx, color);
                }
            }
        });

        // Merge palette removals in serial
        for (auto& [target, c] : removals) {
            nodes[target].palette.erase(c);
        }

        // Remove colored nodes from queue
        nodeq.erase(std::remove_if(nodeq.begin(), nodeq.end(),
            [&](int idx) { return !nodes[idx].conflict; }),
            nodeq.end());

        // Refill palettes
        tbb::parallel_for(size_t(0), nodeq.size(), [&](size_t i) {
            auto& node = nodes[nodeq[i]];
            if (node.palette.size() < 2) {
                node.palette.insert(init_palette + iter);
            }
        });
    }
}


template <typename T>
void map_colors(const std::vector<GraphNode> &nodes, std::vector<std::vector<int>> &colors)
{
    colors.clear();
    if (nodes.empty()) return;

    // Find maximum color index so we can size the vector
    int max_color = -1;
    for (const auto& n : nodes) {
        if (n.color > max_color) {
            max_color = n.color;
        }
    }

    // Preallocate all color buckets
    colors.resize(max_color + 1);

    // Fill in the buckets
    for (const auto& n : nodes) {
        if (n.color >= 0) {
            colors[n.color].push_back(n.index);
        }
    }

    // Remove any empty color buckets
    colors.erase(
        std::remove_if(colors.begin(), colors.end(),
                       [](const std::vector<int>& v) { return v.empty(); }),
        colors.end());

}

#if 0
namespace graphcolor {


	// Interal helper class for setting up the graph and managing colors
	// and node queue (in color_nodes)
	struct GCNode {
		GCNode() : idx(-1), color(-1), conflict(false) {}
		int idx, color; // idx only used for error testing
		bool conflict;
		std::vector<int> neighbors;
		std::unordered_set<int> palette;
	};

	// Creates a directed graph from a symmetric adjacency matrix A (only the upper-triangular is used).
	// Any values in A are considered an edge, no matter the sign/value.
	// So, nodes only see neighbors with a higher index.
	template <typename T>
	static void make_directed_nodes( const SparseMat<T> &A, std::vector<GCNode> &nodes, int stride );

	// Colors the graph created by make_nodes (nodes are modified)
	static void color_nodes( std::vector<GCNode> &nodes );

	// Mapping is color -> list of node indices, e.g.
	//	std::vector<int> color0 = colors[0];
	//	int color0_idx0 = color0[0];
	// where color0_idx0 corresponds to some row in the adjacency matrix
	static void make_map( const std::vector<GCNode> &nodes, std::vector< std::vector<int> > &colors );

	// Colors an adjacency matrix with the above functions
	template <typename T>
	static void color_matrix( const SparseMat<T> &A, std::vector< std::vector<int> > &colors, int stride=1 ){
		std::vector<GCNode> nodes;
		make_directed_nodes( A, nodes, stride );
		color_nodes( nodes );
		make_map( nodes, colors );
	}


} // ns graphcolor

//
// Implementation
//

template <typename T>
static void graphcolor::make_directed_nodes( const SparseMat<T> &A, std::vector<GCNode> &nodes, int stride ){

	// Should also check to make sure it's symmetric but I'll save that for later...
	int dof = A.rows();
	if( dof != A.cols() ){ throw std::runtime_error("**graphcolor::make_nodes Error: Adjacency matrix is not square"); }
	stride = std::max( stride, 1 );
	if( dof % stride != 0 ){ throw std::runtime_error("**graphcolor::make_nodes Error: Bad stride"); }
	const int n_nodes = dof/stride;
	nodes.resize( n_nodes );

	// Create neighbor list
	#pragma omp parallel for schedule(static)
	for( int i=0; i<n_nodes; ++i ){

		GCNode *node = &nodes[i];
		node->idx = i;
		typename Eigen::SparseMatrix<T,Eigen::RowMajor>::InnerIterator it(A,i*stride);

		// Loop until we are past the diagonal
		for( ; it && it.col()/stride <= i; ++it ){}

		// Now add neighbors, avoid adding same neighbor twice
		std::vector<int> already_added = {-1};
		for( ; it; ++it ){
			if( std::abs(it.value()) > T(0) ){
				int idx = it.col()/stride;
				if( already_added.back()==idx ){ continue; }
				already_added.emplace_back(idx);
				node->neighbors.emplace_back(idx);
			}
		}

	} // end loop elements

} // end make nodes


static void graphcolor::color_nodes( std::vector<GCNode> &nodes ){

	int init_palette = 6; // based on observations
	int n_nodes = nodes.size();
	std::vector<int> nodeq( n_nodes ); // node queue
	std::iota(nodeq.begin(), nodeq.end(), 0);

	// Make queue of nodes that need to be processed
	// and initialize their palette.
	#pragma omp parallel for schedule(static)
	for( int i=0; i<n_nodes; ++i ){
		GCNode *node = &nodes[ nodeq[i] ];
		for( int j=0; j<init_palette; ++j ){ node->palette.insert(j); }

		// Make sure graph is directed
		int n_neighbors = node->neighbors.size();
		for( int j=0; j<n_neighbors; ++j ){
			const GCNode *n1 = &nodes[ node->neighbors[j] ];
			if( n1->idx == node->idx ){
				std::cerr << "**graphcolor::color_nodes Error: Self is a neighbor" << std::endl;
				throw std::runtime_error("exiting...");
			}
			if( n1->idx < node->idx ){
				std::cerr << "**graphcolor::color_nodes Error: Not a directed graph" << std::endl;
				std::cerr << "\tnode: " << node->idx << ", neighbor: " << n1->idx << std::endl;
				throw std::runtime_error("exiting...");
			}
		}
	}

	// There is a guarantee on max number of loops based on the max degree
	// of the graph. In my tests I rarely hit that, so I'll just pick some large number.
	const int max_iter = nodes.size();
	for( int rand_iter=0; nodeq.size()>0 && rand_iter < max_iter; ++rand_iter ){

		#pragma omp parallel
		{
			unsigned int tseed = omp_get_thread_num() * time(NULL);

			// Generate a random color
			#pragma omp for
			for( int i=0; i<n_nodes; ++i ){
				GCNode *node = &nodes[ nodeq[i] ];
				// Note about rand_r: posix function, so it may not compile in WIN.
				// Unfortunately there is no std replacement that I know of.
				int c_idx = rand_r(&tseed) % node->palette.size();
				node->color = *std::next( node->palette.begin(), c_idx );
			}

			// Conflict detection
	 		#pragma omp for
			for( int i=0; i<n_nodes; ++i ){
				GCNode *node = &nodes[ nodeq[i] ];
				int curr_c = node->color;
				node->conflict = false;

				// Check neighbors
				int n_neighbors = node->neighbors.size();
				for( int j=0; j<n_neighbors && !node->conflict; ++j ){

					// Hungarian heuristic: node with largest index keeps color
					int cn = nodes[ node->neighbors[j] ].color;
					if( curr_c == cn ){ node->conflict = true; }
				}

			} // end for-loop conflict resolution

		} // end parallel region

		// Remove color from neighbors
		for( int i=0; i<n_nodes; ++i ){
			GCNode *node = &nodes[ nodeq[i] ];
			const int nc = node->color;

			// Look at neighbors. If they are not in conflict,
			// remove their color from your own palette.
			int n_neighbors = node->neighbors.size();
			for( int j=0; j<n_neighbors; ++j ){
				GCNode *n1 = &nodes[ node->neighbors[j] ];
				if( !n1->conflict ){ node->palette.erase(n1->color); }

				// This is the part that is not thread safe:
				if( !node->conflict ){ n1->palette.erase(nc); }
			}

		} // end update neighbor palette

		// Remove colored nodes from the queue
		std::vector<int>::iterator node_iter = nodeq.begin();
		for( ; node_iter != nodeq.end(); ){
			GCNode *node = &nodes[ *node_iter ];
			if( node->conflict ){ node_iter++; } // keep
			else{ node_iter = nodeq.erase(node_iter); } // remove
		}

		// Feed the hungry
		n_nodes = nodeq.size();
		#pragma omp parallel for schedule(static)
		for( int i=0; i<n_nodes; ++i ){
			GCNode *node = &nodes[ nodeq[i] ];
			if( node->palette.size() < 2 ){
				node->palette.insert( init_palette+rand_iter );
			}
		}

	} // end color loop

	// Make sure all nodes are colored
	if( nodeq.size()>0 ){
		throw std::runtime_error("graphcolor::color Error: Nodes remain uncolored");
	}

} // end color


static inline void graphcolor::make_map( const std::vector<GCNode> &nodes, std::vector< std::vector<int> > &colors ){

	// Since the colors are random they are not necessarily 0 to n.
	// So, we'll just remove colors that don't have any nodes but start with a bunch.
	colors.clear();
	for( int i=0; i<14; ++i ){ colors.emplace_back( std::vector<int>() ); }

	int n_nodes = nodes.size();
	for( int i=0; i<n_nodes; ++i ){
		const GCNode *n = &nodes[i];
		// Add colors as needed
		while( n->color >= (int)colors.size() ){ colors.emplace_back( std::vector<int>() ); }
		colors[ n->color ].emplace_back( n->idx );
	}

	// Remove empty colors
	std::vector< std::vector<int> >::iterator it = colors.begin();
	for( ; it != colors.end(); ){ it->size() == 0 ? it = colors.erase(it) : it++; }

} // end make map
#endif

} // end namespace mcl

#endif
