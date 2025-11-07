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

/// @brief Stochastic graph coloring from:
/// "Fast Distributed Algorithms for Brooks-Vizing Colourings" by Grable and Panconesi (2000)
template<typename T>
inline void graph_color(const Eigen::SparseMatrix<T> &adjacency, std::vector<std::vector<int>> &colors);

class GraphColor
{
public:

	struct GraphNode {
		int index = -1;
		int color = -1;
		bool conflict = false;
		std::vector<int> neighbors;
		std::unordered_set<int> palette;
	};

	const int init_palette = 6;
	std::vector<GraphNode> graph; ///< directed graph

	/// @brief Creates an empty graph with n nodes 
	GraphColor(size_t n);

	/// @brief Links nodes i and j.
	void make_union(size_t i, size_t j);

	/// @brief Performs graph coloring.
	void color();

	/// @brief Get colors
	void get_colors(std::vector<std::vector<int>> &colors);

    /// @brief For testing, each node its own color
    void trivial();
};

//
// Implementation
//

GraphColor::GraphColor(size_t n)
{
	graph.resize(n);
	int index = 0;
	for (auto &node : graph) {
		node.index = index++;
	}
}

void GraphColor::make_union(size_t i, size_t j)
{
	if (i == j) {
		return;
	}

	// Hungarian heuristic: node with larger index keeps the color.
	// We'll get that implicitly using a directed graph, i.e.,
	// nodes only see neighbors with larget index.
	if (i > j) {
		std::swap(i, j);
	}

	graph[i].neighbors.emplace_back(j);
}

void GraphColor::color()
{
	if (graph.size() == 0) {
		return;
	}

    std::vector<int> nodeq(graph.size());
    std::iota(nodeq.begin(), nodeq.end(), 0);

    int max_iter = graph.size();
    for (int iter = 0; !nodeq.empty() && iter < max_iter; ++iter) {

        // Assign Random Colors, initialize palette if empty
        tbb::parallel_for(size_t(0), nodeq.size(), [&](size_t i) {
            auto& node = graph[nodeq[i]];
            if (node.palette.empty()) {
                for (int i=0; i<init_palette; ++i) {
                    node.palette.emplace(i);
                }
            }
            static thread_local std::mt19937 gen(std::random_device{}());
            std::uniform_int_distribution<> dist(0, node.palette.size() - 1);
            int offset = dist(gen);
            node.color = *std::next(node.palette.begin(), offset);
        });

        // Conflict Detection
        tbb::parallel_for(size_t(0), nodeq.size(), [&](size_t i) {
            auto& node = graph[nodeq[i]];
            node.conflict = false;
            for (int nbr_idx : node.neighbors) {
                const auto& nbr = graph[nbr_idx];
                if (node.color == nbr.color) {
                    node.conflict = true;
                    break;
                }
            }
        });

        // Each thread collects palette removals locally, then merges
        tbb::concurrent_vector<std::pair<int,int>> removals;
        tbb::parallel_for(size_t(0), nodeq.size(), [&](size_t i) {
            auto& node = graph[nodeq[i]];
            const int color = node.color;
            for (int nbr_idx : node.neighbors) {
                const auto& nbr = graph[nbr_idx];
                // Remove neighbor's color from self palette
				// so we don't choose it in future iterations.
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
            graph[target].palette.erase(c);
        }

        // Remove colored nodes from queue
        nodeq.erase(std::remove_if(nodeq.begin(), nodeq.end(),
            [&](int idx) { return !graph[idx].conflict; }),
            nodeq.end());

        // Refill palettes
        tbb::parallel_for(size_t(0), nodeq.size(), [&](size_t i) {
            auto& node = graph[nodeq[i]];
            if (node.palette.size() < 2) {
                node.palette.insert(init_palette + iter);
            }
        });
	}
}

void GraphColor::get_colors(std::vector<std::vector<int>> &colors)
{
	colors.clear();
    if (graph.empty()) {
		return;
	}

    // Find maximum color index so we can size the vector
    int max_color = -1;
    for (const auto& n : graph) {
        if (n.color > max_color) {
            max_color = n.color;
        }
    }

    // Preallocate all color buckets
    colors.resize(max_color + 1);

    // Fill in the buckets
    for (const auto& n : graph) {
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

void GraphColor::trivial()
{
	for (auto &node : graph) {
        node.palette.clear();
        node.color = node.index;
    }
}

template <typename T>
void graph_color(const Eigen::SparseMatrix<T> &adjacency, std::vector<std::vector<int>> &colors)
{
	if (adjacency.rows() != adjacency.cols()) {
		return;
	}

	GraphColor gc(adjacency.rows());

    // Iterate over outer dimension (columns if column-major, rows if row-major)
    for (int k = 0; k < adjacency.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(adjacency, k); it; ++it) {
            int i = it.row();
            int j = it.col();
			gc.make_union(i, j);
        }
    }

	gc.color();
	gc.get_colors(colors);
}

} // end namespace mcl

#endif
