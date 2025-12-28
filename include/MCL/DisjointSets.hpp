// Copyright Matt Overby 2025.
// Distributed under the MIT License.

#ifndef MCL_GEOM_DISJOINTSETS_HPP
#define MCL_GEOM_DISJOINTSETS_HPP 1

#include <vector>
#include <atomic>
#include <thread>

namespace mcl {

/// @brief Ranked union-find algorithm, supports parallel calls
class DisjointSets {
public:
    /// @brief Constructor
    DisjointSets(int n) : parent(n) {
        rank.resize(n, 0);
        for (int i = 0; i < n; ++i) {
            parent[i].store(i, std::memory_order_relaxed);
        }
    }

    /// @brief Get parent
    int find(int x) {
        int root = x;
        while (root != parent[root].load(std::memory_order_acquire)) {
            root = parent[root].load(std::memory_order_acquire);
        }
        while (x != root) {
            int old_parent = parent[x].load(std::memory_order_acquire);
            parent[x].store(root, std::memory_order_release);
            x = old_parent;
        }
        return root;
    }

    /// @brief Unite nodes x and y
    void make_union(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);        
        if (rootX != rootY) {
            if (rank[rootX] > rank[rootY]) {
                parent[rootY].store(rootX, std::memory_order_release);
            } else if (rank[rootX] < rank[rootY]) {
                parent[rootX].store(rootY, std::memory_order_release);
            } else {
                parent[rootY].store(rootX, std::memory_order_release);
                rank[rootX] += 1;
            }
        }
    }

protected:
    std::vector<std::atomic<int>> parent;  ///< Parent of each element
    std::vector<int> rank;  ///< Rank (or depth) of each tree
};

} // end namespace mcl

#endif // MCL_GEOM_DISJOINTSETS_HPP