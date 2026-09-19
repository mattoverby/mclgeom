// Copyright Matt Overby 2025.
// Distributed under the MIT License.

#ifndef MCL_GEOM_DISJOINTSETS_HPP
#define MCL_GEOM_DISJOINTSETS_HPP 1

#include <atomic>
#include <thread>
#include <utility>
#include <vector>

namespace mcl {

/// @brief Ranked union-find algorithm, supports parallel calls
class DisjointSets
{
  public:
    /// @brief Constructor
    DisjointSets(int n)
        : parent(n)
        , rank(n)
    {
        for (int i = 0; i < n; ++i) {
            parent[i].store(i, std::memory_order_relaxed);
            rank[i].store(0, std::memory_order_relaxed);
        }
    }

    /// @brief Get parent
    int find(int x)
    {
        int root = x;
        while (root != parent[root].load(std::memory_order_acquire)) {
            root = parent[root].load(std::memory_order_acquire);
        }
        while (x != root) {
            int old_parent = parent[x].load(std::memory_order_acquire);
            parent[x].compare_exchange_weak(old_parent, root, std::memory_order_release, std::memory_order_acquire);
            x = old_parent;
        }
        return root;
    }

    /// @brief Unite nodes x and y
    void make_union(int x, int y)
    {
        int rootX = find(x);
        int rootY = find(y);

        while (rootX != rootY) {
            int rankX = rank[rootX].load(std::memory_order_acquire);
            int rankY = rank[rootY].load(std::memory_order_acquire);

            if (rankX < rankY) {
                std::swap(rootX, rootY);
                std::swap(rankX, rankY);
            }

            if (rankX == rankY) {
                int expected = rankX;
                if (!rank[rootX].compare_exchange_weak(
                        expected, rankX + 1, std::memory_order_acq_rel, std::memory_order_acquire)) {
                    rootX = find(x);
                    rootY = find(y);
                    continue;
                }
            }

            int expected_root = rootY;
            if (parent[rootY].compare_exchange_weak(
                    expected_root, rootX, std::memory_order_acq_rel, std::memory_order_acquire)) {
                break;
            }

            rootX = find(x);
            rootY = find(y);
        }
    }

  protected:
    std::vector<std::atomic<int>> parent; ///< Parent of each element
    std::vector<std::atomic<int>> rank;   ///< Rank (or depth) of each tree
};

} // end namespace mcl

#endif // MCL_GEOM_DISJOINTSETS_HPP