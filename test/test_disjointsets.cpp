// Copyright Matt Overby 2025.
// Distributed under the MIT License.
#include <MCL/AssertHandler.hpp>
#include <MCL/DisjointSets.hpp>

#include <iostream>

void
test_disjoint_set(mcl::DisjointSets& ds, int x, int y)
{
    ds.make_union(x, y);
}

int
main(int argc, char* argv[])
{
    (void)(argc);
    (void)(argv);

    const int n = 10;
    mcl::DisjointSets ds(n);

    // Make union in parallel
    std::thread t1(test_disjoint_set, std::ref(ds), 1, 2);
    std::thread t2(test_disjoint_set, std::ref(ds), 3, 4);
    std::thread t3(test_disjoint_set, std::ref(ds), 2, 3);
    std::thread t4(test_disjoint_set, std::ref(ds), 9, 8);

    t1.join();
    t2.join();
    t3.join();
    t4.join();

    for (int i = 0; i < n; ++i) {
        std::cout << "find " << i << ": " << ds.find(i) << std::endl;
    }

    mclAssert(ds.find(0) == 0);
    mclAssert(ds.find(1) == 1);
    mclAssert(ds.find(2) == 1);
    mclAssert(ds.find(3) == 1);
    mclAssert(ds.find(4) == 1);
    mclAssert(ds.find(5) == 5);
    mclAssert(ds.find(6) == 6);
    mclAssert(ds.find(7) == 7);
    mclAssert(ds.find(8) == 9);
    mclAssert(ds.find(9) == 9);

    return EXIT_SUCCESS;
}