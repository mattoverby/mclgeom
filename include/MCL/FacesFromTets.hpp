// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_FACESFROMTETS_HPP
#define MCL_FACESFROMTETS_HPP 1

#include <Eigen/Dense>
#include <set>
#include <vector>

namespace mcl {

// Returns just the 4 faces from a single tet
inline std::vector<Eigen::Vector3i>
faces_from_tet(const Eigen::RowVector4i& t)
{
    using namespace Eigen;
    std::vector<Vector3i> f = {
        Vector3i(t[0], t[1], t[3]), Vector3i(t[0], t[2], t[1]), Vector3i(t[0], t[3], t[2]), Vector3i(t[1], t[2], t[3])
    };
    return f;
}

// Given a tet mesh T, compute surface triangles F.
// True on success
template<typename DerivedT, typename DerivedF>
static inline bool
faces_from_tets(const Eigen::MatrixBase<DerivedT>& T, Eigen::PlainObjectBase<DerivedF>& F)
{
    using namespace Eigen;
    struct FaceKey
    {
        FaceKey()
            : f(Vector3i::Zero())
            , f_sorted(Vector3i::Zero())
        {
        }
        FaceKey(int f0, int f1, int f2)
        {
            f = Vector3i(f0, f1, f2);
            f_sorted = f;
            mcl::sort3(f_sorted[0], f_sorted[1], f_sorted[2]);
        }
        Vector3i f;
        Vector3i f_sorted;
        bool operator<(const FaceKey& other) const
        {
            for (int i = 0; i < 3; ++i) {
                if (f_sorted[i] < other.f_sorted[i]) {
                    return true;
                }
                if (f_sorted[i] > other.f_sorted[i]) {
                    return false;
                }
            }
            return false;
        }
    };

    if (T.rows() == 0 || T.cols() != 4) {
        return false;
    }

    std::set<FaceKey> faces;
    std::set<FaceKey> faces_seen_twice;
    int n_tets = T.rows();
    int total_faces = 0;
    for (int t = 0; t < n_tets; ++t) {
        total_faces += 4;
        int p0 = T(t, 0);
        int p1 = T(t, 1);
        int p2 = T(t, 2);
        int p3 = T(t, 3);
        FaceKey curr_faces[4] = { FaceKey(p0, p1, p3), FaceKey(p0, p2, p1), FaceKey(p0, p3, p2), FaceKey(p1, p2, p3) };

        for (int f = 0; f < 4; ++f) {
            const FaceKey& curr_face = curr_faces[f];
            // We've already seen the face at least twice.
            if (faces_seen_twice.count(curr_face) > 0) {
                continue;
            }
            // Check that we don't already have the face
            typename std::set<FaceKey>::iterator it = faces.find(curr_face);
            if (it != faces.end()) {
                faces.erase(it);
                faces_seen_twice.insert(curr_face);
                continue;
            }
            // Otherwise, add it to faces
            faces.insert(curr_face);
        } // end loop faces

    } // end loop tets

    int nf = faces.size();
    F.resize(nf, 3);
    typename std::set<FaceKey>::const_iterator fit = faces.begin();
    for (int f_idx = 0; fit != faces.end(); ++fit, ++f_idx) {
        F.row(f_idx) = fit->f;
    }

    return true;
}

} // ns mcl

#endif
