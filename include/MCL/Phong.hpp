// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_PHONG_HPP
#define MCL_PHONG_HPP 1

#include <Eigen/Core>

namespace mcl {

class Phong
{
  public:
    enum // presets
    {
        Emerald,
        Jade,
        Obsidian,
        Pearl,
        Ruby,
        Turquoise, // gems
        Brass,
        Bronze,
        Chrome,
        Copper,
        Gold,
        Silver,
        Gunmetal, // metals
        BlackPlastic,
        CyanPlastic,
        GreenPlastic,
        RedPlastic,
        WhitePlastic,
        YellowPlastic, // plastic
        BlackRubber,
        CyanRubber,
        GreenRubber,
        RedRubber,
        WhiteRubber,
        YellowRubber, // rubber
        Cloth,
        Unknown
    };

    Eigen::Vector3f amb, diff, spec;
    float shini; // In constructor, use a value from 0 to 1.

    // These should be typedef'd instead of a creation function.
    static inline Phong get(int preset);

    Phong(Eigen::Vector3f amb_, Eigen::Vector3f diff_, Eigen::Vector3f spec_, float shini_)
        : amb(amb_)
        , diff(diff_)
        , spec(spec_)
        , shini(shini_ * 128.f)
    {
    }
    Phong()
        : amb(1, 0, 0)
        , diff(1, 0, 0)
        , spec(0.3f, 0.3f, 0.3f)
        , shini(32.f)
    {
    }
};

inline Phong
Phong::get(int m)
{
    using namespace Eigen;
    Phong r;
    switch (m) {

        // Gemstones
        case Phong::Emerald:
            r = Phong(Vector3f(0.0215, 0.1745, 0.0215),
                      Vector3f(0.07568, 0.61424, 0.07568),
                      Vector3f(0.633, 0.727811, 0.633),
                      0.6);
            break;
        case Phong::Jade:
            r = Phong(Vector3f(0.135, 0.2225, 0.1575),
                      Vector3f(0.54, 0.89, 0.63),
                      Vector3f(0.316228, 0.316228, 0.316228),
                      0.1);
            break;
        case Phong::Obsidian:
            r = Phong(Vector3f(0.05375, 0.05, 0.06625),
                      Vector3f(0.18275, 0.17, 0.22525),
                      Vector3f(0.332741, 0.328634, 0.346435),
                      0.3);
            break;
        case Phong::Pearl:
            r = Phong(Vector3f(0.25, 0.20725, 0.20725),
                      Vector3f(1.0, 0.829, 0.829),
                      Vector3f(0.296648, 0.296648, 0.296648),
                      0.088);
            break;
        case Phong::Ruby:
            r = Phong(Vector3f(0.1745, 0.01175, 0.01175),
                      Vector3f(0.61424, 0.04136, 0.04136),
                      Vector3f(0.727811, 0.626959, 0.626959),
                      0.6);
            break;
        case Phong::Turquoise:
            r = Phong(Vector3f(0.1, 0.18725, 0.1745),
                      Vector3f(0.396, 0.74151, 0.69102),
                      Vector3f(0.297254, 0.30829, 0.306678),
                      0.1);
            break;

        // Metals
        case Phong::Brass:
            r = Phong(Vector3f(0.329412, 0.223529, 0.027451),
                      Vector3f(0.780392, 0.568627, 0.113725),
                      Vector3f(0.992157, 0.941176, 0.807843),
                      0.21794872);
            break;
        case Phong::Bronze:
            r = Phong(Vector3f(0.2125, 0.1275, 0.054),
                      Vector3f(0.714, 0.4284, 0.18144),
                      Vector3f(0.393548, 0.271906, 0.166721),
                      0.2);
            break;
        case Phong::Chrome:
            r = Phong(Vector3f(0.25, 0.25, 0.25), Vector3f(0.4, 0.4, 0.4), Vector3f(0.774597, 0.774597, 0.774597), 0.6);
            break;
        case Phong::Copper:
            r = Phong(Vector3f(0.19125, 0.0735, 0.0225),
                      Vector3f(0.7038, 0.27048, 0.0828),
                      Vector3f(0.256777, 0.137622, 0.086014),
                      0.6);
            break;
        case Phong::Gold:
            r = Phong(Vector3f(0.24725, 0.1995, 0.0745),
                      Vector3f(0.75164, 0.60648, 0.22648),
                      Vector3f(0.628281, 0.555802, 0.366065),
                      0.4);
            break;
        case Phong::Silver:
            r = Phong(Vector3f(0.19225, 0.19225, 0.19225),
                      Vector3f(0.50754, 0.50754, 0.50754),
                      Vector3f(0.508273, 0.508273, 0.508273),
                      0.4);
            break;
        case Phong::Gunmetal: // I made this one up
            r = Phong(Vector3f(0.1, 0.1, 0.1), Vector3f(0.4, 0.4, 0.4), Vector3f(0.1, 0.1, 0.1), 0.1);
            break;

        // Plastics
        case Phong::BlackPlastic:
            r = Phong(Vector3f(0.0, 0.0, 0.0), Vector3f(0.01, 0.01, 0.01), Vector3f(0.50, 0.50, 0.50), 0.25);
            break;
        case Phong::CyanPlastic:
            r = Phong(Vector3f(0.0, 0.1, 0.06),
                      Vector3f(0.0, 0.50980392, 0.50980392),
                      Vector3f(0.50196078, 0.50196078, 0.50196078),
                      0.25);
            break;
        case Phong::GreenPlastic:
            r = Phong(Vector3f(0.0, 0.0, 0.0), Vector3f(0.1, 0.35, 0.1), Vector3f(0.45, 0.55, 0.45), 0.25);
            break;
        case Phong::RedPlastic:
            r = Phong(Vector3f(0.0, 0.0, 0.0), Vector3f(0.5, 0.0, 0.0), Vector3f(0.7, 0.6, 0.6), 0.25);
            break;
        case Phong::WhitePlastic:
            r = Phong(Vector3f(0.0, 0.0, 0.0), Vector3f(0.55, 0.55, 0.55), Vector3f(0.70, 0.70, 0.70), 0.25);
            break;
        case Phong::YellowPlastic:
            r = Phong(Vector3f(0.0, 0.0, 0.0), Vector3f(0.5, 0.5, 0.0), Vector3f(0.60, 0.60, 0.50), 0.25);
            break;

        // Rubbers
        case Phong::BlackRubber:
            r = Phong(Vector3f(0.02, 0.02, 0.02), Vector3f(0.01, 0.01, 0.01), Vector3f(0.4, 0.4, 0.4), 0.078125);
            break;
        case Phong::CyanRubber:
            r = Phong(Vector3f(0.0, 0.05, 0.05), Vector3f(0.4, 0.5, 0.5), Vector3f(0.04, 0.7, 0.7), 0.078125);
            break;
        case Phong::GreenRubber:
            r = Phong(Vector3f(0.0, 0.05, 0.0), Vector3f(0.4, 0.5, 0.4), Vector3f(0.04, 0.7, 0.04), 0.078125);
            break;
        case Phong::RedRubber:
            r = Phong(Vector3f(0.05, 0.0, 0.0), Vector3f(0.5, 0.4, 0.4), Vector3f(0.7, 0.04, 0.04), 0.078125);
            break;
        case Phong::WhiteRubber:
            r = Phong(Vector3f(0.05, 0.05, 0.05), Vector3f(0.5, 0.5, 0.5), Vector3f(0.7, 0.7, 0.7), 0.078125);
            break;
        case Phong::YellowRubber:
            r = Phong(Vector3f(0.05, 0.05, 0.0), Vector3f(0.5, 0.5, 0.4), Vector3f(0.7, 0.7, 0.04), 0.078125);
            break;
        case Phong::Cloth:
            r = Phong(Vector3f(0.25, 0.20725, 0.20725), Vector3f(1.0, 0.829, 0.829), Vector3f(0., 0., 0.), 0.088);
            break;

        default:
            break;

    } // end switch preset

    // Apply some ambient dampening
    r.amb *= 0.2f;
    return r;

} // end presets

} // end namespace mcl

#endif
