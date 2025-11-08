// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_PHONG_HPP
#define MCL_GEOM_PHONG_HPP 1

#include <Eigen/Dense>

namespace mcl {

class Phong
{
  public:
    Eigen::Vector3f amb, diff, spec;
    float shini; // 0–128

    Phong(const Eigen::Vector3f& amb_, const Eigen::Vector3f& diff_, const Eigen::Vector3f& spec_, float shini_)
        : amb(amb_)
        , diff(diff_)
        , spec(spec_)
        , shini(shini_)
    {
    }

    // Hack to use presets from:
    // https://people.eecs.ku.edu/~jrmiller/Courses/672/InClass/3DLighting/MaterialProperties.html
    Phong(const Eigen::Vector2f& scale,
          const Eigen::Vector3f& amb_,
          const Eigen::Vector3f& diff_,
          const Eigen::Vector3f& spec_,
          float shini_)
        : amb(amb_ * scale[0])
        , diff(diff_)
        , spec(spec_)
        , shini(shini_ * scale[1])
    {
    }

    Phong()
        : amb(1, 0, 0)
        , diff(1, 0, 0)
        , spec(0.3f, 0.3f, 0.3f)
        , shini(32.f)
    {
    }

    // Gemstones
    static const Phong Emerald;
    static const Phong Jade;
    static const Phong Obsidian;
    static const Phong Pearl;
    static const Phong Ruby;
    static const Phong Turquoise;

    // Metals
    static const Phong Brass;
    static const Phong Bronze;
    static const Phong Chrome;
    static const Phong Copper;
    static const Phong Gold;
    static const Phong Silver;
    static const Phong Gunmetal;

    // Plastics
    static const Phong BlackPlastic;
    static const Phong CyanPlastic;
    static const Phong GreenPlastic;
    static const Phong RedPlastic;
    static const Phong WhitePlastic;
    static const Phong YellowPlastic;

    // Rubbers
    static const Phong BlackRubber;
    static const Phong CyanRubber;
    static const Phong GreenRubber;
    static const Phong RedRubber;
    static const Phong WhiteRubber;
    static const Phong YellowRubber;

    // Misc
    static const Phong Cloth;
};

// ---------------------------------------------
//  Definitions of all static constants
// ---------------------------------------------

// Gemstones
inline const Phong Phong::Emerald = Phong({ 0.2f, 128.f },
                                          { 0.0215f, 0.1745f, 0.0215f },
                                          { 0.07568f, 0.61424f, 0.07568f },
                                          { 0.633f, 0.727811f, 0.633f },
                                          0.6f);
inline const Phong Phong::Jade = Phong({ 0.2f, 128.f },
                                       { 0.135f, 0.2225f, 0.1575f },
                                       { 0.54f, 0.89f, 0.63f },
                                       { 0.316228f, 0.316228f, 0.316228f },
                                       0.1f);
inline const Phong Phong::Obsidian = Phong({ 0.2f, 128.f },
                                           { 0.05375f, 0.05f, 0.06625f },
                                           { 0.18275f, 0.17f, 0.22525f },
                                           { 0.332741f, 0.328634f, 0.346435f },
                                           0.3f);
inline const Phong Phong::Pearl = Phong({ 0.2f, 128.f },
                                        { 0.25f, 0.20725f, 0.20725f },
                                        { 1.0f, 0.829f, 0.829f },
                                        { 0.296648f, 0.296648f, 0.296648f },
                                        0.088f);
inline const Phong Phong::Ruby = Phong({ 0.2f, 128.f },
                                       { 0.1745f, 0.01175f, 0.01175f },
                                       { 0.61424f, 0.04136f, 0.04136f },
                                       { 0.727811f, 0.626959f, 0.626959f },
                                       0.6f);
inline const Phong Phong::Turquoise = Phong({ 0.2f, 128.f },
                                            { 0.1f, 0.18725f, 0.1745f },
                                            { 0.396f, 0.74151f, 0.69102f },
                                            { 0.297254f, 0.30829f, 0.306678f },
                                            0.1f);

// Metals
inline const Phong Phong::Brass = Phong({ 0.2f, 128.f },
                                        { 0.329412f, 0.223529f, 0.027451f },
                                        { 0.780392f, 0.568627f, 0.113725f },
                                        { 0.992157f, 0.941176f, 0.807843f },
                                        0.21794872f);
inline const Phong Phong::Bronze = Phong({ 0.2f, 128.f },
                                         { 0.2125f, 0.1275f, 0.054f },
                                         { 0.714f, 0.4284f, 0.18144f },
                                         { 0.393548f, 0.271906f, 0.166721f },
                                         0.2f);
inline const Phong Phong::Chrome =
    Phong({ 0.2f, 128.f }, { 0.25f, 0.25f, 0.25f }, { 0.4f, 0.4f, 0.4f }, { 0.774597f, 0.774597f, 0.774597f }, 0.6f);
inline const Phong Phong::Copper = Phong({ 0.2f, 128.f },
                                         { 0.19125f, 0.0735f, 0.0225f },
                                         { 0.7038f, 0.27048f, 0.0828f },
                                         { 0.256777f, 0.137622f, 0.086014f },
                                         0.6f);
inline const Phong Phong::Gold = Phong({ 0.2f, 128.f },
                                       { 0.24725f, 0.1995f, 0.0745f },
                                       { 0.75164f, 0.60648f, 0.22648f },
                                       { 0.628281f, 0.555802f, 0.366065f },
                                       0.4f);
inline const Phong Phong::Silver = Phong({ 0.2f, 128.f },
                                         { 0.19225f, 0.19225f, 0.19225f },
                                         { 0.50754f, 0.50754f, 0.50754f },
                                         { 0.508273f, 0.508273f, 0.508273f },
                                         0.4f);
inline const Phong Phong::Gunmetal =
    Phong({ 0.2f, 128.f }, { 0.1f, 0.1f, 0.1f }, { 0.4f, 0.4f, 0.4f }, { 0.1f, 0.1f, 0.1f }, 0.1f);

// Plastics
inline const Phong Phong::BlackPlastic =
    Phong({ 0.2f, 128.f }, { 0.0f, 0.0f, 0.0f }, { 0.01f, 0.01f, 0.01f }, { 0.50f, 0.50f, 0.50f }, 0.25f);
inline const Phong Phong::CyanPlastic = Phong({ 0.2f, 128.f },
                                              { 0.0f, 0.1f, 0.06f },
                                              { 0.0f, 0.50980392f, 0.50980392f },
                                              { 0.50196078f, 0.50196078f, 0.50196078f },
                                              0.25f);
inline const Phong Phong::GreenPlastic =
    Phong({ 0.2f, 128.f }, { 0.0f, 0.0f, 0.0f }, { 0.1f, 0.35f, 0.1f }, { 0.45f, 0.55f, 0.45f }, 0.25f);
inline const Phong Phong::RedPlastic =
    Phong({ 0.2f, 128.f }, { 0.0f, 0.0f, 0.0f }, { 0.5f, 0.0f, 0.0f }, { 0.7f, 0.6f, 0.6f }, 0.25f);
inline const Phong Phong::WhitePlastic =
    Phong({ 0.2f, 128.f }, { 0.0f, 0.0f, 0.0f }, { 0.55f, 0.55f, 0.55f }, { 0.70f, 0.70f, 0.70f }, 0.25f);
inline const Phong Phong::YellowPlastic =
    Phong({ 0.2f, 128.f }, { 0.0f, 0.0f, 0.0f }, { 0.5f, 0.5f, 0.0f }, { 0.60f, 0.60f, 0.50f }, 0.25f);

// Rubbers
inline const Phong Phong::BlackRubber =
    Phong({ 0.2f, 128.f }, { 0.02f, 0.02f, 0.02f }, { 0.01f, 0.01f, 0.01f }, { 0.4f, 0.4f, 0.4f }, 0.078125f);
inline const Phong Phong::CyanRubber =
    Phong({ 0.2f, 128.f }, { 0.0f, 0.05f, 0.05f }, { 0.4f, 0.5f, 0.5f }, { 0.04f, 0.7f, 0.7f }, 0.078125f);
inline const Phong Phong::GreenRubber =
    Phong({ 0.2f, 128.f }, { 0.0f, 0.05f, 0.0f }, { 0.4f, 0.5f, 0.4f }, { 0.04f, 0.7f, 0.04f }, 0.078125f);
inline const Phong Phong::RedRubber =
    Phong({ 0.2f, 128.f }, { 0.05f, 0.0f, 0.0f }, { 0.5f, 0.4f, 0.4f }, { 0.7f, 0.04f, 0.04f }, 0.078125f);
inline const Phong Phong::WhiteRubber =
    Phong({ 0.2f, 128.f }, { 0.05f, 0.05f, 0.05f }, { 0.5f, 0.5f, 0.5f }, { 0.7f, 0.7f, 0.7f }, 0.078125f);
inline const Phong Phong::YellowRubber =
    Phong({ 0.2f, 128.f }, { 0.05f, 0.05f, 0.0f }, { 0.5f, 0.5f, 0.4f }, { 0.7f, 0.7f, 0.04f }, 0.078125f);

// Misc
inline const Phong Phong::Cloth =
    Phong({ 0.2f, 128.f }, { 0.25f, 0.20725f, 0.20725f }, { 1.0f, 0.829f, 0.829f }, { 0.f, 0.f, 0.f }, 0.088f);

} // namespace mcl

#endif
