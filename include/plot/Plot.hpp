#ifndef PLOT_HPP
#define PLOT_HPP

#include <plot/Layout.hpp>

namespace plot
{
    constexpr int to_root(Color color)
    {
        switch (color) {
            case Color::Red:     return kRed + 1;
            case Color::Blue:    return kBlue + 1;
            case Color::Green:   return kGreen + 2;
            case Color::Orange:  return kOrange + 7;
            case Color::Violet:  return kViolet + 1;
            case Color::Cyan:    return kCyan + 1;
            case Color::Magenta: return kMagenta + 1;
            case Color::Yellow:  return kYellow + 1;
            case Color::Teal:    return kTeal + 2;
            case Color::Pink:    return kPink + 1;
        }

        return kBlack;
    }

    inline constexpr double CANVAS_WIDTH = 1280;
    inline constexpr double CANVAS_HEIGHT = 720;

    void save(std::string const& path, Layout const& layout);

    template<typename... Panels>
    requires (std::same_as<std::remove_cvref_t<Panels>, Panel> && ...)
    void save(std::string const& path, Panels const& ... panels)
    {
        save(path, Layout::horizontal(panels...));
    }

    template<typename... Drawables>
    void save(std::string const& path, Drawables const& ... drawables)
    {
        save(path, Panel{ drawables... });
    }
}
#endif //PLOT_HPP