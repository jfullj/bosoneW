#ifndef PANEL_HPP
#define PANEL_HPP

#include <plot/Drawable.hpp>
#include <functional>
namespace plot
{  

    struct LegendSettings
    {
        bool enabled = true;

        double x1 = 0.65;
        double y1 = 0.75;
        double x2 = 0.88;
        double y2 = 0.88;
    };

    struct PadSettings
    {
        double left_margin = 0.15,
            right_margin = 0.05,
            top_margin = 0.12,
            bottom_margin = 0.05;
    };

    class Panel
    {
    public:
        template<typename ...Drawables>
        requires (std::derived_from<std::remove_cvref_t<Drawables>, Drawable> && ...)
        explicit Panel(Drawables const& ... drawables) : m_drawables{std::cref(drawables)...} {}
        
        AxisSettings& x_axis();
        AxisSettings& y_axis();
        LegendSettings& legend();
        PadSettings& pad();

        AxisSettings const& x_axis() const;
        AxisSettings const& y_axis() const;
        LegendSettings const& legend() const;
        PadSettings const& pad() const;

        void draw(TPad* pad) const;
        
        static void set_default_x_title(std::string title);
        static void set_default_y_title(std::string title);
    private:
        std::vector<std::reference_wrapper<const Drawable>> m_drawables;

        mutable AxisSettings m_x_axis,
                m_y_axis;

        PadSettings m_pad;

        LegendSettings m_legend;

        inline static std::optional<std::string> s_default_x_title,
                        s_default_y_title;
    };
}

#endif //PANEL_HPP