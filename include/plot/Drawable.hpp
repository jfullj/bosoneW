#ifndef DRAWABLE_HPP
#define DRAWABLE_HPP

#include <string>
#include <optional>
#include <TPad.h>

namespace plot
{
    struct AxisSettings
    {
        std::optional<std::string> title;
        std::optional<double> min;
        std::optional<double> max;
        bool log = false;
        bool grid = false;

        double title_size = 0.05;
        double title_offset = 1.0;

        double label_size = 0.04;
        double label_offset = 0.005;
    };

    enum class Color
    {
        Red,
        Blue,
        Green,
        Orange,
        Violet,
        Cyan,
        Magenta,
        Yellow,
        Teal,
        Pink,
        Count
    };

    struct DrawProps {
        std::optional<Color> color;
        std::string label;

        int line_style = 1;
        int line_width = 2;

        int marker_style = 20;
        double marker_size = 1.0;

        bool show_in_legend = true;
    };

    struct Context
    {
        TPad* pad;
        Color color;

        bool first = false;

        AxisSettings const& x_axis;
        AxisSettings const& y_axis;
    };

    class Drawable
    {
    public:
        DrawProps& props(){ return m_props; }
        DrawProps const& props() const { return m_props; }

        virtual void draw(Context const& context) const = 0;
        virtual TObject* root_object() const = 0;
        
    private:
        DrawProps m_props;
    };
}

#endif //DRAWABLE_HPP