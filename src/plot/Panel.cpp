#include <plot/Panel.hpp>

#include <TLegend.h>

plot::AxisSettings& plot::Panel::x_axis() { return m_x_axis; }
plot::AxisSettings& plot::Panel::y_axis() { return m_y_axis; }
plot::LegendSettings& plot::Panel::legend() { return m_legend; }
plot::PadSettings& plot::Panel::pad() { return m_pad; }

plot::AxisSettings const& plot::Panel::x_axis() const { return m_x_axis; }
plot::AxisSettings const& plot::Panel::y_axis() const { return m_y_axis; }
plot::LegendSettings const& plot::Panel::legend() const { return m_legend; }
plot::PadSettings const& plot::Panel::pad() const { return m_pad; }


void plot::Panel::set_default_x_title(std::string title) { s_default_x_title = std::move(title); }
void plot::Panel::set_default_y_title(std::string title) { s_default_y_title = std::move(title); }

void plot::Panel::draw(TPad* rpad) const
{
    if(!m_x_axis.title.has_value() && s_default_x_title.has_value())
        m_x_axis.title = s_default_x_title.value();
    if(!m_y_axis.title.has_value() && s_default_y_title.has_value())
        m_y_axis.title = s_default_y_title.value();

    rpad->SetGridx(m_x_axis.grid);
    rpad->SetGridy(m_y_axis.grid);
    rpad->SetLogx(m_x_axis.log);
    rpad->SetLogy(m_y_axis.log);

    rpad->SetLeftMargin(pad().left_margin);
    rpad->SetRightMargin(pad().right_margin);
    rpad->SetBottomMargin(pad().bottom_margin);
    rpad->SetTopMargin(pad().top_margin);

    for(std::size_t i{}; i < m_drawables.size(); ++i)
    {
        auto& d{ m_drawables[i].get() };

        Context c{
            .pad = rpad,
            .color = static_cast<Color>(i % static_cast<int>(Color::Count)),
            .first = (i == 0),
            .x_axis = m_x_axis,
            .y_axis = m_y_axis
        };

        d.draw(c);
    }

    if(legend().enabled)
    {
        TLegend pad_legend{
            legend().x1, legend().y1,
            legend().x2, legend().y2
        };

        for(auto& d_ref : m_drawables)
        {
            auto& d{ d_ref.get() };
            if(d.props().show_in_legend)
                pad_legend.AddEntry(d.root_object(), d.props().label.c_str());
        }
        pad_legend.Draw();
    }
}
