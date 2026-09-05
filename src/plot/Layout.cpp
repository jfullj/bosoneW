#include <plot/Layout.hpp>
#include <numeric>
    
plot::Layout::Layout(const Layout& other)
: m_dir{ other.m_dir }
, m_weights{ other.m_weights }
, m_children{ copy_children(other.m_children) } {}

plot::Layout::Layout(Layout&& other) noexcept = default;

plot::Layout& plot::Layout::operator=(const Layout& other)
{
    if(this != &other)
    {
        m_dir = other.m_dir;
        m_weights = other.m_weights;
        m_children = copy_children(other.m_children);
    }

    return *this;
}
plot::Layout& plot::Layout::operator=(Layout&& other) noexcept = default;

plot::Layout::~Layout() = default;
plot::Layout plot::Layout::blank() { return Layout{Direction::None}; }

std::vector<plot::Layout::Child> const& plot::Layout::children() const { return m_children; }
std::vector<double> const& plot::Layout::weights() const { return m_weights; }

std::vector<plot::Rect> plot::Layout::get_children_boundaries(
    Direction dir,
    Rect const& parent,
    std::vector<double> const& weights)
{
    const auto total_weight{ std::accumulate(weights.begin(), weights.end(), 0.0) };
    const auto width{ parent.x2 - parent.x1 };
    const auto height{ parent.y2 - parent.y1 };

    std::vector<Rect> rects;
    if(dir == Direction::Horizontal)
    {   
        Rect r{ 
            .y1 = parent.y1,
            .y2 = parent.y2
        };
        auto x_offset{ parent.x1 };
        for(std::size_t i{}; i < weights.size(); ++i)
        {
            const auto side{ weights[i] / total_weight * width };
            r.x1 = x_offset;
            r.x2 = x_offset + side;
            rects.push_back(r);

            x_offset += side;
        }
    }
    else if(dir == Direction::Vertical)
    {
        Rect r{ 
            .x1 = parent.x1,
            .x2 = parent.x2
        };
        auto y_offset{ parent.y1 };
        for(std::size_t i{}; i < weights.size(); ++i)
        {
            const auto side{ weights[i] / total_weight * height };
            r.y1 = y_offset;
            r.y2 = y_offset + side;
            rects.push_back(r);

            y_offset += side;
        }
    }

    return rects;
}

std::vector<plot::Layout::Child> plot::Layout::copy_children(std::vector<Child> const& children)
{
    std::vector<Child> result;
    result.reserve(children.size());
    for(auto const& child: children)
    {
        result.push_back(
            std::visit([](auto const& c) -> Child {
                if constexpr(std::is_same_v<std::remove_cvref_t<decltype(c)>, Panel>)
                    return c;
                else 
                    return std::make_unique<Layout>(*c);
            }, child)
        );
    }
    return result;
}
