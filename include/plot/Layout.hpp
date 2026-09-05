#ifndef LAYOUT_HPP
#define LAYOUT_HPP

#include <plot/Panel.hpp>
#include <vector>
#include <variant>
#include <memory>
#include <utility>
#include <concepts>
#include <type_traits>
#include <ranges>

namespace plot
{
    struct Rect
    {
        double x1, y1, x2, y2;
    };

    inline constexpr Rect default_rect = Rect{0.0, 0.0, 1.0, 1.0};

    class Layout
    {
    public:
        using Child = std::variant<Panel, std::unique_ptr<Layout>>;

        enum class Direction
        {
            Horizontal,
            Vertical, 
            None
        };

        Layout(const Layout& other);
        Layout(Layout&& other) noexcept;

        Layout& operator=(const Layout& other);
        Layout& operator=(Layout&& other) noexcept;

        ~Layout();

        template<typename... Args>
        static Layout horizontal(Args&& ... args) { return Layout{Direction::Horizontal, std::forward<Args>(args)... }; }
        template<typename... Args>
        static Layout vertical(Args&& ... args) { return Layout{Direction::Vertical, std::forward<Args>(args)... }; }
        static Layout blank();

        std::vector<Child> const& children() const;
        std::vector<double> const& weights() const;

        template<typename Callable>
        requires std::invocable<Callable, Panel const&, Rect const&>
        void visit(Rect const& rect, Callable&& callback) const
        {
            auto boundaries{ get_children_boundaries(m_dir, rect, m_weights) };
            for(auto const& [child, boundary] : std::views::zip(m_children, boundaries))
            {
                std::visit([&](auto const& c)
                {
                    if constexpr(std::is_same_v<std::remove_cvref_t<decltype(c)>, Panel>)
                        callback(c, boundary);
                    else 
                        c->visit(boundary, callback);
                }, child);
            }
        }

    private:

        template<typename... Args>
        requires ((std::same_as<std::remove_cvref_t<Args>, Panel> || std::same_as<std::remove_cvref_t<Args>, Layout>) && ...)
        Layout(Direction dir, std::vector<double> weights, Args&& ... args)
        : m_dir{ dir }
        , m_weights{ std::move(weights) }
        , m_children{ make_child(std::forward<Args>(args))... }{}

        template<typename... Args>
        requires ((std::same_as<std::remove_cvref_t<Args>, Panel> || std::same_as<std::remove_cvref_t<Args>, Layout>) && ...)
        Layout(Direction dir, Args&& ... args)
        : Layout{
            dir,
            std::vector<double>(sizeof...(Args), 1.0),
            std::forward<Args>(args)...
        }{}

        Direction m_dir;
        std::vector<double> m_weights;
        std::vector<Child> m_children;

        static std::vector<Rect> get_children_boundaries(Direction dir, Rect const& parent, std::vector<double> const& weights);

        template<typename T>
        static Child make_child(T&& t)
        {
            if constexpr(std::is_same_v<std::remove_cvref_t<T>, Panel>)
                return t;
            else if constexpr(std::is_same_v<std::remove_cvref_t<T>, Layout>)
                return std::make_unique<Layout>(std::forward<T>(t));
        }

        static std::vector<Child> copy_children(std::vector<Child> const& children);
    };
}

#endif //LAYOUT_HPP