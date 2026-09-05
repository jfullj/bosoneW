#include <plot/Plot.hpp>

#include <TCanvas.h>

void plot::save(std::string const& path, Layout const& layout)
{
    TCanvas canvas{ "", "",
        static_cast<Int_t>(CANVAS_WIDTH),
        static_cast<Int_t>(CANVAS_HEIGHT)
    };

    canvas.cd();

    std::vector<std::unique_ptr<TPad>> pads;
    layout.visit(default_rect, [&](Panel const& panel, Rect const& rect){
        auto pad{ std::make_unique<TPad>("", "", rect.x1, rect.y1, rect.x2, rect.y2) };
        
        pad->cd();
        panel.draw(pad.get());

        pads.push_back(std::move(pad));
    });

    canvas.cd();
    canvas.SaveAs(path.c_str());
}