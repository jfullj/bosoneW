#ifndef MASSSENSITIVITY_HPP
#define MASSSENSITIVITY_HPP

#include <fstream>

namespace Mass
{
    inline const std::size_t CANVAS_WIDTH = 1280;
    inline const std::size_t CANVAS_HEIGHT = 720;

    inline const char* output_file{ DATA_DIR "/results/sigma_mass.txt" };
    inline const char* output_hist_file{ DATA_DIR "/results/template_comparison_plot.png" };
    inline const char* pT_dist_dir{ DATA_DIR "/results/pT_dist.png" };

    void junk(double wmass0, double wmass1, double wwidth, std::size_t event_count);
}

#endif //MASSSENSITIVITY_HPP