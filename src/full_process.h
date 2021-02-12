#include "madata.h"

namespace masb {
    struct full_parameters {
        int k=10;
        Scalar initial_radius=200;
        bool nan_for_initr=false;
        double denoise_preserve=(PI/180.0)*20;
        double denoise_planar=(PI/180.0)*32;
        bool kd_tree_reorder=true;
        double epsilon{};
        double cellsize{};
        double bisec_threshold{};
        double elevation_threshold{};
        double maximum_density{};
        int dimension{};
        bool only_inner{};
        bool squared{};
        bool compute_lfs{};
    };

    void compute_masb_points(full_parameters &input_parameters, ma_data &madata);

    void compute_normals(full_parameters &input_parameters, ma_data &madata);

    void simplify_lfs(full_parameters &input_parameters, ma_data& madata);

}
