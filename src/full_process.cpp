/*
Copyright (c) 2016 Ravi Peters

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <iostream>
#include <fstream>
#include <string>

// cnpy
#include <cnpy/cnpy.h>
// tclap
#include <tclap/CmdLine.h>

// typedefs
#include "madata.h"
#include "full_process.h"


using namespace masb;

// https://stackoverflow.com/questions/5607589/right-way-to-split-an-stdstring-into-a-vectorstring
std::vector<std::string> split(std::string str, const std::string &token) {
    std::vector<std::string> result;
    while (!str.empty()) {
        int index = str.find(token);
        if (index != std::string::npos) {
            result.push_back(str.substr(0, index));
            str = str.substr(index + token.size());
            if (str.empty())result.push_back(str);
        } else {
            result.push_back(str);
            str = "";
        }
    }
    return result;
}

void simplifyData(full_parameters input_parameters, PointList &coords, ma_data madata, float offset[3]) {

    std::vector<int> duplicates;

    for (int i = 0; i < coords.size(); i++) {
        for (int i2 = i; i2 < coords.size(); i2++) {
            if (i != i2 && coords[i][0] == coords[i2][0] && coords[i][1] == coords[i2][1] && coords[i][2] == coords[i2][2]) {
                duplicates.push_back(i2);
                std::cerr << "Found duplicate: " << i << " == " << i2 << std::endl;
            }
        }
    }

    std::cerr << "Found all duplicates! " << duplicates.size() << std::endl;

    std::cerr << "Coords size before erasing: " << coords.size() << std::endl;

    std::reverse(duplicates.begin(), duplicates.end());

    for (auto index : duplicates) {
        coords.erase(coords.begin() + index);
    }

    std::cerr << "Erased duplicates! Coords size now: " << coords.size() << std::endl;

    madata.m = coords.size();

    VectorList normals(madata.m);

    madata.normals = &normals;
    madata.coords = &coords;

    // Stage 1: Compute the normals of all points currently in memory
    compute_normals(input_parameters, madata);

    std::cerr << "Computed normals" << std::endl;

    // Storage space for our results:
    PointList ma_coords(2 * madata.m);

    madata.ma_coords = &ma_coords;
    madata.ma_qidx = new int[2 * madata.m];

    // Stage 2: Compute the MASB
    compute_masb_points(input_parameters, madata);

    std::cerr << "Computed MASB" << std::endl;

    madata.lfs = new float[madata.m];
    madata.mask = new bool[madata.m];

    // Stage 3: Simplify the LFS of the points
    simplify_lfs(input_parameters, madata);

    std::cerr << "Simplified LFS" << std::endl;

    // Stage 4: Release points remaining after simplification to stdout
    for (int i = 0; i < madata.m; i++) {
        if (madata.mask[i]) {
            std::cout << "v " << std::setprecision(9) << coords[i][0] + offset[0] << " " << coords[i][1] + offset[1] << " " << std::setprecision(4) << coords[i][2] + offset[2] << std::endl;
        }
    }

    std::cerr << "Output complete" << std::endl;

    // Free memory
    delete[] madata.mask;
    madata.mask = nullptr;
    delete[] madata.lfs;
    madata.lfs = nullptr;
    delete[] madata.ma_qidx;
    madata.ma_qidx = nullptr;

    // Clear bbox
    madata.bbox = Box();
}

int main(int argc, char **argv) {
    // parse command line arguments
    try {
        TCLAP::CmdLine cmd("Estimates normals using PCA, see also https://github.com/tudelft3d/masbcpp", ' ', "0.1");

        TCLAP::ValueArg<int> kArg("k", "kneighbors", "number of nearest neighbours to use for PCA", false, 10, "int", cmd);

        TCLAP::SwitchArg reorder_kdtreeSwitch("N", "no-kdtree-reorder", "Don't reorder kd-tree points: slower computation but lower memory use", cmd, true);

        TCLAP::ValueArg<double> denoise_preserveArg("d", "preserve", "denoise preserve threshold", false, 20, "double", cmd);
        TCLAP::ValueArg<double> denoise_planarArg("p", "planar", "denoise planar threshold", false, 32, "double", cmd);
        TCLAP::ValueArg<double> initial_radiusArg("r", "radius", "initial ball radius", false, 200, "double", cmd);

        TCLAP::SwitchArg nan_for_initrSwitch("a", "nan", "write nan for points with radius equal to initial radius", cmd, false);

        TCLAP::ValueArg<double> epsilonArg("e", "epsilon", "Control the degree of simplification, higher values mean more simplification. Typical values are in the range [0.01,0.6].", false, 0.1, "double", cmd);
        TCLAP::ValueArg<double> cellsizeArg("c", "cellsize", "Cellsize used during grid-based lfs simplification (in units of your dataset). Large cellsize means faster processing, but potentially more noticable jumps in point density at cell boundaries.", false, 1, "double", cmd);
        TCLAP::ValueArg<double> bisecArg("b", "bisec", "Bisector threshold used to clean the MAT points before LFS computation. With lower values more aggressive cleaning is performed which means more robustness to noise (in the MAT) but also less features will be detected. Typical range [0.1,10] (degrees).", false, 1, "double", cmd);
        TCLAP::ValueArg<double> maxdensArg("m", "max", "Upper bound point density in pts/unit^2", false, 1, "double", cmd);

        TCLAP::ValueArg<double> fake3dArg("f", "fake3d", "Use 2D grid instead of 3D grid, intended for 2.5D datasets (eg. buildings without points only on the roof and not on the walls). In addition this mode will try to detect elevation jumps in the dataset (eg. where there should be a wall) and still try to preserve points around those areas, the value for this parameter is the threshold elevation difference (in units of your dataset) within one gridcell that will be used for the elevation jump detection function.", false, 0.5, "double", cmd);
        TCLAP::SwitchArg squaredSwitch("s", "squared", "Use squared LFS during simplification.", cmd, false);
        TCLAP::SwitchArg nolfsSwitch("n", "no-lfs", "Don't recompute lfs.'", cmd, false);

        cmd.parse(argc, argv);

        full_parameters input_parameters;
        input_parameters.k = kArg.getValue();

        input_parameters.kd_tree_reorder = reorder_kdtreeSwitch.getValue();

        input_parameters.initial_radius = float(initial_radiusArg.getValue());
        input_parameters.denoise_preserve = (PI / 180.0) * denoise_preserveArg.getValue();
        input_parameters.denoise_planar = (PI / 180.0) * denoise_planarArg.getValue();

        input_parameters.nan_for_initr = nan_for_initrSwitch.getValue();
        input_parameters.kd_tree_reorder = reorder_kdtreeSwitch.getValue();

        input_parameters.epsilon = epsilonArg.getValue();
        input_parameters.cellsize = cellsizeArg.getValue();
        input_parameters.bisec_threshold = (bisecArg.getValue() / 180.0) * PI;

        input_parameters.compute_lfs = !nolfsSwitch.getValue();
        input_parameters.elevation_threshold = fake3dArg.getValue();
        input_parameters.maximum_density = maxdensArg.getValue();
        input_parameters.dimension = 3;
        input_parameters.only_inner = true; //innerSwitch.getValue();
        input_parameters.squared = squaredSwitch.getValue();
        if (fake3dArg.isSet())
            input_parameters.dimension = 2;

        bool sprinkling = true;
        float offset[3];

        PointList coords;
        ma_data madata;

        for (std::string line; std::getline(std::cin, line);) {
            std::vector<std::string> splitLine = split(line, " ");

            // All sprinkle points have passed, so can safely start simplifying
            if (splitLine[0] == "#" && splitLine[1] == "endsprinkle") {
                sprinkling = false;

            } else if (splitLine[0] == "b") {
                float minX = std::stof(splitLine[1]);
                float minY = std::stof(splitLine[2]);

                float maxX = std::stof(splitLine[3]);
                float maxY = std::stof(splitLine[4]);

                float minZ = std::stof(splitLine[5]);
                float maxZ = std::stof(splitLine[6]);

                offset[0] = minX + (maxX - minX) / 2;
                offset[1] = minY + (maxY - minY) / 2;
                offset[2] = minZ + (maxZ - minZ) / 2;

                std::cout << line << std::endl;

            } else if (splitLine[0] == "v" && !sprinkling) {

                float x = std::stof(splitLine[1]);
                float y = std::stof(splitLine[2]);
                float z = std::stof(splitLine[3]);

                Point newPoint = Point(x - offset[0], y - offset[1], z - offset[2]);

                if (madata.bbox.isNull()) {
                    madata.bbox = Box(newPoint, newPoint);
                }

                coords.push_back(newPoint);
                madata.bbox.addPoint(newPoint);

            } else if (splitLine[0] == "x") {

                if (!coords.empty()) {
                    simplifyData(input_parameters, coords, madata, offset);

                    // Reset coords and bbox
                    coords.resize(0);
                    madata.bbox = Box();
                }

                std::cout << line << std::endl;

            }  else {
                std::cout << line << std::endl;
            }

        }

        // Process remaining points, if any
        if (!coords.empty()) {
            simplifyData(input_parameters, coords, madata, offset);
        }

    } catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

    return 0;
}
