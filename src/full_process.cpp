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
std::vector<float> split(std::string str, const std::string &token) {
    std::vector<float> result;
    while (!str.empty()) {
        int index = str.find(token);
        if (index != std::string::npos) {
            result.push_back(std::stof(str.substr(0, index)));
            str = str.substr(index + token.size());
            if (str.empty())result.push_back(std::stof(str));
        } else {
            result.push_back(std::stof(str));
            str = "";
        }
    }
    return result;
}

int main(int argc, char **argv) {
    // parse command line arguments
    try {
        TCLAP::CmdLine cmd("Estimates normals using PCA, see also https://github.com/tudelft3d/masbcpp", ' ', "0.1");

//        TCLAP::UnlabeledValueArg<std::string> inputArg("input", "path to input las/laz file", true, "", "input file", cmd);
//        TCLAP::UnlabeledValueArg<std::string> outputArg("output", "path to output laz file.", false, "", "output file", cmd);

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

        TCLAP::ValueArg<int> pointsToProcess("x", "ptp", "", false, 10, "int", cmd);

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

//        std::cout << "Parameters: k=" << input_parameters.k << "\n";
//        std::cout << "Parameters: denoise_preserve=" << denoise_preserveArg.getValue() << ", denoise_planar="
//                  << denoise_planarArg.getValue() << ", initial_radius=" << input_parameters.initial_radius << "\n";

        ma_data madata = {};

        PointList coords;

        int THRES = pointsToProcess.getValue();

        for (std::string line; std::getline(std::cin, line);) {
            std::vector<float> splitLine = split(line, " ");

            // Has x, y, z in line
            if (splitLine.size() == 3) {

                Point newPoint = Point(splitLine[0], splitLine[1], splitLine[2]);

                if (madata.bbox.isNull()) {
                    madata.bbox = Box(newPoint, newPoint);
                }

                coords.push_back(newPoint);
                madata.bbox.addPoint(newPoint);

            } else {
//                madata.m = splitLine[0];
//                coords.resize(madata.m);
//                normals.resize(madata.m);
            }

            // When threshold is reached; process available points and dump result to stdout
            if (coords.size() == THRES) {

                madata.m = THRES;

                VectorList normals(THRES);

                madata.normals = &normals;
                madata.coords = &coords;

//                std::cout << "PROCESSING " << coords.size() << " points" << "\n";

                compute_normals(input_parameters, madata);

                // Storage space for our results:
                PointList ma_coords(2 * madata.m);

//                std::cout << "NORMALS CALCULATED" << "\n";

                madata.coords = &coords;
                madata.normals = &normals;
                madata.ma_coords = &ma_coords;
                madata.ma_qidx = new int[2 * madata.m];

                compute_masb_points(input_parameters, madata);

//                std::cout << "MASB CALCULATED" << "\n";

                madata.lfs = new float[madata.m];
                madata.mask = new bool[madata.m];

//                std::cout << "bbox: " << madata.bbox.min[0] << " " << madata.bbox.min[1] << " " << madata.bbox.min[2] << " "
//                  << madata.bbox.max[0] << " " << madata.bbox.max[1] << " " << madata.bbox.max[2] << std::endl;

                simplify_lfs(input_parameters, madata);

//                std::cout << "SIMPLIFIED!!" << "\n";

                for (int i = 0; i < madata.m; i++) {
                    if (madata.mask[i]) {
                        std::cout << coords[i][0] << " " << coords[i][1] << " " << coords[i][2] << std::endl;
                    }
                }

                coords.resize(0);

            }

        }

//        cnpy::NpyArray coords_npy = cnpy::npy_load(input_coords_path);
//        auto *coords_carray = reinterpret_cast<float *>(coords_npy.data);
//
//        madata.m = coords_npy.shape[0];
//        PointList coords(madata.m);
//
//        for (unsigned int i = 0; i < madata.m; i++) {
//            coords[i] = Point(&coords_carray[i * 3]);
//        }

//        VectorList normals(madata.m);
//        madata.normals = &normals;
//        madata.coords = &coords;

//        // Perform the actual processing
//        compute_normals(input_parameters, madata);
//
//        // Output results
//        auto *normals_carray = new Scalar[madata.m * 3];
//        for (int i = 0; i < normals.size(); i++)
//            for (int j = 0; j < 3; j++)
//                normals_carray[i * 3 + j] = normals[i][j];
//
//        const auto c_size = (unsigned int) normals.size();
//        const unsigned int shape[] = {c_size, 3};
//
//        madata.m = coords_npy.shape[0];
//
//        // Storage space for our results:
//        PointList ma_coords(2 * madata.m);
//
//        madata.coords = &coords;
//        madata.normals = &normals;
//        madata.ma_coords = &ma_coords;
//        madata.ma_qidx = new int[2 * madata.m];
//
//        // Perform the actual processing
//        compute_masb_points(input_parameters, madata);
//
//        // Write out the results for the inside
//        auto *ma_coords_in_carray = new Scalar[madata.m * 3];
//        for (int i = 0; i < madata.m; i++)
//            for (int j = 0; j < 3; j++)
//                ma_coords_in_carray[i * 3 + j] = ma_coords[i][j];
//
//        // Write out the results for the outside
//        auto *ma_coords_out_carray = new Scalar[madata.m * 3];
//        for (int i = 0; i < madata.m; i++)
//            for (int j = 0; j < 3; j++)
//                ma_coords_out_carray[i * 3 + j] = ma_coords[i + madata.m][j];
//
//        madata.m = coords_npy.shape[0];
//        madata.bbox = Box(Point(&coords_carray[0]), Point(&coords_carray[0]));
//        for (int i = 0; i < madata.m; i++) {
//            coords[i] = Point(&coords_carray[i * 3]);
//            madata.bbox.addPoint(coords[i]);
//        }
//        coords_npy.destruct();
//        std::cout << "bbox: " << madata.bbox.min[0] << " " << madata.bbox.min[1] << " " << madata.bbox.min[2] << " "
//                  << madata.bbox.max[0] << " " << madata.bbox.max[1] << " " << madata.bbox.max[2] << std::endl;
//
//        madata.lfs = new float[madata.m];
//
//        madata.coords = &coords; // don't own this memory
//        // madata.normals = &normals;
//        madata.ma_coords = &ma_coords; // don't own this memory
//
//        madata.mask = new bool[madata.m];
//
//        // Perform the actual processing
//        simplify_lfs(input_parameters, madata);
//
//        // count number of remaining points
//        unsigned int cnt = 0;
//        for (int i = 0; i < madata.m; i++)
//            if (madata.mask[i]) cnt++;
//        std::cout << cnt << " out of " << madata.m << " points remaining [" << int(100 * float(cnt) / madata.m) << "%]" << std::endl;
//
//        // Output results
//        cnpy::npy_save(output_filtermask, madata.mask, shape, 1, "w");
//        cnpy::npy_save(output_lfs, madata.lfs, shape, 1, "w");
//
//        // Free memory
//        delete[] madata.mask;
//        madata.mask = nullptr;
//        delete[] madata.lfs;
//        madata.lfs = nullptr;
//        delete[] madata.ma_qidx;
//        madata.ma_qidx = nullptr;


    } catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

    return 0;
}
