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

int *getCell(float x, float y, int cellSize, float minX, float minY) {
    int result[2];
    result[0] = floor((x - minX) / cellSize);
    result[1] = floor((y - minY) / cellSize);
    return result;
}

void simplifyData(int cellX,
                  int cellY,
                  full_parameters input_parameters,
                  std::vector<std::vector<std::vector<Point>>> &gridCells,
                  float offset[3]) {

    ma_data madata = {};
    PointList coords;

    std::vector<int> duplicates;

    for (int i = 0; i < gridCells[cellX][cellY].size(); i++) {
        for (int i2 = i; i2 < gridCells[cellX][cellY].size(); i2++) {
            if (gridCells[cellX][cellY][i][0] == gridCells[cellX][cellY][i2][0]
            && gridCells[cellX][cellY][i][1] == gridCells[cellX][cellY][i2][1]
            && gridCells[cellX][cellY][i][2] == gridCells[cellX][cellY][i2][2]) {
                if (i != i2) {
                    duplicates.push_back(i2);
//                    std::cout << "FOUND A DUPLICATE: " << i << " == " << i2 << std::endl;
                }
            }
        }
    }

    std::reverse(gridCells[cellX][cellY].begin(), gridCells[cellX][cellY].end());

    for (auto index : duplicates) {
        gridCells[cellX][cellY].erase(gridCells[cellX][cellY].begin() + index);
    }

    for (auto point : gridCells[cellX][cellY]) {
        if (madata.bbox.isNull()) {
//            std::cout << "Assigned new BBox" << std::endl;
            madata.bbox = Box(point, point);
        }
        else {
            madata.bbox.addPoint(point);
        }
        coords.push_back(point);
    }

//    std::cout << "Loaded coords: " << coords.size() << " points" << std::endl;

    madata.m = coords.size();

    VectorList normals(madata.m);

    madata.normals = &normals;
    madata.coords = &coords;

    // Stage 1: Compute the normals of all points currently in memory
    compute_normals(input_parameters, madata);

//    std::cout << "Computed normals" << std::endl;

    // Storage space for our results:
    PointList ma_coords(2 * madata.m);

    madata.coords = &coords;
    madata.normals = &normals;
    madata.ma_coords = &ma_coords;
    madata.ma_qidx = new int[2 * madata.m];

    // Stage 2: Compute the MASB
    compute_masb_points(input_parameters, madata);

    std::cout << "Computed MASB" << std::endl;

    madata.lfs = new float[madata.m];
    madata.mask = new bool[madata.m];

    // Stage 3: Simplify the LFS of the points
    simplify_lfs(input_parameters, madata);

    std::cout << "Simplified LFS" << std::endl;

    // Set decimal precision for floating points
    std::cout << std::setprecision(3) << std::fixed;

    // Stage 4: Release points remaining after simplification to stdout
    for (int i = 0; i < madata.m; i++) {
        if (madata.mask[i]) {
            std::cout << "v " << coords[i][0] + offset[0] << " " << coords[i][1] + offset[1] << " " << coords[i][2] + offset[2] << std::endl;
        }
    }

//    std::cout << "Output" << std::endl;

    // Free memory
    delete[] madata.mask;
    madata.mask = nullptr;
    delete[] madata.lfs;
    madata.lfs = nullptr;
    delete[] madata.ma_qidx;
    madata.ma_qidx = nullptr;

    // Delete all points that have been processed in this batch to start fresh
    coords.resize(0);
    gridCells[cellX][cellY].resize(0);

    // Clear bbox
    madata.bbox = Box();
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

        TCLAP::ValueArg<int> pointsToProcess("x", "ptp", "How many points to process per 'batch'", false, 10000, "int", cmd);

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

        int NumPointsToProcessPerBatch = pointsToProcess.getValue();

        bool sprinkling = true;
        float offset[3];

        float minX;
        float minY;
        int cellSize;

        std::vector<std::vector<std::vector<Point>>> gridCells;

        for (std::string line; std::getline(std::cin, line);) {
            std::vector<std::string> splitLine = split(line, " ");

//            std::cout << line << std::endl;

            // All sprinkle points have passed, so can safely start simplifying
            if (splitLine[0] == "#" && splitLine[1] == "endsprinkle") {
                sprinkling = false;

            } else if (splitLine[0] == "s") {
                cellSize = std::stoi(splitLine[1]);

                std::cout << line << std::endl;

            } else if (splitLine[0] == "b") {
                minX = std::stof(splitLine[1]);
                minY = std::stof(splitLine[2]);

                float maxX = std::stof(splitLine[3]);
                float maxY = std::stof(splitLine[4]);

                float minZ = std::stof(splitLine[5]);
                float maxZ = std::stof(splitLine[6]);

                offset[0] = minX + (maxX - minX) / 2;
                offset[1] = minY + (maxY - minY) / 2;
                offset[2] = minZ + (maxZ - minZ) / 2;

                std::cout << line << std::endl;

            } else if (splitLine[0] == "c") {
                int gridSize = std::stoi(splitLine[1]);

                gridCells.resize(gridSize);

                for (auto &cell : gridCells) {
                    cell.resize(gridSize);
                }

                std::cout << line << std::endl;

            } else if (splitLine[0] == "v" && !sprinkling) {

                float x = std::stof(splitLine[1]);
                float y = std::stof(splitLine[2]);
                float z = std::stof(splitLine[3]);

                int * cellId = getCell(x, y, cellSize, minX, minY);

                Point newPoint = Point(x - offset[0], y - offset[1], z - offset[2]);

                gridCells[cellId[0]][cellId[1]].push_back(newPoint);

            } else if (splitLine[0] == "x") {

                int cellX = std::stoi(splitLine[1]);
                int cellY = std::stoi(splitLine[2]);

                if (!gridCells[cellX][cellY].empty()) {
//                    std::cout << "Starting simp" << std::endl;
                    simplifyData(cellX, cellY, input_parameters, gridCells, offset);
//                    std::cout << "Ending simp" << std::endl;

                    // Empty the cell
                    gridCells[cellX][cellY].resize(0);
                }

                std::cout << line << std::endl;

            }  else {
                std::cout << line << std::endl;
            }

        }

        // Process remaining points
        // simplifyData(cellX, cellY, input_parameters, gridCells, coords, madata, offset);

    } catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

    return 0;
}
