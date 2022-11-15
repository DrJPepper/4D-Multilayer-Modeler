#include "./graddesc.h"

GradDesc::GradDesc() {
    angTol = 0.2;
    kAng = 3.0;
    kZ = 40.0;
    maxKappa = 10.0;
    analysisSteps = 1000;
    analysisLines = 100;
    restAng = 0.0;
    isSMFBool = false;
}

void GradDesc::initNodeMatricesSMF() {
    double maxDist = 0.0;
    int shrinkage;
    double rayBox[4] = {uStart, vStart, uLength, vLength};
    auto rayTraceResults = rayTrace(inputFileName, vertices, faces, rayBox, node_count_v+1, node_count_u+1);
    MatrixXd intersectPointsNodes = std::get<2>(rayTraceResults);
    auto xPoints = std::get<0>(rayTraceResults);
    auto yPoints = std::get<1>(rayTraceResults);

    for (int z = 0; z < distanceU.rows(); z++) {
        MatrixXd splinePoints = intersectPointsNodes.block(z * xPoints, 0, xPoints, 3);
        MatrixXd tangents = MatrixXd::Zero(0, 3);
        genSpline(splinePoints, tangents, shrinkage);
        for (int i = 0; i < distanceU.cols(); i++) {
            MatrixXd hermitePoints = intersectPointsNodes.block(z * xPoints + i, 0, 2, 3);
            MatrixXd hermiteTans = tangents.block(i, 0, 2, 3);
            MatrixXd bezierPoints = hermiteToBezier(&hermitePoints, &hermiteTans);
            curvatureU(z, i) = curvatureFromUSigned(&bezierPoints, 0.5);
            distanceU(z, i) = distance(hermitePoints.row(0), hermitePoints.row(1));
            if (distanceU(z, i) > maxDist)
                maxDist = distanceU(z, i);
        }
    }
    for (int z = 0; z < distanceV.cols(); z++) {
        MatrixXd splinePoints = intersectPointsNodes(seqN(z, yPoints, xPoints), all);
        MatrixXd tangents = MatrixXd::Zero(0, 3);
        genSpline(splinePoints, tangents, shrinkage);
        for (int i = 0; i < distanceV.rows(); i++) {
            MatrixXd hermitePoints = intersectPointsNodes(seqN(i * xPoints + z, 2, xPoints), all);
            MatrixXd hermiteTans = tangents.block(i, 0, 2, 3);
            MatrixXd bezierPoints = hermiteToBezier(&hermitePoints, &hermiteTans);
            curvatureV(i, z) = curvatureFromUSigned(&bezierPoints, 0.5);
            distanceV(i, z) = distance(hermitePoints.row(0), hermitePoints.row(1));
            if (distanceV(i, z) > maxDist && distanceV(i, z) < 10.0)
                maxDist = distanceV(i, z);
        }
    }

    spacing = maxDist;
    cout << "Spacing (INM): " << spacing << endl;
    populateSheetOfNodesGeneric();
}

// Populates data structures for a sheet of nodes based on an SMF file
void GradDesc::populateSMFNodes() {
    isSMFBool = true;
    double currDist;
    int shrinkage;
    double rayBox[4] = {uStart, vStart, uLength, vLength};
    auto rayTraceResults = rayTrace(inputFileName, vertices, faces, rayBox, 0, 0);
    intersectPoints = std::get<2>(rayTraceResults);
    auto xPoints = std::get<0>(rayTraceResults);
    auto yPoints = std::get<1>(rayTraceResults);
    const int sampleCount = 10;

    smfDistancesU = MatrixXd::Constant(xPoints - 1, yPoints, -1);
    smfDistancesV = MatrixXd::Constant(xPoints, yPoints - 1, -1);
    smfCurvaturesU = MatrixXd::Constant(xPoints - 1, yPoints, 0);
    smfCurvaturesV = MatrixXd::Constant(xPoints, yPoints - 1, 0);
    curvatureMax = Array2d::Zero(2);
    distanceMax = Array2d::Zero(2);
    distanceMin = Array2d::Constant(numeric_limits<double>::max());

    for (int z = 0; z < yPoints; z++) {
        MatrixXd curvePoints = MatrixXd::Zero(sampleCount * (xPoints - 1), 3);
        MatrixXd splinePoints = intersectPoints.block(z * xPoints, 0, xPoints, 3);
        MatrixXd tangents = MatrixXd::Zero(0, 3), *temp;
        genSpline(splinePoints, tangents, shrinkage);
        currDist = 0.0;
        for (int i = 0; i < xPoints - 1; i++) {
            MatrixXd hermitePoints = intersectPoints.block(z * xPoints + i, 0, 2, 3);
            MatrixXd hermiteTans = tangents.block(i, 0, 2, 3);
            MatrixXd bezierPoints = hermiteToBezier(&hermitePoints, &hermiteTans);
            curvePoints.block(i * sampleCount, 0, sampleCount, 3) =
                QLoop(bezierPoints, sampleCount);
            smfCurvaturesU(i, z) = 0.0;
            for (int w = 0; w < 100; w++) {
                smfCurvaturesU(i, z) += curvatureFromUSigned(&bezierPoints, w / 100.0);
            }
            smfCurvaturesU(i, z) *= 0.01;
            smfDistancesU(i, z) = distance(hermitePoints.row(0), hermitePoints.row(1));
            currDist += smfDistancesU(i, z);
            // Excludes 0 curvature
            if (round(smfCurvaturesU(i, z) * 1e6))
                curvaturesVec.push_back(abs(smfCurvaturesU(i, z)));
            if (smfCurvaturesU(i, z) > curvatureMax(0)) curvatureMax(0) = smfCurvaturesU(i, z);
        }
        if (currDist > distanceMax(0))
            distanceMax(0) = currDist;
        if (currDist < distanceMin(0))
            distanceMin(0) = currDist;
        // If I don't write it this way the matrices get freed automatically
        // after we leave the function's scope
        temp = new MatrixXd;
        (*temp) = curvePoints;
        curvePointsU.push_back(temp);
    }
    for (int z = 0; z < xPoints; z++) {
        MatrixXd curvePoints = MatrixXd::Zero(sampleCount * (yPoints - 1), 3);
        MatrixXd splinePoints = intersectPoints(seqN(z, yPoints, xPoints), all);
        MatrixXd tangents = MatrixXd::Zero(0, 3), *temp;
        genSpline(splinePoints, tangents, shrinkage);
        currDist = 0.0;
        for (int i = 0; i < yPoints - 1; i++) {
            MatrixXd hermitePoints = intersectPoints(seqN(i * xPoints + z, 2, xPoints), all);
            MatrixXd hermiteTans = tangents.block(i, 0, 2, 3);
            MatrixXd bezierPoints = hermiteToBezier(&hermitePoints, &hermiteTans);
            curvePoints.block(i * sampleCount, 0, sampleCount, 3) =
                QLoop(bezierPoints, sampleCount);
            smfDistancesV(z, i) = distance(hermitePoints.row(0), hermitePoints.row(1));
            smfCurvaturesV(z, i) = curvatureFromUSigned(&bezierPoints, 0.5);
            currDist += smfDistancesV(z, i);
            curvaturesVec.push_back(abs(smfDistancesV(z, i)));
            if (smfCurvaturesV(z, i) > curvatureMax(1)) curvatureMax(1) = smfCurvaturesV(z, i);
        }
        if (currDist > distanceMax(1))
            distanceMax(1) = currDist;
        if (currDist < distanceMin(1))
            distanceMin(1) = currDist;
        temp = new MatrixXd;
        (*temp) = curvePoints;
        curvePointsV.push_back(temp);
    }
    analyzeSurfaceGeneric();
    node_u = uRows;
    node_v = vCols;
    uRows += (uRows - 1) * (node_count_u - 1);
    vCols += (vCols - 1) * (node_count_v - 1);
    initVectors();
    initNodeMatricesSMF();
}

// Populates data structures for a sheet of nodes based on a bezier patch
void GradDesc::analyzeSurfaceBezier() {
    int j;
    MatrixXd points = MatrixXd::Zero(0, 3);
    initialize(&points, inputFileName);
    // How many steps to take on a line
    const double step = 1.0 / (analysisSteps - 1);
    // How many lines to read
    const double uStep = 1.0 / (analysisLines - 1);
    double currCurvU, currDistU, currCurvV, currDistV;
    curvatureMax = Array2d::Zero(2);
    distanceMax = Array2d::Zero(2);
    distanceMin = Array2d::Constant(numeric_limits<double>::max());
    for (j = 0; j < analysisLines; j++) {
        currDistU = 0;
        currDistV = 0;
        for (int w = 0; w < analysisSteps; w++) {
            if (w < analysisSteps - 1) {
                currDistU +=
                    distanceInU(&points, w * step, (w + 1) * step, j * uStep);
                currDistV +=
                    distanceInV(&points, j * uStep, w * step, (w + 1) * step);
            }
            currCurvU = abs(curvatureInU(&points, j * uStep, w * step));
            currCurvV = abs(curvatureInV(&points, w * step, j * uStep));
            curvaturesVec.push_back(abs(currCurvU));
            curvaturesVec.push_back(abs(currCurvV));
            if (currCurvU > curvatureMax(0)) curvatureMax(0) = currCurvU;
            if (currCurvV > curvatureMax(1)) curvatureMax(1) = currCurvV;
        }
        if (currDistU > distanceMax(0)) distanceMax(0) = currDistU;
        if (currDistV > distanceMax(1)) distanceMax(1) = currDistV;
        if (currDistU < distanceMin(0)) distanceMin(0) = currDistU;
        if (currDistV < distanceMin(1)) distanceMin(1) = currDistV;
    }
    analyzeSurfaceGeneric();
}

void GradDesc::analyzeSurfaceGeneric() {
    double linearMaxShrink = 1.0 - (distanceMin / distanceMax).minCoeff();
    sort(curvaturesVec.begin(), curvaturesVec.end());
    maxKappa = curvaturesVec[floor((curvaturesVec.size() - 1) * pValue)];
    Array2d mk(1);
    mk << maxKappa, maxKappa;
    Array2d R = 1.0 / mk;
    Array2d T = R * maxShrinkage;
    Array2d nodeCount(2);
    maxZHeight = (T / (layers - 1)).minCoeff();
    minNodes = (distanceMax / 2 / R).ceil();
    aspectRatio = distanceMax(0) / distanceMax(1);
    if (aspectRatio > 1.0) {
        node_count_v = node_count_long;
        node_count_u = ceil(node_count_v / aspectRatio);
    } else {
        node_count_u = node_count_long;
        node_count_v = ceil(node_count_u * aspectRatio);
    }
    if (node_count_u < minNodes(0)) {
        cout << "\033[1;33m[WARN] Node count in U too small, defaulting to "
                "minimum value of "
             << minNodes(0) << "\033[0m\n";
        node_count_u = minNodes(0);
        node_count_v = ceil(node_count_u / aspectRatio);
    }
    if (node_count_v < minNodes(1)) {
        cout << "\033[1;33m[WARN] Node count in V too small, defaulting to "
                "minimum value of "
             << minNodes(1) << "\033[0m\n";
        node_count_v = minNodes(1);
        node_count_u = ceil(node_count_v / aspectRatio);
    }
    zHeight = min(maxZHeight, zHeight);
    if (zHeight == maxZHeight)
        cout << "\033[1;33m[WARN] Z Height too large, defaulting to maximum "
                "value of "
             << maxZHeight << "\033[0m\n";
    nodeCount << node_count_u, node_count_v;
    cout << "Max material shrinkage allowed: " << maxShrinkage << endl;
    cout << "Max in plane shrinkage required: " << linearMaxShrink << endl;
    cout << "Max Curvature in U: " << curvatureMax(0)
         << " (R: " << 1.0 / curvatureMax(0) << ")" << endl;
    cout << "Max Curvature in V: " << curvatureMax(1)
         << " (R: " << 1.0 / curvatureMax(1) << ")" << endl;
    cout << "Max Allowed Curvature (p of " << pValue << "): " << maxKappa
         << " (R: " << 1.0 / maxKappa << ")" << endl;
    cout << "Max Z Height: " << maxZHeight << endl;
    cout << "Z Height: " << maxZHeight << endl;
    cout << "Min Nodes in U: " << minNodes(0) << endl;
    cout << "Min Nodes in V: " << minNodes(1) << endl;
    cout << "U/V Aspect Ratio: " << aspectRatio << endl;
    cout << "Max Distance in U: " << distanceMax(0) << endl;
    cout << "Max Distance in V: " << distanceMax(1) << endl;
    cout << "Min Distance in U: " << distanceMin(0) << endl;
    cout << "Min Distance in V: " << distanceMin(1) << endl;
    cout << "Nodes in U: " << node_count_u << endl;
    cout << "Nodes in V: " << node_count_v << endl;
    // plotCurvHisto(curvUs);
}

void GradDesc::plotCurvHisto(vector<double> &curvs) {
    vector<double> v = curvs;
    sort(v.begin(), v.end());
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();

    vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(),
                   [mean](double x) { return x - mean; });
    double sq_sum =
        std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / v.size());

    cout << "Mean: " << mean << endl;
    cout << "StDev: " << stdev << endl;
    cout << "p of 0.5: " << v[floor(v.size() * 0.5)] << endl;
    cout << "p of 0.75: " << v[floor(v.size() * 0.75)] << endl;
    cout << "p of 0.9: " << v[floor(v.size() * 0.9)] << endl;
    cout << "p of 0.95: " << v[floor(v.size() * 0.95)] << endl;

    Gnuplot gp;
    gp << "n=100\nmin=0.\nmax=" << curvatureMax(0)
       << "\nwidth=(max-min)/n\nhist(x,width)=width*floor(x/width)+width/2.0\n"
       << "set boxwidth width*0.9\n"
       << "set style fill solid 0.5\n";

    // Data will be sent via a temporary file.  These are erased when you call
    // gp.clearTmpfiles() or when gp goes out of scope.  If you pass a filename
    // (e.g. "gp.file1d(pts, 'mydata.dat')"), then the named file will be
    // created and won't be deleted (this is useful when creating a script).
    gp << "plot" << gp.file1d(curvs)
       << "u (hist($1,width)):(1.0) smooth freq w boxes lc rgb\"green\" "
          "notitle\n";
}

// Populates data structures for a sheet of nodes based on a bezier patch
void GradDesc::populateSheetOfNodesAuto() {
    analyzeSurfaceBezier();
    node_u = uRows;
    node_v = vCols;
    uRows += (uRows - 1) * (node_count_u - 1);
    vCols += (vCols - 1) * (node_count_v - 1);
    initVectors();
    populateSheetOfNodesBezierMS();
}

// Populates data structures for a sheet of nodes based on a bezier patch
void GradDesc::populateSheetOfNodesBezierMS() {
    json entry;
    jsOut["list"].push_back(entry);
    points = MatrixXd::Zero(0, 3);
    initialize(&points, inputFileName);
    double step_u = 1.0 / distanceU.cols(), step_v = 1.0 / distanceV.rows();
    double currU, currV;
    double maxDist = 0.0;
    for (int i = 0; i < distanceV.cols(); i++) {
        for (int j = 0; j < distanceU.rows(); j++) {
            if (i < distanceU.cols() && j < distanceU.rows()) {
                distanceU(j, i) = 0.0;
                curvatureU(j, i) = 0.0;
            }
            if (j < distanceV.rows() && i < distanceV.cols()) {
                distanceV(j, i) = 0.0;
                curvatureV(j, i) = 0.0;
            }
            for (int w = 0; w < 100; w++) {
                if (i < distanceU.cols() && j < distanceU.rows()) {
                    currU = i * step_u + step_u * 0.01 * w + step_u * 0.005;
                    curvatureU(j, i) += curvatureInU(&points, currU, j * step_v);
                    if (w < 99)
                        distanceU(j, i) += distanceInU(&points, currU,
                                currU + step_u * 0.01,
                                j * step_v) / (node_u - 1);
                }
                if (j < distanceV.rows() && i < distanceV.cols()) {
                    currV = j * step_v - 0.5 * step_v + step_v * 0.01 * w + step_v * 0.005;
                    saveToJS = w == 50;
                    curvatureV(j, i) += curvatureInV(&points, i * step_u, currV);
                    saveToJS = false;

                    if (w < 99)
                        distanceV(j, i) += distanceInV(&points, i * step_u, currV,
                                currV + step_v * 0.01) / (node_v - 1);
                }
            }
            if (i < distanceU.cols() && j < distanceU.rows()) {
                if (distanceU(j, i) > maxDist) maxDist = distanceU(j, i);
                curvatureU(j, i) *= 0.01;
            }
            if (j < distanceV.rows() && i < distanceV.cols()) {
                if (distanceV(j, i) > maxDist) maxDist = distanceV(j, i);
                curvatureV(j, i) *= 0.01;
            }
        }
    }
    spacing = maxDist;
    cout << "Spacing (MP): " << spacing << endl;
    populateSheetOfNodesGeneric();
}

void GradDesc::populateSheetOfNodesGeneric() {
    double L0X, L0Y;

    layers = 1;

    radius = spacing / node_count_u;

    for (int i = 0; i < positionsVec[0]->rows(); i++) {
        positionsVec[0]->row(i) << (i / vCols) * spacing, (i % vCols) * spacing,
            0;
    }
    for (int i = 0; i < distanceV.cols(); i++) {
        for (int j = 0; j < distanceU.rows(); j++) {
            if (i < distanceU.cols() && j < distanceU.rows()) {
                L0X = distanceU(j, i) * (node_u - 1);
                rest_uVec[0]->block(j * (node_v - 1), i * (node_u - 1),
                                    node_v - 1, node_u - 1) =
                    MatrixXd::Constant(node_v - 1, node_u - 1,
                                       L0X / (node_u - 1));
            }
            if (j < distanceV.rows() && i < distanceV.cols()) {
                L0Y = distanceV(j, i) * (node_v - 1);
                rest_vVec[0]->block(j * (node_v - 1), i * (node_u - 1),
                                    node_v - 1, node_u - 1) =
                    MatrixXd::Constant(node_v - 1, node_u - 1,
                                       L0Y / (node_v - 1));
            }
        }
    }
    stringDeq.push_back("Run in plane optimization");
    void (GradDesc::*optFunc)() = &GradDesc::shapeOptimizeGrid;
    funDeq.push_back(optFunc);
    stringDeq.push_back("Duplicate L0");
    funDeq.push_back(&GradDesc::duplicateL0);
    stringDeq.push_back("Run full optimization");
    funDeq.push_back(optFunc);
}

void GradDesc::duplicateL0() {
    layers = static_cast<int>(positionsVec.size());
    int i, j, n;
    for (n = 1; n < layers; n++) {
        for (i = 0; i < positionsVec[0]->rows(); i++) {
            positionsVec[n]->row(i) << positionsVec[0]->row(i)(0),
                positionsVec[0]->row(i)(1), n * zHeight;
        }
    }
    double L0X, L0Y, maxCurvShrinkage = 0.0;;
    double LnX, RX, RXmin;
    double LnY, RY, RYmin;
    radius = spacing / node_count_u;

    double angleX, angleY, kappaX, kappaY;
    int index;
    for (i = 0; i < distanceV.cols(); i++) {
        for (j = 0; j < distanceU.rows(); j++) {
            if (i < distanceU.cols() && j < distanceU.rows()) {
                L0X = distanceU(j, i) * (node_u - 1);
                kappaX = ((curvatureU(j, i) > 0) - (curvatureU(j, i) < 0)) *
                         min(abs(curvatureU(j, i)), maxKappa);
                RXmin = (layers - 1) * zHeight / maxShrinkage;
                if (kappaX == 0.0) {
                    RX = -1.0;
                    angleX = 0.0;
                } else {
                    RX = 1 / abs(kappaX);
                    RX = RX < RXmin ? RXmin : RX;
                    angleX = ((kappaX > 0) - (kappaX < 0)) *
                             asin(L0X / 2.0 / RX) * 2.0;
                }
            }
            if (j < distanceV.rows() && i < distanceV.cols()) {
                L0Y = distanceV(j, i) * (node_v - 1);
                kappaY = ((curvatureV(j, i) > 0) - (curvatureV(j, i) < 0)) *
                         min(abs(curvatureV(j, i)), maxKappa);
                RYmin = (layers - 1) * zHeight / maxShrinkage;
                if (kappaY == 0.0) {
                    RY = -1.0;
                    angleY = 0.0;
                } else {
                    RY = 1 / abs(kappaY);
                    RY = RY < RYmin ? RYmin : RY;
                    angleY = ((kappaY > 0) - (kappaY < 0)) *
                             asin(L0Y / 2.0 / RY) * 2.0;
                }
            }
            for (n = 0; n < layers; n++) {
                if (i < distanceU.cols() && j < distanceU.rows()) {
                    if (RX < 0.0)
                        LnX = L0X;
                    else {
                        LnX = L0X - 2 * n * zHeight * sin(abs(angleX) / 2);
                        double shrink = (L0X - LnX) / L0X;
                        if (shrink > maxCurvShrinkage)
                            maxCurvShrinkage = shrink;
                    }
                    if (angleX < 0.0) {
                        index = layers - n - 1;
                    } else {
                        index = n;
                    }
                    rest_uVec[index]->block(j * (node_v - 1),
                                            i * (node_u - 1), node_v - 1,
                                            node_u - 1) =
                        MatrixXd::Constant(node_v - 1, node_u - 1,
                                            LnX / (node_u - 1));
                }
                if (j < distanceV.rows() && i < distanceV.cols()) {
                    if (RX < 0.0)
                        LnY = L0Y;
                    else {
                        LnY = L0Y - 2 * n * zHeight * sin(abs(angleY) / 2);
                        double shrink = (L0Y - LnY) / L0Y;
                        if (shrink > maxCurvShrinkage)
                            maxCurvShrinkage = shrink;
                    }
                    if (angleY < 0.0) {
                        index = layers - n - 1;
                    } else {
                        index = n;
                    }
                    rest_vVec[index]->block(j * (node_v - 1),
                                            i * (node_u - 1), node_v - 1,
                                            node_u - 1) =
                        MatrixXd::Constant(node_v - 1, node_u - 1,
                                            LnY / (node_v - 1));
                }
            }
        }
    }
    cout << "Max Curvature Shrinkage: " << maxCurvShrinkage << endl;
}

// Generic input parsing for all configuration types
void GradDesc::populateMatrices(ifstream *inFile) {
    inputFileName = "patchPoints.txt";
    enum InputType { BEZNODESAUTO, BEZNODESMS, SMF, NONE };
    InputType inputType = NONE;
    string line;
    spacing = 0.0;
    kUV = 5.0;
    pValue = 0.9;
    uLength = 0.0;
    vLength = 0.0;
    uStart = 0.0;
    vStart = 0.0;
    const regex grid_regex("^grid:\\s*([0-9]+),\\s*([0-9]+),\\s*([0-9]+)\\s*$");
    const regex ray_start_regex("^ray_start:\\s*([0-9]+\\.[0-9]+),\\s*([0-9]+\\.[0-9]+)\\s*$");
    const regex ray_box_regex("^ray_box_size:\\s*([0-9]+\\.[0-9]+),\\s*([0-9]+\\.[0-9]+)\\s*$");
    const regex spacing_regex("^spacing:\\s*([0-9]+)\\s*$");
    // const regex angular_const("^spring_constants_angular:\\s*([0-9]+)\\s*$");
    const regex radius_regex("^radius:\\s*([0-9]+)\\s*$");
    const regex k_regex("^spring_constant:\\s*(.+)\\s*$");
    const regex z_regex("^z_height:\\s*([0-9]+\\.[0-9]+)\\s*$");
    const regex node_count_u_regex("^node_count(_u)?:\\s*([0-9]+)\\s*$");
    const regex node_count_long_regex(
        "^long_edge_node_count:\\s*([0-9]+)\\s*$");
    const regex node_count_regex("^node_count:\\s*([0-9]+),\\s*([0-9]+)\\s*$");
    const regex type_regex("^type:\\s*(.+)$");
    const regex patch_regex("^patch:\\s*(.+)$");
    const regex smf_regex("^model:\\s*(.+)$");
    const regex k_u_regex("^\\s*spring_constants_u:\\s*$");
    const regex k_v_regex("^\\s*spring_constants_v:\\s*$");
    const regex rest_u_regex("^\\s*rest_lengths_u:\\s*$");
    const regex rest_v_regex("^\\s*rest_lengths_v:\\s*$");
    const regex layer_regex("^layer\\s*\\{$");
    const regex bracket_regex("^\\s*\\}\\s*$");
    const regex shrinkage_regex("^max_shrinkage:\\s*([0-9]+\\.[0-9]+)\\s*$");
    const regex pValue_regex("^curvature_p_value:\\s*([0-9]+\\.[0-9]+)\\s*$");
    const regex granularity_regex(
        "^node_granularity:\\s*([0-9]+\\.[0-9]+)\\s*$");
    smatch sm;
    getline(*inFile, line);
    while (!inFile->eof()) {
        if (line.front() == '#' || !line.compare("")) {
        } else if (regex_search(line, sm, patch_regex)) {
            inputFileName = sm[1];
        } else if (regex_search(line, sm, smf_regex)) {
            inputFileName = sm[1];
        } else if (regex_search(line, sm, type_regex)) {
            if (!sm[1].compare("beznodes_multistep")) {
                inputType = BEZNODESMS;
            } else if (!sm[1].compare("beznodes_auto")) {
                inputType = BEZNODESAUTO;
            } else if (!sm[1].compare("smfnodes")) {
                inputType = SMF;
            } else {
                cout << "Error: invalid type declaration: " << sm[1] << endl;
                exit(1);
            }
        } else if (regex_search(line, sm, grid_regex)) {
            uRows = stoi(sm[1]);
            vCols = stoi(sm[2]);
            layers = stoi(sm[3]);
            if (inputType == BEZNODESMS) {
                node_u = uRows;
                uRows += (uRows - 1) * (node_count_u - 1);
                node_v = vCols;
                vCols += (vCols - 1) * (node_count_v - 1);
                initVectors();
            }
        } else if (regex_search(line, sm, ray_start_regex)) {
            uStart = boost::lexical_cast<double>(sm[1]);
            vStart = boost::lexical_cast<double>(sm[2]);
        } else if (regex_search(line, sm, ray_box_regex)) {
            uLength = boost::lexical_cast<double>(sm[1]);
            vLength = boost::lexical_cast<double>(sm[2]);
        } else if (regex_search(line, sm, spacing_regex)) {
            spacing = static_cast<double>(stoi(sm[1]));
        } else if (regex_search(line, sm, radius_regex)) {
            radius = static_cast<double>(stoi(sm[1]));
        } else if (regex_search(line, sm, k_regex)) {
            k_in = boost::lexical_cast<double>(sm[1]);
            kUV = k_in;
        } else if (regex_search(line, sm, z_regex)) {
            zHeight = static_cast<double>(stod(sm[1]));
        } else if (regex_search(line, sm, node_count_u_regex)) {
            node_count_u = static_cast<int>(stoi(sm[2]));
        } else if (regex_search(line, sm, node_count_long_regex)) {
            node_count_long = static_cast<int>(stoi(sm[1]));
        } else if (regex_search(line, sm, node_count_regex)) {
            node_count_u = static_cast<int>(stoi(sm[1]));
            node_count_v = static_cast<int>(stoi(sm[2]));
        } else if (regex_search(line, sm, shrinkage_regex)) {
            maxShrinkage = static_cast<double>(stod(sm[1]));
        } else if (regex_search(line, sm, pValue_regex)) {
            pValue = static_cast<double>(stod(sm[1]));
        } else if (regex_search(line, sm, granularity_regex)) {
            nodeGranularity = static_cast<double>(stod(sm[1]));
        } else {
            cout << "Error: no matching pattern for input line: " << line
                 << endl;
            exit(1);
        }
        getline(*inFile, line);
    }
    if (inputType == SMF) {
        populateSMFNodes();
    } else if (inputType == BEZNODESMS) {
        populateSheetOfNodesBezierMS();
    } else if (inputType == BEZNODESAUTO) {
        populateSheetOfNodesAuto();
    } else {
        cout << "Error: No input type was specified" << endl;
        exit(1);
    }
}

void GradDesc::initVectors() {
    MatrixXd *temp;
    distanceU = MatrixXd::Constant(uRows, vCols - 1, -1);
    curvatureU = MatrixXd::Constant(uRows, vCols - 1, -1);
    distanceV = MatrixXd::Constant(uRows - 1, vCols, -1);
    curvatureV = MatrixXd::Constant(uRows - 1, vCols, -1);
    for (int i = 0; i < layers; i++) {
        temp = new MatrixXd;
        (*temp) = MatrixXd::Zero(uRows * vCols, 3);
        positionsVec.push_back(temp);
        temp = new MatrixXd;
        (*temp) = MatrixXd::Constant(uRows, vCols - 1, -1);
        rest_uVec.push_back(temp);
        temp = new MatrixXd;
        (*temp) = MatrixXd::Constant(uRows - 1, vCols, -1);
        rest_vVec.push_back(temp);
    }
}

void GradDesc::updateVtkGrid(Scene &scene) {
    vtkSphereSource *currSphere;
    vtkLineSource *currLine;
    int row, col;
    int count = 0, cylCount = 0;
    vtkNew<vtkNamedColors> colors;

    for (int layer = 0; layer < layers; layer++) {
        auto pv = (*positionsVec[layer]);
        for (int i = 0; i < positionsVec[0]->rows(); i++) {
            currSphere = scene.spheres[count];
            count++;
            row = i / vCols;
            col = i % vCols;
            currSphere->SetCenter(pv(i, 0), pv(i, 1), pv(i, 2));
            if (col != vCols - 1) {
                Vector3d color = ratioToRGB(distance(pv.row(i), pv.row(i + 1)) /
                                            (*rest_uVec[layer])(row, col));
                scene.lineActors[cylCount]->GetProperty()->SetColor(
                    color(0), color(1), color(2));
                currLine = scene.lines[cylCount];
                currLine->SetPoint1(pv(i, 0), pv(i, 1), pv(i, 2));
                currLine->SetPoint2(pv(i + 1, 0), pv(i + 1, 1), pv(i + 1, 2));
                cylCount++;
            }
            if (row != uRows - 1) {
                Vector3d color =
                    ratioToRGB(distance(pv.row(i), pv.row(i + vCols)) /
                               (*rest_vVec[layer])(row, col));
                scene.lineActors[cylCount]->GetProperty()->SetColor(
                    color(0), color(1), color(2));
                currLine = scene.lines[cylCount];
                currLine->SetPoint1(pv(i, 0), pv(i, 1), pv(i, 2));
                currLine->SetPoint2(pv(i + vCols, 0), pv(i + vCols, 1),
                                    pv(i + vCols, 2));
                cylCount++;
            }
            if (layer != 0) {
                Vector3d color = ratioToRGB(
                    distance(pv.row(i), positionsVec[layer - 1]->row(i)) /
                    zHeight);
                scene.lineActors[cylCount]->GetProperty()->SetColor(
                    color[0], color[1], color[2]);
                currLine = scene.lines[cylCount];
                currLine->SetPoint1(pv(i, 0), pv(i, 1), pv(i, 2));
                currLine->SetPoint2((*positionsVec[layer - 1])(i, 0),
                                    (*positionsVec[layer - 1])(i, 1),
                                    (*positionsVec[layer - 1])(i, 2));
                cylCount++;
            }
        }
    }
}

double GradDesc::getRestU(int row, int col, int layer) {
    return (*rest_uVec[layer])(row, col);
}

double GradDesc::getRestV(int row, int col, int layer) {
    return (*rest_vVec[layer])(row, col);
}

Array2d GradDesc::getAvgSpringDisp() {
    Vector3d currPt, otherPt, mVec, ldir, rdir, ddir, udir, zddir, zudir;
    Array2d disp;
    disp << 0.0, 0.0;
    Array2i count;
    count << 0, 0;
    double rest_z = zHeight, currRest;
    int row, col, layer, i;
    // Could be handled more flexibly but currently Z stuff is constant for all
    // input types
    bool left, right, up, down, zup, zdown;
    for (layer = 0; layer < layers; layer++) {
        // Loop through all nodes
        for (i = 0; i < positionsVec[layer]->rows(); i++) {
            row = i / vCols;
            col = i % vCols;
            // Boolean checks for what neighbors we have
            left = col != 0;
            right = col != vCols - 1;
            up = row != 0;
            down = row != uRows - 1;
            zup = layer != 0;
            zdown = layer != layers - 1;
            currPt = positionsVec[layer]->row(i);
            mVec << 0.0, 0.0, 0.0;
            // Left
            if (left) {
                otherPt = positionsVec[layer]->row(i - 1);
                ldir = (otherPt - (currPt)).normalized();
                currRest = (*rest_uVec[layer])(row, col - 1);
                disp[0] += abs(distance(otherPt, currPt) - currRest) / currRest;
                count[0]++;
            }
            // Right
            if (right) {
                otherPt = positionsVec[layer]->row(i + 1);
                rdir = (otherPt - (currPt)).normalized();
                currRest = (*rest_uVec[layer])(row, col);
                disp[0] += abs(distance(otherPt, currPt) - currRest) / currRest;
                count[0]++;
            }
            // Up
            if (up) {
                otherPt = positionsVec[layer]->row(i - vCols);
                udir = (otherPt - (currPt)).normalized();
                currRest = (*rest_vVec[layer])(row - 1, col);
                disp[0] += abs(distance(otherPt, currPt) - currRest) / currRest;
                count[0]++;
            }
            // Down
            if (down) {
                otherPt = positionsVec[layer]->row(i + vCols);
                ddir = (otherPt - (currPt)).normalized();
                currRest = (*rest_vVec[layer])(row, col);
                disp[0] += abs(distance(otherPt, currPt) - currRest) / currRest;
                count[0]++;
            }
            // Z Up
            if (zup) {
                otherPt = positionsVec[layer - 1]->row(i);
                zudir = (otherPt - (currPt)).normalized();
                currRest = rest_z;
                disp[0] += abs(distance(otherPt, currPt) - currRest) / currRest;
                count[0]++;
            }
            // Z Down
            if (zdown) {
                otherPt = positionsVec[layer + 1]->row(i);
                zddir = (otherPt - (currPt)).normalized();
                currRest = rest_z;
                disp[0] += abs(distance(otherPt, currPt) - currRest) / currRest;
                count[0]++;
            }
            // Add in all the angle values
            if (zup && left) {
                disp[1] += abs(acos(zudir.dot(ldir)) - H_PI) / H_PI;
                count[1]++;
            }
            if (zup && right) {
                disp[1] += abs(acos(zudir.dot(rdir)) - H_PI) / H_PI;
                count[1]++;
            }
            if (zup && up) {
                disp[1] += abs(acos(zudir.dot(udir)) - H_PI) / H_PI;
                count[1]++;
            }
            if (zup && down) {
                disp[1] += abs(acos(zudir.dot(ddir)) - H_PI) / H_PI;
                count[1]++;
            }
            if (zdown && left) {
                disp[1] += abs(acos(zddir.dot(ldir)) - H_PI) / H_PI;
                count[1]++;
            }
            if (zdown && right) {
                disp[1] += abs(acos(zddir.dot(rdir)) - H_PI) / H_PI;
                count[1]++;
            }
            if (zdown && up) {
                disp[1] += abs(acos(zddir.dot(udir)) - H_PI) / H_PI;
                count[1]++;
            }
            if (zdown && down) {
                disp[1] += abs(acos(zddir.dot(ddir)) - H_PI) / H_PI;
                count[1]++;
            }
            if (up && left) {
                disp[1] += abs(acos(udir.dot(ldir)) - H_PI) / H_PI;
                count[1]++;
            }
            if (down && left) {
                disp[1] += abs(acos(ddir.dot(ldir)) - H_PI) / H_PI;
                count[1]++;
            }
            if (up && right) {
                disp[1] += abs(acos(udir.dot(rdir)) - H_PI) / H_PI;
                count[1]++;
            }
            if (down && right) {
                disp[1] += abs(acos(ddir.dot(rdir)) - H_PI) / H_PI;
                count[1]++;
            }
        }
    }
    disp[0] /= count[0];
    disp[1] /= count[1];
    return disp;
}

double GradDesc::getSysEnergy() {
    Vector3d currPt, otherPt, ldir, rdir, ddir, udir, zddir, zudir;
    double energy;
    double systemEnergy = 0.0;
    double rest_z = zHeight;
    int row, col, layer, i;
    // Could be handled more flexibly but currently Z stuff is constant for all
    // input types
    bool left, right, up, down, zup, zdown;
    for (layer = 0; layer < layers; layer++) {
        // Loop through all nodes
        for (i = 0; i < positionsVec[layer]->rows(); i++) {
            // TODO: Angular values temporarily disabled for this manual energy
            // calculation (still enabled in ShapeOp)
            energy = 0.0;
            row = i / vCols;
            col = i % vCols;
            // Boolean checks for what neighbors we have
            left = col != 0;
            right = col != vCols - 1;
            up = row != 0;
            down = row != uRows - 1;
            zup = layer != 0;
            zdown = layer != layers - 1;
            currPt = positionsVec[layer]->row(i);
            // Left
            if (left) {
                otherPt = positionsVec[layer]->row(i - 1);
                ldir = (otherPt - (currPt)).normalized();
                energy += 0.5 * kUV *
                          pow2(distance(otherPt, currPt) -
                               (*rest_uVec[layer])(row, col - 1));
            }
            // Right
            if (right) {
                otherPt = positionsVec[layer]->row(i + 1);
                rdir = (otherPt - (currPt)).normalized();
                energy += 0.5 * kUV *
                          pow2(distance(otherPt, currPt) -
                               (*rest_uVec[layer])(row, col));
            }
            // Up
            if (up) {
                otherPt = positionsVec[layer]->row(i - vCols);
                udir = (otherPt - (currPt)).normalized();
                energy += 0.5 * kUV *
                          pow2(distance(otherPt, currPt) -
                               (*rest_vVec[layer])(row - 1, col));
            }
            // Down
            if (down) {
                otherPt = positionsVec[layer]->row(i + vCols);
                ddir = (otherPt - (currPt)).normalized();
                energy += 0.5 * kUV *
                          pow2(distance(otherPt, currPt) -
                               (*rest_vVec[layer])(row, col));
            }
            // Z Up
            if (zup) {
                otherPt = positionsVec[layer - 1]->row(i);
                zudir = (otherPt - (currPt)).normalized();
                energy += 0.5 * kZ * pow2(distance(otherPt, currPt) - rest_z);
            }
            // Z Down
            if (zdown) {
                otherPt = positionsVec[layer + 1]->row(i);
                zddir = (otherPt - (currPt)).normalized();
                energy += 0.5 * kZ * pow2(distance(otherPt, currPt) - rest_z);
            }
            // Add in all the angle values
            /*if (zup && left) {
                energy += kAng *
                          (1.0 / (1.0 - pow2(zudir.dot(ldir) - restAng)) - 1.0);
            }
            if (zup && right) {
                energy += kAng *
                          (1.0 / (1.0 - pow2(zudir.dot(rdir) - restAng)) - 1.0);
            }
            if (zup && up) {
                energy += kAng *
                          (1.0 / (1.0 - pow2(zudir.dot(udir) - restAng)) - 1.0);
            }
            if (zup && down) {
                energy += kAng *
                          (1.0 / (1.0 - pow2(zudir.dot(ddir) - restAng)) - 1.0);
            }
            if (zdown && left) {
                energy += kAng *
                          (1.0 / (1.0 - pow2(zddir.dot(ldir) - restAng)) - 1.0);
            }
            if (zdown && right) {
                energy += kAng *
                          (1.0 / (1.0 - pow2(zddir.dot(rdir) - restAng)) - 1.0);
            }
            if (zdown && up) {
                energy += kAng *
                          (1.0 / (1.0 - pow2(zddir.dot(udir) - restAng)) - 1.0);
            }
            if (zdown && down) {
                energy += kAng *
                          (1.0 / (1.0 - pow2(zddir.dot(ddir) - restAng)) - 1.0);
            }
            if (up && left) {
                energy +=
                    kAng * (1.0 / (1.0 - pow2(udir.dot(ldir) - restAng)) - 1.0);
            }
            if (down && left) {
                energy +=
                    kAng * (1.0 / (1.0 - pow2(ddir.dot(ldir) - restAng)) - 1.0);
            }
            if (up && right) {
                energy +=
                    kAng * (1.0 / (1.0 - pow2(udir.dot(rdir) - restAng)) - 1.0);
            }
            if (down && right) {
                energy +=
                    kAng * (1.0 / (1.0 - pow2(ddir.dot(rdir) - restAng)) - 1.0);
            }*/
            systemEnergy += energy;
        }
    }
    cout << "(Linear Only) ";
    return systemEnergy;
}

int GradDesc::getLayers() {
    return layers;
}

int GradDesc::getRows() {
    return uRows;
}

int GradDesc::getCols() {
    return vCols;
}

MatrixXd GradDesc::getVertices() {
    return vertices;
}

MatrixXi GradDesc::getFaces() {
    return faces;
}

double GradDesc::getSpacing() {
    return spacing;
}

void GradDesc::shapeOptimizeGrid() {
    shapeOptimizeGrid(15000);
}

MatrixXd GradDesc::getPVLayer(int layer) {
    return *positionsVec[layer];
}

MatrixXd GradDesc::getMiddleLayer() {
    int index = round((positionsVec.size() - 1) / 2.0);
    return *positionsVec[index];
}

void GradDesc::assessResults() {
    Vector3d vec1, vec2, ptRight, ptDown, ptDiag, center, currPt, pt1, pt2;
    int row, col, i, countU, countV;
    MatrixXd positions = getMiddleLayer();
    int topBottomLayers[2] = {0, positionsVec.size()-1};
    bool right, down;
    double uGoal/*, vGoal*/, uError, distUFin, distVFin, distUGoal,
           distVGoal, theta, kappa, L[2];
    double totalUError = 0.0, totalUErrorLn = 0.0;
    double totalUErrorCurv = 0.0;
    countU = 0;
    countV = 0;
    for (i = 0; i < positions.rows(); i++) {
        row = i / vCols;
        col = i % vCols;
        // Boolean checks for what neighbors we have
        right = col < vCols - 2;
        down = row < uRows - 2;
        currPt = positions.row(i);
        if (down && right) {
            if (row < curvatureU.cols() && col < curvatureU.rows()) {
                distUGoal = distanceU(row, col);
                uGoal = abs(curvatureU(row, col));
            }
            if (col < curvatureV.rows() && row < curvatureV.cols()) {
                distVGoal = distanceV(row, col);
            }

            for (int j = 0; j < 2; j++) {
                MatrixXd l = *positionsVec[topBottomLayers[j]];
                if (row < curvatureU.cols() && col < curvatureU.rows()) {
                    L[j] = distance(l.row(i), l.row(i+1));
                }
            }
            double L0 = max(L[0], L[1]);
            double Ln = min(L[0], L[1]);
            double theta2 = 2.0 * asin(L0 * uGoal / 2.0);
            double LnGoal = L0 - 2.0 * (layers - 1) * zHeight * sin(theta2 / 2.0);

            theta = 2.0 * asin(abs(L[0] - L[1]) / (2 * (layers - 1) * zHeight));
            kappa = 2.0 * sin(theta / 2.0) / L0;
            uError = abs(uGoal - kappa) / kappa;
            totalUErrorCurv += uError;

            if (row < curvatureU.cols() && col < curvatureU.rows()) {
                distUFin = max(L[0], L[1]);
                uError = (distUGoal - distUFin) / distUFin;
                totalUError += uError;
                totalUErrorLn += abs(LnGoal - Ln) / Ln;
                countU++;
            }
            if (col < curvatureV.rows() && row < curvatureV.cols()) {
                distVFin = distance(currPt, positions.row(i+vCols));
                distVFin += distance(positions.row(i+1), positions.row(i+vCols+1));
                distVFin /= 2.0;
                countV++;
            }
        }
    }
    totalUError /= (double)countU;
    totalUErrorLn /= (double)countU;
    totalUErrorCurv /= (double)countU;
    cout << "L0 Linear Error (percent): " << totalUError << endl;
    cout << "Ln Linear Error (percent): " << totalUErrorLn << endl;
}

bool GradDesc::isSMF() {
    return isSMFBool;
}

void GradDesc::shapeOptimizeGrid(int iters) {
    bool left, right, up, down, zup, zdown;
    double currRest;
    double zWeight = 40.0;
    int ptsPerLayer = positionsVec[0]->rows(), layer, row, col, upPt, downPt,
        leftPt, rightPt, zUpPt, zDownPt;
    SO::Matrix3X points = SO::Matrix3X::Zero(3, ptsPerLayer * layers);
    for (layer = 0; layer < layers; layer++)
        points.block(0, layer * ptsPerLayer, 3, ptsPerLayer) =
            positionsVec[layer]->transpose();
    shared_ptr<SO::Solver> s;
    s = make_shared<SO::Solver>();
    s->setPoints(points);
    vector<shared_ptr<SO::Constraint>> constraints;
    for (int i = 0; i < points.cols(); i++) {
        layer = i / ptsPerLayer;
        row = i % ptsPerLayer / vCols;
        col = i % ptsPerLayer % vCols;
        // Boolean checks for what neighbors we have
        left = col != 0;
        right = col != vCols - 1;
        up = row != 0;
        down = row != uRows - 1;
        zup = layer != 0;
        zdown = layer != layers - 1;
        //  Left
        if (left) {
            leftPt = i - 1;
            currRest = (*rest_uVec[layer])(row, col - 1);
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(leftPt);
            auto constraint = make_shared<SO::EdgeStrainConstraint>(
                pts, kUV, points, 1.0, 1.0);
            constraint->setEdgeLength(currRest);
            constraints.push_back(constraint);
        }
        // Right
        if (right) {
            rightPt = i + 1;
            currRest = (*rest_uVec[layer])(row, col);
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(rightPt);
            auto constraint = make_shared<SO::EdgeStrainConstraint>(
                pts, kUV, points, 1.0, 1.0);
            constraint->setEdgeLength(currRest);
            constraints.push_back(constraint);
        }
        // Up
        if (up) {
            upPt = i - vCols;
            currRest = (*rest_vVec[layer])(row - 1, col);
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(upPt);
            auto constraint = make_shared<SO::EdgeStrainConstraint>(
                pts, kUV, points, 1.0, 1.0);
            constraint->setEdgeLength(currRest);
            constraints.push_back(constraint);
        }
        // Down
        if (down) {
            downPt = i + vCols;
            currRest = (*rest_vVec[layer])(row, col);
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(downPt);
            auto constraint = make_shared<SO::EdgeStrainConstraint>(
                pts, kUV, points, 1.0, 1.0);
            constraint->setEdgeLength(currRest);
            constraints.push_back(constraint);
        }
        // Z Up
        if (zup) {
            zUpPt = i - ptsPerLayer;
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(zUpPt);
            auto constraint = make_shared<SO::EdgeStrainConstraint>(
                pts, zWeight, points, 1.0, 1.0);
            constraint->setEdgeLength(zHeight);
            constraints.push_back(constraint);
        }
        // Z Down
        if (zdown) {
            zDownPt = i + ptsPerLayer;
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(zDownPt);
            auto constraint = make_shared<SO::EdgeStrainConstraint>(
                pts, zWeight, points, 1.0, 1.0);
            constraint->setEdgeLength(zHeight);
            constraints.push_back(constraint);
        }

        // Add in all the angle values
        if (zup && left) {
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(zUpPt);
            pts.push_back(leftPt);
            auto constraint = make_shared<SO::AngleConstraint>(
                pts, kAng, points, H_PI - angTol, H_PI + angTol);
            constraints.push_back(constraint);
        }
        if (zup && right) {
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(zUpPt);
            pts.push_back(rightPt);
            auto constraint = make_shared<SO::AngleConstraint>(
                pts, kAng, points, H_PI - angTol, H_PI + angTol);
            constraints.push_back(constraint);
        }
        if (zup && up) {
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(zUpPt);
            pts.push_back(upPt);
            auto constraint = make_shared<SO::AngleConstraint>(
                pts, kAng, points, H_PI - angTol, H_PI + angTol);
            constraints.push_back(constraint);
        }
        if (zup && down) {
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(zUpPt);
            pts.push_back(downPt);
            auto constraint = make_shared<SO::AngleConstraint>(
                pts, kAng, points, H_PI - angTol, H_PI + angTol);
            constraints.push_back(constraint);
        }
        if (zdown && left) {
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(zDownPt);
            pts.push_back(leftPt);
            auto constraint = make_shared<SO::AngleConstraint>(
                pts, kAng, points, H_PI - angTol, H_PI + angTol);
            constraints.push_back(constraint);
        }
        if (zdown && right) {
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(zDownPt);
            pts.push_back(rightPt);
            auto constraint = make_shared<SO::AngleConstraint>(
                pts, kAng, points, H_PI - angTol, H_PI + angTol);
            constraints.push_back(constraint);
        }
        if (zdown && up) {
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(zDownPt);
            pts.push_back(upPt);
            auto constraint = make_shared<SO::AngleConstraint>(
                pts, kAng, points, H_PI - angTol, H_PI + angTol);
            constraints.push_back(constraint);
        }
        if (zdown && down) {
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(zDownPt);
            pts.push_back(downPt);
            auto constraint = make_shared<SO::AngleConstraint>(
                pts, kAng, points, H_PI - angTol, H_PI + angTol);
            constraints.push_back(constraint);
        }
        if (up && left) {
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(upPt);
            pts.push_back(leftPt);
            auto constraint = make_shared<SO::AngleConstraint>(
                pts, kAng, points, H_PI - angTol, H_PI + angTol);
            constraints.push_back(constraint);
        }
        if (down && left) {
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(downPt);
            pts.push_back(leftPt);
            auto constraint = make_shared<SO::AngleConstraint>(
                pts, kAng, points, H_PI - angTol, H_PI + angTol);
            constraints.push_back(constraint);
        }
        if (up && right) {
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(upPt);
            pts.push_back(rightPt);
            auto constraint = make_shared<SO::AngleConstraint>(
                pts, kAng, points, H_PI - angTol, H_PI + angTol);
            constraints.push_back(constraint);
        }
        if (down && right) {
            vector<int> pts;
            pts.push_back(i);
            pts.push_back(downPt);
            pts.push_back(rightPt);
            auto constraint = make_shared<SO::AngleConstraint>(
                pts, kAng, points, H_PI - angTol, H_PI + angTol);
            constraints.push_back(constraint);
        }
    }
    for (auto c : constraints) s->addConstraint(c);
    shared_ptr<SO::ClosenessConstraint> c1;
    vector<int> p1;
    p1.push_back(0);
    c1 = make_shared<SO::ClosenessConstraint>(p1, 50.0, points);
    s->addConstraint(c1);

    s->initialize();
    s->solve(iters);
    for (layer = 0; layer < layers; layer++) {
        (*positionsVec[layer])
            << s->getPoints()
                   .block(0, layer * ptsPerLayer, 3, ptsPerLayer)
                   .transpose();
    }
}
