#pragma once
#include "./bezier.h"
#include "./raytrace.h"
#include "./scene.h"
#include "./spline.h"
#include "./util.h"

class GradDesc {
    public:
        GradDesc();
        void populateMatrices(ifstream *inFile);
        void populateSheetOfNodesBezierMS();
        void populateSheetOfNodesGeneric();
        void duplicateL0();
        void populateSheetOfNodesAuto();
        void populateSMFNodes();
        void analyzeSurfaceGeneric();
        void analyzeSurfaceBezier();
        void updateVtkGrid(Scene &scene);
        void assessResults();
        double getSysEnergy();
        Array2d getAvgSpringDisp();
        int getRows();
        int getCols();
        double getSpacing();
        int getLayers();
        void shapeOptimizeGrid();
        void shapeOptimizeGrid(int iters);
        MatrixXd getVertices();
        MatrixXd getMiddleLayer();
        MatrixXd getPVLayer(int);
        MatrixXi getFaces();
        double getRestU(int row, int col, int layer);
        double getRestV(int row, int col, int layer);
        MatrixXd intersectPoints;
        vector<MatrixXd*> curvePointsU, curvePointsV;
        deque<void (GradDesc::*)()> funDeq;
        deque<string> stringDeq;
        bool isSMF();

    private:
        void plotCurvHisto(vector<double> &curvs);
        void initVectors();
        void initNodeMatricesSMF();
        bool isSMFBool;
        int uRows, node_u, vCols, node_v, node_count_u, node_count_v,
            node_count_long, layers, analysisSteps, analysisLines;
        double kAng, kUV, kZ, zHeight, radius, spacing, maxShrinkage,
            nodeGranularity, maxZHeight, aspectRatio, k_in, maxKappa, angTol,
            restAng, uLength, vLength, uStart, vStart, pValue;
        string inputFileName;
        Array2d curvatureMax, distanceMax, distanceMin, minNodes;
        MatrixXd points, curvatureU, curvatureV, distanceU, distanceV, vertices,
                  smfDistancesU, smfDistancesV, smfCurvaturesU, smfCurvaturesV;
        MatrixXi faces;
        vector<MatrixXd *> positionsVec, rest_uVec, rest_vVec;
        vector<double> curvaturesVec;
};
