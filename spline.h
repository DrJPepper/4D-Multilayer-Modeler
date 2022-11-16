#pragma once

#include "./bezier.h"
#include "./util.h"

// Builds on bezier.h/cpp and extends those functions to the splines generated
// for SMF models
MatrixXd genSpline(MatrixXd &points, MatrixXd &tangents, int &shrinkage);
MatrixXd initialize(MatrixXd &, MatrixXd *, float, int &shrinkage);
void printIvFile(MatrixXd *, MatrixXd *);
MatrixXd Q(MatrixXd *, float u);
MatrixXd tangent(MatrixXd *, float);
MatrixXd genCurveMat(MatrixXd *, int);
MatrixXd hermiteToBezier(MatrixXd *, MatrixXd *, int);
void smoothPoints(MatrixXd &points, int index);
