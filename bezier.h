#pragma once

#include "./util.h"

void initialize(MatrixXd*, string);
MatrixXd Q(MatrixXd *points, float u);
MatrixXd QLoop(MatrixXd &points, int count);
MatrixXd p(MatrixXd*, float, float);
float B(float, int);
int factorial(int);
int kChoosei(int n, int k);
MatrixXd genPatchMat(MatrixXd *, int, int);
Vector3d dQdu(MatrixXd *, float, float);
Vector3d dQdv(MatrixXd *, float, float);
MatrixXd genNormMat(MatrixXd *, int, int);
MatrixXd normal(MatrixXd *, float, float);

Vector3d qPrime(MatrixXd *, float);
Vector3d qDoublePrime(MatrixXd *, float);
double curvature(Vector3f, Vector3f);
double curvatureFromU(MatrixXd *, float);
double curvatureFromUSigned(MatrixXd *, float);
double curvatureInU(MatrixXd *points, float u, float v);
double curvatureInV(MatrixXd *points, float u, float v);
double distanceInU(MatrixXd *points, float ustart, float uend, float v);
double distanceInV(MatrixXd *points, float u, float vstart, float vend);
double curvatureSign(MatrixXd *points, float u);
MatrixXd hermiteToBezier(MatrixXd *points, MatrixXd *tangents);
