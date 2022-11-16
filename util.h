#pragma once

// This file handles all imports and holds a few very generic functions

#include <QtCore>
#include <QtCore/QtGlobal>
#include <QApplication>
#include <QWidget>
#include <QVTKOpenGLNativeWidget.h>

#include <float.h>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <omp.h>
#include <qsurfaceformat.h>
#include <vtkActor.h>
#include <vtkTriangle.h>
#include <vtkActor2D.h>
#include <vtkActorCollection.h>
#include <vtkBoxWidget.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkColor.h>
#include <vtkConeSource.h>
#include <vtkCylinderSource.h>
#include <vtkDataSetMapper.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkLight.h>
#include <vtkLineSource.h>
#include <vtkNamedColors.h>
#include <vtkOBJExporter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkTubeFilter.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>

#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkPropPicker.h>
#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>
#include <vtkRenderWindowInteractor.h>

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <cstddef>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>

#include "doctest.h"
#include "Constraint.h"
#include "Solver.h"
#include "./lib/gnuplot-iostream.h"
#include "lib/json.hpp"

#define PI 3.14159265
#define H_PI 1.57079632
#define MatrixXd Matrix<double,Dynamic,Dynamic,RowMajor>
#define MatrixXi Matrix<int,Dynamic,Dynamic,RowMajor>

using namespace std;
using namespace Eigen;
using json = nlohmann::json;
namespace SO = ShapeOp;

using Eigen::all;

Vector3d ratioToRGB(double ratio);
double distance(MatrixXd one, MatrixXd two);
double vectAngle(Vector3d one, Vector3d two);
double pow2(double x);
json makeEntity(string type, float red, float green, float blue, string description, MatrixXd& position);

// The master JSON logging object. It is global so that it can be added to from any class.
static json jsOut;
extern bool saveToJS;
