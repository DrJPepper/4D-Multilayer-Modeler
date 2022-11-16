#include "spline.h"

// Takes the points that need intersected and returns the stitched spline
MatrixXd genSpline(MatrixXd &points, MatrixXd &tangents, int &shrinkage) {
    int n = 10;
    float tension = 0.0;

    // Run primary functions
    points = initialize(points, &tangents, tension, shrinkage);
    MatrixXd curve = MatrixXd::Zero((points.rows() - 1) * (n + 1), 3);
    for (int i = 0; i < points.rows() - 1; i++) {
        MatrixXd temp = hermiteToBezier(&points, &tangents, i);
        curve.block(i * (n + 1), 0, 1 + n, 3) = genCurveMat(&temp, n);
    }
    return points;
}

// Blends the controls points together from one curve to the next
void smoothPoints(MatrixXd &points, int index) {
    Vector3d zero(0, 0, 0);
    Vector3d p0, p1 = zero, p2 = zero;
    double t, d;
    int counter = index - 1;
    if (index > 1 && index < points.rows()) {
        p0 = points.row(index);
        if (p0 != zero) {
            while (p1 == zero && counter >= 0 && (index - counter) < 7) {
                p1 = points.row(counter);
                counter--;
            }
            while (p2 == zero && counter >= 0 && (index - counter) < 7) {
                p2 = points.row(counter);
                counter--;
            }
            if (p1 != zero && p2 != zero) {
                t = -(p1 - p0).dot(p2 - p1) / pow2((p2 - p1).norm());
                d = sqrt(pow2((p1(0) - p0(0)) + (p2(0) - p1(0)) * t) +
                         pow2((p1(1) - p0(1)) + (p2(1) - p1(1)) * t) +
                         pow2((p1(2) - p0(2)) + (p2(2) - p1(2)) * t));
                if (d >= 1.0) {
                    points.row(index) = zero;
                }
            }
        }
    }
}

// Populates matrix of points
MatrixXd initialize(MatrixXd &points, MatrixXd *tangents, float tension,
                    int &shrinkage) {
    // Store the last tangent in a temp file and first into our real matrix
    MatrixXd zero = MatrixXd::Zero(1, 3);
    shrinkage = 0;
    Vector3d vec, startRow, endRow;
    MatrixXd pointsNew = points;
    int counter = 0, start, end, div;
    bool tempB;
    while (counter < pointsNew.rows() - 1) {
        smoothPoints(pointsNew, counter);
        tempB = pointsNew.row(counter) == zero;
        if (tempB) {
            start = counter;
            while (tempB && counter < pointsNew.rows() - 1) {
                counter++;
                smoothPoints(pointsNew, counter);
                tempB = pointsNew.row(counter) == zero;
            }
            if (counter != pointsNew.rows() - 1) counter--;
            end = counter;
            if (start == 0 && end == pointsNew.rows() - 1) {
                shrinkage = points.rows();
                return zero;
            } else if (start == 0) {
                shrinkage += end + 1;
                counter = -1;
                MatrixXd tempTemp =
                    pointsNew(seq(end + 1, pointsNew.rows() - 1), all);
                pointsNew = tempTemp;
            } else if (end == (pointsNew.rows() - 1)) {
                shrinkage += end - start + 1;
                MatrixXd tempTemp2 = pointsNew(seq(0, start - 1), all);
                pointsNew = tempTemp2;
                counter = pointsNew.rows();
            } else {
                startRow = pointsNew.row(start - 1);
                endRow = pointsNew.row(end + 1);
                vec = endRow - startRow;
                div = end - start + 2;
                for (int q = 1; q <= end - start + 1; q++) {
                    pointsNew.row(start + q - 1) = startRow + vec / div * q;
                }
                counter = end;
            }
        }
        counter++;
    }
    points = pointsNew;
    MatrixXd lastTangent = MatrixXd::Zero(1, 3);
    tangents->conservativeResize(tangents->rows() + 1, tangents->cols());
    tangents->row(tangents->rows() - 1) = points.row(1) - points.row(0);
    tangents->row(tangents->rows() - 1) *= (1 - tension);
    lastTangent = points.row(points.rows() - 1) - points.row(points.rows() - 2);
    lastTangent *= (1 - tension);
    points = pointsNew;

    // TODO: I may be generating an extra point and tan here, but it's working
    // correctly with graddesc so I'm not going to mess with it further for now
    MatrixXd temp = MatrixXd::Zero(2, 3);
    // Calculate tangents
    for (int i = 1; i < points.rows() - 1; i++) {
        tangents->conservativeResize(tangents->rows() + 1, tangents->cols());
        temp.row(0) = points.row(i - 1);
        temp.row(1) = points.row(i + 1);
        tangents->row(tangents->rows() - 1) = tangent(&temp, tension);
    }

    // Move over the final tangent
    tangents->conservativeResize(tangents->rows() + 1, tangents->cols());
    tangents->row(tangents->rows() - 1) = lastTangent;
    return pointsNew;
}

MatrixXd tangent(MatrixXd *inTans, float tension) {
    return (inTans->row(1) - inTans->row(0)) / 2.0 * (1 - tension);
}

// Prints things in open inventor format
void printIvFile(MatrixXd *points, MatrixXd *curve) {
    ofstream outFile;
    outFile.open("./out.iv");
    outFile << fixed;
    outFile << setprecision(6);
    outFile << "#Inventor V2.0 ascii\n\n";

    outFile << "Separator {LightModel {model BASE_COLOR} Material {diffuseColor "
            "1.0 "
            "1.0 1.0}\n";
    outFile << "\tCoordinate3 { 	point [\n";
    for (int i = 0; i < curve->rows(); i++)
        outFile << "\t\t" << curve->row(i) << ",\n";
    outFile << "\t] }\n";
    outFile << "\tIndexedLineSet {coordIndex [\n\t";
    for (int i = 0; i < curve->rows(); i++) outFile << i << ", ";
    outFile << "-1," << endl;
    outFile << "\t] } }\n";
}

MatrixXd hermiteToBezier(MatrixXd *points, MatrixXd *tangents, int row1) {
    MatrixXd result = MatrixXd::Zero(4, 3);
    result.row(0) = points->row(row1);
    result.row(1) = points->row(row1) + (tangents->row(row1) / 3.0);
    result.row(2) = points->row(row1 + 1) - (tangents->row(row1 + 1) / 3.0);
    result.row(3) = points->row(row1 + 1);
    return result;
}

double tempCounter = 0;

// Populates a matrix with values from the Q calculations
MatrixXd genCurveMat(MatrixXd *points, int n) {
    MatrixXd curve = MatrixXd::Zero(n + 1, 3);
    float du = 1.0 / (float)n;
    for (int i = 0; i <= n; i++) {
        curve.row(i) = Q(points, i * du);
        curve(i, 2) = curvatureFromU(points, i * du);
        double qDP = qDoublePrime(points, i * du)(0);
        curve(i, 2) = ((qDP > 0) - (qDP < 0)) * curve(i, 2);
        if (i < n) {
            tempCounter += distance(Q(points, i * du), Q(points, (i + 1) * du));
        }
    }
    return curve;
}
