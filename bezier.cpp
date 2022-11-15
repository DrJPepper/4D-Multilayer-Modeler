#include "./bezier.h"

// Populates matrix of points
void initialize(MatrixXd *points, string inFileName) {
    ifstream inFile(inFileName);
    // Check the file
    if (!inFile.good()) {
        cout << "ERROR: input file invalid or doesn't exist" << endl;
        exit(1);
    }

    float x, y, z;

    for (string line; getline(inFile, line);) {
        points->conservativeResize(points->rows() + 1, points->cols());
        istringstream in(line);
        in >> x >> y >> z;
        points->row(points->rows() - 1) << x, y, z;
    }
}

// Runs Q equation for a whole Bezier curve
MatrixXd QLoop(MatrixXd &points, int count) {
    MatrixXd result = MatrixXd::Zero(count, 3);
    double step = 1.0 / static_cast<double>(count), u;
    for (int w = 0; w < count; w++) {
        MatrixXd acc = MatrixXd::Zero(1, 3);
        u = step * static_cast<double>(w);
        for (int i = 0; i < 4; i++) {
            acc += points.row(i) * kChoosei(3, i) * pow(1 - u, 3 - i) * pow(u, i);
        }
        result.row(w) = acc;
    }
    //cout << result << endl;
    //exit(0);
    return result;
}

// Runs Q equation
MatrixXd Q(MatrixXd *points, float u) {
    MatrixXd acc = MatrixXd::Zero(1, 3);
    for (int i = 0; i < 4; i++) {
        acc += points->row(i) * kChoosei(3, i) * pow(1 - u, 3 - i) * pow(u, i);
    }
    return acc;
}

TEST_CASE("testing the Q function") {
    MatrixXd points(4, 3);
    MatrixXd result(1, 3);
    result << 6, 1.5, 0;
    points << 0, 1.5, 0, 2, 1.5, 0, 4, 1.5, 0, 6, 1.5, 0;
    CHECK(Q(&points, 1.0) == result);
}

// Runs p equation
MatrixXd p(MatrixXd *points, float u, float v) {
    MatrixXd acc = MatrixXd::Zero(1, 3);
    for (int i = 0; i <= 3; i++) {
        for (int j = 0; j <= 3; j++) {
            acc += B(u, i) * B(v, j) * points->row(i + j * 4);
        }
    }
    return acc;
}

// Runs B equation
float B(float u, int i) {
    // n is always 3
    return kChoosei(3, i) * pow(u, i) * pow(1 - u, 3 - i);
}

int factorial(int x) {
    if (x == 0) return 1;
    return x * factorial(x - 1);
}

TEST_CASE("testing the factorial function") {
    CHECK(factorial(0) == 1);
    CHECK(factorial(1) == 1);
    CHECK(factorial(2) == 2);
    CHECK(factorial(3) == 6);
    CHECK(factorial(10) == 3628800);
}

int kChoosei(int k, int i) {
    return factorial(k) / (factorial(i) * factorial(k - i));
}

MatrixXd genPatchMat(MatrixXd *points, int u, int v) {
    MatrixXd patch = MatrixXd::Zero(u * v, 3);
    float du, dv;
    int count = 0;
    for (int i = 0; i < v; i++) {
        for (int j = 0; j < u; j++) {
            du = j / static_cast<float>(u - 1);
            dv = i / static_cast<float>(v - 1);
            patch.row(count) = p(points, du, dv);
            count++;
        }
    }
    return patch;
}

// Populates a matrix with the normal values at each point
MatrixXd genNormMat(MatrixXd *points, int u, int v) {
    MatrixXd normals = MatrixXd::Zero(u * v, 3);
    float du, dv;
    int count = 0;
    for (int i = 0; i < v; i++) {
        for (int j = 0; j < u; j++) {
            du = 1.0 / static_cast<float>(u - 1);
            dv = 1.0 / static_cast<float>(v - 1);
            normals.row(count) = normal(points, j * du, i * dv).normalized();
            count++;
        }
    }
    return normals;
}

// Calculates the partial derivative with respect to u at a point
Vector3d dQdu(MatrixXd *points, float u, float v) {
    MatrixXd curveV = MatrixXd::Zero(4, 3);
    MatrixXd temp = MatrixXd::Zero(4, 3);
    for (int i = 0; i < 4; i++) {
        temp.row(0) = points->row(i);
        temp.row(1) = points->row(i + 4);
        temp.row(2) = points->row(i + 8);
        temp.row(3) = points->row(i + 12);
        curveV.row(i) = Q(&temp, v);
    }
    return -3 * (1 - u) * (1 - u) * (Vector3d)curveV.row(0) +
           (3 * (1 - u) * (1 - u) - 6 * u * (1 - u)) * (Vector3d)curveV.row(1) +
           (6 * u * (1 - u) - 3 * u * u) * (Vector3d)curveV.row(2) +
           3 * u * u * (Vector3d)curveV.row(3);
}

TEST_CASE("testing the dQdu function") {
    MatrixXd points(16, 3);
    Vector3d result;
    result << 6, 0, 0;
    points << 0, 0, 0, 2, 0, 0, 4, 0, 0, 6, 0, 0, 0, 2, 0, 2, 2, 0, 4, 2, 0, 6,
        2, 0, 0, 4, 0, 2, 4, 0, 4, 4, 0, 6, 4, 0, 0, 6, 0, 2, 6, 0, 4, 6, 0, 6,
        6, 0;
    // Question mark on these results?
    CHECK(dQdu(&points, 0.75, 0.25).transpose() == result.transpose());
}

// Calculates the partial derivative with respect to v at a point
Vector3d dQdv(MatrixXd *points, float u, float v) {
    MatrixXd curveU = MatrixXd::Zero(4, 3);
    for (int i = 0; i < 4; i++) {
        MatrixXd temp = points->block(i * 4, 0, 4, 3);
        curveU.row(i) = Q(&temp, u);
    }
    return -3 * (1 - v) * (1 - v) * (Vector3d)curveU.row(0) +
           (3 * (1 - v) * (1 - v) - 6 * v * (1 - v)) * (Vector3d)curveU.row(1) +
           (6 * v * (1 - v) - 3 * v * v) * (Vector3d)curveU.row(2) +
           3 * v * v * (Vector3d)curveU.row(3);
}

// Calculates the normal vector at a point
MatrixXd normal(MatrixXd *points, float u, float v) {
    Vector3d tanU = dQdu(points, u, v);
    Vector3d tanV = dQdv(points, u, v);
    return tanV.transpose().cross(tanU);
}

// Runs Q' equation
Vector3d qPrime(MatrixXd *points, float u) {
    Vector3d acc;
    acc << 0, 0, 0;
    for (int i = 0; i < 3; i++)
        acc += kChoosei(2, i) * pow(u, i) * pow(1 - u, 2 - i) *
               (3 * (Vector3d)(points->row(i + 1) - points->row(i)));
    return acc;
}

TEST_CASE("testing the qPrime function") {
    MatrixXd points(4, 3);
    Vector3d result;
    result << 6, 0, 0;
    points << 0, 1.5, 0, 2, 1.5, 0, 4, 1.5, 0, 6, 1.5, 0;
    // I guess it makes sense this is the same as the dQdu test case
    CHECK(qPrime(&points, 0.75).transpose() == result.transpose());
}

// Runs Q'' equation
Vector3d qDoublePrime(MatrixXd *points, float u) {
    Vector3d acc;
    acc << 0, 0, 0;
    for (int i = 0; i < 2; i++)
        acc += pow(u, i) * pow(1 - u, 1 - i) * 6 *
               (Vector3d)(points->row(i + 2) - 2 * points->row(i + 1) +
                          points->row(i));
    return acc;
}

TEST_CASE("testing the qDoublePrime function") {
    MatrixXd points(4, 3);
    Vector3d result;
    result << 0, 0, 0;
    points << 0, 1.5, 0, 2, 1.5, 0, 4, 1.5, 0, 6, 1.5, 0;
    CHECK(qDoublePrime(&points, 0.75).transpose() == result.transpose());
}

double curvature(Vector3d prime, Vector3d doublePrime) {
    return (prime.cross(doublePrime)).norm() / pow(prime.norm(), 3);
}

double curvatureFromUSigned(MatrixXd *points, float u) {
    Vector3d prime = qPrime(points, u);
    Vector3d doublePrime = qDoublePrime(points, u);
    return ((doublePrime(0) > 0) - (doublePrime(0) < 0)) * curvature(prime, doublePrime);
}

double curvatureFromU(MatrixXd *points, float u) {
    Vector3d prime = qPrime(points, u);
    Vector3d doublePrime = qDoublePrime(points, u);
    return curvature(prime, doublePrime);
}

MatrixXd hermiteToBezier(MatrixXd *points, MatrixXd *tangents) {
    MatrixXd result = MatrixXd::Zero(4, 3);
    result.row(0) = points->row(0);
    result.row(1) = points->row(0) + (tangents->row(0) / 3.0);
    result.row(2) = points->row(1) - (tangents->row(1) / 3.0);
    result.row(3) = points->row(1);
    return result;
}

double distanceInV(MatrixXd *points, float u, float vstart, float vend) {
    MatrixXd one = p(points, u, vstart);
    MatrixXd two = p(points, u, vend);
    return distance(one, two);
}

double distanceInU(MatrixXd *points, float ustart, float uend, float v) {
    MatrixXd one = p(points, ustart, v);
    MatrixXd two = p(points, uend, v);
    return distance(one, two);
}

TEST_CASE("testing the distanceInU function") {
    MatrixXd points(16, 3);
    points << 0, 0, 0, 2, 0, 0, 4, 0, 0, 6, 0, 0, 0, 2, 0, 2, 2, 0, 4, 2, 0, 6,
        2, 0, 0, 4, 0, 2, 4, 0, 4, 4, 0, 6, 4, 0, 0, 6, 0, 2, 6, 0, 4, 6, 0, 6,
        6, 0;
    CHECK(distanceInU(&points, 0.0, 0.5, 0.5) == 3.0);
}

/*
 * points: 16 Bezier patch control points
 * u, v: u and v values
 * inU: a flag to signal which direction we want it run in
 */
double curvatureSign(MatrixXd *points, float u, float v, bool inU) {
    double result;
    Vector3d N, qDP;
    MatrixXd bezPtsU, bezPtsV;
    bool check = u < 0.1 && saveToJS;
    auto ps = p(points, u, v);
    if (check) {
        auto entity = makeEntity("point", 1.0, 1.0, 1.0, fmt::format("Point at {}, {}", u, v), ps);
        jsOut["list"][0]["entities"].push_back(entity);
    }
    // Generate 4 Bezier curve points for V
    MatrixXd hermitePts = MatrixXd::Zero(2, 3);
    MatrixXd hermiteTans = MatrixXd::Zero(2, 3);
    hermitePts.row(0) = p(points, u, 0.0);
    hermitePts.row(1) = p(points, u, 1.0);
    hermiteTans.row(0) = dQdv(points, u, 0.0);
    hermiteTans.row(1) = dQdv(points, u, 1.0);
    bezPtsV = hermiteToBezier(&hermitePts, &hermiteTans);
    // Generate 4 Bezier curve points for U
    hermitePts = MatrixXd::Zero(2, 3);
    hermiteTans = MatrixXd::Zero(2, 3);
    hermitePts.row(0) = p(points, 0.0, v); 
    hermitePts.row(1) = p(points, 1.0, v); 
    hermiteTans.row(0) = dQdu(points, 0.0, v);
    hermiteTans.row(1) = dQdu(points, 1.0, v);
    bezPtsU = hermiteToBezier(&hermitePts, &hermiteTans);
    // Calculate N and Q (which should be equivalent to P)
    // double prime in the requested direction
    N = qPrime(&bezPtsU, u).cross(qPrime(&bezPtsV, v));
    if (inU) {
        qDP = qDoublePrime(&bezPtsU, u);
    } else {
        qDP = qDoublePrime(&bezPtsV, v);
    }
    // Get the final dot product result
    result = N.dot(qDP);
    if (check) {
        double index = 0;
        MatrixXd psVec(2, 3);
        psVec.row(0) = ps;
        psVec.row(1) = N;
        jsOut["list"][index]["entities"].push_back(makeEntity("vector", 1.0, 0.0, 0.0, "N vector", psVec));
        psVec.row(1) = qDP;
        jsOut["list"][index]["entities"].push_back(makeEntity("vector", 0.0, 1.0, 1.0, "Q double prime vector", psVec));
        psVec.row(1) = qPrime(&bezPtsU, u);
        jsOut["list"][index]["entities"].push_back(makeEntity("vector", 0.0, 1.0, 0.0, "Q Prime in u", psVec));
        psVec.row(1) = qPrime(&bezPtsV, v);
        jsOut["list"][index]["entities"].push_back(makeEntity("vector", 0.0, 1.0, 0.0, "Q Prime in v", psVec));
        if (!jsOut["list"][index].contains("description"))
            jsOut["list"][index]["description"] = "N dot qDP results:";
        jsOut["list"][index]["description"] = fmt::format("{}\n  at {},{}: {}", jsOut["list"][index]["description"], u, v, result);
        ofstream o("out.json");
        o << jsOut << endl;
        o.close();
    }
    // Return the sign of the result as a floating point
    return static_cast<double>(-1.0 * ((result > 0) - (result < 0)));
}

double curvatureInV(MatrixXd *points, float u, float v) {
    double curvature, qDP;
    MatrixXd bezPts;
    MatrixXd hermitePts = MatrixXd::Zero(2, 3);
    MatrixXd hermiteTans = MatrixXd::Zero(2, 3);
    hermitePts.row(0) = p(points, u, 0.0);
    hermitePts.row(1) = p(points, u, 1.0);
    hermiteTans.row(0) = dQdv(points, u, 0.0);
    hermiteTans.row(1) = dQdv(points, u, 1.0);
    bezPts = hermiteToBezier(&hermitePts, &hermiteTans);
    curvature = curvatureFromU(&bezPts, v);
    qDP = qDoublePrime(&bezPts, v)(2);
    return curvatureSign(points, u, v, false) * curvature;
}

double curvatureInU(MatrixXd *points, float u, float v) {
    double curvature, qDP;
    MatrixXd bezPts;
    MatrixXd hermitePts = MatrixXd::Zero(2, 3);
    MatrixXd hermiteTans = MatrixXd::Zero(2, 3);
    hermitePts.row(0) = p(points, 0.0, v);
    hermitePts.row(1) = p(points, 1.0, v);
    hermiteTans.row(0) = dQdu(points, 0.0, v);
    hermiteTans.row(1) = dQdu(points, 1.0, v);
    bezPts = hermiteToBezier(&hermitePts, &hermiteTans);
    curvature = curvatureFromU(&bezPts, u);
    qDP = qDoublePrime(&bezPts, u)(2);
    return curvatureSign(points, u, v, true) * curvature;
}

TEST_CASE("testing the curvatureInU function") {
    MatrixXd points(16, 3);
    points << 0, 0, 0, 2, 0, 0, 4, 0, 0, 6, 0, 0, 0, 2, 0, 2, 2, 0, 4, 2, 0, 6,
        2, 0, 0, 4, 0, 2, 4, 0, 4, 4, 0, 6, 4, 0, 0, 6, 0, 2, 6, 0, 4, 6, 0, 6,
        6, 0;
    CHECK(curvatureInU(&points, 0.5, 0.5) == 0.0);
}
