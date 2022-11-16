#include "util.h"

bool saveToJS = false;

// Converts a ratio (in the case of this project, a spring's rest length
// vs its current length) to RGB values moving from red to green to blue
Vector3d ratioToRGB(double ratio) {
    double cap = 0.2, green, other;
    double minCap = 1.0 - cap;
    double maxCap = 1.0 + cap;
    Vector3d result = Vector3d::Zero(3);

    if (ratio > maxCap)
        result << 255.0, 0.0, 0.0;
    else if (ratio < minCap)
        result << 0.0, 0.0, 255.0;
    else if (ratio > 1.0) {
        other = (ratio - 1.0) / cap * 255.0;
        green = 255 - other;
        result << other, green, 0.0;
    } else if (ratio <= 1.0) {
        other = (1.0 - ratio) / cap * 255.0;
        green = 255 - other;
        result << 0.0, green, other;
    } else {
        cout << "NaN in ratioToRGB" << endl;
    }

    return result / 255.0;
}

double distance(MatrixXd one, MatrixXd two) {
    return sqrt(pow2(two(0) - one(0)) + pow2(two(1) - one(1)) +
                pow2(two(2) - one(2)));
}

double vectAngle(Vector3d one, Vector3d two) {
    return acos((one.dot(two)) / (one.norm() * two.norm()));
}

double pow2(double x) {
    return x * x;
}

// A generic function that returns a JSON object in a form compatible with my
// primitives visualizer program. This is used for debugging only.
json makeEntity(string type, float red, float green, float blue, string description, MatrixXd& position) {
    json entity = {
        {"type", type},
        {"color", {red, green, blue}},
        {"description", description},
    };
    if (!type.compare("point"))
        entity["position"] = {position(0), position(1), position(2)};
    else if (!type.compare("vector") || !type.compare("line"))
        entity["position"] = {position(0, 0), position(0, 1), position(0, 2), position(1, 0), position(1, 1), position(1, 2)};
    else
        entity = {"Errored", fmt::format("Invalid type {}", type)};

    return entity;
}
