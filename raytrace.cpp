#define PI 3.14159265
// NOTE: Don't change this, multithreading not fully set up for these functions
#define NUM_THREADS 1

#include "./raytrace.h"

std::tuple<int, int, MatrixXd&> rayTrace(string inputFileName, MatrixXd &vertices, MatrixXi &faces, double *rayBox, int xPoints, int yPoints) {
    Settings *setts;

    ifstream inFile(inputFileName);
    if (!inFile.good()) {
        cout << "ERROR: SMF file invalid or doesn't exist\n";
        exit(1);
    }
    setts = new Settings(&inFile, rayBox);
    if (xPoints && xPoints) {
        setts->x = xPoints;
        setts->y = yPoints;
    }

    // Run primary functions
    MatrixXd pixelMap = MatrixXd::Zero(setts->x * setts->y, 3);
    intersectForAllPixels(objects, setts, &pixelMap);
    vertices = reinterpret_cast<SMFModel*>(objects[0]->obj)->vertices;
    faces = reinterpret_cast<SMFModel*>(objects[0]->obj)->facesStored;
    MatrixXd *intPoints = new MatrixXd;
    (*intPoints) = tBuffer.block(0, 1, tBuffer.rows(), 3);
    return {setts->x, setts->y, std::ref(*intPoints)};
}

// Loops over all pixels
void intersectForAllPixels(vector<ACGBVH *> objects, Settings *setts,
                           MatrixXd *result) {
    Vector3d djk = setts->Zv;
    int j, k, index;
    int counter = 0;
    int size = static_cast<int>(objects.size());
    double uStep, vStep;
    if (!setts->uLength || !setts->vLength) {
        uStep = objects[0]->obj->maxVec(1) / setts->x;
        vStep = objects[0]->obj->maxVec(2) / setts->y;
    } else {
        uStep = setts->uLength / setts->x;
        vStep = setts->vLength / setts->y;
    }

    // t, x, y, z (of intersection), x, y, z (of normal)
    tBuffer = MatrixXd::Zero(result->rows(), 7);
    // Stores the front object for each pixel since tBuffer is a double matrix
    // and it seems like a bad idea to store indices as doubles
    tObjsBuffer = VectorXi::Zero(result->rows());
    for (int i = 0; i < tBuffer.rows(); i++) tObjsBuffer(i) = -1;
    for (int i = 0; i < size; i++) {
#pragma omp parallel for num_threads(NUM_THREADS) private(j, k, index)
        for (j = 0; j < setts->x; j++) {
            for (k = 0; k < setts->y; k++) {
                index = k * setts->x + j;
                setts->cameraLoc << -50.0, setts->uStart + uStep * j, setts->vStart + vStep * k;
                objects[i]->intersectPixel(djk.normalized(), setts, index, i);
            }
#pragma omp critical
            {
                // Keep a running counter for every row processed/to process
                counter++;
            }
        }
    }
}

// Performs the Phong shading calculations
void shadePixels(vector<ACGBVH *> objects, Settings *setts, MatrixXd *result) {
    Intersectable *currObj;
    Vector3d r, n, intersect, L, IL, v;
    double nDotL, rDotV;
    for (int i = 0; i < result->rows(); i++) {
        if (tObjsBuffer(i) != -1) {
            for (int j = 0; j < lights.rows(); j++) {
                currObj = objects[tObjsBuffer(i)]->obj;
                n = tBuffer.block(i, 4, 1, 3).transpose();
                intersect = tBuffer.block(i, 1, 1, 3).transpose();
                L = (lights.block(j, 0, 1, 3).transpose() - intersect)
                        .normalized();
                v = (setts->cameraLoc - intersect).normalized();
                IL = lights.block(j, 3, 1, 3).transpose();
                nDotL = max(0.0, static_cast<double>(n.dot(L)));
                r = ((2 * nDotL * n) - L).normalized();
                rDotV = max(0.0, static_cast<double>(r.dot(v)));
                result->row(i) +=
                    (IL * currObj->kd * nDotL).cwiseProduct(currObj->color);
                result->row(i) += (IL * currObj->ks * pow(rDotV, currObj->n))
                                      .cwiseProduct(currObj->color);
                result->row(i) +=
                    currObj->color.cwiseProduct(setts->ambient) * currObj->ka;
            }
        }
    }
}

// Parse the input file and generate resulting settings and objects
Settings::Settings(ifstream *inFile, double *rayBox) {
    const regex thetaR("^theta:\\s*(-?[0-9]+\\.[0-9]+)\\s*$");
    const regex dR("^d:\\s*(-?[0-9]+\\.[0-9]+)\\s*$");
    const regex cameraLocR(
        "^camera_loc:\\s*(-?[0-9]+\\.[0-9]+),\\s*(-?[0-9]+\\.["
        "0-9]+),\\s*(-?[0-9]+\\.[0-9]+)\\s*$");
    const regex ZvR(
        "^Zv:\\s*(-?[0-9]+\\.[0-9]+),\\s*(-?[0-9]+\\.[0-9]+),\\s*(-?["
        "0-9]+\\.[0-9]+)\\s*$");
    const regex VupR(
        "^Vup:\\s*(-?[0-9]+\\.[0-9]+),\\s*(-?[0-9]+\\.[0-9]+),\\s*(-"
        "?[0-9]+\\.[0-9]+)\\s*$");
    const regex ambientR(
        "^ambient:\\s*(-?[0-9]+\\.[0-9]+),\\s*(-?[0-9]+\\.[0-9]+"
        "),\\s*(-?[0-9]+\\.[0-9]+)\\s*$");
    const regex lightR(
        "^light:\\s*(-?[0-9]+\\.[0-9]+),\\s*(-?[0-9]+\\.[0-9]+),"
        "\\s*(-?[0-9]+\\.[0-9]+),\\s*(-?[0-9]+\\.[0-9]+),\\s*(-?["
        "0-9]+\\.[0-9]+),\\s*(-?[0-9]+\\.[0-9]+)\\s*$");
    const regex xR("^image_x:\\s*([0-9]+)\\s*$");
    const regex yR("^image_y:\\s*([0-9]+)\\s*$");
    const regex meshR("^mesh:\\s*(.+)\\s*$");
    const regex circleR("^circle:\\s*(-?[0-9]+\\.[0-9]+)\\s*$");

    int maxObjCount = 200;

    string line;
    smatch sm;
    // Default values for some things if not specified in input file
    d = 3.0;
    x = 512;
    y = 512;
    uStart = rayBox[0];
    vStart = rayBox[1];
    uLength = rayBox[2];
    vLength = rayBox[3];
    cameraLoc << 3.0, 0.0, 0.0;
    Zv << -1.0, 0.0, 0.0;
    Vup << 0.0, 1.0, 0.0;
    lights = MatrixXd::Zero(0, 6);
    getline(*inFile, line);
    while (!inFile->eof()) {
        if (line.front() == '#' || !line.compare("")) {
        } else if (regex_search(line, sm, thetaR)) {
            theta = stod(sm[1]) * PI / 180.0;
        } else if (regex_search(line, sm, dR)) {
            d = stod(sm[1]);
        } else if (regex_search(line, sm, cameraLocR)) {
            cameraLoc = Vector3d(stod(sm[1]), stod(sm[2]), stod(sm[3]));
        } else if (regex_search(line, sm, ZvR)) {
            Zv = Vector3d(stod(sm[1]), stod(sm[2]), stod(sm[3]));
        } else if (regex_search(line, sm, ambientR)) {
            ambient = Vector3d(stod(sm[1]), stod(sm[2]), stod(sm[3]));
        } else if (regex_search(line, sm, VupR)) {
            Vup = Vector3d(stod(sm[1]), stod(sm[2]), stod(sm[3]));
        } else if (regex_search(line, sm, lightR)) {
            lights.conservativeResize(lights.rows() + 1, lights.cols());
            lights.row(lights.rows() - 1) << stod(sm[1]), stod(sm[2]),
                stod(sm[3]), stod(sm[4]), stod(sm[5]), stod(sm[6]);
        } else if (regex_search(line, sm, xR)) {
            x = stoi(sm[1]);
        } else if (regex_search(line, sm, yR)) {
            y = stoi(sm[1]);
            // meshR and circleR result in creation of objects
            // smf files are parsed in the SMFModel constructor
        } else if (regex_search(line, sm, meshR)) {
            ifstream currMesh(sm[1]);
            if (!currMesh.good()) {
                cout << "ERROR: input mesh file " << sm[1]
                     << " invalid or doesn't exist\n";
                exit(1);
            }
            Intersectable *temp = new SMFModel(&currMesh, inFile);
            ACGBVH *newObj = new ACGBVH(temp, maxObjCount);
            objects.push_back(newObj);
        } else if (regex_search(line, sm, circleR)) {
            Intersectable *temp = new Sphere(stod(sm[1]), inFile);
            ACGBVH *newObj = new ACGBVH(temp, maxObjCount);
            objects.push_back(newObj);
        } else {
            cout << "Ignoring input line: " << line << endl;
        }
        getline(*inFile, line);
    }
    Xv = Zv.cross(Vup);
    Yv = Xv.cross(Zv);
    Xv.normalize();
    Yv.normalize();
    Zv.normalize();
    h = d * tan(theta / 2.0);
    Sj = 2 * h;
    Sk = Sj * (static_cast<double>(y) / static_cast<double>(x));
    P00 = cameraLoc + d * Zv - (Sj / 2) * Xv + (Sk / 2) * Yv;
}

Sphere::Sphere(double radius, ifstream *inSetts) {
    rs = radius;

    type = SPHERE;

    initSetts(inSetts);

    Xc = translate(0);
    Yc = translate(1);
    Zc = translate(2);
    center = Vector3d(Xc, Yc, Zc);
}

// Both check for intersection and populate relevant data stores
// The boolean return value is not currently being checked against although I'm
// leaving it in case it becomes relevant again
bool Sphere::intersectPixel(Vector3d Rd, Settings *setts, int pixelNum,
                            int objNum) {
    Vector3d Ro = setts->cameraLoc;
    Vector3d tempVec, tempNorm;
    double A = pow(Rd(0), 2) + pow(Rd(1), 2) + pow(Rd(2), 2);
    double B = 2 * (Rd(0) * (Ro(0) - Xc) + Rd(1) * (Ro(1) - Yc) +
                    Rd(2) * (Ro(2) - Zc));
    double C = pow(Ro(0) - Xc, 2) + pow(Ro(1) - Yc, 2) + pow(Ro(2) - Zc, 2) -
               pow(rs, 2);
    double disc = pow(B, 2) - 4 * A * C;
    if (disc >= 0) {
        double t0 = (-B - sqrt(disc)) / 2.0;
        double t1 = (-B + sqrt(disc)) / 2.0;
        double t;
        if (t1 >= 0.0 && t1 <= t0)
            t = t1;
        else
            t = t0;
        if (tObjsBuffer(pixelNum) == -1 || t <= tBuffer(pixelNum, 0)) {
            tempVec = Ro + (Rd * t);
            tempNorm = Vector3d((tempVec(0) - Xc) / rs, (tempVec(1) - Yc) / rs,
                                (tempVec(2) - Zc) / rs)
                           .normalized();
#pragma omp critical
            {
                tBuffer.row(pixelNum) << t, tempVec(0), tempVec(1), tempVec(2),
                    tempNorm(0), tempNorm(1), tempNorm(2);
                tObjsBuffer(pixelNum) = objNum;
            }
            return true;
        }
    }
    return false;
}

// Parses "Intersectable" specific input file lines
void Intersectable::initSetts(ifstream *inSetts) {
    const regex translateR(
        "^\\s+translate:\\s*(-?[0-9]+\\.[0-9]+),\\s*(-?[0-9]+"
        "\\.[0-9]+),\\s*(-?[0-9]+\\.[0-9]+)\\s*$");
    const regex rotateR(
        "^\\s+rotate:\\s*(-?[0-9]+\\.[0-9]+),\\s*(-?[0-9]+\\.[0-"
        "9]+),\\s*(-?[0-9]+\\.[0-9]+)\\s*$");
    const regex colorR(
        "^\\s+color:\\s*(-?[0-9]+\\.[0-9]+),\\s*(-?[0-9]+\\.[0-9]+"
        "),\\s*(-?[0-9]+\\.[0-9]+)\\s*$");
    const regex scaleR("^\\s+scale:\\s*([0-9]+\\.[0-9]+)\\s*$");
    const regex kdR("^\\s+kd:\\s*([0-9]+\\.[0-9]+)\\s*$");
    const regex ksR("^\\s+ks:\\s*([0-9]+\\.[0-9]+)\\s*$");
    const regex kaR("^\\s+ka:\\s*([0-9]+\\.[0-9]+)\\s*$");
    const regex nR("^\\s+n:\\s*([0-9]+\\.[0-9]+)\\s*$");
    const regex indentedR("^\\s+");
    bool doneSetts = false;
    smatch sm;
    string line;

    getline(*inSetts, line);
    while (!doneSetts) {
        if (line.front() == '#') {
        } else if (regex_search(line, sm, translateR)) {
            translate = Vector3d(stod(sm[1]), stod(sm[2]), stod(sm[3]));
        } else if (regex_search(line, sm, rotateR)) {
            rotate =
                Vector3d(stod(sm[1]), stod(sm[2]), stod(sm[3])) * PI / 180.0;
        } else if (regex_search(line, sm, colorR)) {
            color = Vector3d(stod(sm[1]), stod(sm[2]), stod(sm[3]));
        } else if (regex_search(line, sm, scaleR)) {
            scale = stod(sm[1]);
        } else if (regex_search(line, sm, kdR)) {
            kd = stod(sm[1]);
        } else if (regex_search(line, sm, ksR)) {
            ks = stod(sm[1]);
        } else if (regex_search(line, sm, kaR)) {
            ka = stod(sm[1]);
        } else if (regex_search(line, sm, nR)) {
            n = stod(sm[1]);
        } else if (regex_search(line, sm, indentedR)) {
            cout << "Ignoring input line (in Intersectable parsing): " << line
                 << endl;
        } else if (!line.compare("")) {
            doneSetts = true;
        } else {
            cout << "Ignoring input line: " << line << endl;
        }
        if (!doneSetts) getline(*inSetts, line);
    }
}

// smf files are parsed in this constructor
SMFModel::SMFModel(ifstream *inSMF, ifstream *inSetts) {
    const regex vertexR("^v\\s*(.*)$");
    const regex faceR("^f\\s*(.*)$");
    smatch sm;
    string line;
    double temp1, temp2, temp3;
    MatrixXd tempRot = MatrixXd::Zero(4, 4);
    Vector4d tempPt;
    Vector3d a, b, c, tempNorm;

    type = SMF;
    initSetts(inSetts);

    vertices = MatrixXd::Zero(0, 3);
    facesStored = MatrixXi::Zero(0, 3);

    minVec << DBL_MAX, DBL_MAX, DBL_MAX;
    maxVec << DBL_MIN, DBL_MIN, DBL_MIN;

    getline(*inSMF, line);

    while (!inSMF->eof()) {
        if (line.front() == '#' || !line.compare("")) {
        } else if (regex_search(line, sm, vertexR)) {
            vertices.conservativeResize(vertices.rows() + 1, NoChange);
            stringstream ss(sm[1]);
            ss >> temp1;
            ss >> temp2;
            ss >> temp3;
            vertices.row(vertices.rows() - 1) << temp1, temp2, temp3;
        } else if (regex_search(line, sm, faceR)) {
            facesStored.conservativeResize(facesStored.rows() + 1, NoChange);
            stringstream ss(sm[1]);
            ss >> temp1;
            ss >> temp2;
            ss >> temp3;
            facesStored.row(facesStored.rows() - 1) << temp1, temp2, temp3;
        }
        getline(*inSMF, line);
    }
    int i;
    // Transformations
#pragma omp parallel for num_threads(NUM_THREADS) \
    firstprivate(tempPt, tempRot) private(i)
    for (i = 0; i < vertices.rows(); i++) {
        tempPt << vertices(i, 0), vertices(i, 1), vertices(i, 2), 1;
        tempRot << 1, 0, 0, 0, 0, cos(rotate(0)), -sin(rotate(0)), 0, 0,
            sin(rotate(0)), cos(rotate(0)), 0, 0, 0, 0, 1;
        tempPt = tempRot * tempPt;
        tempRot << cos(rotate(1)), 0, sin(rotate(1)), 0, 0, 1, 0, 0,
            -sin(rotate(1)), 0, cos(rotate(1)), 0, 0, 0, 0, 1;
        tempPt = tempRot * tempPt;
        tempRot << cos(rotate(2)), -sin(rotate(2)), 0, 0, sin(rotate(2)),
            cos(rotate(2)), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
        tempPt = tempRot * tempPt;
        for (int l = 0; l < 3; l++) {
            if (tempPt(l) < minVec(l)) minVec(l) = tempPt(l);
            if (tempPt(l) > maxVec(l)) maxVec(l) = tempPt(l);
        }
#pragma omp critical
        {
            vertices.row(i) = tempPt.head(3);
            vertices.row(i) += translate;
        }
    }
    vertices = vertices.rowwise() - minVec.transpose();
    maxVec -= minVec;
    minVec -= minVec;
    // Normal averaging calculation
    normals = MatrixXd::Zero(vertices.rows(), 3);
    for (i = 0; i < facesStored.rows(); i++) {
        a = vertices.row(facesStored(i, 0) - 1);
        b = vertices.row(facesStored(i, 1) - 1);
        c = vertices.row(facesStored(i, 2) - 1);
        tempNorm = (b - a).cross(c - a).normalized();
        for (int j = 0; j < 3; j++)
            normals.row(facesStored(i, j) - 1) += tempNorm;
    }
    for (i = 0; i < normals.rows(); i++)
        normals.row(i) = normals.row(i).normalized();
}

bool SMFModel::intersectPixel(Vector3d ray, Settings *setts, int pixelNum,
                              int objNum) {
    return intersectPixel(ray, setts, &facesStored, pixelNum, objNum);
}

// Intersects pixel with triangles of an SMF model and updates data stores
// accordingly
bool SMFModel::intersectPixel(Vector3d ray, Settings *setts, MatrixXi *facesIn,
                              int pixelNum, int objNum) {
    MatrixXi faces = *facesIn;
    Matrix3d temp;
    Vector3d tempVec, tempNorm, normA, normB, normC;
    MatrixXd calc = MatrixXd::Zero(3, 4);
    Vector3d R = setts->cameraLoc;
    Vector3d a, b, c;
    double Adet, beta, gamma, t;
    bool updated = false;
    for (int i = 0; i < faces.rows(); i++) {
        a = vertices.row(faces.row(i)(0) - 1);
        b = vertices.row(faces.row(i)(1) - 1);
        c = vertices.row(faces.row(i)(2) - 1);
        // Cache these calculations for small performance gain
        calc << a(0) - b(0), a(0) - c(0), a(0) - R(0), ray(0), a(1) - b(1),
            a(1) - c(1), a(1) - R(1), ray(1), a(2) - b(2), a(2) - c(2),
            a(2) - R(2), ray(2);
        temp = calc(all, {0, 1, 3});
        Adet = temp.determinant();
        temp = calc(all, {2, 1, 3});
        beta = temp.determinant() / Adet;
        if (beta >= 0) {
            temp = calc(all, {0, 2, 3});
            gamma = temp.determinant() / Adet;
            if (gamma >= 0 && beta + gamma <= 1) {
                temp = calc(all, {0, 1, 2});
                t = temp.determinant() / Adet;
                if (t >= 0) {
                    if (tObjsBuffer(pixelNum) == -1 ||
                        t <= tBuffer(pixelNum, 0)) {
                        updated = true;
                        tempVec = (1 - beta - gamma) * a + beta * b + gamma * c;
                        tempNorm = tempVec;
                        tempNorm =
                            ((1 - beta - gamma) * normals.row(faces(i, 0) - 1) +
                             beta * normals.row(faces(i, 1) - 1) +
                             gamma * normals.row(faces(i, 2) - 1))
                                .normalized();
#pragma omp critical
                        {
                            tBuffer.row(pixelNum) << t, tempVec(0), tempVec(1),
                                tempVec(2), tempNorm(0), tempNorm(1),
                                tempNorm(2);
                            tObjsBuffer(pixelNum) = objNum;
                        }
                    }
                }
            }
        }
    }
    return updated;
}

Box::Box(Sphere *sphere) {
    minVec = sphere->center.array() + sphere->rs;
    maxVec = sphere->center.array() - sphere->rs;
    type = BOX;
}

Box::Box(SMFModel *smf) {
    minVec = Vector3d(smf->vertices.col(0).minCoeff(),
                      smf->vertices.col(1).minCoeff(),
                      smf->vertices.col(2).minCoeff());
    maxVec = Vector3d(smf->vertices.col(0).maxCoeff(),
                      smf->vertices.col(1).maxCoeff(),
                      smf->vertices.col(2).maxCoeff());
    type = BOX;
}

Box::Box(SubSMF *subsmf) {
    MatrixXi faces = subsmf->faces;
    // Just defaulting these to something
    Vector3d currRow = subsmf->smf->vertices.row(faces(0, 0) - 1);
    minVec = currRow.replicate(1, 1);
    maxVec = currRow.replicate(1, 1);
    // Loop through faces in sub smf only
    for (int i = 0; i < faces.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            currRow = subsmf->smf->vertices.row(faces(i, j) - 1);
            for (int k = 0; k < 3; k++) {
                if (currRow(k) < minVec(k))
                    minVec(k) = currRow(k);
                else if (currRow(k) > maxVec(k))
                    maxVec(k) = currRow(k);
            }
        }
    }
    type = BOX;
}

// Perform box intersection algorithm from the slides
bool Box::intersectPixel(Vector3d Rd, Settings *setts, int pixelNum,
                         int objNum) {
    double temp;
    Vector3d t1 = (minVec - setts->cameraLoc).array() / Rd.array();
    Vector3d t2 = (maxVec - setts->cameraLoc).array() / Rd.array();
    for (int i = 0; i < 3; i++) {
        if (t1(i) > t2(i)) {
            temp = t1(i);
            t1(i) = t2(i);
            t2(i) = temp;
        }
    }
    double tnear = t1.maxCoeff();
    double tfar = t2.minCoeff();
    return (tnear <= tfar) && (tfar >= 0);
}

// Sets up the tree data structure for an object's bounding volume hierarchy
ACGBVH::ACGBVH(Intersectable *inObj, int inMaxObjCount) {
    type = HIERARCH;
    obj = inObj;
    maxObjCount = inMaxObjCount;
    if (obj->type == SPHERE) {
        boundingVol = new Box(reinterpret_cast<Sphere *>(inObj));
        isLeaf = true;
    } else if (obj->type == SMF) {
        boundingVol = new Box(reinterpret_cast<SMFModel *>(inObj));
        if ((reinterpret_cast<SMFModel *>(inObj))->facesStored.rows() > maxObjCount) {
            isLeaf = false;
            createChildren();
        } else {
            isLeaf = true;
        }
    } else if (obj->type == SUBSMF) {
        boundingVol = new Box(reinterpret_cast<SubSMF *>(inObj));
        if ((reinterpret_cast<SubSMF *>(inObj))->faces.rows() > maxObjCount) {
            isLeaf = false;
            createChildren();
        } else {
            isLeaf = true;
        }
    } else {
        cout << "Error: invalid shape type '" << obj->type
             << "' passed to ACGBVH\n";
        exit(1);
    }
}

// Defers to its children or object stored
bool ACGBVH::intersectPixel(Vector3d ray, Settings *setts, int pixelNum,
                            int objNum) {
    if (boundingVol->intersectPixel(ray, setts, pixelNum, objNum)) {
        if (isLeaf) {
            return obj->intersectPixel(ray, setts, pixelNum, objNum);
        } else {
            children[0]->intersectPixel(ray, setts, pixelNum, objNum);
            children[1]->intersectPixel(ray, setts, pixelNum, objNum);
            return true;
        }
    }
    return false;
}

// Does the volume division
void ACGBVH::createChildren() {
    Array3d diff = boundingVol->maxVec - boundingVol->minVec;
    Array3d currRow;
    double split;
    int axis;
    SMFModel *smf;
    MatrixXi newFaces1 = MatrixXi::Zero(0, 3);
    MatrixXi newFaces2 = MatrixXi::Zero(0, 3);
    MatrixXi currFaces;
    if (diff(0) >= diff(1) && diff(0) >= diff(2))
        axis = 0;
    else if (diff(1) > diff(0) && diff(1) >= diff(2))
        axis = 1;
    else
        axis = 2;
    split = boundingVol->minVec(axis) + (diff(axis) / 2.0);
    if (obj->type == SUBSMF) {
        smf = (reinterpret_cast<SubSMF *>(obj))->smf;
        currFaces = (reinterpret_cast<SubSMF *>(obj))->faces.replicate(1, 1);
    } else {
        smf = reinterpret_cast<SMFModel *>(obj);
        currFaces = smf->facesStored.replicate(1, 1);
    }
    for (int i = 0; i < currFaces.rows(); i++) {
        currRow << smf->vertices(currFaces(i, 0) - 1, axis),
            smf->vertices(currFaces(i, 1) - 1, axis),
            smf->vertices(currFaces(i, 2) - 1, axis);
        if ((currRow > split).any()) {
            newFaces1.conservativeResize(newFaces1.rows() + 1,
                                         newFaces1.cols());
            newFaces1.row(newFaces1.rows() - 1) = currFaces.row(i);
        } else {
            newFaces2.conservativeResize(newFaces2.rows() + 1,
                                         newFaces2.cols());
            newFaces2.row(newFaces2.rows() - 1) = currFaces.row(i);
        }
    }
    SubSMF *subsmf1 = new SubSMF(smf, newFaces1);
    SubSMF *subsmf2 = new SubSMF(smf, newFaces2);
    ACGBVH *temp1 = new ACGBVH(subsmf1, maxObjCount);
    children.push_back(temp1);
    ACGBVH *temp2 = new ACGBVH(subsmf2, maxObjCount);
    children.push_back(temp2);
}

SubSMF::SubSMF(SMFModel *smfIn, MatrixXi facesIn) {
    type = SUBSMF;
    smf = smfIn;
    faces = facesIn;
}

bool SubSMF::intersectPixel(Vector3d ray, Settings *setts, int pixelNum,
                            int objNum) {
    return smf->intersectPixel(ray, setts, &faces, pixelNum, objNum);
}
