#pragma once
#include "util.h"

class Settings;
class Intersectable;
class Box;
class Sphere;
class SMFModel;
class ACGBVH;
class SubSMF;

enum Shape { SMF, SUBSMF, SPHERE, BOX, HIERARCH };

// Class containing settings for image
class Settings {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        int x, y;
        double theta, d, h, Sj, Sk, uStart, vStart, uLength, vLength;
        Vector3d cameraLoc, Zv, Vup, Xv, Yv, P00, ambient;

        Settings(ifstream *inFile, double*);
};

// Interface for objects in scene
class Intersectable {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        virtual bool intersectPixel(Vector3d ray, Settings *setts, int pixelNum,
                                    int objNum) = 0;
        double scale = 1.0;
        double kd = 0.5;
        double ks = 0.5;
        double ka = 0.5;
        double n = 10.0;
        Shape type;
        Vector3d translate = Vector3d(0.0, 0.0, 0.0);
        Vector3d rotate = Vector3d(0.0, 0.0, 0.0);
        Vector3d color = Vector3d(0.1, 0.1, 0.1);
        Vector3d minVec, maxVec;

    protected:
        void initSetts(ifstream *inSetts);
};

// Bounding box implementation
class Box : public Intersectable {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Box(SubSMF *subsmf);
        Box(Sphere *sphere);
        Box(SMFModel *smf);
        bool intersectPixel(Vector3d ray, Settings *setts, int pixelNum,
                            int objNum);
        Vector3d maxVec, minVec;
};

// Sphere for intersecting
class Sphere : public Intersectable {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Sphere(double radius, ifstream *inSetts);

        bool intersectPixel(Vector3d Rd, Settings *setts, int pixelNum,
                            int objNum);
        Vector3d center;
        double rs;

    private:
        double Xc, Yc, Zc;
};

// Triangle mesh class for intersecting
class SMFModel : public Intersectable {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        SMFModel(ifstream *inSMF, ifstream *inSetts);
        bool intersectPixel(Vector3d ray, Settings *setts, int pixelNum,
                            int objNum);
        bool intersectPixel(Vector3d ray, Settings *setts, MatrixXi *facesIn,
                            int pixelNum, int objNum);
        MatrixXd vertices;
        MatrixXi facesStored;
        // Vector3d minVec, maxVec;

    private:
        MatrixXd normals;
};

// Workaround class to only store some faces in a node of the ACGBVH tree
// without copying everything stored in the SMFModel
class SubSMF : public Intersectable {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        SubSMF(SMFModel *smfIn, MatrixXi facesIn);
        bool intersectPixel(Vector3d ray, Settings *setts, int pixelNum,
                            int objNum);
        SMFModel *smf;
        MatrixXi faces;
};

// Class to set up tree structure and hold other objects involved
class ACGBVH : public Intersectable {
    public:
        ACGBVH(Intersectable *inObj, int inMaxObjCount);
        Intersectable *obj;
        bool intersectPixel(Vector3d ray, Settings *setts, int pixelNum,
                            int objNum);

    private:
        int maxObjCount;
        bool isLeaf;
        Box *boundingVol;
        vector<ACGBVH *> children;
        void createChildren();
};

static vector<ACGBVH *> objects;
static MatrixXd lights;
static MatrixXd tBuffer;
static VectorXi tObjsBuffer;

void printPNGFile(MatrixXd pixelMap, Settings *setts);
void intersectForAllPixels(vector<ACGBVH *> objects, Settings *setts,
                           MatrixXd *result);
void shadePixels(vector<ACGBVH *> objects, Settings *setts, MatrixXd *result);
std::tuple<int, int, MatrixXd&> rayTrace(string inputFileName, MatrixXd&, MatrixXi&, double*, int, int);
