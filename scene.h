#pragma once

#include "util.h"

// Basically just a struct to hold the VTK actors and souurces
class Scene {
    public:
        Scene();
        vector<vtkSphereSource *> spheres;
        vector<vtkLineSource *> lines;
        vector<vtkActor *> sphereActors;
        vector<vtkActor *> lineActors;
};
