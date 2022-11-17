# 4D Multilayer Modeler

![Screenshot egg chair](https://www.cs.drexel.edu/~jcp353/images/egg_chair.png "Screenshot egg chair")

This repo contains the code from my paper [Design and Simulation of Multi-Tiered 4D Printed Objects](http://cad-journal.net/files/vol_20/CAD_20(3)_2023_489-506.pdf). See the paper for details of its capabilities.

## Dependencies
- Build tools: CMake, Make and GCC
- Libraries: Qt5, VTK and OpenMP must be installed on your system. The other dependencies are included in the `lib` directory. CMake may request more dependencies when you go to build the `Makefile`, I've found what it asks for to be somewhat inconsistent however, so you should just install those one at a time following any CMake error messages.

## Compatibility
I have successfully built and ran this on Arch Linux and Ubuntu, using only packages found in the official repositories. Despite some effort, I have not been able to get a version working on Mac so this is currently Linux only.

## Building

    mkdir build
    cd build
    cmake ..
    make

I have not seen this when running similar Qt/VTK projects on Linux, but I have needed to pass `-DQt5_ROOT=...` to CMake on Mac before, so if CMake cannot find Qt5 give that argument a try.

## Running
Make sure you are in the `build` directory, then run

    ./4d_modeler -f ../input/INPUT_FILE

## Operation

![Screenshot face tessellation](https://www.cs.drexel.edu/~jcp353/images/face_tess.png "Screenshot face tessellation")

- The `Continue` button moves the program to the next step in the process.
- The `Center` button recenters the model within the window.
- The `Print` button outputs the model as a `.stp` file.
- Standard mouse and scroll wheel controls can be used to pan and zoom the model.

## Performance Caveat
The way I set up VTK has each individual edge and vertex as distinct VTK "actors." This is the most intuitive way to set up VTK programs, works great for small scenes, and allows for entity picking which I used during debugging. However, I discovered that once you get a few hundred/thousand actors in a scene, the rendering performance drops off considerably. The alternative is to use `vtkGlyph3DMapper`'s and/or other similar ways of bundling multiple primitives into one actor. In my WIP follow up to this project I've switched to that configuration, however I don't intend to backport the change here. So, just be warned that panning and zooming can be laggy for larger models.
