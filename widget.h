#pragma once

#include "./graddesc.h"
#include "./util.h"
#include "./scene.h"
#include "./raytrace.h"

#ifndef WIDGET_H
    #define WIDGET_H

namespace Ui {
    class Widget;
}

struct VertexData {
    double u, v;
    int row, col;
    MatrixXd pt1;
};

// Handles everything related to the Qt window and rendering the scene
class Widget : public QWidget {
        Q_OBJECT
    public:
        explicit Widget(QWidget *parent = 0, GradDesc * = nullptr);
        std::unordered_map<vtkActor*, VertexData*> vertexDataHash;
        GradDesc *gd;
        Ui::Widget *ui;
        ~Widget();

    private Q_SLOTS:
        void on_continueButton_clicked();
        void on_centerButton_clicked();
        void on_printButton_clicked();
        void runShapeOp();

    private:
        void addLayers();
        void triangulateGrid();
        void renderSMF();
        int layersDrawn;
        bool hasAddedLayers, triangulate, fullyDone;
        Scene scene;
        vtkSmartPointer<vtkRenderer> backRenderer;
        vtkSmartPointer<vtkOBJExporter> exporter;
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
};

void clickCallbackFunction(vtkObject* caller, long unsigned int eventId,
        void* clientData, void* callData);

#endif  // WIDGET_H
