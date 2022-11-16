#include "./widget.h"
#include "./ui_widget.h"

Widget::Widget(QWidget *parent, GradDesc *gdIn) : QWidget(parent), ui(new Ui::Widget) {
    ui->setupUi(this);

    gd = gdIn;
    scene = Scene();

    backRenderer = vtkSmartPointer<vtkRenderer>::New();
    vtkNew<vtkNamedColors> colors;

    layersDrawn = 0;
    hasAddedLayers = false;
    triangulate = true;
    fullyDone = false;

    // Set up all the VTK boilerplate
    backRenderer->SetBackground(colors->GetColor3d("Black").GetData());

    vtkSmartPointer<vtkGenericOpenGLRenderWindow> backWindow =
        vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    backWindow->AddRenderer(backRenderer);
    ui->qvtkWidget->setRenderWindow(backWindow);

    renderWindowInteractor = backWindow->GetInteractor();
    vtkNew<vtkCallbackCommand> clickCallback;
    clickCallback->SetCallback(clickCallbackFunction);
    clickCallback->SetClientData(this);
    renderWindowInteractor->AddObserver(vtkCommand::LeftButtonPressEvent,
            clickCallback);

    exporter = vtkSmartPointer<vtkOBJExporter>::New();
    exporter->SetActiveRenderer(backRenderer);
    exporter->SetRenderWindow(backWindow);
    const char *objFileName = "model_export";
    exporter->SetFilePrefix(objFileName);
    if (gd->isSMF()) {
        renderSMF();
    } else {
        addLayers();
        runShapeOp();
    }
}

// These are global for convenience with the clickCallbackFunction since it
// cannot be a class instance method
static vtkSmartPointer<vtkActor> lastActor = nullptr;
static double lastColor[3];

// Handles click inputs from the user
void clickCallbackFunction(vtkObject* caller,
        long unsigned int eventId,
        void* clientData,
        void* vtkNotUsed(callData))
{
    vtkNew<vtkNamedColors> colors;
    auto interactor = reinterpret_cast<vtkRenderWindowInteractor*>(caller);
    auto widget = reinterpret_cast<Widget*>(clientData);
    auto gd = widget->gd;
    int* clickPos = interactor->GetEventPosition();
    // Pick from this location.
    vtkNew<vtkPropPicker> picker;
    picker->Pick(clickPos[0], clickPos[1], 0, interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());


    auto pickedActor = picker->GetActor();
    if (pickedActor == nullptr)
    {
        std::cout << "No actor picked." << std::endl;
    }
    else
    {
        double* pos = picker->GetPickPosition();
        auto lookup = widget->vertexDataHash.find(pickedActor);
        if (lookup != widget->vertexDataHash.end()) {
            auto vd = lookup->second;
            stringstream ss;
            ss << vd->pt1;
            std::string output = fmt::format(
                    "Picked actor: {}\nrow,col: ({}, {}) | u,v: ({:.2f}, {:.2f})\nPick Position: ({:.2f}, {:.2f}, {:.2f})\nInput position: {}\n",
                    static_cast<void*>(picker->GetActor()), vd->row, vd->col, vd->u, vd->v,
                    pos[0], pos[1], pos[2], ss.str());
            widget->ui->pickInfo->setPlainText(
                    QString::fromStdString(output));
        }
        double temp[3];
        temp[0] = pickedActor->GetProperty()->GetColor()[0];
        temp[1] = pickedActor->GetProperty()->GetColor()[1];
        temp[2] = pickedActor->GetProperty()->GetColor()[2];
        pickedActor->GetProperty()->SetColor(
                colors->GetColor3d("Yellow").GetData());
        if (lastActor != nullptr) {
            lastActor->GetProperty()->SetColor(lastColor);
        }
        lastActor = pickedActor;
        lastColor[0] = temp[0];
        lastColor[1] = temp[1];
        lastColor[2] = temp[2];
    }
}

// Renders the VTK scene for the first ray tracing step when the input is an SMF
void Widget::renderSMF() {

    vtkNew<vtkNamedColors> colors;
    MatrixXd verts = gd->getVertices();
    MatrixXi faces = gd->getFaces();
    faces = faces - MatrixXi::Ones(faces.rows(), faces.cols());
    vtkNew<vtkSphereSource> sphereSource;
    auto centroid = verts.colwise().sum() / verts.rows();
    sphereSource->SetCenter(centroid(0), centroid(1), centroid(2));
    sphereSource->SetRadius(1.0);
    // Make the surface smooth.
    sphereSource->SetPhiResolution(6);
    sphereSource->SetThetaResolution(6);

    vtkNew<vtkPolyDataMapper> mapper1;
    mapper1->SetInputConnection(sphereSource->GetOutputPort());

    vtkNew<vtkActor> actor1;
    actor1->SetMapper(mapper1);
    actor1->GetProperty()->SetColor(
            colors->GetColor3d("Red").GetData());
    scene.sphereActors.push_back(actor1.GetPointer());
    backRenderer->AddActor(actor1);

    for (auto pt : gd->intersectPoints.rowwise()) {
        vtkNew<vtkSphereSource> sphereSource1;
        sphereSource1->SetCenter(pt(0), pt(1), pt(2));
        sphereSource1->SetRadius(0.07);
        // Make the surface smooth.
        sphereSource1->SetPhiResolution(6);
        sphereSource1->SetThetaResolution(6);
        vtkNew<vtkPolyDataMapper> mapper2;
        mapper2->SetInputConnection(sphereSource1->GetOutputPort());
        vtkNew<vtkActor> actor2;
        actor2->SetMapper(mapper2);
        actor2->GetProperty()->SetColor(
                colors->GetColor3d("Green").GetData());
        scene.sphereActors.push_back(actor2.GetPointer());
        backRenderer->AddActor(actor2);
    }

    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();

    for (auto pt : verts.rowwise()) {
        points->InsertNextPoint(pt(0), pt(1), pt(2));
    }

    vtkNew<vtkCellArray> triangles;

    for (auto t : faces.rowwise()) {
        vtkSmartPointer<vtkTriangle> triangle =
            vtkSmartPointer<vtkTriangle>::New();

        triangle->GetPointIds()->SetId(0, t(0));

        triangle->GetPointIds()->SetId(1, t(1));

        triangle->GetPointIds()->SetId(2, t(2));

        triangles->InsertNextCell(triangle);
    }

    // Create a polydata object
    vtkSmartPointer<vtkPolyData> polyData =
        vtkSmartPointer<vtkPolyData>::New();

    // Add the geometry and topology to the polydata
    polyData->SetPoints(points);
    polyData->SetPolys(triangles);
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(polyData);

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(
            colors->GetColor3d("Cornsilk").GetData());
    scene.sphereActors.push_back(actor.GetPointer());
    backRenderer->AddActor(actor);

    // Splines

    // Something of a typing monster but it works
    std::vector<std::tuple<std::vector<MatrixXd*>, const char*>> directions;
    directions.push_back(std::make_tuple(gd->curvePointsU, "Tomato"));
    directions.push_back(std::make_tuple(gd->curvePointsV, "Blue"));
    for (auto curvePointsX : directions) {
        for (MatrixXd *curvePoints : std::get<0>(curvePointsX)) {
            vtkNew<vtkPoints> vtkPoints;
            for (auto point : curvePoints->rowwise())
                vtkPoints->InsertNextPoint(point(0), point(1), point(2));
            vtkNew<vtkPolyLine> polyLine;
            int ptCount = curvePoints->rows();
            polyLine->GetPointIds()->SetNumberOfIds(ptCount);
            for (int i = 0; i < ptCount; i++)
                polyLine->GetPointIds()->SetId(i, i);
            vtkNew<vtkCellArray> cells;
            cells->InsertNextCell(polyLine);
            // Create a polydata to store everything in
            vtkNew<vtkPolyData> polyDataLine;

            // Add the points to the dataset
            polyDataLine->SetPoints(vtkPoints);

            // Add the lines to the dataset
            polyDataLine->SetLines(cells);

            vtkNew<vtkTubeFilter> tubes;
            tubes->SetInputData(polyDataLine);
            tubes->SetRadius(0.04);
            tubes->SetNumberOfSides(10);

            // Set up actor and mapper
            vtkNew<vtkPolyDataMapper> mapperLine;
            mapperLine->SetInputConnection(tubes->GetOutputPort());

            vtkNew<vtkActor> actorLine;
            actorLine->SetMapper(mapperLine);
            actorLine->GetProperty()->SetColor(colors->GetColor3d(std::get<1>(curvePointsX)).GetData());
            backRenderer->AddActor(actorLine);
        }
    }

    backRenderer->GetActiveCamera()->SetPosition(-500.0, centroid(1), centroid(2));
    backRenderer->ResetCamera();
    backRenderer->GetActiveCamera()->SetFocalPoint(centroid(0), centroid(1), centroid(2));
}

// Adds more layers to the grid scene after in plane optimization has terminated
void Widget::addLayers() {
    int row, col, cols = gd->getCols(), rows = gd->getRows();
    double spacing = gd->getSpacing();
    vtkNew<vtkNamedColors> colors;

    for (int layer = layersDrawn; layer < gd->getLayers(); layer++) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                row = i;
                col = j;
                // Create a sphere
                vtkNew<vtkSphereSource> sphereSource;
                scene.spheres.push_back(sphereSource.GetPointer());
                sphereSource->SetCenter(i * spacing, j * spacing,
                        layer * spacing);
                sphereSource->SetRadius(spacing * 0.075);
                // Make the surface smooth.
                sphereSource->SetPhiResolution(6);
                sphereSource->SetThetaResolution(6);

                vtkNew<vtkPolyDataMapper> mapper;
                mapper->SetInputConnection(sphereSource->GetOutputPort());

                vtkNew<vtkActor> actor;
                actor->SetMapper(mapper);
                actor->GetProperty()->SetColor(
                        colors->GetColor3d("Cornsilk").GetData());
                scene.sphereActors.push_back(actor.GetPointer());
                backRenderer->AddActor(actor);
                VertexData *vd = new VertexData;
                vd->row = row;
                vd->col = col;
                vd->u = static_cast<double>(row) / static_cast<double>(rows - 1);
                vd->v = static_cast<double>(col) / static_cast<double>(cols - 1);
                vertexDataHash[actor.GetPointer()] = vd;
                MatrixXd pt1(1,3), pt2(1,3);
                pt1 << i * spacing, j * spacing, layer * spacing;
                vd->pt1 = pt1;
                // Right
                if (col != cols - 1) {
                    vtkNew<vtkLineSource> lineSource;
                    scene.lines.push_back(lineSource.GetPointer());
                    lineSource->SetPoint1(pt1(0), pt1(1), pt1(2));
                    pt2 << (i + 1) * spacing, j * spacing,
                            layer * spacing;
                    lineSource->SetPoint2(pt2(0), pt2(1), pt2(2));
                    lineSource->SetResolution(6);
                    lineSource->Update();
                    vtkNew<vtkTubeFilter> tubeFilter;
                    tubeFilter->SetInputConnection(lineSource->GetOutputPort());
                    tubeFilter->SetNumberOfSides(8);
                    tubeFilter->SetRadius(spacing * 0.025);
                    tubeFilter->Update();
                    vtkNew<vtkDataSetMapper> mapper;
                    mapper->SetInputConnection(tubeFilter->GetOutputPort());
                    mapper->ScalarVisibilityOff();

                    vtkNew<vtkActor> actor;
                    actor->SetMapper(mapper);

                    if (layer == 0) {
                        Vector3d color = ratioToRGB(distance(pt1, pt2) / gd->getRestU(row, col, layer));
                        actor->GetProperty()->SetColor(color(0), color(1), color(2));
                    }
                    scene.lineActors.push_back(actor.GetPointer());

                    backRenderer->AddActor(actor);
                }
                // Down in plane
                if (row != rows - 1) {
                    vtkNew<vtkLineSource> lineSource;
                    scene.lines.push_back(lineSource.GetPointer());
                    lineSource->SetPoint1(pt1(0), pt1(1), pt1(2));
                    pt2 << i * spacing, (j + 1) * spacing,
                            layer * spacing;
                    lineSource->SetPoint2(pt2(0), pt2(1), pt2(2));
                    lineSource->SetResolution(6);
                    lineSource->Update();
                    vtkNew<vtkTubeFilter> tubeFilter;
                    tubeFilter->SetInputConnection(lineSource->GetOutputPort());
                    tubeFilter->SetRadius(spacing * 0.025);
                    tubeFilter->SetNumberOfSides(8);
                    tubeFilter->Update();
                    vtkNew<vtkDataSetMapper> mapper;
                    mapper->SetInputConnection(tubeFilter->GetOutputPort());
                    mapper->ScalarVisibilityOff();

                    vtkNew<vtkActor> actor;
                    actor->SetMapper(mapper);
                    if (layer == 0) {
                        Vector3d color = ratioToRGB(distance(pt1, pt2) / gd->getRestV(row, col, layer));
                        actor->GetProperty()->SetColor(color(0), color(1), color(2));
                    }
                    scene.lineActors.push_back(actor.GetPointer());

                    backRenderer->AddActor(actor);
                }
                // Down in z
                if (layer != 0) {
                    vtkNew<vtkLineSource> lineSource;
                    scene.lines.push_back(lineSource.GetPointer());
                    lineSource->SetPoint1(pt1(0), pt1(1), pt1(2));
                    pt2 << i * spacing, j * spacing,
                            (layer - 1) * spacing;
                    lineSource->SetPoint2(pt2(0), pt2(1), pt2(2));
                    lineSource->SetResolution(6);
                    lineSource->Update();
                    vtkNew<vtkTubeFilter> tubeFilter;
                    tubeFilter->SetInputConnection(lineSource->GetOutputPort());
                    tubeFilter->SetRadius(spacing * 0.025);
                    tubeFilter->SetNumberOfSides(8);
                    tubeFilter->Update();
                    vtkNew<vtkDataSetMapper> mapper;
                    mapper->SetInputConnection(tubeFilter->GetOutputPort());
                    mapper->ScalarVisibilityOff();

                    vtkNew<vtkActor> actor;
                    actor->SetMapper(mapper);
                    scene.lineActors.push_back(actor.GetPointer());

                    backRenderer->AddActor(actor);
                }
            }
        }
    }
    layersDrawn = gd->getLayers();
}

// Calls ShapeOp
void Widget::runShapeOp() {
    cout << "Initial System Energy: " << gd->getSysEnergy() << endl;
    Vector2d disps = gd->getAvgSpringDisp();
    cout << "Initial Linear Error: " << disps[0] << endl;
    cout << "Initial Angular Error: " << disps[1] << endl;
    if (gd->stringDeq.size()) {
        ui->stepTextBrowser->setPlainText(
                QString::fromStdString(gd->stringDeq.front()));
    } else {
        gd->shapeOptimizeGrid();
        cout << "Final System Energy: " << gd->getSysEnergy() << endl;
        disps = gd->getAvgSpringDisp();
        cout << "Final Linear Error: " << disps[0] << endl;
        cout << "Final Angular Error: " << disps[1] << endl;
        ui->stepTextBrowser->setPlainText(
                QString::fromStdString("No remaining operations"));
    }
    gd->updateVtkGrid(scene);
    backRenderer->GetRenderWindow()->Render();
}

Widget::~Widget() {
    delete ui;
}

// The very final step of the GUI process where the middle layer of the grid is
// tessellated for easier comparison to the input patch
void Widget::triangulateGrid() {
    auto actors = backRenderer->GetActors();
    actors->InitTraversal();
    while (actors->GetNumberOfItems())
        backRenderer->RemoveActor(actors->GetNextItem());
    vtkNew<vtkNamedColors> colors;
    MatrixXd positions = gd->getMiddleLayer();
    
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();

    int rows = gd->getRows();
    int cols = gd->getCols();
    int count = 0;
    for (int i = 0; i < positions.rows(); i++) {
        int row = i / cols;
        int col = i % cols;
        if (row < rows - 1 && col < cols - 1) {
            Vector3d pt1 = positions.row(i);
            Vector3d pt2 = positions.row(i+1);
            Vector3d pt3 = positions.row(i+cols);
            Vector3d pt4 = positions.row(i+1+cols);
            points->InsertNextPoint(pt1(0), pt1(1), pt1(2));
            points->InsertNextPoint(pt4(0), pt4(1), pt4(2));
            points->InsertNextPoint(pt2(0), pt2(1), pt2(2));
            points->InsertNextPoint(pt1(0), pt1(1), pt1(2));
            points->InsertNextPoint(pt3(0), pt3(1), pt3(2));
            points->InsertNextPoint(pt4(0), pt4(1), pt4(2));
            count += 6;
        }
    }
    vtkNew<vtkCellArray> triangles;

    for (int i = 0; i < count; i += 3) {
        vtkSmartPointer<vtkTriangle> triangle =
            vtkSmartPointer<vtkTriangle>::New();

        triangle->GetPointIds()->SetId(0, i);

        triangle->GetPointIds()->SetId(1, i+1);

        triangle->GetPointIds()->SetId(2, i+2);

        triangles->InsertNextCell(triangle);
    }

    // Create a polydata object
    vtkSmartPointer<vtkPolyData> polyData =
        vtkSmartPointer<vtkPolyData>::New();

    // Add the geometry and topology to the polydata
    polyData->SetPoints(points);
    polyData->SetPolys(triangles);
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(polyData);

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(
            colors->GetColor3d("Cornsilk").GetData());
    scene.sphereActors.push_back(actor.GetPointer());
    backRenderer->AddActor(actor);
    backRenderer->GetRenderWindow()->Render();
}

// Calls from the queue of function pointers when the continue button is clicked
void Widget::on_continueButton_clicked() {
    string isFinal;
    if (gd->isSMF() && !hasAddedLayers) {
        auto actors = backRenderer->GetActors();
        actors->InitTraversal();
        while (actors->GetNumberOfItems())
            backRenderer->RemoveActor(actors->GetNextItem());
        addLayers();
        gd->updateVtkGrid(scene);
        backRenderer->GetActiveCamera()->SetPosition(0, 0, 500.0);
        backRenderer->ResetCamera();
        backRenderer->GetRenderWindow()->Render();
        hasAddedLayers = true;
    } else if (gd->funDeq.size()) {
        // TODO: Following line doesn't work, not sure why
        // ui->stepTextBrowser->setPlainText(
        // QString::fromStdString("Running..."));
        (gd->*gd->funDeq.front())();
        gd->funDeq.pop_front();
        cout << "============\n"
            << "Energies after '" << gd->stringDeq.front() << "'\n";
        gd->stringDeq.pop_front();
        if (!gd->stringDeq.size())
            isFinal = "Final ";
        else
            isFinal = "";
        cout << isFinal << "System Energy: " << gd->getSysEnergy() << endl;
        Vector2d disps = gd->getAvgSpringDisp();
        cout << isFinal << "Linear Error: " << disps[0] << endl;
        cout << isFinal << "Angular Error: " << disps[1] << endl;
        addLayers();
        gd->updateVtkGrid(scene);
        backRenderer->GetRenderWindow()->Render();
    }
    if (gd->stringDeq.size()) {
        ui->stepTextBrowser->setPlainText(
                QString::fromStdString(gd->stringDeq.front()));
    } else if (triangulate) {
        ui->stepTextBrowser->setPlainText(
                QString::fromStdString("Triangulate grid"));
        triangulate = false;
    } else if (!fullyDone) {
        triangulateGrid();
        gd->assessResults();
        ui->stepTextBrowser->setPlainText(
                QString::fromStdString("No remaining operations"));
        fullyDone = true;
    }
}

void Widget::on_printButton_clicked() {
    exporter->Update();
}

void Widget::on_centerButton_clicked() {
    backRenderer->ResetCamera();
    backRenderer->GetRenderWindow()->Render();
}
