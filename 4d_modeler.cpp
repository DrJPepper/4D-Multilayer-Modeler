#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"
#include "graddesc.h"
#include "util.h"
#include "widget.h"

GradDesc *gd;

int main(int argc, char *argv[]) {
    doctest::Context ctx;
    ctx.setOption("abort-after", 5);  // default - stop after 5 failed asserts
    ctx.applyCommandLine(argc, argv); // apply command line - argc / argv
    ctx.setOption("no-breaks", true); // override - don't break in the debugger
    int res = ctx.run();              // run test cases unless with --no-run
    if(ctx.shouldExit())              // query flags (and --exit) rely on this
        return res;

    string inputFileName = "input.txt";
    // Parse arguments
    if (argc > 1) {
        if (strcmp("-f", argv[1]) != 0) {
            cout << "ERROR: invalid argument, expected only '-f filename'\n";
            return 1;
        }
        inputFileName = argv[2];
    }

    ifstream inFile(inputFileName);
    if (!inFile.good()) {
        cout << "ERROR: input file invalid or doesn't exist\n";
        return 1;
    }

    // Needed to ensure appropriate OpenGL context is created for VTK rendering.
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());

    QApplication a(argc, argv);
    gd = new GradDesc();
    gd->populateMatrices(&inFile);
    Widget *w = new Widget(0, gd);
    w->show();
    return a.exec();
}
