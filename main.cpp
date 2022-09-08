#include <sys/stat.h>
#include <string>
#include "headers/Graph.h"
#include "headers/Model.h"
#include <chrono>

int main(int argc, const char *argv[])
{
    if (argc < 4)
    {
        cout << "./MaxService graph.txt param.txt result.txt" << endl;
        return 0;
    }
    else
    {
        auto *graph = new Graph(argv[1], argv[2], argv[3]);
        // graph->showGraph();
        // getchar();
        graph->MVE(argv[3], "prep.txt");
        graph->finishPreprocessing(argv[3], true, false);
        auto *model = new Model(graph);
        model->initialize();
        model->initModel();
        model->solve("3600");
        model->writeSolution(argv[3], 0);
    }
    return 0;
}
