//
// Created by carlos on 06/03/19.
//

#ifndef MRP_MODEL_H
#define MRP_MODEL_H

#include <iostream>
#include <vector>
#include "string"
#include <iomanip>
#include <bits/ios_base.h>
#include <algorithm>
#include <fstream>
#include <gurobi_c++.h>
#include "Graph.h"
#include <boost/algorithm/string.hpp>

using namespace std;

class Model
{

public:
  Graph *graph;
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  vector<vector<GRBVar>> y;
  vector<GRBVar> z, lambda, xi;
  GRBVar delta_min;

  void objectiveFunction();

  void allNodesAttended();

  void calcLambdaXi();

  void limDelayAndJitter();

  void limVariation();

  Model(Graph *graph);

  void initialize();

  void initModel();

  void solve(string timeLimit);

  bool checkSolution();

  void writeSolution(string instance, int preprocessingTime);
};

#endif // MRP_MODEL_H
