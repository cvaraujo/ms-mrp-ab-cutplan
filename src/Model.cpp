//
// Created by Carlos on 06/03/19.
// Updated by Carlos on 19/08/22.
//

#include <chrono>
#include <cstdio>
#include "../headers/Model.h"

using namespace lemon;

class cyclecallback : public GRBCallback
{
public:
  double lastiter, lastnode;
  int numvars;
  vector<vector<GRBVar>> y;
  typedef ListGraph G;
  typedef G::Edge Edge;
  typedef G::EdgeIt EdgeIt;
  typedef G::Node Node;
  typedef G::EdgeMap<double> LengthMap;
  typedef G::NodeMap<bool> BoolNodeMap;

  Graph *graph;

  cyclecallback(Graph *xgraph, int xnumvars, vector<vector<GRBVar>> xy)
  {
    lastiter = lastnode = 0;
    numvars = xnumvars;
    y = xy;
    graph = xgraph;
  }

protected:
  void callback()
  {
    try
    {
      if (where == GRB_CB_MIPNODE)
      {
        double nodecnt = getDoubleInfo(GRB_CB_MIP_NODCNT);
        int mipStatus = getIntInfo(GRB_CB_MIPNODE_STATUS);
        lastnode++;
        if (mipStatus == GRB_OPTIMAL)
        {
          int i, n = graph->getN();
          double cutValue;
          G g;
          vector<Edge> setEdges;
          vector<Node> setNodes = vector<Node>(n);

          for (int i = 0; i < n; i++)
          {
            setNodes[i] = g.addNode();
          }
          LengthMap map(g);
          LengthMap length(g);

          for (int i = 0; i < n; i++)
          {
            for (auto arc : graph->arcs[i])
            {
              if (getNodeRel(y[i][arc->getD()]) > 0)
              {
                setEdges.push_back(g.addEdge(setNodes[i], setNodes[arc->getD()]));
                length[setEdges[setEdges.size() - 1]] = getNodeRel(y[i][arc->getD()]);
              }
            }
          }

          GomoryHu<G, LengthMap> gh(g, length);
          gh.run();
          BoolNodeMap bm(g);

          for (auto k : graph->terminals)
          {
            cutValue = gh.minCutMap(setNodes[graph->getRoot()], setNodes[k], bm);
            if (cutValue < 1)
            {
              GRBLinExpr expr;
              for (i = 0; i < n; i++)
              {
                for (auto *arc : graph->arcs[i])
                {
                  if (!bm[setNodes[arc->getD()]])
                  {
                    expr += y[i][arc->getD()];
                  }
                }
              }
              addCut(expr >= 1);
            }
          }
          setEdges.clear();
          setNodes.clear();
        }
      }
    }
    catch (GRBException e)
    {
      cout << "Error number: " << e.getErrorCode() << endl;
      cout << e.getMessage() << endl;
    }
    catch (...)
    {
      cout << "Error during callback" << endl;
    }
  }
};

Model::Model(Graph *graph)
{
  if (graph != nullptr)
  {
    this->graph = graph;
  }
  else
    exit(EXIT_FAILURE);
}

void Model::initialize()
{
  int o, d, n = graph->getN(), m = graph->getM();
  try
  {

    env.set("LogFile", "MS_mip.log");
    env.start();

    y = vector<vector<GRBVar>>(n, vector<GRBVar>(n));
    z = vector<GRBVar>(n);
    lambda = vector<GRBVar>(n);
    xi = vector<GRBVar>(n);

    char name[40];
    for (o = 0; o < n; o++)
    {
      for (auto *arc : graph->arcs[o])
      {
        d = arc->getD();
        sprintf(name, "y_%d_%d", o, d);
        y[o][d] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
      }
    }

    for (auto i : graph->terminals)
    {
      sprintf(name, "z_%d", i);
      z[i] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
    }

    for (int i = 0; i < n; i++)
    {
      sprintf(name, "lambda_%d", i);
      lambda[i] = model.addVar(0.0, graph->getParamDelay() + 1, 0.0, GRB_CONTINUOUS, name);
      sprintf(name, "xi_%d", i);
      xi[i] = model.addVar(0.0, graph->getParamJitter() + 1, 0.0, GRB_CONTINUOUS, name);
    }

    delta_min = model.addVar(1.0, graph->getParamDelay() - graph->getParamVariation(), 0.0, GRB_CONTINUOUS, "delta_min");

    model.update();
    cout << "Create variables" << endl;
  }
  catch (GRBException &ex)
  {
    cout << ex.getMessage() << endl;
    cout << ex.getErrorCode() << endl;
    exit(EXIT_FAILURE);
  }
}

void Model::initModel()
{
  cout << "Begin the model creation" << endl;
  objectiveFunction();
  allNodesAttended();
  calcLambdaXi();
  limDelayAndJitter();
  limVariation();
  cout << "All done!" << endl;
}

void Model::objectiveFunction()
{
  GRBLinExpr objective;
  for (auto k : graph->terminals)
    objective += z[k];
  model.setObjective(objective, GRB_MINIMIZE);
  model.update();
  cout << "Objective Function was added successfully!" << endl;
}

void Model::allNodesAttended()
{
  model.addConstr(y[graph->getRoot()][0] == 1);
  for (auto j : graph->DuS)
  {
    if (graph->removed[j])
    {
      model.addConstr(y[0][j] == 1, "in_arcs_" + to_string(j));
      model.addConstr(lambda[j] == graph->getParamDelay() + 1, "lambda_" + to_string(0) + "_" + to_string(j));
      model.addConstr(xi[j] == graph->getParamJitter() + 1, "xi_" + to_string(0) + "_" + to_string(j));
    }
    else
    {
      GRBLinExpr inArcs;
      for (int i = 0; i < graph->getN(); i++)
      {
        for (auto *arc : graph->arcs[i])
        {
          if (arc->getD() == j)
          {
            inArcs += y[i][j];
          }
        }
      }
      model.addConstr(inArcs == 1, "in_arcs_" + to_string(j));
    }
  }
  model.update();
  cout << "All nodes are inserted" << endl;
}

void Model::calcLambdaXi()
{
  int j, bigM;
  for (int i = 0; i < graph->getN(); i++)
  {
    for (auto *arc : graph->arcs[i])
    {
      j = arc->getD();
      bigM = graph->getParamDelay() + 1;
      model.addConstr(lambda[j] >= (lambda[i] + (arc->getDelay() * y[i][j])) - (bigM * (1 - y[i][j])), "lambda_" + to_string(i) + "_" + to_string(j));
      model.addConstr(lambda[j] <= (lambda[i] + arc->getDelay()) + (bigM * (1 - y[i][j])), "lambda_minus_" + to_string(i) + "_" + to_string(j));

      bigM = graph->getParamJitter() + 1;
      model.addConstr(xi[j] >= xi[i] + (arc->getJitter() * y[i][j]) - (bigM * (1 - y[i][j])), "xi_" + to_string(i) + "_" + to_string(j));
      // model.addConstr(xi[j] <= (xi[i] + arc->getJitter()) + (bigM * (1 - y[i][j])), "xi_minus_" + to_string(i) + "_" + to_string(j));
    }
  }
  model.update();
  cout << "Computing lambda and xi values" << endl;
}

void Model::limDelayAndJitter()
{
  for (auto k : graph->terminals)
  {
    model.addConstr(lambda[k] <= graph->getParamDelay() + (graph->getBigMDelay() - graph->getParamDelay()) * z[k], "delay_limit_" + to_string(k));
    model.addConstr(xi[k] <= graph->getParamJitter() + (graph->getBigMJitter() - graph->getParamJitter()) * z[k], "jitter_limit_" + to_string(k));
  }
  model.update();
  cout << "Delay and Jitter limits" << endl;
}

void Model::limVariation()
{
  model.addConstr(delta_min <= graph->getParamDelay() - graph->getParamVariation());
  model.addConstr(delta_min >= 1);

  for (auto k : graph->terminals)
  {
    model.addConstr(lambda[k] >= delta_min - (graph->getParamDelay()) * z[k], "min_delay_variation_" + to_string(k));
    model.addConstr(lambda[k] <= (delta_min + graph->getParamVariation()) + graph->getParamDelay() * z[k], "max_delay_variation_" + to_string(k));
  }

  // int o, d, bigMK, bigML;
  // for (auto k : graph->terminals)
  // {
  //   for (auto l : graph->terminals)
  //   {
  //     if (k != l)
  //     {
  //       bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
  //       bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
  //       model.addConstr(lambda[k] - lambda[l] <= graph->getParamVariation() + (bigMK * z[k]) + (bigML * z[l]), "limit_of_variation_between_pairs_" + to_string(k) + "_" + to_string(l));
  //     }
  //   }
  // }
  model.update();
  cout << "Delay variation limits" << endl;
}

void Model::solve(string timeLimit)
{
  try
  {
    model.set("TimeLimit", timeLimit);
    // model.set(GRB_DoubleParam_Heuristics, 0.0);
    cyclecallback cb = cyclecallback(graph, 0, y);
    model.update();
    model.setCallback(&cb);
    // model.computeIIS();
    // model.set("OutputFlag", "0");
    model.write("modelo.lp");
    model.optimize();
  }
  catch (GRBException &ex)
  {
    cout << ex.getMessage() << endl;
  }
}

bool Model::checkSolution()
{
  int i, n = graph->getN();
  int root = graph->getRoot();
  vector<int> pred = vector<int>(n);

  // Get the tree
  for (i = 0; i < n; i++)
    for (auto *arc : graph->arcs[i])
      if (y[i][arc->getD()].get(GRB_DoubleAttr_X) > 0.5)
        pred[arc->getD()] = i;

  // Check the path for each terminal
  for (auto k : graph->terminals)
  {
    double delay_path = 0, jitter_path = 0;
    int act_node = k;

    while (act_node != root)
    {
      delay_path += graph->getDelay(pred[act_node], act_node);
      jitter_path += graph->getJitter(pred[act_node], act_node);

      act_node = pred[act_node];
    }

    if ((delay_path > graph->getParamDelay() || jitter_path > graph->getParamJitter()) && z[k].get(GRB_DoubleAttr_X) <= 0)
    {
      cout << "ERROR" << endl;
      return false;
    }
  }

  return true;
}

void Model::writeSolution(string instance, int preprocessingTime)
{
  try
  {
    if (!checkSolution())
    {
      return;
    }
    ofstream output;
    output.open(instance, ofstream::app);
    output << "Prep. Time: " << preprocessingTime << endl;
    double ub = model.get(GRB_DoubleAttr_ObjVal), lb = model.get(GRB_DoubleAttr_ObjBound);
    output << "UB: " << ub << endl;
    output << "LB: " << lb << endl;
    if (ub != 0)
      output << "gap: " << (ub - lb) / ub << endl;

    output << "N. Nodes: " << model.get(GRB_DoubleAttr_NodeCount) << endl;
    output << "Runtime: " << model.get(GRB_DoubleAttr_Runtime) << endl;

    output << "----- Solution -----" << endl;
    for (int i = 0; i < graph->getN(); i++)
    {
      if (!graph->removed[i])
        for (auto *arc : graph->arcs[i])
        {
          if (y[i][arc->getD()].get(GRB_DoubleAttr_X) > 0.1)
          {
            output << i << " - " << arc->getD() << endl;
          }
        }
    }
    output << graph->getParamDelay() << ", " << graph->getParamJitter() << ", " << graph->getParamVariation() << ", " << delta_min.get(GRB_DoubleAttr_X) << endl;
    for (auto i : graph->terminals)
      output << i << " = " << lambda[i].get(GRB_DoubleAttr_X) << ", " << xi[i].get(GRB_DoubleAttr_X) << endl;
    output.close();
  }
  catch (GRBException &ex)
  {
    cout << ex.getMessage() << endl;
  }
}
