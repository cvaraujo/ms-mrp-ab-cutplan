//
// Created by carlos on 06/03/19.
//

#include <chrono>
#include "../headers/Model.h"

class cyclecallback: public GRBCallback {
public:
  double lastiter, lastnode;
  int numvars;
  vector<vector<GRBVar>> y;
  
  Graph *graph;

  cyclecallback(Graph *xgraph, int xnumvars, vector<vector<GRBVar>> xy){
    lastiter = lastnode = -1;
    numvars = xnumvars;
    y = xy;
    graph = xgraph;
  }

protected:
  void callback() {
    try {
      if (where == GRB_CB_MIPNODE) {
	cout << "*** New node ***" << endl;
	double nodecnt = getDoubleInfo(GRB_CB_MIP_NODCNT);
	cout << "Total of nodes: " << nodecnt << endl;
	
	int n = graph->getN();
        	
	for (int i = 0; i < graph->getN(); i++) {
	  for (auto arc : graph->arcs[i]) {
	    cout << "y[" << i << "][" << arc->getD() << "]" << " = " << getNodeRel(y[i][arc->getD()]) << endl;
	  }
	}
	cout << "Acho que deu Certo" << endl;
	getchar();
      }
    } catch(GRBException e) {
      cout << "Error number: " << e.getErrorCode() << endl;
      cout << e.getMessage() << endl;
    } catch (...) {
      cout << "Error during callback" << endl;
    }
  }
  
};

Model::Model(Graph *graph) {
  if (graph != nullptr) {
    this->graph = graph;
  } else exit(EXIT_FAILURE);
}

void Model::initialize() {
  int o, d, n = graph->getN(), m = graph->getM();
  try {

    env.set("LogFile", "MS_mip.log");
    env.start();
    
    y = vector<vector<GRBVar>>(n, vector<GRBVar>(n));
    z = vector<GRBVar>(n);
    lambda = vector<GRBVar>(n);
    xi = vector<GRBVar>(n);
    
    char name[30];
    for (o = 0; o < n; o++) {
      for (auto *arc : graph->arcs[o]) {
	d = arc->getD();
	sprintf(name, "y_%d_%d", o, d);
	y[o][d] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
      }
    }
    
    for (auto i : graph->terminals) {
      sprintf(name, "z_%d", i);
      z[i] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
    }
    
    for (int i = 0; i < n; i++) {
      sprintf(name, "lambda_%d", i);
      lambda[i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
      sprintf(name, "xi_%d", i);    
      xi[i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
    }
    
    model.update();
    cout << "Create variables" << endl;
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
    cout << ex.getErrorCode() << endl;
    exit(EXIT_FAILURE);
  }
}

void Model::initModel() {
  cout << "Begin the model creation" << endl;
  objectiveFunction();
  allNodesAttended(); calcLambdaXi();
  limDelayAndJitter(); limVariation();
  cout << "All done!" << endl;
}

void Model::objectiveFunction() {
  GRBLinExpr objective;
  for (auto k : graph->terminals) objective += z[k];
  model.setObjective(objective, GRB_MINIMIZE);
  model.update();
  cout << "Objective Function was added successfully!" << endl;
}

void Model::allNodesAttended() {
  model.addConstr(y[graph->getRoot()][0] == 1);
  for (auto j : graph->DuS) {
    if (graph->removed[j]) {
      model.addConstr(y[0][j] == 1, "in_arcs_" + to_string(j));
      model.addConstr(lambda[j] == graph->getParamDelay() + 1, "lambda_" + to_string(0) + "_" + to_string(j));
      model.addConstr(xi[j] == graph->getParamJitter() + 1, "xi_" + to_string(0) + "_" + to_string(j));
    } else {
      GRBLinExpr inArcs;
      for (int i = 0; i < graph->getN(); i++) {
	for (auto *arc : graph->arcs[i]) {
	  if (arc->getD() == j) {
	    inArcs += y[i][j];
	  }
	}
      }
      model.addConstr(inArcs == 1, "in_arcs_" + to_string(j));   
    }
  }
  //     for (int j = 0; i < graph->getN(); j++) {
  //         if (graph->removed[j]) {
  //             continue;
  //             cout << "Aqui 1" << endl;
  //             
  //             cout << "Aqui 2" << endl;
  //         } else if (j != graph->getRoot()) {
  //             GRBLinExpr inArcs;
  //             for (auto *arc : graph->arcs[j]) {
  //                 inArcs += y[arc->getD()][j];
  //             }
  //             inArcs += y[0][j];
  //             model.addConstr(inArcs == 1, "in_arcs_" + to_string(j));
  //         }
  //     }
  model.update();
  cout << "All nodes are inserted" << endl;
}

void Model::calcLambdaXi() {
  int j, bigM;
  for (int i = 0; i < graph->getN(); i++) {
      for (auto *arc : graph->arcs[i]) {
	j = arc->getD();
	bigM = graph->getBigMDelay();
            
	model.addConstr(lambda[j] >= (lambda[i] + (arc->getDelay() * y[i][j])) - (bigM * (1 - y[i][j])), "lambda_" + to_string(i) + "_" + to_string(j));
	model.addConstr(lambda[j] <= (lambda[i] + arc->getDelay()) + (bigM * (1 - y[i][j])), "lambda_minus_" + to_string(i) + "_" + to_string(j));
        
	bigM = graph->getBigMJitter();
	model.addConstr(xi[j] >= xi[i] + (arc->getJitter() * y[i][j]) - (bigM * (1 - y[i][j])), "xi_" + to_string(i) + "_" + to_string(j));
	model.addConstr(xi[j] <= (xi[i] + arc->getJitter()) + (bigM * (1 - y[i][j])), "xi_minus_" + to_string(i) + "_" + to_string(j));
      }
  }
  model.update();
  cout << "Computing lambda and xi values" << endl;
}

void Model::limDelayAndJitter() {
  for (auto k : graph->terminals) {
    model.addConstr(lambda[k] <= graph->getParamDelay() + (graph->getBigMDelay() - graph->getParamDelay()) * z[k], "delay_limit_" + to_string(k));
        
    model.addConstr(xi[k] <= graph->getParamJitter() + (graph->getBigMJitter() - graph->getParamJitter()) * z[k], "jitter_limit_" + to_string(k));
        
  }
  model.update();
  cout << "Delay and Jitter limits" << endl;
}

void Model::limVariation() {
  int o, d, bigMK, bigML;
  for (auto k : graph->terminals) {
    for (auto l : graph->terminals) {
      if (k != l) {
	bigMK = graph->getBigMDelay() - min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
	bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
          
	model.addConstr(lambda[k] - lambda[l] <= graph->getParamVariation() + (bigMK * z[k]) + (bigML * z[l]), "limit_of_variation_between_pairs_" + to_string(k) + "_" + to_string(l));
      }
    }
  }
  model.update();
  cout << "Delay variation limits" << endl;
}

void Model::solve(string timeLimit) {
  try {
    model.set("TimeLimit", timeLimit);
    model.set(GRB_DoubleParam_Heuristics, 0.0);
    cyclecallback cb = cyclecallback(graph, 0, y);
    model.update();
    model.setCallback(&cb);
    //model.computeIIS();
    //model.set("OutputFlag", "0");
    model.write("modelo.lp");
    model.optimize();
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
  }
}

void Model::writeSolution(string instance, int preprocessingTime) {
  try {
    ofstream output;
    output.open(instance, ofstream::app);
    output << "Prep. Time: " << preprocessingTime << endl;
    double ub = model.get(GRB_DoubleAttr_ObjVal), lb = model.get(GRB_DoubleAttr_ObjBound);
    output << "UB: " << ub << endl;
    output << "LB: " << lb << endl;
    if (ub != 0) output << "gap: " << (ub - lb) / ub << endl;
    
        
    output << "N. Nodes: " << model.get(GRB_DoubleAttr_NodeCount) << endl;
    output << "Runtime: " << model.get(GRB_DoubleAttr_Runtime) << endl;

    output << "----- Solution -----" << endl;
    for (int i = 0; i < graph->getN(); i++) {
      if (!graph->removed[i])
	for (auto *arc : graph->arcs[i]) {
	  if (y[i][arc->getD()].get(GRB_DoubleAttr_X) > 0.1) {
	    output << i << " - " << arc->getD() << endl;
	  }
	}
    }
    output << graph->getParamDelay() << ", " << graph->getParamJitter() << ", " << graph->getParamVariation() << endl;
    for (auto i : graph->terminals)
      output << i << " = " << lambda[i].get(GRB_DoubleAttr_X) << ", " << xi[i].get(GRB_DoubleAttr_X)<< endl;
    output.close();
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
  }

}
