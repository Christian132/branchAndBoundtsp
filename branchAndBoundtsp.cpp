// Auther: Christian Markmoeller
// Course: CECS 570
// Project title: Traveling sales man problem, using branch and bound
// File: branchAndBoundtsp.cpp

// This program finds the shortest round trip given an undirected graph. It uses the branch
// and bound method to solve this problem. It prunes solutions when they are out of the scope
// of being a better soltion (a better solution in this case is a shorter round trip to all nodes or simply premature cycles)
// The program always assumes the start node is 0
// The graph in this program is represented by a matrix
#include "branchAndBoundtsp.h"

// Function adds an edge in both direction implying UNDIRECTED graph
std::vector<std::vector<int>> add_edge(std::vector<std::vector<int> > graph, int startVertex, int endVertex, int edgeWeight) {
    graph[startVertex][endVertex] = edgeWeight;
    graph[endVertex][startVertex] = edgeWeight;
    return graph;
}

// Delete an edge from being used
std::vector<std::vector<int>> del_edge(std::vector<std::vector<int> > graph, int startVertex, int endVertex) {
    graph[startVertex][endVertex] = 0;
    graph[endVertex][startVertex] = 0;
    return graph;
}

// Delete the remaining edges that cannot be used (if a vertex have two edges in use already)
std::vector<std::vector<int>> del_leftover_edges(std::vector<std::vector<int> > graph, int startVertex, int endVertex1, int endVertex2) {
    for (int i = 0; i<graph[startVertex].size(); i++){
        if(i != endVertex1 && i != endVertex2){
            graph[startVertex][i] = 0;
            graph[i][startVertex] = 0;
        }
    }
    return graph;
}

// Function assumes at least two edges for a vertex. Returns the lower bound for a given graph. It is the initial cost to take the route
// Number is found by taking the two lowest edges from each vertex, summing them up and divide it my two
double find_lower_bound(std::vector<std::vector<int>> graph){
    int firstLowestEdge = std::numeric_limits<int>::max(); // Max value of int, to give an upper bound
    int secondLowestEdge = std::numeric_limits<int>::max(); // Max value of int, to give an upper bound
    int numOfVert = graph.size();
    // Set excluded edges to 0
    int lowerBoundCost = 0;
    bool isSet = false;  // Value to check wheather the value has been set in "firstLowestEdge" varible
    for(int i = 0; i < numOfVert; i++){
        for(int j = 0; j<numOfVert; j++){
            if (graph[i][j] <= firstLowestEdge && graph[i][j] != 0){
                if (firstLowestEdge < secondLowestEdge){
                    secondLowestEdge = firstLowestEdge;
                }
                firstLowestEdge = graph[i][j];
                isSet = true;
            }
            if (graph[i][j] < secondLowestEdge && graph[i][j] != 0 && !isSet){
                secondLowestEdge = graph[i][j];
            }
            isSet = false;
        }
        // There is not two edges leaving the node
        if(firstLowestEdge == std::numeric_limits<int>::max() || secondLowestEdge == std::numeric_limits<int>::max()){
            lowerBoundCost = 0;
            return lowerBoundCost;
        }
        else{
          lowerBoundCost = lowerBoundCost + firstLowestEdge + secondLowestEdge;
        }
        firstLowestEdge = std::numeric_limits<int>::max(); // Reset the value
        secondLowestEdge = std::numeric_limits<int>::max(); // Reset the value
    }
    return lowerBoundCost/2.0;
}

// Print the graph
void print_graph(std::vector<std::vector<int>> graph){
    for(int i = 0; i<graph.size(); i++){
        for(int j = 0; j<graph.size();j++){
            std::cout << graph[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

// Visit an edge - meaning setting the visit vector to visit the nodes (assuming undirected graph)
std::vector<std::vector<int>> visit_edge(std::vector<std::vector<int>> edgeVisitList, int startVertex, int endVertex){
    edgeVisitList[startVertex].push_back(endVertex);
    edgeVisitList[endVertex].push_back(startVertex);
    return edgeVisitList;
}

// Check for a finished path - meaning if EVERY vertex has to visited edges (condition to be done)
bool check_for_finished_path(std::vector<std::vector<int>> usedEdges){
    std::vector<bool> vertVisitList(usedEdges.size(),false);
    int start = 0;
    int next = 0;
    int prev = -1; // make sure to not include a valid edge
    for(int i = 0; i < usedEdges.size(); i++){
        if (usedEdges[i].size() != 2){
            return false;
        }
    }
    return true;
}

// Check if there is a premature circle in the graph. Returns true for a premature circle
bool check_for_circle(std::vector<std::vector<int>> usedEdges){
  std::vector<bool> vertVisitList(usedEdges.size(),false);
  int start = 0;
  int next = 0;
  int prev = -1; // make sure to not include a valid edge
  // Check for circle
  for(int i = 0; i < vertVisitList.size(); i++){
    // There is no circle if there is no edge leaving
    if(usedEdges[start].size()>1){
      next = usedEdges[start][0];
      if (next == prev){
        next = usedEdges[start][1];
      }
      // There is a circle!
      if (vertVisitList[next] == true){
          return true;
      }
      else{
          vertVisitList[next] = true;
      }
      prev = start;
      start = next;
    }
    else{
      break;
    }
  }
  return false;
}

// Print the path through the graph that takes the least costly route
void print_finished_path(std::vector<std::vector<int>> finishedPath){
  // Starting in vertex "0"
  int start = 0;
  int next = 0;
  int prev = 0;
  std::cout << start << " - ";
  for(int i = 0; i < finishedPath.size(); i++){
      next = finishedPath[start][0];
      if (next == prev){
        next = finishedPath[start][1];
      }
      std::cout << next;
      if (i != finishedPath.size() - 1){
        std::cout << " - ";
      }
      prev = start;
      start = next;
  }
  std::cout << std::endl;
}

// Return the final cost of the finished route through
int final_path_cost(std::vector<std::vector<int>> finishedPath, std::vector<std::vector<int>> graph){
    // Starting in vertex "0"
  int start = 0;
  int next = 0;
  int prev = 0;
  int cost = 0;
  for(int i = 0; i < finishedPath.size(); i++){
      next = finishedPath[start][0];
      if (next == prev)
          next = finishedPath[start][1];
      cost = cost + graph[start][next];
      prev = start;
      start = next;
  }
  return cost;
}

// Initialize the priority queue to exclude and include the first found edge. Returns a prioroty queue
std::priority_queue<pathDescription, std::vector<pathDescription>, CompPathDes> initPrioQueue(std::vector<std::vector<int>> graph){
  int numOfVert = graph.size();
  int initVert = 0;
  // Finding the initial destination edge from vertex "0"
  for (int i = 0; i < numOfVert; i++){
      if (graph[0][i] > 0){
          initVert = i;
          break;
      }
  }
  std::priority_queue<pathDescription, std::vector<pathDescription>, CompPathDes> prioQueue;
  pathDescription tempPathDesc; // Temporary variable to hold the top element in the queue
  tempPathDesc.usedEdges.resize(numOfVert); // Initialize usedEdges variable

  // Do not include the edge initialization
  std::vector<std::vector<int>> initGraphEx = graph;
  initGraphEx = del_edge(initGraphEx, 0,initVert);
  tempPathDesc.cost = find_lower_bound(initGraphEx);
  tempPathDesc.graph = initGraphEx;
  prioQueue.push(tempPathDesc);

  // Include the edge initialization
  std::vector<std::vector<int>> initGraphVi = graph;
  tempPathDesc.cost = find_lower_bound(initGraphVi);
  tempPathDesc.graph = graph;
  tempPathDesc.usedEdges = visit_edge(tempPathDesc.usedEdges,0,initVert);
  prioQueue.push(tempPathDesc);
  return prioQueue;
}

// This function will try to find the shortest path in a graph given. It will use the branch-and-bound method,
// meaning it will ignore possible solutions that are not a good solution
void find_shortest_path(std::vector<std::vector<int>> graph){
    auto start = std::chrono::system_clock::now();// Timer to see how long all execution took
    int numOfVert = graph.size();
    pathDescription tempPathDesc; // Temporary variable to hold the top element in the queue
    std::priority_queue<pathDescription, std::vector<pathDescription>, CompPathDes> prioQueue;
    prioQueue = initPrioQueue(graph);
    bool foundEdge; // boolean to check wheather an edge has been considered
    double shortestPath = std::numeric_limits<double>::max(); // Shortest path is initialized as "infinity"
    std::vector<std::vector<int>> finishedPath; // Final cycle through the graph
    int numberOfZeroes = 0; // Variable to check number of zeroes in a row in the matrix
    int numberOfChecks = 0; // How many solutions has been considered

    // While loop to find the shortest path. Keeps looping as long as it has a potential in the queue
    while(!prioQueue.empty()){
        tempPathDesc = prioQueue.top();
        prioQueue.pop();
        // Only loop if the cost is less than an already found solution
        if (tempPathDesc.cost < shortestPath && tempPathDesc.cost > 0){
            numberOfChecks += 1;
            foundEdge = false; // Used to break the for loops. "Found" an edge to check
            pathDescription tempPathDescVis;
            pathDescription tempPathDescEx;
            for(int i = 0; i<numOfVert; i++){
                for(int j = 0; j<numOfVert; j++){
                  // Rules to prune a solution. Consider only if the node has not been visited by the start node and if the node has less then two edges to it and if there is no premature circle
                    if (tempPathDesc.graph[i][j] > 0 && !std::binary_search(tempPathDesc.usedEdges[i].begin(), tempPathDesc.usedEdges[i].end(), j) && tempPathDesc.usedEdges[i].size() <= 2 && !check_for_circle(tempPathDesc.usedEdges)){
                        tempPathDescVis = tempPathDesc; // Temp to calculate the Visited edge
                        tempPathDescEx = tempPathDesc; // Temp to calculate the Excluded edge

                        // Exclude the edge in question from the excluded set
                        tempPathDescEx.graph = del_edge(tempPathDesc.graph, i,j);
                        // Check to see if only two edges remains from in graph - if so visit them
                        numberOfZeroes = std::count(tempPathDescEx.graph[i].begin(),tempPathDescEx.graph[i].end(), 0);
                        if (tempPathDescEx.graph.size() - numberOfZeroes == 2 && tempPathDescEx.usedEdges[i].size() < 2){
                            for(int k = 0;k < tempPathDesc.graph.size();k++){
                                if (tempPathDescEx.graph[i][k] != 0 && i != k){
                                    if (tempPathDescEx.usedEdges[i].size() == 0){
                                        tempPathDescEx.usedEdges = visit_edge(tempPathDescEx.usedEdges,i,k);
                                    }
                                    else if (tempPathDescEx.usedEdges[i].size() == 1 && tempPathDescEx.usedEdges[i][0] != k){
                                        tempPathDescEx.usedEdges = visit_edge(tempPathDescEx.usedEdges,i,k);
                                    }
                                }
                            }
                        }
                        tempPathDescEx.cost = find_lower_bound(tempPathDescEx.graph);
                        // Only push new potential solution if it's temporary solution is less than a known solution
                        if(tempPathDescEx.cost <= shortestPath){
                            prioQueue.push(tempPathDescEx);
                        }
                        // Visiting the edge in question
                        tempPathDescVis.usedEdges = visit_edge(tempPathDescVis.usedEdges,i,j);
                        // Exclude rest of edges if we have two edges or more
                        if (tempPathDescVis.usedEdges[i].size() == 2){
                            tempPathDescVis.graph = del_leftover_edges(tempPathDescVis.graph,i,tempPathDescVis.usedEdges[i][0], tempPathDescVis.usedEdges[i][1]);
                        }
                        // Trying to add too many edges - do nothing and discard solution
                        else if (tempPathDescVis.usedEdges[i].size() > 2 ){
                            foundEdge = true;
                            break;
                        }
                        // Include that edge
                        tempPathDescVis.cost = find_lower_bound(tempPathDescVis.graph);
                        prioQueue.push(tempPathDescVis); // Include the visited edge

                        // Check if a path is finished and if it is less than already known solutions for the visited
                        if (check_for_finished_path(tempPathDescVis.usedEdges)){
                            if(tempPathDescVis.cost < shortestPath ){
                                shortestPath = final_path_cost(tempPathDescVis.usedEdges, graph);
                                finishedPath = tempPathDescVis.usedEdges;
                            }
                        }
                        // Check if a path is finished and if it is less than already known solutions for the excluded
                        if (check_for_finished_path(tempPathDescEx.usedEdges)){
                            if(tempPathDescEx.cost < shortestPath ){
                                shortestPath = final_path_cost(tempPathDescEx.usedEdges, graph);
                                finishedPath = tempPathDescEx.usedEdges;
                            }
                        }
                        foundEdge = true;
                        break;
                    }
                }
                // If an edge has been considered, find a new edge
                if(foundEdge){
                    break;
                }
            }
        }
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end-start;
    print_graph(graph);
    std::cout << "\n\n\n\nFinal solution after time :"  << diff.count() << std::endl;
    std::cout << "The shortest path is: " << shortestPath << std::endl;
    std::cout << "Did " << numberOfChecks << " checks in the tree to find the optimal path!" << std::endl;
    print_finished_path(finishedPath);
}

int main(){
  std::vector<std::vector<int>> bigGraph;
  int numVertBigGraph = 15; // Number of verticies in a graph
  bigGraph.resize(numVertBigGraph);
  for(int i = 0; i<numVertBigGraph;i++){
      bigGraph[i].resize(numVertBigGraph);
  }
  int randVal = 0;
  // Create the graph with edges of random cost with a given bound
  for(int i = 0; i<numVertBigGraph;i++){
      for(int j = 0; j<numVertBigGraph;j++){
          randVal = rand()%15; // Mod with the max cost of an edge
          if (i == j)
              bigGraph[i][j] = 0;
          else
              bigGraph = add_edge(bigGraph, i, j, randVal);
      }
  }
  find_shortest_path(bigGraph);
  return 0;
}
