// Auther: Christian Markmoeller
// Course: CECS 570
// Project title: Traveling sales man problem, using branch and bound
// File: branchAndBoundtsp.h

// This program finds the shortest round trip given an undirected graph. It uses the branch
// and bound method to solve this problem. It prunes solutions when they are out of the scope
// of being a better soltion (a better solution in this case is a shorter round trip to all nodes)

#ifndef BRANCHANDBOUNDTSP_H
#define BRANCHANDBOUNDTSP_H

#include <queue>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <stdlib.h>
#include <chrono>



// Struct to hold the current graph, visited edges and the cost
struct pathDescription {
    double cost;
    std::vector<std::vector<int>> graph;
    std::vector<std::vector<int>> usedEdges;
};

// A function to help the priority queue find the lowest element in the pathDescription
struct CompPathDes{
    bool operator()(const pathDescription& a, const pathDescription& b){
        return a.cost > b.cost;
    }
};

// Function to help compare values in the pathDesription struct for the prio queue
bool pathDesriptionCompare(pathDescription compVal1, pathDescription compVal2);
// Function adds an edge in both direction implying UNDIRECTED graph
std::vector<std::vector<int>> add_edge(std::vector<std::vector<int> > graph, int startVertex, int endVertex, int edgeWeight);

// Delete an edge from being used
std::vector<std::vector<int>> del_edge(std::vector<std::vector<int> > graph, int startVertex, int endVertex);
// Delete the remaining edges that cannot be used (if a vertex have two edges in use already)
std::vector<std::vector<int>> del_leftover_edges(std::vector<std::vector<int> > graph, int startVertex, int endVertex1, int endVertex2);
// Function assumes at least two edges for a vertex. Returns the lower bound for a given graph
double find_lower_bound(std::vector<std::vector<int>> graph);

// Print the given graph
void print_graph(std::vector<std::vector<int>> graph);

// Visit an edge - meaning setting the visiting vector to visit the node (assuming undirected graph)
std::vector<std::vector<int>> visit_edge(std::vector<std::vector<int>> edgeVisitList, int startVertex, int endVertex);

// Check for a finished path - meaning if EVERY vertex has to visited edges (condition to be done)
bool check_for_finished_path(std::vector<std::vector<int>> usedEdges);


void print_finished_path(std::vector<std::vector<int>> finishedPath);

int final_path_cost(std::vector<std::vector<int>> finishedPath, std::vector<std::vector<int>> graph);

void find_shortest_path(std::vector<std::vector<int>> graph);

// Check if there is a circle in the visited set
bool check_for_circle(std::vector<std::vector<int>> usedEdges);


#endif
