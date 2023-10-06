#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include <map>


using namespace std;

template<typename VertexT, typename WeightT>//vertex latitude and longitude and weightT is distance between two vectors
class graph {
private:
    map<VertexT, map<VertexT, WeightT>> ListofNodes;
    int edgeSize;



    //
    // We are using adjacency matrix implementation, where rows
    // are the starting vertex and cols are the ending vertex.
    // We keep track of the vertices in the Vertices vector,
    // where the vertex's position in the vector --- 0, 1, 2,
    // 3, 4, 5, ... --- denotes the row in the adjacency matrix
    // where their edges are found.  Example: if vertex "ORD" is
    // in position 1 of the Vertices vector, then row 1 of
    // AdjMatrix are the edges that start at "ORD" and lead to
    // other vertices.
    //

    //
    // _LookupVertex
    //
    // Finds the vertex in the Vertices vector and returns it's
    // index position if found, otherwise returns -1.
    //
    int _LookupVertex(VertexT v) const {
        auto it = ListofNodes.find(v);
        if (it != ListofNodes.end()) {
            return distance(ListofNodes.begin(), it);
        }
        else{
            return -1;
        }

    }

public:
    //
    // constructor:
    //
    // Constructs an empty graph where n is the max # of vertices
    // you expect the graph to contain.
    //
    // NOTE: the graph is implemented using an adjacency matrix.
    // If n exceeds the dimensions of this matrix, an exception
    // will be thrown to let you know that this implementation
    // will not suffice.
    //
    graph() {
        edgeSize=0;
    }


    //
    // NumVertices
    //
    // Returns the # of vertices currently in the graph.
    //
    int NumVertices() const {
        return ListofNodes.size();
    }

    //
    // NumEdges
    //
    // Returns the # of edges currently in the graph.
    //
    int NumEdges() const {
        return edgeSize;
    }

    //
    // addVertex
    //
    // Adds the vertex v to the graph if there's room, and if so
    // returns true.  If the graph is full, or the vertex already
    // exists in the graph, then false is returned.
    //
    bool addVertex(VertexT v) {
        // Check if the vertex already exists
        if (ListofNodes.find(v) == ListofNodes.end()) {
            // If it doesn't exist, add a new entry with an empty map as its value
            ListofNodes[v] = map<VertexT, WeightT>();
            // Return true to indicate that the vertex was added
            return true;
        }
        else{
            // If the vertex already exists, return false
            return false;
        }



    }

    //
    // addEdge
    //
    // Adds the edge (from, to, weight) to the graph, and returns
    // true.  If the vertices do not exist or for some reason the
    // graph is full, false is returned.
    //
    // NOTE: if the edge already exists, the existing edge weight
    // is overwritten with the new edge weight.
    //
    bool addEdge(VertexT from, VertexT to, WeightT weight) {
        // Check if the vertices exist
        if (ListofNodes.find(from) == ListofNodes.end() || ListofNodes.find(to) == ListofNodes.end()) {
            return false;
        }

        // Check if the edge already exists
        auto &edges = ListofNodes[from];
        auto edge_it = edges.find(to);

        if (edge_it == edges.end()) {
            // Edge doesn't exist, add the edge
            ListofNodes[from][to] = weight;
            edgeSize += 1;
            return true;
        } else {
            ListofNodes[from][to] = weight;
            return true;
        }
    }

    //
    // getWeight
    //
    // Returns the weight associated with a given edge.  If
    // the edge exists, the weight is returned via the reference
    // parameter and true is returned.  If the edge does not
    // exist, the weight parameter is unchanged and false is
    // returned.
    //
    bool getWeight(VertexT from, VertexT to, WeightT &weight) const {
        // Check if the vertices exist
        auto from_it = ListofNodes.find(from);
        auto to_it = ListofNodes.find(to);
        if (from_it == ListofNodes.end()) {
            return false;
        }
        if(to_it == ListofNodes.end()){
            return false;
        }

        // Search for the edge in the adjacency list of 'from' vertex
        auto edge_it = from_it->second.find(to);
        if (edge_it != from_it->second.end()) {
            weight = edge_it->second;
            return true;
        }
        else{
            return false;
        }


    }

    //
    // neighbors
    //
    // Returns a set containing the neighbors of v, i.e. all
    // vertices that can be reached from v along one edge.
    // Since a set is returned, the neighbors are returned in
    // sorted order; use foreach to iterate through the set.
    //
    set<VertexT> neighbors(VertexT v) const {
        // Find the vertex in the adjacency list
        auto it = ListofNodes.find(v);
        set<VertexT> result;
        if (it != ListofNodes.end()) {
            // Iterate through the neighbors and add them to the result set
            for (const auto& neighbor : it->second) {
                result.insert(neighbor.first);
            }
        }

        return result;
    }

    //
    // getVertices
    //
    // Returns a vector containing all the vertices currently in
    // the list.
    //
    vector<VertexT> getVertices() const {
        vector<VertexT> result;
        result.reserve(ListofNodes.size());
        for (const auto& entry : ListofNodes) {
            result.push_back(entry.first);
        }

        return result;
    }

    //
    // dump
    //
    // Dumps the internal state of the graph for debugging purposes.
    //
    // Example:
    //    graph<string,int>  G(26);
    //    ...
    //    G.dump(cout);  // dump to console
    //
    void dump(ostream &output) const {
        output << "***************************************************" << endl;
        output << "********************* LIST ***********************" << endl;
        for (const auto &vertexPair : ListofNodes) {
            const VertexT &vertex = vertexPair.first;
            const map<VertexT, WeightT> &edges = vertexPair.second;

            output << "Vertex: " << vertex << "\n";
            output << "Connected to:\n";

            for (const auto &edgePair : edges) {
                const VertexT &connectedVertex = edgePair.first;
                const WeightT &weight = edgePair.second;

                output << "  - " << connectedVertex << " with weight " << weight << "\n";
            }
            output << endl;
        }
        output << "**************************************************" << endl;
    }
};
