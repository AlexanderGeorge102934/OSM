#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <queue>
#include <functional>
#include <limits>

#include "tinyxml2.h"
#include "dist.h"
#include "graph.h"
#include "osm.h"


using namespace std;
using namespace tinyxml2;
const double INF = numeric_limits<double>::max();

BuildingInfo searchBuilding(const vector<BuildingInfo> &Buildings, const string &buildingName) {
    for (const auto &building: Buildings) {
        if (building.Fullname.find(buildingName) != string::npos ||
            building.Abbrev.find(buildingName) != string::npos) {
            return building;
        }
    }
    return BuildingInfo(); // Return an empty building if not found
}


long long findNearestFootwayNode(const Coordinates &point, const vector<FootwayInfo> &Footways,
                                 const map<long long, Coordinates> &Nodes) {
    long long nearestNode = -1;
    double minDist = INF;

    for (const auto &footway: Footways) {
        for (const auto &node: footway.Nodes) {
            double distance = distBetween2Points(point.Lat, point.Lon, Nodes.at(node).Lat, Nodes.at(node).Lon);
            if (distance < minDist) {
                minDist = distance;
                nearestNode = node;
            }
        }
    }

    return nearestNode;
}


void dijkstraShortestPath(
        graph<long long, double>& G,
        long long startV,
        map<long long, pair<double, long long>>& result) {

    vector<long long> vertices = G.getVertices();

    for (const auto& v : vertices) {
        result[v].first = INF;
        result[v].second = 0;
    }
    result[startV].first = 0;

    // A priority queue with pairs in the form (distance, vertex)
    priority_queue<pair<double, long long>,
            vector<pair<double, long long>>,
            greater<pair<double, long long>>> unvisitedQueue;

    set<long long> visited; // Keep track of visited nodes

    for (const auto& v : vertices) {
        unvisitedQueue.push({result[v].first, v});
    }

    while (!unvisitedQueue.empty()) {
        long long currentV = unvisitedQueue.top().second;
        unvisitedQueue.pop();

        if (visited.count(currentV)) {
            continue; // Skip if the node has already been visited
        }
        visited.insert(currentV);

        for (const auto& adjV : G.neighbors(currentV)) {
            double edgeWeight;
            G.getWeight(currentV, adjV, edgeWeight);
            double alternativePathDistance = result[currentV].first + edgeWeight;

            if (alternativePathDistance < result[adjV].first) {
                result[adjV].first = alternativePathDistance;
                result[adjV].second = currentV;
                unvisitedQueue.push({result[adjV].first, adjV}); // Update the queue with the new distance
            }
        }
    }
}


string buildPath(const map<long long, pair<double, long long>> &result, long long startNode, long long endNode) {//FIXME
    if (result.at(endNode).first == INF) {
        return "No path found.";
    }

    string path;
    long long currentNode = endNode;
    while (currentNode != startNode) {
        path = "->" + to_string(currentNode) + path;
        currentNode = result.at(currentNode).second;
    }
    path = to_string(currentNode) + path;
    return path;
}

void printInformation(const map<long long, Coordinates> &Nodes,
                      const BuildingInfo &p1_building,
                      const BuildingInfo &p2_building,
                      const BuildingInfo &destination_building,
                      const Coordinates &p1_coordinates,
                      const Coordinates &p2_coordinates,
                      long long p1_nearest_node,
                      long long p2_nearest_node,
                      long long destination_nearest_node,
                      double p1_distance,
                      double p2_distance,
                      const string &p1_path,
                      const string &p2_path) {

    cout << "Person 1's point:" << endl;
    cout << p1_building.Fullname << endl;
    cout << "(" << p1_coordinates.Lat << ", " << p1_coordinates.Lon << ")" << endl;

    cout << "Person 2's point:" << endl;
    cout << p2_building.Fullname << endl;
    cout << "(" << p2_coordinates.Lat << ", " << p2_coordinates.Lon << ")" << endl;

    cout << "Destination Building:" << endl;
    cout << destination_building.Fullname << endl;
    cout << "(" << destination_building.Coords.Lat << ", " << destination_building.Coords.Lon << ")" << endl;

    cout << "Nearest P1 node:" << endl;
    cout << p1_nearest_node << endl;
    cout << "(" << Nodes.at(p1_nearest_node).Lat << ", " << Nodes.at(p1_nearest_node).Lon << ")" << endl;

    cout << "Nearest P2 node:" << endl;
    cout << p2_nearest_node << endl;
    cout << "(" << Nodes.at(p2_nearest_node).Lat << ", " << Nodes.at(p2_nearest_node).Lon << ")" << endl;

    cout << "Nearest destination node:" << endl;
    cout << destination_nearest_node << endl;
    cout << "(" << Nodes.at(destination_nearest_node).Lat << ", " << Nodes.at(destination_nearest_node).Lon << ")" << endl;

    cout << "Person 1's distance to dest: " << p1_distance << " miles" << endl;
    cout << "Path: " << p1_path << endl;

    cout << "Person 2's distance to dest: " << p2_distance << " miles" << endl;
    cout << "Path: " << p2_path << endl;
}


//
// Implement your standard application here
//
void application(
        map<long long, Coordinates> &Nodes, vector<FootwayInfo> &Footways,
        vector<BuildingInfo> &Buildings, graph<long long, double> &G) {
    string person1Building, person2Building;

    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);

    while (person1Building != "#") {
        cout << "Enter person 2's building (partial name or abbreviation)> ";
        getline(cin, person2Building);

        // Lookup buildings and find nearest start nodes
        BuildingInfo building1 = searchBuilding(Buildings, person1Building);
        BuildingInfo building2 = searchBuilding(Buildings, person2Building);

        if (building1.Fullname.empty()) {
            cout << "Person 1's building not found" << endl;
        }

        if (building2.Fullname.empty()) {
            cout << "Person 2's building not found" << endl;
        }

        if (building1.Fullname.empty() || building2.Fullname.empty()) {
            cout << "Enter person 1's building (partial name or abbreviation), or #> ";
            getline(cin, person1Building);
            continue;
        }

        // Find the center of the two buildings
        Coordinates center = centerBetween2Points(building1.Coords.Lat, building1.Coords.Lon, building2.Coords.Lat,
                                                  building2.Coords.Lon);

        // Locate the building closest to the center
        BuildingInfo destBuilding;
        double minDist = INF;
        for (const auto &building: Buildings) {
            double dist = distBetween2Points(center.Lat, center.Lon, building.Coords.Lat, building.Coords.Lon);
            if (dist < minDist) {
                minDist = dist;
                destBuilding = building;
            }
        }


        //Find Nearest node
        long long nearestP1Node = findNearestFootwayNode(building1.Coords, Footways, Nodes);
        long long nearestP2Node = findNearestFootwayNode(building2.Coords, Footways, Nodes);
        long long nearestDestNode = findNearestFootwayNode(destBuilding.Coords, Footways, Nodes);


        map<long long, pair<double, long long>> result1, result2;

        //dijkstra calculations //FIXME
        dijkstraShortestPath(G, nearestP1Node, result1);
        dijkstraShortestPath(G, nearestP2Node, result2);

        double distance1 = result1[nearestDestNode].first;
        double distance2 = result2[nearestDestNode].first;

        string path1 = buildPath(result1, nearestP1Node, nearestDestNode);
        string path2 = buildPath(result2, nearestP2Node, nearestDestNode);

        printInformation(Nodes,building1, building2, destBuilding, building1.Coords, building2.Coords, nearestP1Node,
                         nearestP2Node, nearestDestNode, distance1, distance2, path1, path2);
        cout << endl;
        cout << "Enter person 1's building (partial name or abbreviation), or #> ";
        getline(cin, person1Building);
    }
}

int main() {
    graph<long long, double> G;

    // maps a Node ID to it's coordinates (lat, lon)
    map<long long, Coordinates> Nodes;
    // info about each footway, in no particular order
    vector<FootwayInfo> Footways;
    // info about each building, in no particular order
    vector<BuildingInfo> Buildings;
    XMLDocument xmldoc;

    cout << "** Navigating UIC open street map **" << endl;
    cout << endl;
    cout << std::setprecision(8);

    string def_filename = "map.osm";
    string filename;

    cout << "Enter map filename> ";
    getline(cin, filename);

    if (filename == "") {
        filename = def_filename;
    }

    //
    // Load XML-based map file
    //
    if (!LoadOpenStreetMap(filename, xmldoc)) {
        cout << "**Error: unable to load open street map." << endl;
        cout << endl;
        return 0;
    }

    //
    // Read the nodes, which are the various known positions on the map:
    //
    int nodeCount = ReadMapNodes(xmldoc, Nodes);

    //
    // Read the footways, which are the walking paths:
    //
    int footwayCount = ReadFootways(xmldoc, Footways);

    //
    // Read the university buildings:
    //
    int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

    //
    // Stats
    //
    assert(nodeCount == (int) Nodes.size());
    assert(footwayCount == (int) Footways.size());
    assert(buildingCount == (int) Buildings.size());

    cout << endl;
    cout << "# of nodes: " << Nodes.size() << endl;
    cout << "# of footways: " << Footways.size() << endl;
    cout << "# of buildings: " << Buildings.size() << endl;

    //Build Graph
    for (pair<long long, Coordinates> node: Nodes) {
        G.addVertex(node.first);
    }
    for (const auto &footway: Footways) {
        for (size_t i = 0; i < footway.Nodes.size() - 1; ++i) {
            long long node1 = footway.Nodes[i];
            long long node2 = footway.Nodes[i + 1];
            double distance1 = distBetween2Points(Nodes[node1].Lat, Nodes[node1].Lon, Nodes[node2].Lat,
                                                  Nodes[node2].Lon);
            G.addEdge(node1, node2, distance1);
            double distance2 = distBetween2Points(Nodes[node2].Lat, Nodes[node2].Lon, Nodes[node1].Lat,
                                                  Nodes[node1].Lon);
            G.addEdge(node2, node1, distance2);
        }
    }

    cout << "# of vertices: " << G.NumVertices() << endl;
    cout << "# of edges: " << G.NumEdges() << endl;
    cout << endl;

    // Execute Application
    application(Nodes, Footways, Buildings, G);

    //
    // done:
    //
    cout << "** Done **" << endl;
    return 0;
}
