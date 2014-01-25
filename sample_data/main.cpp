#include "../hclustering.cpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <cstring>
#include <cassert>
#include <cstdlib>

static unsigned numClusters = 1;
static unsigned dataSize = 0;

struct Point {
    float x, y;
    Point(float x, float y) { this->x = x; this->y = y; };
    Point(void) { x = 0.0; y = 0.0; };
};

float distPoint (Point &p1, Point &p2)
{
    return (((p1.x - p2.x)*(p1.x - p2.x)) +
	    ((p1.y - p2.y)*(p1.y - p2.y)));
}

static void printForGnuPlot(Point *data, unsigned *result)
{
    // This will create numClusters files in /tmp/hclustering/
    // with each file given only a number (clusterId) and points
    // belonging to that clustered entered in the file.
    // Also a plot.gp file is generated which can be sourced into
    // gnuplot to plot the final graph.

    std::map<unsigned,unsigned> clusterIdMap;
    unsigned idCounter = 1, i;
    std::ofstream files[numClusters+1];
    char filename[20];

    // remove existing files if any.
    for (i = 1; i <= numClusters; i++) {
	sprintf(filename, "/tmp/hclustering/%d.txt", i);
	files[i].open(filename, std::ios::out | std::ios::trunc);
    }
    // files[0] will be the gnuplot script.
    files[0].open("/tmp/hclustering/plot.gp", std::ios::out | std::ios::trunc);

    for (i = 0; i < dataSize; i++) {
	if (clusterIdMap[result[i]] == 0)
	    clusterIdMap[result[i]] = idCounter++;
	files[clusterIdMap[result[i]]] << data[i].x << " " <<
	      data[i].y << std::endl;
    }

    assert(idCounter == (numClusters+1));

    // write out the gnuplot script
    files[0] << "set nokey\n";
    files[0] << "plot ";
    for (i = 1; i <= numClusters; i++) {
	files[0] << "     '/tmp/hclustering/" << i << ".txt'";
	if (i < numClusters)
	    files[0] << ", \\\n";
    }
    files[0] << "\npause -1 \"Hit any key to exit\"\n";

    // Close all files
    for (i = 0; i <= numClusters; i++)
	files[i].close();
}

int main(int argc, char *argv[])
{
    std::ifstream file;
    float x, y;
    std::vector<Point> data;
    std::vector<unsigned> result;

    if (argc != 3) {
	std::cerr << "Usage: " << argv[0] << " filename.txt numClusters";
    }

    file.open(argv[1], std::ios::in);
    numClusters = atoi(argv[2]);

    while (file >> x >> y) {
	data.push_back(Point(x, y));
    }
    dataSize = data.size();

    if (numClusters <= 1 || numClusters > dataSize) {
	std::cerr << "numClusters must be greater than 1 and smaller than data size\n";
    }

    /*
    for (unsigned i = 0; i < dataSize; i++)
	std::cout << data[i].x << "\t" <<  data[i].y << "\n";
    exit(1);
    */

    // Initialize the HClustering object.
    HClustering<Point> 
	hc(distPoint, numClusters, HClustering<Point>::LINK_MAXIMUM);
    // Cluster the data and get results.
    hc.cluster(&data[0], dataSize);
    result.resize(dataSize);
    hc.getResult(&result[0]);

    printForGnuPlot(&data[0], &result[0]);

    return 0;
}
