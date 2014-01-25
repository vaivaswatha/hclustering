#include <vector>
#include <cassert>
#include <cstdlib>
#include "hclustering.h"

#define INF -1.0
#define IS_REP(x) (clusters[(x)].rep == clusters[(x)].id)

// Currently not caching distances at all. distFunc() will
// be called everytime we need a distance.
template <class Pt>
void HClustering<Pt>::cacheDistances(Pt *data, unsigned nElm)
{
    // NOP, just initialize getDistance.
    getDistance(0, 0, data);
}

template <class Pt>
float HClustering<Pt>::getDistance(unsigned i, unsigned j, Pt *p_data)
{
    static Pt *data = NULL;

    if (data == NULL) {
	// data must be initialized before getDistance is called with p_data being NULL.
	assert(p_data != NULL);
	data = p_data;
    }
    return distFunc(data[i], data[j]);
}

template <class Pt>
float HClustering<Pt>::computeLinkage(unsigned c1, unsigned c2)
{
    unsigned i, j;
    float retVal = INF;

    assert(IS_REP(c1) && IS_REP(c2));

    // for looping convinience, we add the rep of a node to
    // its "members" list and pop it back later
    clusters[c1].members.push_back(c1);
    clusters[c2].members.push_back(c2);

    switch (linkageType) {
    case LINK_MINIMUM:
    {
	float minDistanceSeen = INF;
	for (i = 0; i < clusters[c1].members.size(); i++) {
	    for (j = 0; j < clusters[c2].members.size(); j++) {
		float distance = getDistance(clusters[c1].members[i],  
					     clusters[c2].members[j]);
		if (minDistanceSeen == INF ||
		    minDistanceSeen > distance)
		    {
			minDistanceSeen = distance;
		    }
	    }
	}
	retVal = minDistanceSeen;
    }
    break;
    case LINK_MAXIMUM:
    {
	// abusing INF here.
	float maxDistanceSeen = INF;
	for (i = 0; i < clusters[c1].members.size(); i++) {
	    for (j = 0; j < clusters[c2].members.size(); j++) {
		float distance = getDistance(clusters[c1].members[i],  
					     clusters[c2].members[j]);
		if (maxDistanceSeen == INF ||
		    maxDistanceSeen < distance)
		{
		    maxDistanceSeen = distance;
		}
	    }
	}
	retVal = maxDistanceSeen;
    }
    break;
    case LINK_AVERAGE:
	assert(0 && "LINK_AVERAGE not yet supported");
	break;
    }

    // pop back the elements artificially inserted earlier.
    clusters[c1].members.pop_back();
    clusters[c2].members.pop_back();

    assert(retVal != INF);
    return retVal;
}

template <class Pt>
void HClustering<Pt>::cluster(Pt *data, unsigned nElm)
{
    unsigned i, j, k, n;
    typedef std::pair<unsigned,float> shortestDistanceRow;
    std::vector<shortestDistanceRow> 
	shortestDistance(nElm, shortestDistanceRow(nElm, INF));

    // Initially there are nElm clusters, so have these many
    // TreeNodes and link them to the input data.
    clusters.resize(nElm);
    for (i = 0; i < nElm; i++) {
	clusters[i].id = i;
	clusters[i].rep = i;
    }

    cacheDistances(data, nElm);

    // compute the shortest distance array.
    for (i = 0; i < nElm; i++) {
	for (j = 0; j < nElm; j++) {
	    float linkageIJ;
	    // No need to check for cluster rep here since
	    // at this point everybody is his own rep.
	    assert(IS_REP(i));
	    if (i == j)
		continue;
	    linkageIJ = computeLinkage(i, j);
	    if (shortestDistance[i].second == INF ||
		linkageIJ < shortestDistance[i].second)
	    {
		shortestDistance[i].first = j;
		shortestDistance[i].second = linkageIJ;
	    }
	}
    }

    // Run the clustering loop till there are only numClusters remaining.
    for (n = nElm; n > numClusters; n--) {
	float shortestDistanceSeen = INF;
	unsigned shortestDistanceFrom = nElm;

	// TODO: Keep a list of elements that represent a cluster.
	// TODO: Maybe some kind of sorted list? I don't know
	for (i = 0; i < nElm; i++) {
	    // If this element is not the rep for the cluster it is in:
	    if (!IS_REP(i))
		continue;
	    if (shortestDistanceSeen == INF ||
		shortestDistanceSeen > shortestDistance[i].second)
	    {
		shortestDistanceSeen = shortestDistance[i].second;
		shortestDistanceFrom = i;
	    }
	}

	assert (shortestDistanceSeen != INF && shortestDistanceFrom < nElm);
	i = shortestDistanceFrom;
	j = shortestDistance[i].first;
	assert(i != j);
	// i and j must both represent their clusters already
	// (since we only look at reps during during each iteration).
	assert(IS_REP(i) && IS_REP(j));

	// If the shortest distance was from i to j, then j to i
	// must also be the shortest distance.
	assert(shortestDistance[j].second == shortestDistanceSeen);

	if (clusters[i].members.size() < clusters[j].members.size()) {
	    // let 'j' represent the smaller cluster.
	    unsigned temp = i;
	    i = j;
	    j = temp;
	}
	// We now need to combine the clusters i and j,
	// So make 'j' belong to 'i' as it is smaller.
	clusters[i].members.push_back(j);
	clusters[i].members.insert(clusters[i].members.end(), 
				   clusters[j].members.begin(),
				   clusters[j].members.end());
	clusters[j].rep = i;
	{
	    // Just to free up memory used for the members array.
	    std::vector<unsigned> tmp;
	    clusters[j].members.swap(tmp);
	}
	// The shortest distance from i needs recomputation.
	shortestDistance[i].first = nElm;
	shortestDistance[i].second = INF;
	
	// compute linkage between clusters (w.r.t i).
	for (k = 0; k < nElm; k++) {
	    float distance;

	    // ignore non rep data. also ignore i itself.
	    if (!IS_REP(k) || k == i)
		continue;
	    distance = computeLinkage(i, k);
	    assert(distance != INF);
	    // Now update the shortestDistance array.
	    if ((shortestDistance[i].second == INF ||
		 shortestDistance[i].second > distance))
	    {
		shortestDistance[i].second = distance;
		shortestDistance[i].first = k;
	    }
	    // for 'k', the shortestDistance array may be
	    // pointing to a non-rep now, or the distance
	    // to cluster 'i' may change (as cluster 'i' has
	    // more members now - not always LINK_MINIMUM).
	    // TODO: Maybe this can be done on-demand, i.e.
	    // when we see (i, j) distance with j being non-rep.
	    if (!IS_REP(shortestDistance[k].first) ||
		shortestDistance[k].first == i) 
	    {
		assert(clusters[shortestDistance[k].first].rep ==
		       clusters[i].id);
		if (linkageType == LINK_MINIMUM) {
		    // In case of LINK_MINIMUM the nearest
		    // cluster would now be 'i' (since 'j'
		    // now belongs to 'i').
		    assert(shortestDistance[k].second == distance);
		    shortestDistance[k].first = i;
		} else {
		    // computeLinkage(i, k) cannot be smaller than
		    // the known shortest distance for node k
		    assert(distance >= shortestDistance[k].second);
		    // reset shortestDistance for 'k' since its now invalid.
		    shortestDistance[k].second = INF;
		    shortestDistance[k].first = nElm;
		    // Will now need to look through ALL
		    // clusters for a new nearest linkage.
		    for (j = 0; j < nElm; j++) {
			float linkageKJ;
			if (j == k || !IS_REP(j))
			    continue;
			linkageKJ = computeLinkage(k, j);
			if (shortestDistance[k].second == INF ||
			    linkageKJ < shortestDistance[k].second)
			{
				shortestDistance[k].first = j;
				shortestDistance[k].second = linkageKJ;
			}
		    }
		}
		// At this point, the shortest distance from k must be to a rep.
		assert(IS_REP(shortestDistance[k].first));
	    }
	}
    }
}

template <class Pt>
unsigned HClustering<Pt>::getRepRoot(unsigned i)
{
    while (!IS_REP(i)) {
	// bad programming to modify 'i'.
	i = clusters[i].rep;
    }

    return i;
}

// result[] is assumed to be of size nElm
template <class Pt>
void HClustering<Pt>::getResult(unsigned result[])
{
    unsigned i;

    for (i = 0; i < clusters.size(); i++) {
	result[i] = getRepRoot(i);
    }
}

#ifdef TEST_CODE

#include <ctime>
#include <iostream>
#include <fstream>
#include <cstring>
#include <random>
#include <map>

#define SIZE 1000
#define NUM_CLUSTERS 50

// #define GNUPLOT

struct Point {
    float x, y;
};

float distPoint (Point &p1, Point &p2)
{
    return (((p1.x - p2.x)*(p1.x - p2.x)) +
	    ((p1.y - p2.y)*(p1.y - p2.y)));
}

static void printForGnuPlot(Point *data, unsigned *result)
{
    // This will create NUM_CLUSTERS files in /tmp/hclustering/
    // with each file given only a number (clusterId) and points
    // belonging to that clustered entered in the file.
    // Also a plot.gp file is generated which can be sourced into
    // gnuplot to plot the final graph.

    std::map<unsigned,unsigned> clusterIdMap;
    unsigned idCounter = 1, i;
    std::ofstream files[NUM_CLUSTERS+1];
    char filename[20];

    // remove existing files if any.
    for (i = 1; i <= NUM_CLUSTERS; i++) {
	sprintf(filename, "/tmp/hclustering/%d.txt", i);
	files[i].open(filename, std::ios::out | std::ios::trunc);
    }
    // files[0] will be the gnuplot script.
    files[0].open("/tmp/hclustering/plot.gp", std::ios::out | std::ios::trunc);

    for (i = 0; i < SIZE; i++) {
	if (clusterIdMap[result[i]] == 0)
	    clusterIdMap[result[i]] = idCounter++;
	files[clusterIdMap[result[i]]] << data[i].x << " " <<
	      data[i].y << std::endl;
    }

    assert(idCounter == (NUM_CLUSTERS+1));

    // write out the gnuplot script
    files[0] << "set nokey\n";
    files[0] << "plot ";
    for (i = 1; i <= NUM_CLUSTERS; i++) {
	files[0] << "     '/tmp/hclustering/" << i << ".txt'";
	if (i < NUM_CLUSTERS)
	    files[0] << ", \\\n";
    }
    files[0] << "\npause -1 \"Hit any key to exit\"\n";

    // Close all files
    for (i = 0; i <= NUM_CLUSTERS; i++)
	files[i].close();
}

int main (void)
{
    unsigned i, result[SIZE];
    Point data[SIZE];
    // Initialize the HClustering object.
    HClustering<Point> 
	hc(distPoint, NUM_CLUSTERS, HClustering<Point>::LINK_MAXIMUM);

    // Some stuff to generate a uniform distribution 
    // of reals in the given range (Using 1.0-25.0 now).
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(1.0f, 255.0f);
    // generate random points.
    for (i = 0; i < SIZE; i++) {
	data[i].x = dist(gen);
	data[i].y = dist(gen);
    }

    // Cluster the data and get results.
    hc.cluster(data, SIZE);
    hc.getResult(result);

#ifdef GNUPLOT
        // print results for gnuplot.
    printForGnuPlot(data, result);
#else
    // Print the data along with the clusters they are assigned to.
    for (i = 0; i < SIZE; i++) {
	std::cout << data[i].x << " " << data[i].y <<
	    " " << result[i] << std::endl;
    }
#endif // GNUPLOT

    return 0;
}

#endif // TEST_CODE
