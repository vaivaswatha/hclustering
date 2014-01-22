#include <vector>
#include "hclustering.h"

template <class Pt>
struct TreeNode
{
    Pt *el;
    TreeNode *root;
    std::vector<TreeNode *> members;
    // assert ((root == this && members.size() > 0) ||
    //         (root != this && members.size() == 0))

    TreeNode(void) { 
	el = NULL;
	root = NULL;
    };
};

unsigned make1D(unsigned i, unsigned j, unsigned nElm)
{
    return ((nElm * i) + j);
}

#define INF -1.0

template <class Pt>
void HClustering<Pt>::cluster(Pt *data, unsigned nElm)
{
    unsigned i, j, n;
    std::vector<TreeNode<Pt> > clusters;
    std::vector<float> distanceMatrix;
    std::vector<std::pair<unsigned, float> > shortestDistance (nElm+1, INF);

    // Initially there are nElm clusters, so have these many
    // TreeNodes and link them to the input data.
    clusters.resize(nElm);
    for (i = 0; i < nElm; i++) {
	clusters[i].el = &data[i];
    }

    // TODO: This can be a triangural matrix to save space.
    distanceMatrix.resize(nElm * nElm);

    // fill up the distance matrix.
    for (i = 0; i < nElm; i++)
	for (j = 0; j < nElm; j++)
	    distanceMatrix[make1D(i, j, nElm)] = distFunc(*data[i], *data[j]);

    // compute the shortest distance array.
    for (i = 0; i < nElm; i++) {
	for (j = i+1; j < nElm; j++) {
	    if (shortestDistance[i].second == INF ||
		distanceMatrix[make1D(i, j, nElm)] < shortestDistance[i].second)
	    {
		shortestDistance[i].first = j;
		shortestDistance[i].second = distanceMatrix[make1D(i, j, nElm)];
		// symmetry!
		shortestDistance[j].first = i;
		shortestDistance[j].second = shortestDistance[i].second;
	    }
	}
    }

    // Run the clustering loop till there are only numClusters remaining.
    for (n = nElm; n >= numClusters; n--) {
	
    }
}

#ifdef TEST_CODE
int main (void)
{
    

}
#endif // TEST_CODE
