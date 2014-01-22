#include <vector>
#include <cassert>
#include "hclustering.h"

#define INF -1.0
#define IS_REP(x) (clusters[(x)].rep == &clusters[(x)])

inline unsigned make1D(unsigned i, unsigned j, unsigned nElm)
{
    return ((nElm * i) + j);
}

template <class Pt>
float HClustering<Pt>::computeLinkage(unsigned c1, unsigned c2)
{
    unsigned i, j;
    float retVal = INF;

    assert(IS_REP(c1) && IS_REP(c2));

    // for looping convinience, we add the rep of a node to
    // its "members" list and pop it back later
    clusters[c1].members.push_back(&clusters[c1]);
    clusters[c2].members.push_back(&clusters[c2]);

    switch (linkageType) {
    case LINK_MINIMUM:
    {
	float minDistanceSeen = INF;
	for (i = 0; i < clusters[c1].members.size(); i++) {
	    for (j = 0; j < clusters[c2].members.size(); j++) {
		float distance = distanceMatrix[make1D(clusters[c1].members[i], 
						       clusters[c2].members[j], 
						       clusters.size())];
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
		float distance = distanceMatrix[make1D(clusters[c1].members[i], 
					        clusters[c2].members[j], 
						clusters.size())];
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
    std::vector<std::pair<unsigned, float> > shortestDistance (nElm, INF);

    // Initially there are nElm clusters, so have these many
    // TreeNodes and link them to the input data.
    clusters.resize(nElm);
    for (i = 0; i < nElm; i++) {
	clusters[i].el = &data[i];
	clusters[i].id = i;
    }

    // TODO: This can be a triangural matrix to save space.
    distanceMatrix.resize(nElm * nElm);

    // fill up the distance matrix.
    for (i = 0; i < nElm; i++)
	for (j = 0; j < nElm; j++)
	    distanceMatrix[make1D(i, j, nElm)] = distFunc(*data[i], *data[j]);

    // compute the shortest distance array.
    for (i = 0; i < nElm; i++) {
	for (j = 0; j < nElm; j++) {
	    float linkageIJ;
	    // No need to check for cluster rep here since
	    // at this everybody is his own rep.
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
		shortestDistanceFrom = shortestDistance[i].first;
	    }
	}

	assert (shortestDistanceSeen != INF && shortestDistanceFrom < nElm);
	i = shortestDistanceFrom;
	j = shortestDistance[i].first;
	assert(i != j);
	// i and j must both represent their clusters already
	// (since we only look at reps during during each iteration).
	assert(IS_REP(i) && IS_REP(j));

	if (clusters[i].members.size() < clusters[j].members.size()) {
	    // let 'j' represent the smaller cluster.
	    unsigned temp = i;
	    i = j;
	    j = temp;
	}
	// We now need to combine the clusters i and j,
	// So make 'j' belong to 'i' as it is smaller.
	clusters[i].members.insert(clusters[i].members.end(), 
				   clusters[j].members.begin(),
				   clusters[j].members.end());
	clusters[j].rep = &clusters[i];
	{
	    // Just to free up memory used for the members array.
	    std::vector<TreeNode *> tmp;
	    clusters[j].members.swap(tmp);
	}
	
	// compute linkage between clusters (w.r.t i).
	for (k = 0; k < nElm; k++) {
	    float distance;

	    // ignore non rep data. also ignore i itself.
	    if (!IS_REP(k) || k == i)
		continue;
	    distance = computeLinkage(i, k);
	    assert(distance != INF);
	    // Now update the shortestDistance array.
	    if (shortestDistance[i].second > distance) {
		shortestDistance[i].second = distance;
		shortestDistance[i].first = k;
	    }
	    // for 'k', the shortestDistance array may be
	    // pointing to a non-rep now. Correct that.
	    if (!IS_REP(shortestDistance[k].first)) {
		// This must have just happened.
		assert(clusters[shortestDistance[k].first].rep ==
		       &clusters[i]);
		if (linkageType == LINK_MINIMUM) {
		    // In case of LINK_MINIMUM the nearest
		    // cluster would now be 'i' (since 'j'
		    // now belongs to 'i').
		    assert(shortestDistance[k].second == distance);
		    shortestDistance[k].first = i;
		} else if (distance > shortestDistance[k].second) {
		    // NOTE: This can never happen in LINK_MINIMUM.
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
		// computeLinkage(i, k) cannot be smaller than
		// the known shortest distance for node k
		assert(distance >= shortestDistance[k].second);
		// At this point, the shortest distance must be to a rep.
		assert(IS_REP(clusters[k].first));
	    }
	}
    }
}

#ifdef TEST_CODE
int main (void)
{
    
}
#endif // TEST_CODE
