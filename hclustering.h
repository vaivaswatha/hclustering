#ifndef _HCLUSTERING_
#define _HCLUSTERING_

#include <cstddef>
#include <climits>
#include <vector>

template <class Pt>
class HClustering
{
 public:
    typedef float (*DistFunc)(Pt &el1, Pt &el2);
    enum LinkageType {
	LINK_MINIMUM,
	LINK_MAXIMUM,
	LINK_AVERAGE,
    };

    // p_distFunc is a pointer to the distance computation
    // function. It takes in two <Pt> elements and returns
    // a float. This function is expected to be symmetric.
    // p_numClusters is the number of clusters that need to
    // be formed by this algorithm.
    // p_linkageType is an optional parameter that defines
    // how to compute linkage/distance between two clusters.
    HClustering(DistFunc p_distFunc, unsigned p_numClusters, 
		LinkageType p_linkageType = LINK_MAXIMUM) 
    {
	distFunc = p_distFunc;
	numClusters = p_numClusters;
	linkageType = p_linkageType;
    }

    void cluster(Pt *data, unsigned nElm);
    // getResult will assign a cluster id (not necessarily
    // contiguous) to each slot in result[], assuming they
    // correspond to each slot in data[] above.
    // NOTE: result[] is assumed to be of size nElm
    void getResult(unsigned result[]);

 private:

    struct TreeNode
    {
	// Can be used to index "clusters[]"
	unsigned id;
	// Immediate rep. Follow through for root.
	unsigned rep;
	// If this node is rep, who are other members?
	// This includes all members in the sub-tree, 
	// not just the immediate children. Makes it
	// easy to merge clusters this way.
	std::vector<unsigned> members;
	// assert ((rep == id) ||
	//         (rep != id && members.size() == 0))

	TreeNode(void) { 
	    // not really needed.
	    id = UINT_MAX;
	    // Initially every node is a cluster by itself.
	    rep = UINT_MAX;
	};
    };

    // User defined distance function.
    DistFunc distFunc;
    // Number of clusters to be formed ultimately.
    unsigned numClusters;
    LinkageType linkageType;
    std::vector<TreeNode> clusters;
    void cacheDistances(Pt *data, unsigned nElm);
    float getDistance(unsigned i, unsigned j, Pt *p_data = NULL);

    HClustering(void) { ; };
    // follow through the rep field till a rep TreeNode.
    unsigned getRepRoot(unsigned i);
    // compute linkage between clusters i and k.
    float computeLinkage(unsigned i, unsigned k);
};

// References:
//   http://elki.dbs.ifi.lmu.de/wiki/Tutorial/HierarchicalClustering
//   http://home.deib.polimi.it/matteucc/Clustering/tutorial_html/hierarchical.html
//   http://nlp.stanford.edu/IR-book/html/htmledition/single-link-and-complete-link-clustering-1.html

#endif // _HCLUSTERING_
