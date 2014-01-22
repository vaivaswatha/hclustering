#ifndef _HCLUSTERING_
#define _HCLUSTERING_

#include <cstddef>
#include <climits>

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

    HClustering(DistFunc p_distFunc, unsigned p_numClusters, 
		LinkageType p_linkageType = LINK_MINIMUM) 
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
	Pt *el;
	unsigned id;
	unsigned rep;
	std::vector<unsigned> members;
	// assert ((rep == id && members.size() > 0) ||
	//         (rep != id && members.size() == 0))

	TreeNode(void) { 
	    el = NULL;
	    // Initially every node is a cluster by itself.
	    rep = UINT_MAX;
	};
    };

    DistFunc distFunc;
    unsigned numClusters;
    LinkageType linkageType;
    std::vector<TreeNode> clusters;
    std::vector<float> distanceMatrix;

    HClustering(void) { ; };
    // compute linkage between clusters i and k.
    float computeLinkage(unsigned i, unsigned k);
};

#endif // _HCLUSTERING_
