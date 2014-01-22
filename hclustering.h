#ifndef _HCLUSTERING_
#define _HCLUSTERING_

#include <cstddef>

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

 private:

    struct TreeNode
    {
	Pt *el;
	unsigned id;
	TreeNode *root;
	std::vector<TreeNode *> members;
	// assert ((root == this && members.size() > 0) ||
	//         (root != this && members.size() == 0))

	TreeNode(void) { 
	    el = NULL;
	    // Initially every node is a cluster by itself.
	    root = this;
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
