#ifndef _HCLUSTERING_
#define _HCLUSTERING_

#include <cstddef>

template <class Pt>
class HClustering
{
 public:
    typedef float (*DistFunc)(Pt &el1, Pt &el2);

    HClustering(DistFunc p_distFunc, unsigned p_numClusters) {
	distFunc = p_distFunc;
	numClusters = p_numClusters;
    }

    void cluster(Pt *data, unsigned nElm);

 private:
    DistFunc distFunc;
    unsigned numClusters;
    HClustering(void) { ; };
};

#endif // _HCLUSTERING_
