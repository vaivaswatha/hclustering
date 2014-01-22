#ifndef _HCLUSTERING_
#define _HCLUSTERING_

#include <cstddef>

template <class Pt>
class HClustering
{
 public:
    typedef float (*DistFunc)(Pt &el1, Pt &el2);

 private:
    DistFunc distFunc;
    HClustering(void) { ; };

 public:
    HClustering(DistFunc p_distFunc) {
	distFunc = p_distFunc;
    }

    void cluster(Pt *data, unsigned nElm);
};

#endif // _HCLUSTERING_
