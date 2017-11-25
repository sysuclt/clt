#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>

class Tokamak {
	public:
	void readin() {
		freopen("data/Bp_nova.dat","r",stdin);
		//注意是从根目录下开始的
	}
};