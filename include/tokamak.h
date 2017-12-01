#include <petscsys.h>
#define int PetscInt
#define float PetscScalar

class Tokamak {
private:
	int rank, size;
	int r_size, t_size, p_size;
	float ***b_r, ***b_t, ***b_p;
	float ***v_r, ***v_t, ***v_p;
	float ***pres, ***dens;
	float ***toco_r, ***toco_p;
	float ****gather;
public:
	int step;
	Tokamak();//构造函数
	void readin();//从文件中读入磁场、速度、压强、密度的环坐标数据
	void switch2Grid();//四个环坐标变四个网格
	void calculateDeltaTime();//计算隐格式时间步长
	void makeMatrixVec();//由四个网格组装成一个矩阵和右端项
	void calculate();//计算出下一个时间步的结果
	void scatterVec();//将计算结果还原为四个网格
	void switch2Cyc();//转换为环坐标
	void extraPolation();//向外插值
	void print();
};