#include <petscsys.h>
#define int PetscInt
#define float PetscScalar

class Tokamak {
private:
	int rank, size;

	//以下是环坐标系下变量
	//toco is short for toroidal coordinates环坐标
	int r_size, t_size, p_size;
	float ***b_r, ***b_t, ***b_p;
	//b means 磁场
	float ***v_r, ***v_t, ***v_p;
	//v means 速度
	float ***pres_toco, ***dens_toco;
	//pres is short for pressure压强
	//dens is short for density密度
	float ***toco_r, ***toco_p;/*NO DATA!!!!!!!*/
	//坐标上每个点实际的r、p的值
	float ****toco_gather;

	//以下是柱坐标系下变量
	//grid means 网格
	int x_size, y_size, z_size;
	float ***b_x, ***b_y, ***b_z;
	float ***v_x, ***v_y, ***v_z;
	float ***pres_grid, ***dens_grid;
	float ***grid_x, ***grid_z;/*NO DATA!!!!!!!*/
	float ****grid_gather;

public:
	int step;
	Tokamak();//构造函数
	void readin();//从文件中读入磁场、速度、压强、密度的环坐标数据
	void toco2grid();//四个环坐标变四个网格
	void calculateDeltaTime();//计算隐格式时间步长
	void makeMatrixVec();//由四个网格组装成一个矩阵和右端项
	void calculate();//计算出下一个时间步的结果
	void scatterVec();//将计算结果还原为四个网格
	void grid2toco();//转换为环坐标
	void extraPolation();//向外插值
	void print();
};