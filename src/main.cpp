#include <petscsys.h>
#include "tokamak.h"

int main(int argc, char* argv[]) {
	char help[] = "Find any useful things in readme.md\n";
	//PetscOptionsSetValue(NULL,"-log_view","");
	//PetscOptionsSetValue(NULL,"--with-scalar-type","real");
	PetscInitialize(&argc, &argv, (char*)0, help);
	PetscOptionsSetFromOptions(NULL);
	//运行时用 -help 就能看到这个help
	PetscInt size, rank;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	Tokamak tokamak;
	tokamak.readin();//从文件中读入磁场、速度、压强、密度的环坐标数据
	for (int i = 0; i < tokamak.step; i++) {
		tokamak.switch2Grid();//四个环坐标变四个网格
		tokamak.calculateDeltaTime();//计算隐格式时间步长
		tokamak.makeMatrixVec();//由四个网格组装成一个矩阵和右端项
		tokamak.calculate();//计算出下一个时间步的结果
		tokamak.scatterVec();//将计算结果还原为四个网格
		tokamak.switch2Cyc();//转换为环坐标
		tokamak.extraPolation();//向外插值
	}
	tokamak.print();
	PetscFinalize();
	return 0;
}