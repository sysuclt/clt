#include "tokamak.h"
#include "shepard.h"
#include <cstdio>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#define int PetscInt
#define float PetscScalar

Tokamak::Tokamak() {//后面记得加入许多接口配置
	//系统初始化
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	step = 100;

	//初始化环坐标
	r_size = 111;
	t_size = 200;
	p_size = 8;
	toco_gather = new float***[10];
	//申请动态内存
	if (rank == 0) {
		for (int i = 0; i < 10; i++) {
			toco_gather[i] = new float**[r_size];
			for (int j = 0; j < r_size; j++) {
				toco_gather[i][j] = new float*[t_size];
				for (int k = 0; k < t_size; k++)
					toco_gather[i][j][k] = new float[p_size];
			}
		}
		b_r			= toco_gather[0];
		b_t			= toco_gather[1];
		b_p			= toco_gather[2];
		v_r			= toco_gather[3];
		v_t			= toco_gather[4];
		v_p			= toco_gather[5];
		pres_toco	= toco_gather[6];
		dens_toco	= toco_gather[7];
		toco_r		= toco_gather[8];
		toco_p		= toco_gather[9];
	}
	else
		for (int i = 0; i < 10; i++)
			toco_gather[i] = NULL;

	//初始化柱坐标
	x_size = 256;
	y_size = 8;
	z_size = 256;
	grid_gather = new float***[8];
	//申请动态内存
	if (rank == 0) {
		for (int i = 0; i < 10; i++) {
			grid_gather[i] = new float**[x_size];
			for (int j = 0; j < x_size; j++) {
				grid_gather[i][j] = new float*[y_size];
				for (int k = 0; k < y_size; k++)
					grid_gather[i][j][k] = new float[z_size];
			}
		}
		b_x			= grid_gather[0];
		b_y			= grid_gather[1];
		b_z			= grid_gather[2];
		v_x			= grid_gather[3];
		v_y			= grid_gather[4];
		v_z			= grid_gather[5];
		pres_grid	= grid_gather[6];
		dens_grid	= grid_gather[7];
		grid_x		= toco_gather[8];
		grid_z		= toco_gather[9];
	}
	else
		for (int i = 0; i < 8; i++)
			grid_gather[i] = NULL;
}

void Tokamak::readin() {
	//从文件中读入磁场、速度、压强、密度的环坐标数据
	if (rank != 0)
		return;
	char data_file_path[] = "./data/";
	char *data_file_name[] = {"Br_nova.dat", "Bt_nova.dat", "Bp_nova.dat", \
							  "vr_nova.dat", "vt_nova.dat", "vp_nova.dat", \
							  "P_nova.dat", "rho_nova.dat"};
	for (int i = 0; i < 8; i++) {
		//获得文件完整路径名并打开文件
		char *file_path_name = new char[strlen(data_file_path) + strlen(data_file_name[i])+1];
		file_path_name[0] = 0;
		strcat(file_path_name, data_file_path);
		strcat(file_path_name, data_file_name[i]);
		//printf("%s\n", file_path_name);
		FILE *fpr = fopen(file_path_name, "r");
		//开始读取文件
		for (int j = 0; j < r_size; j++) {
			for (int k = 0; k < t_size; k++) {
				for (int u = 0; u < p_size; u++) {
					fscanf(fpr, "%lf", &toco_gather[i][j][k][u]);
					//printf("%d %d %d %d : %e\n", i, j, k, u, toco_gather[i][j][k][u]);
				}
			}
		}
	}
}

void Tokamak::toco2grid(){
	//四个环坐标变四个网格
	if (rank != 0)
		return;
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < p_size; j++) {
			double *x = new double[r_size*t_size];
			double *y = new double[r_size*t_size];
			double *f = new double[r_size*t_size];
			for (int k = 0; k < r_size; k++)
				for (int u = 0; u < t_size; u++) {
					x[k*t_size+u] = toco_r[k][u][j];/*顺序可能有问题！！！*/
					y[k*t_size+u] = toco_p[k][u][j];/*顺序可能有问题！！！*/
					f[k*t_size+u] = toco_gather[i][k][u][j];
				}
			//transformer.Set_Missing();//这是啥！！！
			CShepard2d transformer;
			transformer.Interpolate(x, y, f, r_size*t_size, 33, 89);/*参数调整！！！*/
			for (int k = 0; k < x_size; k++)
				for (int u = 0; u < z_size; u++)
					transformer.GetValue(grid_x[k][j][u], grid_z[k][j][u], grid_gather[i][k][j][u]);
		}
	}
}

void Tokamak::calculateDeltaTime() {
	//计算隐格式时间步长
}

void Tokamak::makeMatrixVec() {
	//由四个网格组装成一个矩阵和右端项
}

void Tokamak::calculate() {
	//计算出下一个时间步的结果
}

void Tokamak::scatterVec() {
	//将计算结果还原为四个网格
}

void Tokamak::grid2toco() {
	//转换为环坐标
}

void Tokamak::extraPolation() {
	//向外插值
}

void Tokamak::print() {
	//输出
}