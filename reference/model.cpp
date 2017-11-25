/*
	如果出现了不能写文件的操作，请在命令行执行以下命令：
	echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
*/

#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <cmath>

char help[] = "解一个偏微分方程\n\
-△u	= cos(3x)sin(πy)	(x,y) ∈ G=(0,π)x(0,1)\n\
u(x,0) = u(x,1) = 0, 0 <= x <= 1\n\
u(0,y) =  sin(πy)/(9+π²)\n\
u(π,y) = -sin(πy)/(9+π²)\n\n";

PetscScalar x_start = 0, x_end = M_PI;
//x的取值范围
PetscScalar y_start = 0, y_end = 1;
//y的取值范围
PetscScalar h = 0.002;//这里设置为0.002需要8核运行1分钟左右，但设置为0.001就可能会在运行完之前就被操作系统击杀[手动捂脸]
//网格的边长为0.01
PetscInt x_size = (x_end-x_start)/h+1;
//x轴方向上的网格数
PetscInt y_size = (y_end-y_start)/h+1;
//y轴方向上的网格数

PetscScalar f(PetscScalar x, PetscScalar y) {
	return cos(3*x)*sin(M_PI*y);
}

enum matType{ AIJ, DENSE };
class EMat {//封装的矩阵
	private:
	PetscErrorCode ierr;
	PetscErrorCode init(matType type, PetscInt rows, PetscInt columns) {
		this -> globalRows = rows;
		this -> globalColumns = columns;
		ierr = MatCreate(PETSC_COMM_WORLD, &self);CHKERRQ(ierr);
		ierr = MatSetSizes(self, PETSC_DECIDE, PETSC_DECIDE, rows, columns);CHKERRQ(ierr);
		switch (type) {
			case AIJ:
				ierr = MatSetType(self, MATMPIAIJ);CHKERRQ(ierr);
				ierr = MatMPIAIJSetPreallocation(self, 5, NULL, 4, NULL);CHKERRQ(ierr);
				break;
			case DENSE:
				ierr = MatSetType(self, MATMPIDENSE);CHKERRQ(ierr);
				ierr = MatMPIDenseSetPreallocation(self, NULL);CHKERRQ(ierr);
				break;
		}
		ierr = MatGetOwnershipRange(self, &rowStart, &rowEnd);CHKERRQ(ierr);
	}
	public:
	Mat self;
	PetscInt globalRows, globalColumns;
	PetscInt rowStart, rowEnd;
	EMat(matType type, PetscInt rows, PetscInt columns) {
		init(type, rows, columns);
	}
	PetscErrorCode setValue(PetscInt row, PetscInt column, PetscScalar value) {
		ierr = MatSetValue(self, row, column, value, INSERT_VALUES);CHKERRQ(ierr);
	}
	PetscErrorCode startAsm() {
		ierr = MatAssemblyBegin(self, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}
	PetscErrorCode endAsm() {
		ierr = MatAssemblyEnd(self, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}
	PetscErrorCode print() {
		ierr = MatView(self, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	}
	PetscErrorCode setName(char *name) {
		PetscObjectSetName((PetscObject)self, name);
	}
};

class EVec {//封装的向量
	private:
	PetscErrorCode ierr;
	PetscScalar *array;
	PetscErrorCode init(PetscInt size) {
		this->size = size;
		ierr = VecCreate(PETSC_COMM_WORLD, &self);CHKERRQ(ierr);
		ierr = VecSetType(self, VECMPI);CHKERRQ(ierr);
		ierr = VecSetSizes(self, PETSC_DECIDE, size);CHKERRQ(ierr);
		ierr = VecGetOwnershipRange(self, &start, &end);CHKERRQ(ierr);
		ierr = VecSet(self, 0);
	}
	public:
	Vec self;
	PetscInt size, start, end;
	EVec(PetscInt size) {
		init(size);
	}
	PetscErrorCode setValue(PetscInt position, PetscScalar value) {
		ierr = VecSetValue(self, position, value, ADD_VALUES);CHKERRQ(ierr);
	}
	PetscErrorCode startAsm() {
		ierr = VecAssemblyBegin(self);CHKERRQ(ierr);
	}
	PetscErrorCode endAsm() {
		ierr = VecAssemblyEnd(self);CHKERRQ(ierr);
	}
	PetscErrorCode print() {
		ierr = VecView(self, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	}
	PetscErrorCode setName(char *name) {
		PetscObjectSetName((PetscObject)self, name);
	}
	PetscErrorCode startGetValues() {
		ierr = VecGetArray(self, &array);CHKERRQ(ierr);
	}
	PetscErrorCode endGetValues() {
		ierr = VecRestoreArray(self, &array);CHKERRQ(ierr);
	}
	PetscScalar& operator[](int i) {
		return array[i-start];
	}
};

class EKsp {//封装的KSP解法器
	private:
	PetscErrorCode ierr;
	public:
	EKsp(){}
	PetscInt calculate(EMat &A, EVec &x, EVec &b) {
		KSP ksp;
		ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
		ierr = KSPSetOperators(ksp, A.self, A.self);CHKERRQ(ierr);
		ierr = KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT,PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);
		ierr = KSPSolve(ksp, b.self, x.self);CHKERRQ(ierr);
		PetscInt iterationTimes;
		ierr = KSPGetIterationNumber(ksp, &iterationTimes);CHKERRQ(ierr);
		return iterationTimes;
	}
};

class MatlabPrinter {//封装的用于输出到MATLAB（ASCII形式）的打印机
	private:
	PetscViewer printer;
	void init(PetscObject obj) {
		char *obj_name, *full_name;
		PetscObjectGetName(obj, &obj_name);
		full_name = (char *)malloc(sizeof(char)*(strlen(obj_name)+10));
		strcpy(full_name, obj_name);
		strcpy(full_name+strlen(obj_name), "_output.m");
		//printf("%s\n", full_name);
		PetscViewerASCIIOpen(PETSC_COMM_WORLD, full_name, &printer);
		PetscViewerPushFormat(printer, PETSC_VIEWER_ASCII_MATLAB);
	}
	public:
	MatlabPrinter(EMat &A) {
		init((PetscObject)A.self);
		MatView(A.self, printer);
	}
	MatlabPrinter(EVec &b) {
		init((PetscObject)b.self);
		VecView(b.self, printer);
	}
};
 
void row2position(PetscInt i, PetscInt &x, PetscInt &y) {//向量映射到矩阵
	x = i%x_size;
	y = i/x_size;
}

PetscInt position2row(PetscInt x, PetscInt y) {//矩阵映射到向量
	return y*x_size+x;
}

bool isBoundary(PetscInt x, PetscInt y) {//矩阵边界判别
	if (x == 0 || x == x_size-1)
		return true;
	if (y == 0 || y == y_size-1)
		return true;
	return false;
}

int main(int argc, char* argv[]) {
	//初始化MPI环境
	PetscInitialize(&argc, &argv, (char*)0, help);
	PetscInt size, rank;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	PetscErrorCode ierr;

	//系数矩阵和右端项
	EMat A(AIJ, x_size*y_size, x_size*y_size);
	EVec b(x_size*y_size);

	//设置系数矩阵
	for (int i = A.rowStart; i < A.rowEnd; i++) {
		PetscInt xi, yi;
		A.setValue(i,i,1);
		row2position(i, xi, yi);
		if (isBoundary(xi,yi) == false) {
			A.setValue(i, position2row(xi-1,yi), -0.25);
			A.setValue(i, position2row(xi+1,yi), -0.25);
			A.setValue(i, position2row(xi,yi-1), -0.25);
			A.setValue(i, position2row(xi,yi+1), -0.25);
		}
	}
	A.startAsm();

	//设置右端项
	const PetscScalar denominator = 9+M_PI*M_PI;
	for (int i = b.start; i < b.end; i++) {
		PetscInt xi, yi;
		row2position(i, xi, yi);
		PetscScalar x = (PetscScalar)xi*h, y = (PetscScalar)yi*h;
		if (isBoundary(xi,yi)) {
			if (xi == 0)
				b.setValue(i, sin(M_PI*y)/denominator);
			else if (xi == x_size-1)
				b.setValue(i, -sin(M_PI*y)/denominator);
		}
		else  b.setValue(i, f(x, y));
	}
	b.startAsm();
	b.endAsm();
	A.endAsm();
	
	//解方程
	EVec x(x_size*y_size);
	(new EKsp())->calculate(A, x, b);

	//向量映射到矩阵
	x.startGetValues();
	EMat ans(DENSE, y_size, x_size);
	for (int i = x.start; i < x.end; i++) {
		PetscInt xi, yi;
		row2position(i, xi, yi);
		ans.setValue(yi, xi, x[i]);
	}
	x.endGetValues();
	ans.startAsm();
	ans.endAsm();

	//打印矩阵
	ans.setName("ans");
	new MatlabPrinter(ans);

	return 0;
}
