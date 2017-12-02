#万能makefile，适用于所有单源文件编译
#使用方法如下：
#比如你的代码文件名为：myCode.cpp 或 myCode.c
#那么在命令行中输入 make myCode.exe 即可
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test
PETSC_ALL_LIB := ${PETSC_MAT_LIB} ${PETSC_VEC_LIB} ${PETSC_KSP_LIB}

%.exe : %.o
	-${CLINKER} $^ -o $@  ${PETSC_ALL_LIB}
	${RM} %.o