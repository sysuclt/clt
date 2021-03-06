# 说明文档

## 更新日志

### 2017/12/08 22:00
- 研发了FEI（Formula & Equation Interperter）
- 尝试将FEI应用于之前写的解泊松方程的程序中，但是出现了未知bug

### 2017/12/01 22:00
- 对函数名做出了一定调整
- 分离了`tokamak.h`和`tokamak.cpp`
- 完成了数据读入函数，但部分数据缺失
- 在网上找到了别人已经写好了的把环坐标变换乘网格的[c++程序](https://github.com/UoA-eResearch/saga-gis/blob/master/saga-gis/src/modules/grid/grid_gridding/Shepard.cpp)
- 完成了从环坐标变换到网格的程序（直接移植了别人封装好的代码）`Tokamak::toco2grid()`

#### 有待解决的问题
- 缺少环坐标每个点的实际坐标值的数据（包括第三维）
- 缺少网格坐标每个点的实际坐标的定义（包括第三维）
- 忘了网格坐标的规模，现在暂时随手写成了256×256×8
- 考证环坐标变换成网格程序的正确性，顺便质疑一下只用二维面的数据没问题么
- 用了别人开源代码要怎么声明别人的版权呀
- 环坐标变换成网格需要参数输入，需要如何确定
- 暂时还没写自定义接口
- 怎么都没研究明白用[`PetscOptionsSetValue()`](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsSetValue.html)配置[`PetscInt`](https://www.mcs.anl.gov/petsc/petsc-3.7/docs/manualpages/Sys/PetscInt.html#PetscInt)位数的方法

### 2017/11/15 18:00
- 程序框架基本没有问题，加了一点小改动
- 原来的数据处理函数应该没问题
- 了解了环坐标以及柱坐标体系，详见后文坐标体系
- qs2grd并不是完整的转换函数，确认了新的转换函数
    - `map_nova`环坐标到柱坐标转换函数
        - `qshep2`构造差值函数
        - `qs2grd`使用插值函数进行坐标映射
    - `map_xz2st`柱坐标到环坐标转换函数
- 计划使用并行方法

#### 接下来的工作任务
- 设计数据结构
- 确定程序框架细节
- 翻译所有转换函数

### 2017/11/25 10:00
- 确定了项目文档结构
- 完整地搭建了开发环境
- 搭建了方便快捷的开发环境
- 在GitHub上建立了开发小组，将项目在GitHub上同步，请把GitHub账户发给我拉你们进小组

#### 有待确认的问题
- 怀疑原来的处理数据函数有问题，确认是否需要重新实现
- 确认程序流程是否有问题（见`src/main.cpp`）
- 了解环坐标体系，确认网格坐标范围
- qs2grd 还要计算导数？表示眉头一皱，发现事情并不简单

## 文件结构
- `makefile`可以编译整个工程的makefile
- `readme.md`有问题，就看这里
- `toBuildSys/`一个可以快去去到本计算机sublime的build system文件夹的快捷方式
- `src/`源文件，如main.cpp、myClass.cpp等
    - `main.cpp`主程序
- `reference/`一些参考代码
    - `model.cpp`Petsc程序参考模板（不用管）
    - `old.cpp`原来的tokmark-2.cpp，在这里换了个名（不用管）
- `build/`编译过程中产生的链接文件，（不用管）
- `bin/`编译生成的可执行文件
- `data/`输入数据
    - `Bp_nova.dat`
    - ……等，还不清楚这是干嘛的，等知道了再写
- `output/`输出数据，这里暂时没有东西

## 开发环境（实验室电脑）
- 使用`~/fgn/sublime.sh`可打开**能使用中文**的sublime，在任意文件下使用`Ctrl + b`就可以编译运行整个工程
- 使用`atom`编写`readme.md`，点击工具栏`Packages -> Markdown Previewer`即可实时预览
- 在`atom`下用`ctrl+9`可以打开`github`
- 当使用远程ssh连接时，推荐使用XShell+XFtp（官网有6.0beta版，无需注册与付费）

## 坐标系统
#### 环坐标系
- `t` 在截面圆上的极角，同极坐标系下的θ
- `r` 在截面圆上的半径，同极坐标系下的ρ
- `p` 环形托卡马克的第几个截面，截面是个圆
- 环坐标系在项目中都以`toco`表示，是圆环坐标系TOroidal COordinates的缩写
- 注意这里的环坐标系与维基百科上的圆环坐标系体系不同

#### 柱坐标系
- 这里我们先把托卡马克想象成放在水平桌子上的一个甜甜圈
- `x` 把一根筷子竖直立在桌子上，穿过甜甜圈的中心，r是点到筷子的最短距离
- `y` 环形托卡马克的第几个截面
- `z` 用刀把甜甜圈砍成上下一样大的两半，z是点到刀切面的距离（有正负）
- 柱坐标系在项目中都以`grid`表示，意为网格
