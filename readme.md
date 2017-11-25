# 说明文档

## 更新日志

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
- 当使用远程ssh连接时，推荐使用XShell+XFtp（官网有6.0beta版，无需注册与付费）