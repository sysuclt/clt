# 说明文档
### 文件结构
- `makefile`可以编译整个工程的makefile
- `readme.md`有问题，就看这里
- `toBuildSys/`一个可以快去去到本计算机sublime的build system文件夹的快捷方式
- `src/`源文件，如main.cpp、myClass.cpp等
    - `main.cpp`主程序
- ‘reference/’一些参考代码
    - `model.cpp`Petsc程序参考模板（不用管）
    - `old.cpp`原来的tokmark-2.cpp，在这里换了个名（不用管）
- `build/`编译过程中产生的链接文件，（不用管）
- `bin/`编译生成的可执行文件
- `data/`输入数据
    - `Bp_nova.dat`
    - ……等，还不清楚这是干嘛的，等知道了再写
- `output/`输出数据，这里暂时没有东西

### 开发环境（实验室电脑）
- 使用`~/fgn/sublime.sh`可打开**能使用中文**的sublime，在任意文件下使用`Ctrl + b`就可以编译运行整个工程
- 使用`atom`编写`readme.md`，点击工具栏`Packages -> Markdown Previewer`即可实时预览

### 有待确认的问题
- 确认程序流程是否有问题
- 了解环坐标体系，缺失网格坐标范围
- qs2grd 还要计算导数？这个函数的功能十分迷，看不懂参数
