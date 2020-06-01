# GreedyCut

This is a C++ implement of the following paper:

[Tianyu Zhu](https://Mizuki7fan.github.io), [Chunyang Ye](https://chunyangye.github.io), [Shuangming Chai](https://https://kfckfckf.github.io), and [Xiao-Ming Fu](http://staff.ustc.edu.cn/~fuxm).
Greedy cut construction for parameterizations.
Computer Graphics Forum (Eurographics), 2020.

The code is written by Tianyu Zhu using Microsoft Visual Studio 2019 with OpenMP, and has passed the test on Windows 10, Ubuntu 18.04 and macOS 10.14.3.

## External Libraries
---
All the required libraries are the latest versions.

* [Eigen](http://eigen.tuxfamily.org/)
* [OpenMesh](https://www.openmesh.org/)
* [PARDISO](https://www.pardiso-project.org/) (academic license) or [MKL Pardiso Solver](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html)
* OpenMP

## Clone Repository
---
```
git clone https://github.com/Mizuki7fan/GreedyCut
```

## Build and Run
---
### build
```
cd GreedyCut
```
**Edit CmakeLists.txt** to set the solver and configuration path.

```
mkdir build && cd build
cmake ..
```

### run
---
```
GreedyCut.exe <mesh> opt.txt
```
Arguments in opt.txt:
* PointSampling_method: the method to find second vertex.
* PointFinding_method: the metric of distance in PointFinding stage.
* BanArea_Method: weather to connect local maximizer before finding second cut.
* BanArea_Dn: the radius of forbidden region.
* BanArea_Alpha: the threshold in dual cut strategy.
* BanArea_Metric: the metric of distance of the radius of forbidden region.
* BanArea_ShrinkRate: the shrink rate of radius of forbidden region.
* Influence_Threshold: the scope of influence threshold.
* Distortion_Threshold: the distortion threshold of filtering stage.
* Trimming_Rate: the length threshold of add auxiliary point stage.
* Max_AddCount: the maximum of adding auxiliary point.

### example
```
GreedyCut.exe alien.obj opt.txt
```

## Output Files
* landmark.txt: the landmark find by algorithm.
* cut.txt: the cut find by algorithm
* cuted_mesh.obj: the open mesh cuted by the cut.