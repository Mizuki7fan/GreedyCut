#include <iostream>
#include "Auxiliary/OptionReader.h"
#include "MeshDefinition/MeshDefinition.h"
#include "MeshDefinition/MeshCache.h"
#include "PointSampling/PointSampling.h"
#include "MeshCut/MeshCut.h"
//MinGW和MSVC都需要包含这两个
#ifdef _WIN64
#include <direct.h>
#include <io.h>
#endif
//MINGW需要包含这些
#ifdef __MINGW64__
#define OM_STATIC_BUILD 1
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>
#endif


//主函数
int main(int argc, char* argv[])
{
#ifdef __MINGW64__
	std::cout << "The Complier is MinGW" << std::endl;	
	OpenMesh::IO::_OBJReader_();
	OpenMesh::IO::_OBJWriter_();
#elif _MSC_VER
	std::cout << "The Complier is MSVC" << std::endl;
#endif // __GNUC

	//网格名称需要单独读取，为了批量化的测试
	std::string ModelPath = argv[1];
	//读取配置文件
	Option opt(argv[2],ModelPath);
	std::cout << "start processing [" << opt.ModelName << "] ......" << std::endl;

	//创建保存运行结果的文件夹
#ifdef _WIN64
	if (_access(opt.OutputDir.c_str(),0)==-1)
		int i = _mkdir(opt.OutputDir.c_str());
#endif // _WIN64

	//读网格，新建备份类
	Mesh mesh;
	OpenMesh::IO::read_mesh(mesh, argv[1]);
	if (mesh.n_vertices() == 0)
	{
		std::cerr << "网格顶点数为0" << std::endl;
		return 0;
	}
	MeshCache MC(mesh);
	//std::vector<int> test_lmk = { 1,567,345,8762,10001,55,6,5506,2365,6988,2359,11110,1234 };
	//Algorithm::Dijkstra(MC, test_lmk);
	//MC.updataCapacity();
	//return 1;
	std::cout << "Find Cut 1" << std::endl;
	//生成第一条割缝，根据Option中的值来确定
	//三种最远点采样的方法
	switch (opt.FirstcutFirstPoint)
	{
	case 1:
		std::cout << "FirstcutMethod:ReadDistance" << std::endl;
	case 2:
		std::cout << "FirstcutMethod:GeodesicDistance" << std::endl;
	case 3:
		std::cout << "FirstcutMethod:LargestEdgeLength" << std::endl;
	default:
		std::cerr << "FirstcutMethod:Wrong" << std::endl;
		return -1;
		break;
	}
	PointSampling ps(MC, 0, 10, opt.FirstcutMethod);
	std::vector<int> SamplePoints= ps.ComputeSamples();
	//初始的一刀
	MeshCut MCut(mesh, MC);
	//传入采样点
	MCut.SetCondition(SamplePoints);
	//算Cut
	MCut.CalcCut();
	//








}