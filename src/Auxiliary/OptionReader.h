#pragma once
#include <string>
#include <fstream>

class Option
{
public:
	//根据网格路径和option文件的名字来获得相应的路径名字
	Option(std::string,std::string);
	
	//输出目录
	std::string ModelName;
	std::string ModelPath;

	std::string OutputDir = "../Output/";
	std::string SampleMethod = "Dijkstra";
	std::string SampleFirstPoint = "Random";
	std::string MeshcutOutput = "Yes";
	std::string KPNewtonOutput = "Yes";
	std::string VertexPriorityMetric = "RealDis";
	int Dn = 10;
	std::string BanAreaMethod = "NonConnect";


	double filtering_rate;
	double trimming_rate;

};
