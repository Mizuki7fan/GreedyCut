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
	std::string OutputDir = "Output/";
	int SampleMethod = 0;
	int SampleFirstPoint = 0;
	int MeshcutOutput = 1;
	int KPNewtonOutput = 0;

	double filtering_rate;
	double trimming_rate;

};
