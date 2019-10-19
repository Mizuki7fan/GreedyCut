#pragma once
#include <string>
#include <fstream>

class Option
{
public:
	//��������·����option�ļ��������������Ӧ��·������
	Option(std::string,std::string);
	
	//���Ŀ¼
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
