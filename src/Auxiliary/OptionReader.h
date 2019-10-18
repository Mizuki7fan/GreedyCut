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
	int FirstcutMethod = 0;
	int FirstcutFirstPoint = 0;
	double filtering_rate;
	double trimming_rate;
};
