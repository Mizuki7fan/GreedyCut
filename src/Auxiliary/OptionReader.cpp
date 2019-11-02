#include "OptionReader.h"

Option::Option(std::string s,std::string ModelPath)
:ModelPath(ModelPath)
{
	//获取网格名称
	ModelName = ModelPath.substr(ModelPath.find_last_of('\\') + 1, ModelPath.length());
	
	std::ifstream input(s);
	std::string line,key,value;
	int position;
	while (std::getline(input, line))
	{
		if (line[0] == '#')
			continue;
		position = line.find("=");
		if (position == line.npos)
			continue;
		key = line.substr(0, position);
		value = line.substr(position + 1, line.length());
		if (key == "OutputDir")
			OutputDir = value + ModelName.substr(0, ModelName.length() - 4);
		else if (key == "SampleMethod")
			SampleMethod = value;
		else if (key == "SampleFirstPoint")
			SampleFirstPoint = value;
		else if (key == "MeshcutOutput")
			MeshcutOutput = value;
		else if (key == "KPNewtonResultOutput")
			KPNewtonOutput = value;
		else if (key == "VertexPriorityMetric")
			VertexPriorityMetric = value;
		else if (key == "Dn")
			Dn = std::stoi(value);
		else if (key == "BanAreaMethod")
			BanAreaMethod = value;

	}
}