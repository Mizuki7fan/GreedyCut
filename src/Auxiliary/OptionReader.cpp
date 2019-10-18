#include "OptionReader.h"

Option::Option(std::string s,std::string ModelPath)
:ModelPath(ModelPath)
{
	//��ȡ��������
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
			OutputDir = value + ModelName;
		else if (key == "FirstcutMethod")
			FirstcutMethod = std::stoi(value);
		else if (key == "FirstcutFirstPoint")
			FirstcutFirstPoint = std::stoi(value);
	}
}