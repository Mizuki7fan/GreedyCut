#pragma once
#include "../MeshDefinition/MeshDefinition.h"
#include "../MeshDefinition/MeshCache.h"
#include "../Auxiliary/Algorithm.h"

class AAP
{
public:
	AAP(Mesh & mesh,MeshCache& MC,std::vector<int>& landmark);
	void GenGraph();
	void Op();
	void Greed1_Op();
	std::vector<int> GetResult();
	void Set(int TR,int AddCount,int ParaCount);
	~AAP();

private:
	Mesh& mesh;
	MeshCache& MC;
	std::vector<int> landmark;
	double TrimmingRate=0.01;
	std::priority_queue<PathInfo> PresavedPath;
	//std::vector<PathInfo> PresavedPath;
	double TreeLength;
	int MaxAddCount = 30;
	int ParaCount = 6;
	bool greed1_ok = false;

};