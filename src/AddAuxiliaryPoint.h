#pragma once
#include "MeshCache.h"
#include "Auxiliary.h"

class AAP
{
public:
	AAP(Mesh& mesh, MeshCache& MC, std::vector<int>& landmark);
	void Set(double TrimmingRate, int MaxAddCount);
	void Run();
	std::vector<int> GetResult() { return landmark; }
	~AAP();

private:

	void GenGraph();
	void Greed1_Op();

	Mesh& mesh;
	MeshCache& MC;
	std::vector<int> landmark;
	double TrimmingRate = 0.01;
	std::priority_queue<PathInfo> PresavedPath;
	double TreeLength;
	int MaxAddCount = 30;
	bool greed1_ok = false;
};
