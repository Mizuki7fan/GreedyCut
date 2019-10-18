#pragma once
#include "../MeshDefinition/MeshDefinition.h"
#include "../MeshDefinition/MeshCache.h"
#include "../Auxiliary/Algorithm.h"
#include <omp.h>
#include <iostream>



class GreedCut
{
public:
	GreedCut(Mesh& mesh, MeshCache& MC, std::vector<int>& landmark,double trimming_rate);
	void GenGraph();
	double CalcTreeLength(std::priority_queue<PathInfo>GenTreePath, std::vector<int> k);
	void OP();
	void Greed1_Op();
	double CalcExtraLandmark(int p);
	void GetResult(std::vector<int>& k) { k = landmark; };
private:
	Mesh& mesh;
	MeshCache& MC;
	std::vector<int> landmark;
	std::vector<int> back_landmark;
	std::priority_queue<PathInfo> PresavedPath;
	double tree_length;
	bool greed1_ok;
	double trimming_rate=0.01;
};