#pragma once
#include "../MeshDefinition/MeshDefinition.h"
#include "../MeshDefinition/MeshCache.h"
#include "../Auxiliary/Algorithm.h"

class GAP
{
public:
	GAP(const Mesh& mesh,MeshCache& MCache);
	void Set(const std::vector<std::pair<int, double>>& FP, std::string PriorityMetric, int FixThreshold, int ForbiddenRadius);
	//计算将所有固定点链接而得到的初始切割
	void GenFirstCut();
	void ClassifyFeaturePoints(int threshold);

private:
	const Mesh& mesh;
	MeshCache MCache;
	//点的优先级的度量
	std::string VertexPriorityMetric = "Neighbourhood";
	int FixThreshold = 20;
	int GAPForBiddenRadius = 5;


	//备选的特征点
	std::vector<std::pair<int, double>> FeaturePoints;
	std::vector<int> FixPoints;
	std::vector<int> CanditatePoints;
};