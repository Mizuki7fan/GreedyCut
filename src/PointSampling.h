#pragma once
#include "MeshCache.h"
#include "Auxiliary.h"

class PointSampling
{
public:
	PointSampling(MeshCache& MC);
	~PointSampling();
	void Set(std::string method_id);
	void ComputeSamples(std::vector<int>& Result);
	void ComputeSamplesFromSelectedRegion(std::string method, std::vector<int>& Region, std::vector<int>& SampleVertices, std::vector<int>& SampleEdges);

private:
	double PointSetRealDistance(const int& vid, const std::vector<int>& pointset) const;
	double PointSetPathLength(const int& vid, const std::vector<int>& pointset) const;

	MeshCache& MCache;
	std::string method;
	int first_point;
	int n_samples = 2;
};
