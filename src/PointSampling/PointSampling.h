#pragma once
#include "../MeshDefinition/MeshCache.h"
#include "../Auxiliary/Algorithm.h"
#include <random>


class PointSampling
{
public:
	PointSampling(MeshCache & m);
	~PointSampling(void);
	void Set(std::string method,std::string FirstPoint,int n_samples);
	void ComputeSamples(std::vector<int>& Result);
	void ComputeSamplesFromSelectedRegion(std::string method, std::vector<int>& Region, std::vector<int>& SampleVertices,std::vector<int>& SampleEdges);

private:
	// Compute distance from a point to a set

	double PointSetRealDistance(const int & vid, const std::vector<int>& pointset) const;
	double PointSetPathLength(const int& vid, const std::vector<int>& pointset) const;


	MeshCache & MCache; // const reference to the mesh
	std::string method;//1-RealDistance 2-GeodesticDistance 3-LargestEdgeLength
	int n_samples;
	int FirstPoint;

};

