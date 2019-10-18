#pragma once
#include "../MeshDefinition/MeshCache.h"
#include "../Auxiliary/Algorithm.h"


class PointSampling
{
public:
	PointSampling(MeshCache & m,int FirstPoint,int n_samples,int method_id);
	~PointSampling(void);

	std::vector<int> ComputeSamples();

private:
	// Compute distance from a point to a set

	double PointSetRealDistance(const int & vid, const std::vector<int>& pointset) const;
	double PointSetPathLength(const int& vid, const std::vector<int>& pointset) const;


	MeshCache & MCache; // const reference to the mesh
	int method_id;//1-RealDistance 2-GeodesticDistance 3-LargestEdgeLength
	int n_samples;
	int FirstPoint;

};

