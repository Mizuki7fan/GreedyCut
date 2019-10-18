#pragma once
#include "..//MeshDefinition/MeshDefinition.h"
#include "..//MeshDefinition/MeshCache.h"

class TraverseVertexFinding
{
public:
	TraverseVertexFinding(const Mesh& mesh,const Mesh& paraed_mesh,MeshCache& MC);
	~TraverseVertexFinding();
	
	void ComputeDistortion();
	void PrepareComputeDistortion(void);
	void CalcVertexDistortion();
	void FindLandmarkPoints(std::vector<std::pair<int, int>>&);

private:
	const Mesh& origin_mesh;
	const Mesh& paraed_mesh;
	MeshCache& MC;
	std::vector<double> facearea, fpx1, fpx2, fpy2; // for distortion
	std::vector<double> vertex_distortion;
	std::vector<double> facedistortion; // updated in ComputeDistortion()
	std::vector<std::pair<int, double>> featurepoints; // result feature points
	int exe_id;
	std::string m_name;
};