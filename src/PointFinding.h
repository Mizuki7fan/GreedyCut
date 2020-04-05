#pragma once
#include "MeshCache.h"
//用于从网格中选取顶点
class PointFinding
{
public:
	PointFinding(const Mesh&,const Mesh&,MeshCache& MCache);
	~PointFinding();
	void Set(std::string metric);
	//找结果的点
	void Find(std::vector<std::pair<int,double>>& result);
	void FindLocalMaximizer();
	std::vector<int> GetLocalMaximizer();

private:
	void PrepareComputeDistortion();
	void ComputeFaceDistortion();
	void ComputeVertexDistortion();
	void FindByRealDis(std::vector<std::pair<int,double>>&);
	void FindByNeighbourhood(std::vector<std::pair<int, double>>&);

	//输入：原始网格，参数化后网格，若干设置
	const Mesh& OriMesh;
	const Mesh& ParaedMesh;
	MeshCache& MC;
	std::string Metric;

	std::vector<double> facearea, fpx1, fpx2, fpy2; // for distortion
	std::vector<double> vertex_distortion;
	std::vector<double> facedistortion; // updated in ComputeDistortion()
	
	std::vector<int> LocalMaximizer;


};