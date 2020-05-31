#pragma once
#include "MeshCache.h"
#include "Auxiliary.h"

class MeshCut
{
public:
	MeshCut(Mesh& mesh, MeshCache& MCache);
	~MeshCut();
	void Set(const std::vector<int>& lmk, const std::vector<int>& initseam = std::vector<int>());
	void Connect();
	void MakeSeam();
	Mesh GetCutedMesh() const { return cuted_mesh; }
	std::vector<int>& GetCutvertex() { return cutVertex; }
	std::vector<int>& GetCutedge() { return cutEdge; }
	void SetBanCondition(const std::vector<int>& banB, const std::vector<int>& banE, const std::string BanMethod);
	void CalcBanArea(int Dn, std::string Metric, double alpha, double decay);
	std::vector<int>& GetMaxConnectedRegion() { return  MaxConnectedRegion; }
	void SetCut(std::vector<int>& cV, std::vector<int> cE);
	void GetCorrespondence(std::vector<int>& he2idx, std::vector<int>& idx2meshvid);

private:
	const Mesh& oriMesh;
	MeshCache& MC;
	Mesh cuted_mesh;
	std::vector<int> lmk;
	std::vector<int> initseam;
	std::vector<int> cutVertex;
	std::vector<int> cutEdge;

	std::vector<int> ban_vertex;
	std::vector<int> ban_edge;

	std::vector<int> he_to_idx; // to vertex index of a halfedge
	std::vector<int> idx_to_mesh_idx; // vertex indices correspondence

	std::vector<int> MaxConnectedRegion;
};