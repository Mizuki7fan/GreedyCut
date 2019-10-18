#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "..//MeshDefinition/MeshDefinition.h"
#include "..//MeshDefinition/MeshCache.h"
#include "../Auxiliary/Algorithm.h"

#include <vector>
#include <iostream>

class MeshCut
{
public:
	//�����и��࣬������ԭʼ�����Լ���Cache
	MeshCut(Mesh& mesh, MeshCache& MCache);
	//����landmark������еĸ���
	void SetCondition(const std::vector<int>& lmk,const std::vector<int>& initseam=std::vector<int>());
	//��������ĸ��
	void Connect();
	//Cut_to_Seam
	void MakeSeam();
	//�����п��������
	Mesh GetCutedMesh() const { return cuted_mesh; };
	std::vector<int>& GetCutvertex() { return cutVertex; }
	std::vector<int>& GetCutEdge() { return cutEdge; }
	void CalcForbiddenArea();
	void get_correspondence(std::vector<int>& he2idx, std::vector<int>& idx2meshvid);

private:
	const Mesh& orimesh;
	MeshCache& MCache;
	Mesh cuted_mesh;
	std::vector<int> lmk;
	std::vector<int> initseam;
	std::vector<int> cutVertex;
	std::vector<int> cutEdge;
	std::vector<int> ban_vertex;
	std::vector<int> ban_edge;

	std::vector<int> he_to_idx; // to vertex index of a halfedge
	std::vector<int> idx_to_mesh_idx; // vertex indices correspondence
}; 