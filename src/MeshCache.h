#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <string>
#include <queue>

struct MeshTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;
	VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	EdgeAttributes(OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);
};

typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> Mesh;

struct node
{
	int id;
	double dis;
	node(int id, double d)
	{
		this->id = id;
		this->dis = d;
	}
	bool operator<(const node& rhs) const { return dis > rhs.dis; }
};
struct PathInfo
{
	int s_p, e_p;
	double length;
	std::vector<int> path;
	PathInfo(int a, int b, double c)
	{
		s_p = a;
		e_p = b;
		length = c;
	}
	bool operator<(const PathInfo& rhs) const
	{
		return length > rhs.length;
	}
};

class MeshTools
{
public:
	static double Area(const Mesh& mesh);
};

class MeshCache
{
public:
	MeshCache(Mesh &mesh);
	~MeshCache();

	int NVertices=-1;
	int NEdges;
	int NFaces;
	std::vector<double> Vx, Vy, Vz;
	std::vector<std::vector<double>> Vd;
	std::vector<std::vector<int>> DijkstraIsVisited;
	std::vector<std::priority_queue<node>> DijkstraCache;
	std::vector<std::vector<int>> VVp;

	std::vector<std::vector<int>> VV;
	std::vector<std::vector<int>> VE;
	std::vector<double> EL;
	std::vector<std::vector<int>> EV;
	std::vector<std::map<int, int>> VVE;
	double AVG_EL;

	std::vector<std::vector<int>> Neighbour;
	std::vector<int> Max_Neighbour;

	std::vector<double> FA;
	std::vector<std::vector<int>> FV;
	double AVG_FA;
	double Lbb;


};