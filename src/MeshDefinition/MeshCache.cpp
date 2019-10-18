#include "MeshCache.h"

//����MeshCache��
MeshCache::MeshCache(Mesh& mesh)
{	
	n_vertices = mesh.n_vertices();
	n_edges = mesh.n_edges();
	dijkstra_cache.resize(n_vertices);
	dijkstra_isvisited.resize(n_vertices);
	V_VP.resize(n_vertices);
	V_D.resize(n_vertices);

	vv.resize(n_vertices, std::vector<int>());
	ve.resize(n_vertices, std::vector<int>());
	vve.resize(n_vertices);
	Vx.resize(n_vertices);
	Vy.resize(n_vertices);
	Vz.resize(n_vertices);
	for (const auto& v : mesh.vertices())
	{
		Mesh::Point p = mesh.point(v);
		Vx[v.idx()] = p.data()[0];
		Vy[v.idx()] = p.data()[1];
		Vz[v.idx()] = p.data()[2];

		for (const auto& viheh : mesh.vih_range(v))
		{
			Mesh::VertexHandle v1 = mesh.from_vertex_handle(viheh);
			Mesh::EdgeHandle e = mesh.edge_handle(viheh);
			vv[v.idx()].push_back(v1.idx());
			ve[v.idx()].push_back(e.idx());
			vve[v.idx()][v1.idx()] = e.idx();
		}
	}
	el.resize(n_edges);
	ev.resize(n_edges, std::vector<int>(2));
	for (const auto& e : mesh.edges())
	{
		Mesh::HalfedgeHandle he = mesh.halfedge_handle(e, 0);
		ev[e.idx()][0] = mesh.to_vertex_handle(he).idx();
		ev[e.idx()][1] = mesh.from_vertex_handle(he).idx();
		el[e.idx()] = mesh.calc_edge_length(e);
	}
}

MeshCache::~MeshCache()
{
}
void MeshCache::updataCapacity()
{
	capacity = 0;
	for (int i = 0; i < n_vertices; i++)
	{
		capacity += dijkstra_cache[i].size() * sizeof(node);
		capacity += dijkstra_isvisited[i].size() * sizeof(int);
		capacity += V_VP[i].size() * sizeof(int);
		capacity += V_D[i].size() * sizeof(double);
	}
	std::cout << "capacity : " <<capacity<< std::endl;
}
/*
void MeshCache::UpdateNeighbourInfo(int k)
{
	Neighbour.resize(n_vertices, std::vector<std::vector<int>>(k + 1, std::vector<int>()));
	std::vector<std::vector<int>> v_v_dis(n_vertices, std::vector<int>(n_vertices, -1));
#pragma omp parallel for num_threads(16)
	for (int i = 0; i < n_vertices; i++)
	{
		Neighbour[i][0].push_back(i);
		Neighbour[i][1] = vv[i];
		v_v_dis[i][i] = 0;
		for (auto a : vv[i])
			v_v_dis[i][a] = 1;
	}
	for (int u = 2; u <= k; u++)
	{
#pragma omp parallel for num_threads(16)
		for (int i = 0; i < n_vertices; i++)
		{
			for (auto a : Neighbour[i][u - 1])
			{
				for (auto b : vv[a])
				{
					if (v_v_dis[i][b]==-1)
					{
						v_v_dis[i][b] = u;
						Neighbour[i][u].push_back(b);
					}
				}
			}
		}
	}
}
*/