#include "MeshCache.h"

//构造MeshCache类
MeshCache::MeshCache(Mesh& mesh)
{	
	n_vertices = mesh.n_vertices();
	n_edges = mesh.n_edges();
	n_faces = mesh.n_faces();
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
	avg_el = 0;
	for (const auto& e : mesh.edges())
	{
		Mesh::HalfedgeHandle he = mesh.halfedge_handle(e, 0);
		ev[e.idx()][0] = mesh.to_vertex_handle(he).idx();
		ev[e.idx()][1] = mesh.from_vertex_handle(he).idx();
		el[e.idx()] = mesh.calc_edge_length(e);
		avg_el += el[e.idx()];
	}
	avg_el /= n_edges;
	//预存邻域信息
	Neighbour.resize(n_vertices);
	//记录顶点已经求取的最大邻域值
	Max_Neighbour.resize(n_vertices, 0);

	fa.resize(n_faces);
	fv.resize(n_faces, std::vector<int>(3));
	avg_fa = 0;
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		const auto& f = mesh.face_handle(i);
		auto heh = mesh.halfedge_handle(f);
		const auto& p0 = mesh.point(mesh.from_vertex_handle(heh));
		const auto& p1 = mesh.point(mesh.to_vertex_handle(heh));
		const auto& p2 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)));
		fa[i] = 0.5 * ((p1 - p0) % (p2 - p0)).norm();
		avg_fa += fa[i];
		Mesh::FaceVertexIter fvh = mesh.fv_begin(f);
		fv[i][0] = fvh->idx();
		fvh++;
		fv[i][1] = fvh->idx();
		fvh++;
		fv[i][2] = fvh->idx();
	}

	avg_fa /= n_faces;


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