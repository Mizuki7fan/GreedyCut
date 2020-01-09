#include "MeshCut.h"
#include <fstream>
#include <string>

MeshCut::MeshCut(Mesh &mesh, MeshCache &MCache) : orimesh(mesh), MCache(MCache)
{
}

void MeshCut::SetCondition(const std::vector<int> &lmk, const std::vector<int> &initseam)
{
	this->lmk = lmk;
	this->initseam = initseam;
}

void MeshCut::SetBanCondition(const std::vector<int> &banV, const std::vector<int> &banE, const std::string BanMethod)
{
	//设置Ban选的相关条件,根据BanMethod进行设置
	if (BanMethod == "NonConnect")
	{
		//不进行连接，直接将所有点作为输入
		std::set<int> ban;
		for (auto a : banV)
			ban.insert(a);
		for (auto a : banE)
			ban.insert(a);
		for (auto a : ban)
			ban_vertex.push_back(a);
	}
	else if (BanMethod == "Connect")
	{
		std::vector<int> lmk = banV;
		std::vector<int> cutVertex = banE;
		std::vector<int> cutEdge;

		int add_count = lmk.size();
		int min_sp, min_ep;
		double min_dis;
		std::vector<int> considered(MCache.n_vertices, 0);
		std::vector<int> is_marked(MCache.n_vertices, 0);
		for (int u = 0; u < cutVertex.size(); u++)
			is_marked[cutVertex[u]] = 1;
		while (add_count!=0)
		{
			int min_sp = -1, min_ep = -1;
			double min_dis = DBL_MAX;
			for (int u = 0; u < lmk.size(); u++)
			{
				if (considered[lmk[u]] == 1)
					continue;
				int s_p = lmk[u];
				if (MCache.V_D[s_p].size() == 0)
				{
					MCache.V_D[s_p].resize(MCache.n_vertices, DBL_MAX);
					continue;
				}
				for (int s = 0; s < cutVertex.size(); s++)
				{
					if (MCache.V_D[s_p][cutVertex[s]] < min_dis)
					{
						min_sp = s_p;
						min_ep = cutVertex[s];
						min_dis = MCache.V_D[s_p][cutVertex[s]];
					}
				}
			}
			if (min_sp != -1)
			{
				while (1)
				{
					int tmp = min_ep;
					min_ep = MCache.V_VP[min_sp][min_ep];
					if (min_ep == tmp)
						break;
					cutVertex.push_back(min_ep);
					cutEdge.push_back(MCache.vve[tmp][min_ep]);
				}
				considered[min_sp] = 1;
				add_count--;
				continue;
			}
			else
			{
				for (int i = 0; i < lmk.size(); i++)
				{
					if (considered[lmk[i]] == 1)
						continue;
					std::vector<int>& is_visited = MCache.dijkstra_isvisited[lmk[i]];
					std::vector<double>& distance = MCache.V_D[lmk[i]];
					std::priority_queue<node>& que = MCache.dijkstra_cache[lmk[i]];
					std::vector<int>& v_p = MCache.V_VP[lmk[i]];
					if (is_visited.size() == 0)
					{
						is_visited.resize(MCache.n_vertices, 0);
						distance.resize(MCache.n_vertices, DBL_MAX);
						distance[lmk[i]] = 0;
						que.push(node(lmk[i], 0));
						v_p.resize(MCache.n_vertices, -1);
						v_p[lmk[i]] = lmk[i];
					}
					if (que.empty())
					{
						que.push(node(lmk[i], 0));
					}
					while (!que.empty())
					{
						node tmp = que.top();
						que.pop();
						if (is_visited[tmp.id] == 1)
							continue;
						if (is_marked[tmp.id] == 1)
						{
							break;
						}
						is_visited[tmp.id] = 1;
						for (int s = 0; s < MCache.vv[tmp.id].size(); s++)
						{
							int v_id = MCache.vv[tmp.id][s];
							int e_id = MCache.ve[tmp.id][s];
							if (distance[tmp.id] + MCache.el[e_id] < distance[v_id])
							{
								distance[v_id] = distance[tmp.id] + MCache.el[e_id];
								que.push(node(v_id, distance[v_id]));
								v_p[v_id] = tmp.id;
							}
						}
					}
				}
			}
		}
		ban_vertex = cutVertex;
	}
	std::ofstream ban("banV.txt");
	for (auto a : ban_vertex)
		ban << a << std::endl;
	ban.close();
	
}

void MeshCut::SetCut(std::vector<int> cV, std::vector<int> cE)
{
	this->cutVertex = cV;
	this->cutEdge = cE;
}

void MeshCut::CalcBanArea(int Dn, std::string metric, double Alpha,double DecayRate)
{ //根据现有的区域来算Dn
	//第二个是优先级的度量
	std::vector<int> BanV(MCache.n_vertices, 0);
	std::vector<int> BanE(MCache.n_edges, 0);
	std::vector<int> BanF(MCache.n_faces, 0);

	while (1)
	{
		BanV.clear(); BanV.resize(MCache.n_vertices, 0);
		BanE.clear(); BanE.resize(MCache.n_edges, 0);
		BanF.clear(); BanF.resize(MCache.n_faces, 0);
		if (metric == "RealDis")
		{ //表示算禁止区域是按照实际的网格距离来计算,这里可以并行
			double range = Dn * MCache.avg_el;
			double total_length = MCache.avg_el * MCache.n_vertices;
#pragma omp parallel for
			for (int u = 0; u < ban_vertex.size(); u++)
			{
				int s_p = ban_vertex[u];
				std::vector<double> &distance = MCache.V_D[s_p];
				if (distance.size() != 0)
				{
					//如果距离数组不为0则表示之前求过，那么
					double max_dis = DBL_MIN;
					for (int vid = 0; vid < distance.size(); vid++)
					{
						if (distance[vid] > total_length)
							continue;
						if (distance[vid] > max_dis)
							max_dis = distance[vid];
						if (distance[vid] < range)
							BanV[vid] = 1;
					}
					if (max_dis > range)
						continue;
				}
				std::vector<int> &is_visited = MCache.dijkstra_isvisited[s_p];
				std::priority_queue<node> &que = MCache.dijkstra_cache[s_p];
				std::vector<int> &v_p = MCache.V_VP[s_p];
				if (distance.size() == 0)
				{
					is_visited.resize(MCache.n_vertices, 0);
					distance.resize(MCache.n_vertices, DBL_MAX);
					distance[s_p] = 0;
					que.push(node(s_p, 0));
					v_p.resize(MCache.n_vertices, -1);
					v_p[s_p] = s_p;
				}
				while (1)
				{
					node tmp = que.top();
					que.pop();
					if (is_visited[tmp.id])
						continue;
					is_visited[tmp.id] = 1;
					for (int u = 0; u < MCache.vv[tmp.id].size(); u++)
					{
						int vid = MCache.vv[tmp.id][u];
						int eid = MCache.ve[tmp.id][u];
						if (distance[tmp.id] + MCache.el[eid] < distance[vid])
						{
							distance[vid] = distance[tmp.id] + MCache.el[eid];
							que.push(node(vid, distance[vid]));
							v_p[vid] = tmp.id;
						}
					}
					if (distance[tmp.id] > range)
						break;
					else
					{
						BanV[tmp.id] = 1;
					}
				}
			}
		}
		else if (metric == "Neighbourhood")
		{
#pragma omp parallel for
			for (int u = 0; u < ban_vertex.size(); u++)
			{
				int a = ban_vertex[u];
				Algorithm::UpdateNeighbourhood(MCache, Dn, a);
				for (int i = 0; i < MCache.n_vertices; i++)
				{
					if (MCache.Neighbour[a][i] <= Dn)
						BanV[i] = 1;
				}
			}
		}
		for (int i = 0; i < MCache.n_edges; i++)
		{
			int v1 = MCache.ev[i][0];
			int v2 = MCache.ev[i][1];
			if (BanV[v1] == 1 || BanV[v2] == 1)
				BanE[i] = 1;
		}
		std::vector<int> v_camp = BanV;
		int max_region = 2;
		std::vector<double> weight(MCache.n_edges, DBL_MAX);
		for (int i = 0; i < MCache.n_edges; i++)
		{
			if (BanE[i] == 0)
				weight[i] = MCache.el[i];
		}
		for (int i = 0; i < MCache.n_vertices; i++)
		{
			if (v_camp[i] != 0)
				continue;
			std::vector<double> d(MCache.n_vertices, DBL_MAX);
			std::vector<int> v_i(MCache.n_vertices, 0);
			Algorithm::Dijkstra_with_restrict(MCache, i, weight, v_i, d);
			double total_length = MCache.avg_el * MCache.n_edges;
			for (int u = 0; u < d.size(); u++)
			{
				if (d[u] < total_length)
					v_camp[u] = max_region;
			}
			max_region++;
		}
		std::vector<double> area_region(max_region, 0);
		//遍历所有的面，统计各个区域的面积
		for (int i = 0; i < MCache.n_faces; i++)
		{
			int v1 = MCache.fv[i][0];
			int v2 = MCache.fv[i][1];
			int v3 = MCache.fv[i][2];
			//同一个面的3个点不可能属于2个不同的阵营
			int f_region = 1;
			if (v_camp[v1] != 1)
				f_region = v_camp[v1];
			if (v_camp[v2] != 1)
				f_region = v_camp[v2];
			if (v_camp[v3] != 1)
				f_region = v_camp[v3];
			area_region[f_region] += MCache.fa[i];
		}
		//area_region记录了每个阵营的面积，选取最大的判断是否要回退
		int max_region_id = -1;
		double max_region_area = DBL_MIN;
		for (int i = 2; i < area_region.size(); i++)
		{
			if (area_region[i] > max_region_area)
			{
				max_region_area = area_region[i];
				max_region_id = i;
			}
		}
		if (max_region_area < Alpha * MCache.avg_fa * MCache.n_faces)
		{
			Dn = std::floor(DecayRate* Dn);
			std::cout << "最大连通区域，回退" << std::endl;
		}
		else
		{
			MaxConnectedRegion.clear();
			MaxConnectedRegion.resize(MCache.n_vertices);
			for (int vid = 0; vid < v_camp.size(); vid++)
			{
				if (v_camp[vid] == max_region_id)
					MaxConnectedRegion[vid] = 1;
				else
					MaxConnectedRegion[vid] = 0;
			}
			break;
		}
	}

}

void MeshCut::Connect()
{
	//如果初始没有割缝，则直接用Kruskal算法算割缝
	if (initseam.size() == 0)
	{
		//求lmk数组中两两点的距离，获得完全图
		Algorithm::Dijkstra(MCache, lmk);
		//之后根据求出的完全图
		Algorithm::Kruskal(MCache, lmk, cutVertex, cutEdge);
	}
}

void MeshCut::MakeSeam()
{
	cuted_mesh.clean();
	std::vector<int> candidate_cut_valence(MCache.n_vertices, 0);
	std::vector<int> candidate_seam_edge(MCache.n_edges, 0);
	std::vector<int> candidate_seam_vertex(MCache.n_vertices, 0);
	for (auto a : cutEdge)
		candidate_seam_edge[a] = 1;
	for (auto a : cutVertex)
		candidate_seam_vertex[a] = 1;
	for (int i = 0; i < cutEdge.size(); i++)
	{
		Mesh::HalfedgeHandle he = orimesh.halfedge_handle(orimesh.edge_handle(cutEdge[i]), 0);
		candidate_cut_valence[orimesh.to_vertex_handle(he).idx()]++;
		candidate_cut_valence[orimesh.from_vertex_handle(he).idx()]++;
	}
	he_to_idx.resize(2 * MCache.n_edges, 0);
	idx_to_mesh_idx.clear();
	Mesh::HalfedgeHandle h_begin;
	for (int i = 0; i < candidate_seam_edge.size(); i++)
	{
		if (candidate_seam_edge[i])
		{
			h_begin = orimesh.halfedge_handle(orimesh.edge_handle(i), 0);
			break;
		}
	}
	auto h_iter = h_begin;
	int uv_idx = 0;
	do
	{
		he_to_idx[h_iter.idx()] = uv_idx;
		idx_to_mesh_idx.push_back(orimesh.to_vertex_handle(h_iter).idx());
		h_iter = orimesh.next_halfedge_handle(h_iter);
		while (!candidate_seam_edge[orimesh.edge_handle(h_iter).idx()])
		{
			h_iter = orimesh.opposite_halfedge_handle(h_iter);
			he_to_idx[h_iter.idx()] = uv_idx;
			h_iter = orimesh.next_halfedge_handle(h_iter);
		}
		uv_idx++;
	} while (h_iter != h_begin);
	for (const auto &vh : orimesh.vertices())
	{
		if (candidate_cut_valence[vh.idx()] > 0)
			continue;
		for (const auto &viheh : orimesh.vih_range(vh))
		{
			he_to_idx[viheh.idx()] = uv_idx;
		}
		idx_to_mesh_idx.push_back(vh.idx());
		uv_idx++;
	}
	for (int i = 0; i < uv_idx; i++)
	{
		cuted_mesh.add_vertex(orimesh.point(orimesh.vertex_handle(idx_to_mesh_idx[i])));
	}
	for (const auto &fh : orimesh.faces())
	{
		std::vector<Mesh::VertexHandle> face_vhandles;
		for (const auto &fheh : orimesh.fh_range(fh))
		{
			face_vhandles.push_back(cuted_mesh.vertex_handle(he_to_idx[fheh.idx()]));
		}
		cuted_mesh.add_face(face_vhandles);
	}
}

void MeshCut::get_correspondence(std::vector<int> &he2idx, std::vector<int> &idx2meshvid)
{
	he2idx = this->he_to_idx;
	idx2meshvid = this->idx_to_mesh_idx;
}
