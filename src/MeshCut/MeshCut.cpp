#include "MeshCut.h"
#include <string>

MeshCut::MeshCut(Mesh &mesh, MeshCache &MCache) : 
	orimesh(mesh), MCache(MCache)
{
}

void MeshCut::SetCondition(const std::vector<int> &lmk, const std::vector<int> &initseam)
{
	this->lmk = lmk;
	this->initseam = initseam;
}

void MeshCut::CalcCut()
{
	//如果初始没有割缝，则直接用Kruskal算法算割缝
	if (initseam.size() == 0)
	{
		Algorithm::Dijkstra(MCache, lmk);
		Algorithm::Kruskal(MCache, lmk, cutVertex, cutEdge);
	}
	else
	{
		cutVertex = initseam;
		int add_count = lmk.size();
		int min_sp, min_ep;
		double min_dis;
		//??lmk??????
		std::vector<int> considered(MCache.n_vertices, 0);
		//???????????
		std::vector<int> is_marked(MCache.n_vertices, 0);
		for (int u = 0; u < cutVertex.size(); u++)
		{
			is_marked[cutVertex[u]] = 1;
		}
		while (add_count != 0)
		{
			int min_sp = -1;
			int min_ep = -1;
			double min_dis = DBL_MAX;
			//???landmark?????cutVertex????????????????????????cutVertex?

			for (int u = 0; u < lmk.size(); u++)
			{
				if (considered[lmk[u]] == 1)
					continue;
				int s_p = lmk[u];
				//????????????????
				if (MCache.V_D[s_p].size() == 0)
				{
					MCache.V_D[s_p].resize(MCache.n_vertices, DBL_MAX);
					continue;
				}
				//??????????
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
			{ //???????????path?
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
			else //???????????????????????????????????
			{
				for (int i = 0; i < lmk.size(); i++)
				{
					if (considered[lmk[i]] == 1)
						continue;
					std::vector<int> &is_visited = MCache.dijkstra_isvisited[lmk[i]];
					std::vector<double> &distance = MCache.V_D[lmk[i]];
					std::priority_queue<node> &que = MCache.dijkstra_cache[lmk[i]];
					std::vector<int> &v_p = MCache.V_VP[lmk[i]];
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

void MeshCut::CalcForbiddenArea()
{
	bool pass = false;
	int k = 10;
	std::vector<int> v_camp(MCache.n_vertices, -2);
	std::vector<int> visited(MCache.n_vertices, 0);
	int max_camp = 0;
	int max_v_camp = 0;
	int max_camp_count = -1;

	while (!pass)
	{
		v_camp.clear();
		v_camp.resize(MCache.n_vertices, -2);
		visited.clear();
		visited.resize(MCache.n_vertices, 0);
		ban_vertex.clear();
		ban_edge.clear();
		ban_vertex.resize(MCache.n_vertices, 0);
		ban_edge.resize(MCache.n_edges, 0);

		for (auto a : cutVertex)
			ban_vertex[a] = 1;
		for (int u = 0; u < k-1; u++)
		{
			std::vector<int> tmp=ban_vertex;
			for (int i = 0; i < MCache.n_vertices; i++)
			{
				//如果点i被ban了，那么就把i的1邻域都给ban掉
				if (ban_vertex[i] == 1)
				{
					for (auto a : MCache.vv[i])
						tmp[a] = 1;
				}
			}
			ban_vertex=tmp;
		}
/*旧版
		for (int u = 0; u < cutVertex.size(); u++)
		{
			for (int i = 0; i < k; i++)
			{
				for (int j = 0; j < MCache.Neighbour[cutVertex[u]][i].size(); j++)
					ban_vertex[MCache.Neighbour[cutVertex[u]][i][j]] = 1;
			}
		}
		*/
		for (int i = 0; i < MCache.n_edges; i++)
		{
			if (ban_vertex[MCache.ev[i][0]] == 1 && ban_vertex[MCache.ev[i][1]] == 1)
				ban_edge[i] = 1;
		}
		std::vector<int> result_vertex, result_edge;
		for (int i = 0; i < ban_vertex.size(); i++)
		{
			if (ban_vertex[i] == 1)
				result_vertex.push_back(i);
		}
		for (int i = 0; i < ban_edge.size(); i++)
		{
			if (ban_edge[i] == 1)
				result_edge.push_back(i);
		}
		max_camp = 0;
		for (int i = 0; i < MCache.n_vertices; i++)
		{
			if (ban_vertex[i] == 1)
			{
				v_camp[i] = -1;
				continue;
			}
			if (v_camp[i] == -2)
			{
				std::queue<int> q;
				q.push(i);
				while (!q.empty())
				{
					int tmp = q.front();
					q.pop();
					v_camp[tmp] = max_camp;
					for (int u = 0; u < MCache.vv[tmp].size(); u++)
					{
						int k = MCache.vv[tmp][u];
						if (ban_vertex[k] == 1 || visited[k] == 1)
							continue;
						else
						{
							visited[k] = 1;
							q.push(k);
						}
					}
				}
				max_camp++;
			}
		}
		std::vector<int> camp_count(max_camp, 0);
		for (int i = 0; i < MCache.n_vertices; i++)
		{
			if (v_camp[i] >= 0)
				camp_count[v_camp[i]]++;
		}
		max_v_camp = -1;
		max_camp_count = -1;
		for (int i = 0; i < camp_count.size(); i++)
		{
			if (camp_count[i] > max_camp_count)
			{
				max_v_camp = i;
				max_camp_count = camp_count[i];
			}
		}
		//?????
		float a = float(max_camp_count) / float(orimesh.n_vertices());
		if (a >= 0.05)
		{
			pass = true;
		}
		else
			k = k / 2;
	}
	int s_p = -1;
	for (int i = 0; i < MCache.n_vertices; i++)
	{
		if (v_camp[i] == max_v_camp)
		{
			s_p = i;
			break;
		}
	}
	std::vector<double> weight(MCache.n_edges, -1);
	double total_length = 0;
	for (int i = 0; i < MCache.n_edges; i++)
	{
		total_length += MCache.el[i];
		if (ban_edge[i] == 1)
			weight[i] = DBL_MAX;
		else
			weight[i] = MCache.el[i];
	}
	//???????????????
	cutVertex.clear();
	cutEdge.clear();
	std::vector<double> d(MCache.n_vertices, DBL_MAX);
	std::vector<int> v_p(MCache.n_vertices, -1);
	Algorithm::Dijkstra_with_restrict(MCache, s_p, weight, v_p, d);
	std::vector<int> path;
	double max_dis = -1;
	int max_id = -1;
	for (int i = 0; i < MCache.n_vertices; i++)
	{
		if (v_camp[i] == max_v_camp && d[i] > max_dis)
		{
			max_id = i;
			max_dis = d[i];
		}
	} //?????
	int e_p = max_id;
	Algorithm::FindPath(v_p, e_p, path);
	cutVertex.push_back(path[0]);
	for (int i = 1; i < path.size(); i++)
	{
		cutVertex.push_back(path[i]);
		cutEdge.push_back(MCache.vve[path[i - 1]][path[i]]);
	}
	initseam.clear();
}

void MeshCut::get_correspondence(std::vector<int> &he2idx, std::vector<int> &idx2meshvid)
{
	he2idx = this->he_to_idx;
	idx2meshvid = this->idx_to_mesh_idx;
}
