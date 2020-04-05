#include "MeshCut.h"

MeshCut::MeshCut(Mesh &mesh,MeshCache& MCache):oriMesh(mesh),MC(MCache)
{
}

MeshCut::~MeshCut()
{
}

void MeshCut::Set(const std::vector<int>& lmk, const std::vector<int>& initseam)
{
	this->lmk = lmk;
	this->initseam = initseam;
}

void MeshCut::Connect()
{
	Algorithm::Dijkstra_group(MC, lmk);
	Algorithm::Kruskal(MC, lmk, cutVertex, cutEdge);
}

void MeshCut::MakeSeam()
{
	cuted_mesh.clean();
	std::vector<int> candidate_cut_valence(MC.NVertices, 0);
	std::vector<int> candidate_seam_edge(MC.NEdges, 0);
	std::vector<int> candidate_seam_vertex(MC.NVertices, 0);
	for (auto a : cutEdge)
		candidate_seam_edge[a] = 1;
	for (auto a : cutVertex)
		candidate_seam_vertex[a] = 1;
	for (int i = 0; i < cutEdge.size(); i++)
	{
		Mesh::HalfedgeHandle he = oriMesh.halfedge_handle(oriMesh.edge_handle(cutEdge[i]), 0);
		candidate_cut_valence[oriMesh.to_vertex_handle(he).idx()]++;
		candidate_cut_valence[oriMesh.from_vertex_handle(he).idx()]++;
	}
	he_to_idx.resize(2 * MC.NEdges, 0);
	idx_to_mesh_idx.clear();
	Mesh::HalfedgeHandle h_begin;
	for (int i = 0; i < candidate_seam_edge.size(); i++)
	{
		if (candidate_seam_edge[i])
		{
			h_begin = oriMesh.halfedge_handle(oriMesh.edge_handle(i), 0);
			break;
		}
	}
	auto h_iter = h_begin;
	int uv_idx = 0;
	do
	{
		he_to_idx[h_iter.idx()] = uv_idx;
		idx_to_mesh_idx.push_back(oriMesh.to_vertex_handle(h_iter).idx());
		h_iter = oriMesh.next_halfedge_handle(h_iter);
		while (!candidate_seam_edge[oriMesh.edge_handle(h_iter).idx()])
		{
			h_iter = oriMesh.opposite_halfedge_handle(h_iter);
			he_to_idx[h_iter.idx()] = uv_idx;
			h_iter = oriMesh.next_halfedge_handle(h_iter);
		}
		uv_idx++;
	} while (h_iter != h_begin);
	for (const auto& vh : oriMesh.vertices())
	{
		if (candidate_cut_valence[vh.idx()] > 0)
			continue;
		for (const auto& viheh : oriMesh.vih_range(vh))
		{
			he_to_idx[viheh.idx()] = uv_idx;
		}
		idx_to_mesh_idx.push_back(vh.idx());
		uv_idx++;
	}
	for (int i = 0; i < uv_idx; i++)
	{
		cuted_mesh.add_vertex(oriMesh.point(oriMesh.vertex_handle(idx_to_mesh_idx[i])));
	}
	for (const auto& fh : oriMesh.faces())
	{
		std::vector<Mesh::VertexHandle> face_vhandles;
		for (const auto& fheh : oriMesh.fh_range(fh))
		{
			face_vhandles.push_back(cuted_mesh.vertex_handle(he_to_idx[fheh.idx()]));
		}
		cuted_mesh.add_face(face_vhandles);
	}
}

void MeshCut::SetBanCondition(const std::vector<int>& banV, const std::vector<int>& banE, const std::string BanMethod)
{
	if (BanMethod == "NonConnect")
	{
		std::set<int> ban;
		for (auto a : banV)
			ban.insert(a);
		for (auto a : banE)
			ban.insert(a);
		for (auto a : ban)
			ban_vertex.push_back(a);
	}
}

	void MeshCut::CalcBanArea(int Dn, std::string Metric, double Alpha, double shrinkRate)
	{
		std::vector<int> BanV(MC.NVertices,0);
		std::vector<int> BanE(MC.NEdges, 0);
		std::vector<int> BanF(MC.NVertices, 0);
		while (1)
		{
			BanV.clear(); BanV.resize(MC.NVertices, 0);
			BanE.clear(); BanE.resize(MC.NEdges, 0);
			BanF.clear(); BanF.resize(MC.NFaces, 0);
			if (Metric == "Dijkstra")
			{
				double range = Dn * MC.AVG_EL;
				double total_length = MC.AVG_EL * MC.NVertices;
				#pragma omp parallel for
				for (int u = 0; u < ban_vertex.size(); u++)
				{
					int s_p = ban_vertex[u];
					std::vector<double>& distance = MC.Vd[s_p];
					if (distance.size() != 0)
					{
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
					std::vector<int>& is_visited = MC.DijkstraIsVisited[s_p];

					std::priority_queue<node>& que = MC.DijkstraCache[s_p];

					std::vector<int>& v_p = MC.VVp[s_p];
					if (distance.size() == 0)
					{
						is_visited.resize(MC.NVertices, 0);
						distance.resize(MC.NVertices, DBL_MAX);
						distance[s_p] = 0;
						que.push(node(s_p, 0));
						v_p.resize(MC.NVertices, -1);
						v_p[s_p] = s_p;
					}
					while (1)
					{
						node tmp = que.top();
						que.pop();
						if (is_visited[tmp.id])
							continue;
						is_visited[tmp.id] = 1;
						for (int u = 0; u < MC.VV[tmp.id].size(); u++)
						{
							int vid = MC.VV[tmp.id][u];
							int eid = MC.VE[tmp.id][u];
							if (distance[tmp.id] + MC.EL[eid] < distance[vid])
							{
								distance[vid] = distance[tmp.id] + MC.EL[eid];
								que.push(node(vid, distance[vid]));
								v_p[vid] = tmp.id;
							}
						}
						if (distance[tmp.id] > range)
						{
							break;
						}
						else
							BanV[tmp.id] = 1;
					}
				}
			}
			else if (Metric == "Neighbourhood")
			{
#pragma omp parallel for
				for (int u = 0; u < ban_vertex.size(); u++)
				{
					int a = ban_vertex[u];
					Algorithm::UpdateNeighbourhood(MC, Dn, a);
					for (int i = 0; i < MC.NVertices; i++)
						if (MC.Neighbour[a][i] <= Dn)
							BanV[i] = 1;
				}
			}
			for (int i = 0; i < MC.NEdges; i++)
			{
				int v1 = MC.EV[i][0];
				int v2 = MC.EV[i][1];
				if (BanV[v1] == 1 || BanV[v2] == 1)
					BanE[i] = 1;
			}
			std::vector<int> v_camp = BanV;
			int max_region = 2;
			std::vector<double> weight(MC.NEdges, DBL_MAX);
			for (int i = 0; i < MC.NEdges; i++)
				if (BanE[i] == 0)
					weight[i] = MC.EL[i];
			for (int i = 0; i < MC.NVertices; i++)
			{
				if (v_camp[i] != 0)
					continue;
				std::vector<double> d(MC.NVertices, DBL_MAX);
				std::vector<int> v_i(MC.NVertices, 0);
				Algorithm::Dijkstra_with_restrict(MC, i, weight, v_i, d);
				double total_length = MC.AVG_EL * MC.NEdges;
				for (int u = 0; u < d.size(); u++)
				{
					if (d[u] < total_length)
						v_camp[u] = max_region;
				}
				max_region++;
			}
			std::vector<double> area_region(max_region, 0);
			for (int i = 0; i < MC.NFaces; i++)
			{
				int v1 = MC.FV[i][0];
				int v2 = MC.FV[i][1];
				int v3 = MC.FV[i][2];
				int f_region = 1;
				if (v_camp[v1] != 1)
					f_region = v_camp[v1];
				if (v_camp[v2] != 1)
					f_region = v_camp[v2];
				if (v_camp[v3] != 1)
					f_region = v_camp[v3];
				area_region[f_region] += MC.FA[i];
			}
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
			if (max_region_area < Alpha * MC.AVG_FA * MC.NFaces)
			{
				Dn = std::floor(shrinkRate * Dn);
				printf("Shrink in Mcut2");
			}
			else
			{
				MaxConnectedRegion.clear();
				MaxConnectedRegion.resize(MC.NVertices);
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
		std::ofstream tmp("BanV.txt");
		tmp << "VERTICES" << std::endl;
		for (auto a : BanV)
			tmp << a << std::endl;
		tmp.close();
}

void MeshCut::SetCut(std::vector<int>& cV, std::vector<int> cE)
{
	this->cutVertex = cV;
	this->cutEdge = cE;
}

void MeshCut::GetCorrespondence(std::vector<int>& he2idx, std::vector<int>& idx2meshvid)
{
	he2idx = this->he_to_idx;
	idx2meshvid = this->idx_to_mesh_idx;
}
