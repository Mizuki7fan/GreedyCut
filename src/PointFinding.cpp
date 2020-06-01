#include "PointFinding.h"
#include <random>


PointFinding::PointFinding(const Mesh& m, const Mesh& pm, MeshCache& MCache) :OriMesh(m), ParaedMesh(pm), MC(MCache)
{

	PrepareComputeDistortion();
	ComputeFaceDistortion();
	ComputeVertexDistortion();
}

PointFinding::~PointFinding()
{
}

void PointFinding::Set(std::string metric)
{
	this->Metric = metric;
}

void PointFinding::Find(std::vector<std::pair<int, double>>& result)
{
	if (Metric == "Dijkstra")
		FindByRealDis(result);
	else if (Metric == "Neighbourhood")
		FindByNeighbourhood(result);
}

void PointFinding::FindLocalMaximizer()
{
	LocalMaximizer.resize(MC.NVertices, 0);
	for (int i = 0; i < MC.NVertices; i++)
	{
		if (LocalMaximizer[i] == 1)
			continue;
		for (auto a : MC.VV[i])
		{
			if (vertex_distortion[a] < vertex_distortion[i])
				LocalMaximizer[a] = 1;
			else if (vertex_distortion[i] < vertex_distortion[a])
				LocalMaximizer[i] = 1;
		}
	}
}

std::vector<int> PointFinding::GetLocalMaximizer()
{
	std::vector<int> Result;
	for (int i = 0; i < MC.NVertices; i++)
		if (LocalMaximizer[i] == 0)
			Result.push_back(i);
	return Result;
}


void PointFinding::PrepareComputeDistortion()
{
	facearea.resize(OriMesh.n_faces());
	fpx1.resize(OriMesh.n_faces());
	fpx2.resize(OriMesh.n_faces());
	fpy2.resize(OriMesh.n_faces());
	for (const auto& fh : OriMesh.faces())
	{
		auto heh = OriMesh.halfedge_handle(fh);
		const auto& p0 = OriMesh.point(OriMesh.from_vertex_handle(heh));
		const auto& p1 = OriMesh.point(OriMesh.to_vertex_handle(heh));
		const auto& p2 = OriMesh.point(OriMesh.to_vertex_handle(OriMesh.next_halfedge_handle(heh)));
		auto np = (p1 - p0) % (p2 - p0);
		double area = np.norm() * 0.5;
		np.normalize();
		fpx1[fh.idx()] = (p1 - p0).norm();
		auto pe1 = (p1 - p0).normalized();
		auto pe2 = (np % pe1).normalized();;
		fpx2[fh.idx()] = (p2 - p0) | pe1;
		fpy2[fh.idx()] = (p2 - p0) | pe2;
		facearea[fh.idx()] = area;
	}

}

void PointFinding::ComputeFaceDistortion()
{
	double alpha = 0.5;
	if (OriMesh.n_faces() != ParaedMesh.n_faces())
	{
		return;
	}
	facedistortion.resize(OriMesh.n_faces());
	for (const auto& fh : ParaedMesh.faces())
	{
		auto fid = fh.idx();
		auto heh = ParaedMesh.halfedge_handle(fh);
		const auto& p0 = ParaedMesh.point(ParaedMesh.from_vertex_handle(heh));
		const auto& p1 = ParaedMesh.point(ParaedMesh.to_vertex_handle(heh));
		const auto& p2 = ParaedMesh.point(ParaedMesh.to_vertex_handle(ParaedMesh.next_halfedge_handle(heh)));
		auto np = (p1 - p0) % (p2 - p0);
		double area = np.norm() * 0.5;
		np.normalize();
		auto qx1 = (p1 - p0).norm();
		auto pe1 = (p1 - p0).normalized();
		auto pe2 = (np % pe1).normalized();;
		auto qx2 = (p2 - p0) | pe1;
		auto qy2 = (p2 - p0) | pe2;
		double a = qx1 / fpx1[fid];
		double b = qy2 / fpy2[fid];
		double det = a * b;
		if (det > 0.0)
		{
			double c = ((-qx1 * fpx2[fid] + qx2 * fpx1[fid]) / fpx1[fid] / fpy2[fid]);
			facedistortion[fid] = alpha * (a * a + b * b + c * c) / det * 0.5 + (1 - alpha) * (det + 1.0 / det) * 0.5;
			//facedistortion[fid] = (det + 1.0 / det) * 0.5;
			//facedistortion[fid] = 1.0 / det;
		}
		else
		{
			facedistortion[fid] = std::numeric_limits<double>::infinity();
		}
	}
}

void PointFinding::ComputeVertexDistortion()
{
	vertex_distortion.resize(OriMesh.n_vertices());
	for (int i = 0; i < vertex_distortion.size(); i++)
	{
		Mesh::VertexHandle v = OriMesh.vertex_handle(i);
		double v_dis = 0;
		int v_valence = 0;
		for (const auto& vfh : OriMesh.vf_range(v))
		{
			v_dis += facedistortion[vfh.idx()];
			++v_valence;
		}
		vertex_distortion[i] = v_dis / v_valence;
	}
}

void PointFinding::FindByRealDis(std::vector<std::pair<int, double>>& result)
{
	std::vector<double> Priority(MC.NVertices, 0);

	for (int i = 0; i < LocalMaximizer.size(); i++)
	{
		if (LocalMaximizer[i] == 1)
		{
			continue;
		}
		double v_dis = vertex_distortion[i];
		std::vector<int>& is_visited = MC.DijkstraIsVisited[i];
		std::vector<double>& distance = MC.Vd[i];
		std::priority_queue<node>& que = MC.DijkstraCache[i];
		std::vector<int>& v_p = MC.VVp[i];
		if (is_visited.size() == 0)
		{
			is_visited.resize(MC.NVertices, 0);
			distance.resize(MC.NVertices, DBL_MAX);
			distance[i] = 0;
			que.push(node(i, 0));
			v_p.resize(MC.NVertices, -1);
			v_p[i] = i;
		}
		else
		{
			double level = MC.AVG_EL * MC.NEdges;
			for (int j = 0; j < MC.Vd[i].size(); j++)
			{
				if (is_visited[j] == 0)
					continue;
				else if (v_dis < vertex_distortion[j])
				{
					if (level > MC.Vd[i][j])
						level = MC.Vd[i][j];
				}
			}
			if (level < MC.AVG_EL * MC.NEdges)
			{
				Priority[i] = level;
				continue;
			}
		}
		while (1)
		{
			if (que.empty())
			{
				Priority[i] = DBL_MAX;
				break;
			}
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
			if (vertex_distortion[tmp.id] > v_dis)
			{
				if (distance[tmp.id] > MC.AVG_EL)
				{
					Priority[i] = distance[tmp.id];
				}
				break;
			}
		}
	}
	for (int i = 0; i < Priority.size(); i++)
		if (Priority[i] >= MC.AVG_EL)
			result.push_back({ i,Priority[i] });
}

void PointFinding::FindByNeighbourhood(std::vector<std::pair<int, double>>& FeaturePoints)
{
	FeaturePoints.clear();
	std::vector<int> Priority(MC.NVertices, 0);
	for (int i = 0; i < OriMesh.n_vertices(); i++)
	{
		if (LocalMaximizer[i] == 1)
			continue;
		double v_dis = vertex_distortion[i];
		std::vector<int> is_visited(MC.NVertices, 0);
		is_visited[i] = 1;
		std::vector<int> now_competitor(0), next_competitor(0);
		now_competitor = MC.VV[i];
		for (int u = 0; u < now_competitor.size(); u++)
		{
			is_visited[now_competitor[u]] = 1;
		}
		int level = 1;
		while (1)
		{
			if (level == 50)
			{
				Priority[i] = 50;
				break;
			}
			next_competitor.clear();
			for (auto a : now_competitor)
			{
				for (auto b : MC.VV[a])
				{
					if (!is_visited[b])
					{
						next_competitor.push_back(b);

					}
				}
			}
			std::sort(next_competitor.begin(), next_competitor.end());

			next_competitor.erase(std::unique(next_competitor.begin(), next_competitor.end()), next_competitor.end());
			now_competitor = next_competitor;
			int mudeki = 1;
			for (auto a : now_competitor)
			{
				if (v_dis < vertex_distortion[a])
				{
					mudeki = 0;
					break;
				}
				else
				{
					is_visited[a] = 1;
				}
			}
			if (mudeki == 0)
			{
				Priority[i] = level;
				break;
			}
			else
			{
				level++;
			}
		}
	}
	for (int i = 0; i < Priority.size(); i++)
		if (Priority[i] >= 1)
			FeaturePoints.push_back({ i,Priority[i] });
}