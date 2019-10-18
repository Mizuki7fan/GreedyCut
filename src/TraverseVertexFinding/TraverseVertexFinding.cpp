#include "TraverseVertexFinding.h"
#include <iostream>

TraverseVertexFinding::TraverseVertexFinding(const Mesh& mesh, const Mesh& paraed_mesh, MeshCache& MC)
	:origin_mesh(mesh),paraed_mesh(paraed_mesh),MC(MC)
{
	PrepareComputeDistortion();
	ComputeDistortion();
	CalcVertexDistortion();

}

TraverseVertexFinding::~TraverseVertexFinding()
{
}

void TraverseVertexFinding::ComputeDistortion()
{
	double alpha = 0.5;
	if (origin_mesh.n_faces() != paraed_mesh.n_faces())
	{
		return;
	}
	facedistortion.resize(origin_mesh.n_faces());
	for (const auto& fh : paraed_mesh.faces())
	{
		auto fid = fh.idx();
		auto heh = paraed_mesh.halfedge_handle(fh);
		const auto& p0 = paraed_mesh.point(paraed_mesh.from_vertex_handle(heh));
		const auto& p1 = paraed_mesh.point(paraed_mesh.to_vertex_handle(heh));
		const auto& p2 = paraed_mesh.point(paraed_mesh.to_vertex_handle(paraed_mesh.next_halfedge_handle(heh)));
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

void TraverseVertexFinding::PrepareComputeDistortion(void)
{
	facearea.resize(origin_mesh.n_faces());
	fpx1.resize(origin_mesh.n_faces());
	fpx2.resize(origin_mesh.n_faces());
	fpy2.resize(origin_mesh.n_faces());
	for (const auto& fh : origin_mesh.faces())
	{
		auto heh = origin_mesh.halfedge_handle(fh);
		const auto& p0 = origin_mesh.point(origin_mesh.from_vertex_handle(heh));
		const auto& p1 = origin_mesh.point(origin_mesh.to_vertex_handle(heh));
		const auto& p2 = origin_mesh.point(origin_mesh.to_vertex_handle(origin_mesh.next_halfedge_handle(heh)));
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

void TraverseVertexFinding::CalcVertexDistortion()
{
	//根据面片扭曲计算顶点扭曲
	vertex_distortion.resize(MC.n_vertices);
	for (int i = 0; i < vertex_distortion.size(); i++)
	{
		Mesh::VertexHandle v = origin_mesh.vertex_handle(i);
		double v_dis = 0;
		int v_valence = 0;
		for (const auto& vfh : origin_mesh.vf_range(v))
		{
			v_dis += facedistortion[vfh.idx()];
			++v_valence;
		}
		vertex_distortion[i] = v_dis / v_valence;
	}
}

void TraverseVertexFinding::FindLandmarkPoints(std::vector<std::pair<int, int>>& FeaturePoints)
{

	std::vector<int> result(MC.n_vertices,0);
	int v_nei_range = 1;
	int k_inf = 50;
	if (k_inf > MC.Neighbour[0].size())
		k_inf = MC.Neighbour[0].size();
	std::vector<int> is_failed(MC.n_vertices, 0);
	for (int i = 0; i < MC.n_vertices; i++)
	{
		//如果该点之前被其相邻的点领先过，则直接删除该点
		if (is_failed[i] == 1)
			continue;
		double v_dis = vertex_distortion[i];
		std::vector<int> is_visited(MC.n_vertices, 0);
		is_visited[i] = 1;
		//等级最高为20
		//初始的对手为1邻域，如果能打倒所有对手，则升级并且根据当前对手查找下一个对手
		std::vector<int> now_competitor(0), next_competitor(0);
		now_competitor = MC.vv[i];
		//先打1邻域的对手，如果扑街则直接跳过，否则让1邻域的所有对手扑街
		for (int u = 0; u < now_competitor.size(); u++)
		{
			double tmp = vertex_distortion[now_competitor[u]];
			if (v_dis > tmp)
			{
				//对手扑街
				is_failed[now_competitor[u]] = 1;
				is_visited[now_competitor[u]] = 1;
			}
			else
			{
				is_failed[i] = 1;
			}
		}
		if (is_failed[i] == 1)
			continue;
		int level = 1;

		while (1)
		{
			if (level == 50)
			{
				result[i] = 50;
				break;
			}
			next_competitor.clear();
			for (auto a : now_competitor)
			{
				for (auto b: MC.vv[a])
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
				{//我怎么会输？
					mudeki = 0;
					break;
				}
				else
				{//彩笔，还有脸再来？
					is_visited[a] = 1;
				}
			}
			//挑战结束
			if (mudeki == 0)
			{
				result[i] = level;
				break;
			}
			else
			{
				level++;
			}
		}
	}
	for (int i = 0; i < result.size(); i++)
	{
		if (result[i] >= 1)
			FeaturePoints.push_back({ i,result[i] });
	}
/*
	for (int i = 0; i < MC.n_vertices; i++)
	{
		v_nei_range = 1;
		int flag = 0;
		double v_dis = vertex_distortion[i];
		while (v_nei_range < k_inf)
		{
			for (int j = 0; j < MC.Neighbour[i][v_nei_range].size(); j++)
			{
				int p = MC.Neighbour[i][v_nei_range][j];
				//比较两个点的扭曲值大小
				if (v_dis < vertex_distortion[p])
				{
					flag = 1;
				}
				else
				{
					if (v_nei_range - 1 < result[p])
						result[p] = v_nei_range - 1;
				}

			}
			if (flag == 0)
			{
				v_nei_range++;
			}
			else
			{
				break;
			}
		}
		result[i] = v_nei_range - 1;
	}
	int min_range = 1;
	for (int i = 0; i < result.size(); i++)
	{
		if (result[i] >= min_range)
			FeaturePoints.push_back({ i,result[i] });
	}
	*/
	std::cout << "Find Over" << std::endl;
}

