#include "Algorithm.h"

Algorithm::Algorithm()
{
}

Algorithm::~Algorithm()
{
}

void Algorithm::Dijkstra(MeshCache &MC, std::vector<int> &lmk)
{
#pragma omp parallel for
	for (int i = 0; i < lmk.size(); i++)
	{
		std::vector<int> is_lmk(MC.n_vertices, 0);
		for (int i = 0; i < lmk.size(); i++)
			is_lmk[lmk[i]] = 1;
		int count = lmk.size();
		int s_p = lmk[i];
		std::vector<int> &is_visited = MC.dijkstra_isvisited[s_p];
		std::vector<double> &distance = MC.V_D[s_p];
		std::priority_queue<node> &que = MC.dijkstra_cache[s_p];
		std::vector<int> &v_p = MC.V_VP[s_p];
		if (is_visited.size() == 0)
		{
			is_visited.resize(MC.n_vertices, 0);
			distance.resize(MC.n_vertices, DBL_MAX);
			distance[s_p] = 0;
			que.push(node(s_p, 0));
			v_p.resize(MC.n_vertices, -1);
			v_p[s_p] = s_p;
		}
		for (int i = 0; i < lmk.size(); i++)
		{
			if (is_visited[lmk[i]] != 0)
				count--;
		}
		while (count)
		{
			node tmp = que.top();
			que.pop();
			if (is_visited[tmp.id] == 1)
				continue;
			if (is_lmk[tmp.id] == 1)
			{
				count--;
				is_lmk[tmp.id] = 0;
			}
			is_visited[tmp.id] = 1;
			for (int u = 0; u < MC.vv[tmp.id].size(); u++)
			{
				int v_id = MC.vv[tmp.id][u];
				int e_id = MC.ve[tmp.id][u];
				if (distance[tmp.id] + MC.el[e_id] < distance[v_id])
				{
					distance[v_id] = distance[tmp.id] + MC.el[e_id];
					que.push(node(v_id, distance[v_id]));
					v_p[v_id] = tmp.id;
				}
			}
		}
	}
}

void Algorithm::Kruskal(MeshCache &MCache, std::vector<int> &lmk, std::vector<int> &cutvertex, std::vector<int> &cutedge)
{
	int nv = lmk.size();
	std::vector<int> spanning_tree_current(nv);
	std::priority_queue<PathInfo> EdgeInfo;
	for (int i = 0; i < lmk.size(); i++)
		for (int j = i + 1; j < lmk.size(); j++)
		{
			EdgeInfo.push(PathInfo(lmk[i], lmk[j], MCache.V_D[lmk[i]][lmk[j]]));
		}
	std::vector<int> edge;
	std::set<int> vertex;
	for (int i = 0; i < nv; ++i)
	{
		spanning_tree_current[i] = i;
	}
	int index = 0;
	int j = 0;
	while (index < nv - 1)
	{
		PathInfo tmp = EdgeInfo.top();
		EdgeInfo.pop();
		int m_id = -1, n_id = -1;
		for (int u = 0; u < lmk.size(); u++)
		{
			if (lmk[u] == tmp.s_p)
				m_id = u;
			if (lmk[u] == tmp.e_p)
				n_id = u;
		}
		m_id = spanning_tree_current[m_id];
		n_id = spanning_tree_current[n_id];
		if (m_id < n_id)
			std::swap(m_id, n_id);
		if (m_id != n_id)
		{ //根据两个点找中间的边
			std::vector<int> path;
			Algorithm::FindPath(MCache.V_VP[tmp.s_p], tmp.e_p, path);
			vertex.insert(path[0]);
			for (int i = 0; i < path.size() - 1; i++)
			{
				vertex.insert(path[i + 1]);
				edge.push_back(MCache.vve[path[i]][path[i + 1]]);
			}
			index++;
			for (int i = 0; i < nv; ++i)
			{
				if (spanning_tree_current[i] == n_id)
					spanning_tree_current[i] = m_id;
			}
		}
	}
	std::set<int> e;
	for (auto a : edge)
		e.insert(a);
	for (auto a : e)
	{
		cutedge.push_back(a);
	}
	for (auto a : vertex)
		cutvertex.push_back(a);
}

void Algorithm::Kruskal(std::vector<int>& lmk,std::priority_queue<PathInfo> que,std::vector<PathInfo>& Res)
{
	//统计有多少个顶点？
	int nv = lmk.size();
	std::vector<int> spanning_tree_current(nv);
	std::set<int> reV, reE;
	//根据存储好的边来算生成树
	std::vector<int> edge;
	std::set<int> vertex;
	for (int i = 0; i < nv; ++i)
	{
		spanning_tree_current[i] = i;
	}
	int index = 0;
	int j = 0;
	while (index < nv - 1)
	{
		PathInfo tmp = que.top();
		que.pop();
		int m_id = -1, n_id = -1;
		for (int u = 0; u <nv; u++)
		{
			if (lmk[u] == tmp.s_p)
				m_id = u;
			if (lmk[u] == tmp.e_p)
				n_id = u;
		}
		m_id = spanning_tree_current[m_id];
		n_id = spanning_tree_current[n_id];
		if (m_id < n_id)
			std::swap(m_id, n_id);
		if (m_id != n_id)
		{ //根据两个点找中间的边
			Res.push_back(tmp);
			index++;
			for (int i = 0; i < nv; ++i)
			{
				if (spanning_tree_current[i] == n_id)
					spanning_tree_current[i] = m_id;
			}
		}
	}
}

void Algorithm::FindPath(std::vector<int> &v_p, int e_p, std::vector<int> &path)
{
	path.clear();
	path.push_back(e_p);
	do
	{
		e_p = v_p[e_p];
		path.push_back(e_p);
	} while (e_p != v_p[e_p]);
}

void Algorithm::Dijkstra_with_restrict(MeshCache &MCache, int s_p, std::vector<double> &weight, std::vector<int> &v_p, std::vector<double> &d)
{
	double total_length = MCache.avg_el * MCache.n_edges;
	std::vector<int> is_visited(MCache.n_vertices, 0);
	d[s_p] = 0;
	is_visited[s_p] = 0;
	v_p[s_p] = s_p;
	std::priority_queue<node> que;
	que.push(node(s_p, 0));
	while (!que.empty())
	{
		node tmp = que.top();
		que.pop();
		if (is_visited[tmp.id]==1)
			continue;
		is_visited[tmp.id]=1;
		for (int u = 0; u < MCache.vv[tmp.id].size(); u++)
		{
			int vid = MCache.vv[tmp.id][u];
			int eid = MCache.ve[tmp.id][u];
			if (d[tmp.id] + weight[eid] < d[vid])
			{
				d[vid] = d[tmp.id] + weight[eid];
				que.push(node(vid, d[vid]));
				v_p[vid] = tmp.id;
			}
		}
	}
}

void Algorithm::Dijkstra_with_nearest2(MeshCache &MC, int s_p, std::vector<int> &is_target, std::vector<int> &path)
{
	int e_p = -1;
	std::vector<double> distance = MC.V_D[s_p];
	std::vector<int> is_visited = MC.dijkstra_isvisited[s_p];
	std::priority_queue<node> que = MC.dijkstra_cache[s_p];
	std::vector<int> v_p = MC.V_VP[s_p];
	if (distance.size() != 0)
	{
		double min_dis = DBL_MAX;
		int min_id = -1;
		for (int i = 0; i < is_target.size(); i++)
		{
			if (is_target[i] == 1)
			{
				if (is_visited[i] == 0)
					continue;
				if (distance[i] < min_dis)
				{
					min_dis = distance[i];
					min_id = i;
				}
			}
		}
		if (min_id != -1)
		{
			e_p = min_id;
		}
	}
	if (e_p == -1)
	{

		if (is_visited.size() == 0)
		{
			is_visited.resize(MC.n_vertices, 0);
			distance.resize(MC.n_vertices, DBL_MAX);
			distance[s_p] = 0;
			que.push(node(s_p, 0));
			v_p.resize(MC.n_vertices, -1);
			v_p[s_p] = s_p;
		}
		while (!que.empty())
		{
			node tmp = que.top();
			que.pop();
			if (is_visited[tmp.id] == 1)
				continue;
			if (is_target[tmp.id] == 1)
			{
				e_p = tmp.id;
				break;
			}
			is_visited[tmp.id] = 1;
			for (int u = 0; u < MC.vv[tmp.id].size(); u++)
			{
				int v_id = MC.vv[tmp.id][u];
				int e_id = MC.ve[tmp.id][u];
				if (distance[tmp.id] + MC.el[e_id] < distance[v_id])
				{
					distance[v_id] = distance[tmp.id] + MC.el[e_id];
					que.push(node(v_id, distance[v_id]));
					v_p[v_id] = tmp.id;
				}
			}
		}
	}

	FindPath(v_p, e_p, path);
}

void Algorithm::Dijkstra_all(MeshCache &MC, int s_p)
{
	//计算k点到网格上所有点的距离
	std::vector<double> &distance = MC.V_D[s_p];
	std::vector<int> &is_visited = MC.dijkstra_isvisited[s_p];
	std::priority_queue<node> &que = MC.dijkstra_cache[s_p];
	std::vector<int> &v_p = MC.V_VP[s_p];
	if (distance.size() == 0)
	{
		is_visited.resize(MC.n_vertices, 0);
		distance.resize(MC.n_vertices, DBL_MAX);
		distance[s_p] = 0;
		que.push(node(s_p, 0));
		v_p.resize(MC.n_vertices, -1);
		v_p[s_p] = s_p;
	}
	while (!que.empty())
	{
		node tmp = que.top();
		que.pop();
		if (is_visited[tmp.id])
			continue;
		is_visited[tmp.id] = 1;
		for (int u = 0; u < MC.vv[tmp.id].size(); u++)
		{
			int vid = MC.vv[tmp.id][u];
			int eid = MC.ve[tmp.id][u];
			if (distance[tmp.id] + MC.el[eid] < distance[vid])
			{
				distance[vid] = distance[tmp.id] + MC.el[eid];
				que.push(node(vid, distance[vid]));
				v_p[vid] = tmp.id;
			}
		}
	}
}
//标记存储顶点v的k邻域
void Algorithm::UpdateNeighbourhood(MeshCache &MC, int k, int v)
{
	if (k <= MC.Max_Neighbour[v])
	{
		std::cout << "Memory Hit" << std::endl;
		return;

	}
	if (k == 1)
	{
		if (MC.Neighbour[v].size() == 0)
		{
			MC.Neighbour[v].resize(MC.n_vertices, INT_MAX);
			MC.Neighbour[v][v] = 0;
			for (auto a : MC.vv[v])
				MC.Neighbour[v][a] = 1;
			MC.Max_Neighbour[v] = 1;
		}
	}
	else
	{
		UpdateNeighbourhood(MC, k - 1, v);
		for (int u = 0; u < MC.n_vertices; u++)
		{
			if (MC.Neighbour[v][u] == k - 1)
			{
				for (int s = 0; s < MC.vv[u].size(); s++)
				{
					if (MC.Neighbour[v][MC.vv[u][s]] > k)
						MC.Neighbour[v][MC.vv[u][s]] = k;
				}
			}
		}
		MC.Max_Neighbour[v] = k;
	}
}
