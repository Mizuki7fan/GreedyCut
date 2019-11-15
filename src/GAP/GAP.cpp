#include "GAP.h"

GAP::GAP(const Mesh& mesh, MeshCache& MCache):mesh(mesh),MCache(MCache)
{

}

void GAP::Set(const std::vector<std::pair<int, double>>& FP,std::string PriorityMetric,int FixThreshold,int ForbiddenRadius)
{
	FeaturePoints = FP;
	this->VertexPriorityMetric = PriorityMetric;
	this->FixThreshold = FixThreshold;
	this->GAPForBiddenRadius = ForbiddenRadius;
}

void GAP::GenFirstCut()
{
	//根据度量进行lmk的分级
	//生成初始的Cut，要先求所有lmk点的5邻域或者5倍平均边长,根据配置文件中设置的度量优先级的标准，先求好禁止区域
	//设置好
	ClassifyFeaturePoints(FixThreshold);
	std::vector<std::vector<int>> ForbiddenArea;
	std::ofstream cand("cand.txt");
	cand << "VERTICES" << std::endl;
	for (auto a : CanditatePoints)
		cand << a << std::endl;
	cand.close();


	if (VertexPriorityMetric == "Neighbourhood")
	{
		//设置禁止区域的点

		//记录每个候选点的最大邻域数
		std::vector<int> RadiusInfo(CanditatePoints.size(), GAPForBiddenRadius);
		//遍历网格的每个顶点
		for (int i = 0; i < CanditatePoints.size(); i++)
		{
			Algorithm::UpdateNeighbourhood(MCache, RadiusInfo[i], CanditatePoints[i]);
			for (int u = 0; u < FixPoints.size(); u++)
			{
				if (MCache.Neighbour[CanditatePoints[i]][FixPoints[u]] < RadiusInfo[i])
				{
					RadiusInfo[i] = MCache.Neighbour[CanditatePoints[i]][FixPoints[u]];
				}
			}
		}

		//#pragma omp parallel for
		int m = 0;
		while (1)
		{
			ForbiddenArea.clear();
			ForbiddenArea.resize(MCache.n_vertices, std::vector<int>());
			for (int i = 0; i < CanditatePoints.size(); i++)
			{
				int v = CanditatePoints[i];
			//	Algorithm::UpdateNeighbourhood(MCache, RadiusInfo[i], v);
				for (int u = 0; u < MCache.n_vertices; u++)
				{
					if (MCache.Neighbour[v][u] < RadiusInfo[i])
					{
						ForbiddenArea[u].push_back(i);
					}
				}
			}
			bool close = true;
			std::vector<int> flag(CanditatePoints.size(), 0);
			for (int i = 0; i < ForbiddenArea.size(); i++)
			{
				if (ForbiddenArea[i].size() >= 2)
				{
					for (auto a : ForbiddenArea[i])
					{
						close = false;
						flag[a] = 1;
					}
				}
			}
			if (close)
				break;
			for (int i = 0; i < flag.size(); i++)
			{
				RadiusInfo[i] -= flag[i];
			}
			m++;

		}
		for (auto a : RadiusInfo)
			a--;
		ForbiddenArea.clear();
		ForbiddenArea.resize(MCache.n_vertices, std::vector<int>());
		for (int i = 0; i < CanditatePoints.size(); i++)
		{
			int v = CanditatePoints[i];
			Algorithm::UpdateNeighbourhood(MCache, RadiusInfo[i], v);
			for (int u = 0; u < MCache.n_vertices; u++)
			{
				if (MCache.Neighbour[v][u] < RadiusInfo[i])
				{
					ForbiddenArea[u].push_back(i);
				}
			}
		}
		std::ofstream out("debug_fbdregion.txt");
		out << "VERTICES" << std::endl;
		for (int i = 0; i < ForbiddenArea.size(); i++)
		{
			if (ForbiddenArea[i].size() != 0)
				out << i << std::endl;
		}
		out.close();
		std::vector<double> weight(MCache.n_edges, -1);
		for (int i = 0; i < MCache.n_edges; i++)
			weight[i] = MCache.el[i];
		//边赋一个正无穷的权
		for (int i = 0; i < ForbiddenArea.size(); i++)
		{
			if (ForbiddenArea[i].size() != 0)
			{
				for (int u = 0; u < MCache.ve[i].size(); u++)
				{
					weight[MCache.ve[i][u]] = DBL_MAX;
				}
			}
		}
		std::priority_queue<PathInfo> que;
		for (int i = 0; i < FixPoints.size(); i++)
		{//给定所有路径获得其生成树

			std::vector<int> v_p(MCache.n_vertices,-1);
			std::vector<double> d(MCache.n_vertices,DBL_MAX);
			Algorithm::Dijkstra_with_restrict(MCache, FixPoints[i], weight, v_p, d);
			for (int j = i+1; j < FixPoints.size(); j++)
			{
				std::vector<int> path;
				int s_p = FixPoints[i];
				int e_p = FixPoints[j];
				double length = 0;
				while (e_p!=s_p)
				{
					path.push_back(e_p);
					length +=MCache.el[MCache.vve[e_p][v_p[e_p]]];
					e_p = v_p[e_p];
				}
				path.push_back(e_p);

				PathInfo tmp(FixPoints[i],FixPoints[j], length);
				tmp.path = path;
				que.push(tmp);
			}
		}
	}
	else if (VertexPriorityMetric == "RealDis")
	{

	}
	
}

void GAP::ClassifyFeaturePoints(int threshold)
{
	FixPoints.clear();
	CanditatePoints.clear();
	//按照优先级进行排序,匿名函数
	std::sort(FeaturePoints.begin(), FeaturePoints.end(), [&](std::pair<int, double>& x, std::pair<int, double>& y) {return x.second > y.second; });
	for (auto a : FeaturePoints)
	{
		if (a.second > threshold)
			FixPoints.push_back(a.first);
		else
			CanditatePoints.push_back(a.first);
	}
}
