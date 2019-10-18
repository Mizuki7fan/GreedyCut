#include "PointSampling.h"


PointSampling::PointSampling(MeshCache& m,int FirstPoint,int n_samples,int method_id)
	:MCache(m),n_samples(n_samples),method_id(method_id),FirstPoint(FirstPoint)
{

}

PointSampling::~PointSampling()
{
}

std::vector<int> PointSampling::ComputeSamples()
{
	std::vector<int> selectpts;
	selectpts.push_back(FirstPoint);
	if (method_id == 0)
	{
		//找网格中实际距离到FirstPoint最远的点
		for (int i = 0; i < n_samples - 1; i++)
		{
			int farthestid = 0;
			double dist = 0.0;
			for (int j = 0; j < MCache.n_vertices; j++)
			{
				//找到点集selectpts所有点最短距离最大的点
				double disttemp = PointSetRealDistance(j, selectpts);
				if (disttemp > dist)
				{
					dist = disttemp;
					farthestid = j;
				}
			}
			selectpts.push_back(farthestid);
		}
	}
	//测地距离
//	else if (method_id==1)
	//两点间的路径
	else if (method_id == 2)
	{
		//对于slectpts的点集，遍历这些点到网格上所有的点的距离
		Algorithm::Dijkstra_all(MCache, selectpts[0]);
		for (int i = 0; i < n_samples - 1; i++)
		{
			Algorithm::Dijkstra_all(MCache, selectpts[selectpts.size()-1]);
			int farthestid = 0;
			double dist = 0.0;
			for (int j = 0; j < MCache.n_vertices; j++)
			{
				//找到点集selectpts所有点最短距离最大的点
				double disttemp = PointSetPathLength(j, selectpts);
				if (disttemp > dist)
				{
					dist = disttemp;
					farthestid = j;
				}
			}
			selectpts.push_back(farthestid);
		}
	}
	return selectpts;
}

//输入一个集合和一个点，获得这个点到这个集合所有点的最短路径的最大值
double PointSampling::PointSetRealDistance(const int & vid, const std::vector<int> & pointset) const
{
	double dist = DBL_MAX;
	for (const auto & id : pointset)
	{
		double dx = MCache.Vx[vid] - MCache.Vx[id];
		double dy = MCache.Vy[vid] - MCache.Vy[id];
		double dz = MCache.Vz[vid] - MCache.Vz[id];
		double disttemp = std::sqrt(dx * dx + dy * dy + dz * dz);
		dist = disttemp < dist ? disttemp : dist;
	}
	return dist;
}

//输入一个集合和一个点，获得这个点到这个集合所有点的最短路径的最大值
double PointSampling::PointSetPathLength(const int& vid, const std::vector<int>& pointset) const
{
	double dist = DBL_MAX;
	for (const auto& id : pointset)
	{
		double disttemp = MCache.V_D[id][vid];
		dist = disttemp < dist ? disttemp : dist;
	}
	return dist;
}
