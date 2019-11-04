#include "PointSampling.h"


PointSampling::PointSampling(MeshCache& m)
	:MCache(m)
{

}

PointSampling::~PointSampling()
{
}

void PointSampling::Set(std::string method_id, std::string FirstPoint, int n_samples)
{
	this->method = method_id;
	if (method == "Dijkstra") { std::cout << "SampleMethod: Dijkstra" << std::endl; }
	else if (method == "GeodesicDistance") {	std::cout << "SampleMethod: GeodesicDistance" << std::endl;}
	switch (method == "RealDistance") { std::cout << "SampleMethod: RealDistance" << std::endl; }

	if (FirstPoint == "Random") { 
		std::cout << "SampleFirstPoint: Random" << std::endl; 
		std::default_random_engine e;
		std::uniform_int_distribution<unsigned> u(0,MCache.n_vertices);
		this->FirstPoint = u(e);
		std::cout << "Generate Random FirstPoint: " << this->FirstPoint << std::endl;
	}
	else
	{
		std::cout << "SampleFirstPoint: " << FirstPoint << std::endl;
		this->FirstPoint = std::stoi(FirstPoint);
	}
	this->n_samples = n_samples;
	std::cout << "SampleCount: " << this->n_samples << std::endl;
}

void PointSampling::ComputeSamples(std::vector<int>& Result)
{
	std::vector<int> selectpts;
	selectpts.push_back(FirstPoint);
	if (method == "RealDistance")
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
	else if (method == "Dijkstra")
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
	Result= selectpts;
}

void PointSampling::ComputeSamplesFromSelectedRegion(std::string method, std::vector<int>& Region, std::vector<int>& SampleVertices, std::vector<int>& SampleEdges)
{
	//输入采样方法和采样区域，输出采样得到的边
	//寻找起点：
	int s_p = -1;
	for (int i = 0; i < Region.size(); i++)
	{
		if (Region[i] == 0)
			continue;
		for (int j = 0; j < MCache.vv[i].size(); j++)
		{
			if (Region[MCache.vv[i][j]] == 0)
			{
				s_p = i;
				break;
			}
		}
	}
	std::vector<double> weight(MCache.n_edges, DBL_MAX);
	for (int i = 0; i < MCache.n_edges; i++)
	{
		int v1 = MCache.ev[i][0];
		int v2 = MCache.ev[i][1];
		if (Region[v1] == 1 && Region[v2] == 1)
		{
			weight[i] = MCache.el[i];
		}
	}
	std::vector<int> v_p(MCache.n_vertices, 0); std::vector<double> d(MCache.n_edges, DBL_MAX);
	Algorithm::Dijkstra_with_restrict(MCache, s_p, weight, v_p, d);
	//根据采样方法决定：
	int e_p = -1;
	if (method == "RealDistance")
	{
		//只找1个点
		int farthestid = 0;
		double dist = 0.0;
		for (int j = 0; j < MCache.n_vertices; j++)
		{
			if (Region[j] == 0)
				continue;
			double dx = MCache.Vx[j] - MCache.Vx[s_p];
			double dy = MCache.Vy[j] - MCache.Vy[s_p];
			double dz = MCache.Vz[j] - MCache.Vz[s_p];
			double disttemp = std::sqrt(dx * dx + dy * dy + dz * dz);
			if (disttemp > dist)
			{
				dist = disttemp;
				farthestid = j;
			}
		}
		e_p = farthestid;
	}
	else if (method == "Dijkstra")
	{
		double dist = 0.0;
		double total_l = MCache.avg_el * MCache.n_vertices;
		for (int i = 0; i < MCache.n_vertices; i++)
		{
			if (d[i] > total_l)
				continue;
			if (d[i] > dist)
			{
				dist = d[i];
				e_p = i;
			}
		}
	}
	//获取最远点e_p；
	std::vector<int> EdgeV;
	SampleEdges.clear();
	//获取路径
	Algorithm::FindPath(v_p, e_p, EdgeV);
	for (int i = 1; i < EdgeV.size(); i++)
		SampleEdges.push_back(MCache.vve[EdgeV[i-1]][EdgeV[i]]);
	std::ofstream cut2("cut2.txt");
	cut2 << "VERTICES" << std::endl;
	for (auto a : EdgeV)
		cut2 << a << std::endl;
	cut2 << "EDGES" << std::endl;
	for (auto a : SampleEdges)
		cut2 << a << std::endl;
	cut2.close();
	SampleVertices = EdgeV;





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
