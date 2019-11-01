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
