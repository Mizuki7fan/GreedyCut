#include "PointSampling.h"
#include <random>

PointSampling::PointSampling(MeshCache& MC)
	:MCache(MC)
{
}

PointSampling::~PointSampling()
{
}

void PointSampling::Set(std::string method_id)
{
	this->method = method_id;
	if (method == "Dijkstra") { printf("SampleMethod: Dijkstra\n"); }
	else if (method == "GeodesicDistance") { printf("SampleMethod: GeodesicDistance\n"); }
	else if (method == "RealDistance") { printf("SampleMethod: RealDistance\n"); }
	std::default_random_engine e(time(NULL));
	std::uniform_int_distribution<unsigned> u(0, MCache.NVertices);
	this->first_point = u(e);
	printf("%s%d\n", "Generate Random FirstPoint:", this->first_point);
}

void PointSampling::ComputeSamples(std::vector<int>& Result)
{
	std::vector<int> selectpts;
	selectpts.push_back(this->first_point);
	if (method == "RealDistance")
	{
		for (int i = 0; i < n_samples - 1; i++)
		{
			int farthestid = 0;
			double dist = 0.0;
			for (int j = 0; j < MCache.NVertices; j++)
			{
				
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
	else if (method == "Dijkstra")
	{
		Algorithm::Dijkstra_all(MCache, selectpts[0]);
		for (int i = 0; i < n_samples - 1; i++)
		{
			Algorithm::Dijkstra_all(MCache, selectpts[selectpts.size() - 1]);
			int farthestid = 0;
			double dist = 0.0;
			for (int j = 0; j < MCache.NVertices; j++)
			{
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
	Result = selectpts;
}

void PointSampling::ComputeSamplesFromSelectedRegion(std::string method, std::vector<int>& Region, std::vector<int>& SampleVertices, std::vector<int>& SampleEdges)
{
	int s_p = -1;
	for (int i = 0; i < Region.size(); i++)
	{
		if (Region[i] == 0)
			continue;
		for (int j = 0; j < MCache.VV[i].size(); j++)
		{
			if (Region[MCache.VV[i][j]] == 0)
			{
				s_p = i;
				break;
			}
		}
	}
	std::vector<double> weight(MCache.NEdges, DBL_MAX);
	for (int i = 0; i < MCache.NEdges; i++)
	{
		int v1 = MCache.EV[i][0];
		int v2 = MCache.EV[i][1];
		if (Region[v1] == 1 && Region[v2] == 1)
		{
			weight[i] = MCache.EL[i];
		}
	}
	std::vector<int> v_p(MCache.NVertices, 0); std::vector<double> d(MCache.NEdges, DBL_MAX);
	Algorithm::Dijkstra_with_restrict(MCache, s_p, weight, v_p, d);
	int e_p = -1;
	if (method == "RealDistance")
	{
		//Âè?Êâ?1‰∏?ÁÇ?
		int farthestid = 0;
		double dist = 0.0;
		for (int j = 0; j < MCache.NVertices; j++)
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
		double total_l = MCache.AVG_EL * MCache.NVertices;
		for (int i = 0; i < MCache.NVertices; i++)
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
	std::vector<int> EdgeV;
	SampleEdges.clear();
	Algorithm::FindPath(v_p, e_p, EdgeV);
	for (int i = 1; i < EdgeV.size(); i++)
		SampleEdges.push_back(MCache.VVE[EdgeV[i - 1]][EdgeV[i]]);
	SampleVertices = EdgeV;
}

double PointSampling::PointSetRealDistance(const int& vid, const std::vector<int>& pointset) const
{
	double dist = DBL_MAX;
	for (const auto& id : pointset)
	{
		double dx = MCache.Vx[vid] - MCache.Vx[id];
		double dy = MCache.Vy[vid] - MCache.Vy[id];
		double dz = MCache.Vz[vid] - MCache.Vz[id];
		double disttemp = std::sqrt(dx * dx + dy * dy + dz * dz);
		dist = disttemp < dist ? disttemp : dist;
	}
	return dist;
}

double PointSampling::PointSetPathLength(const int& vid, const std::vector<int>& pointset) const
{
	double dist = DBL_MAX;
	for (const auto& id : pointset)
	{
		double disttemp = MCache.Vd[id][vid];
		dist = disttemp < dist ? disttemp : dist;
	}
	return dist;
}