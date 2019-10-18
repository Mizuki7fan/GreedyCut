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
		//��������ʵ�ʾ��뵽FirstPoint��Զ�ĵ�
		for (int i = 0; i < n_samples - 1; i++)
		{
			int farthestid = 0;
			double dist = 0.0;
			for (int j = 0; j < MCache.n_vertices; j++)
			{
				//�ҵ��㼯selectpts���е���̾������ĵ�
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
	//��ؾ���
//	else if (method_id==1)
	//������·��
	else if (method_id == 2)
	{
		//����slectpts�ĵ㼯��������Щ�㵽���������еĵ�ľ���
		Algorithm::Dijkstra_all(MCache, selectpts[0]);
		for (int i = 0; i < n_samples - 1; i++)
		{
			Algorithm::Dijkstra_all(MCache, selectpts[selectpts.size()-1]);
			int farthestid = 0;
			double dist = 0.0;
			for (int j = 0; j < MCache.n_vertices; j++)
			{
				//�ҵ��㼯selectpts���е���̾������ĵ�
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

//����һ�����Ϻ�һ���㣬�������㵽����������е�����·�������ֵ
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

//����һ�����Ϻ�һ���㣬�������㵽����������е�����·�������ֵ
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
