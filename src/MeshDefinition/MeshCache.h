#pragma once
#include "MeshDefinition.h"
#include <iostream>
#include <queue>

class MeshCache
{
public:
	MeshCache(Mesh& mesh);
	~MeshCache();
	//����one-ring 
	std::vector<std::vector<int>> vv;
	//�����ڱ�
	std::vector<std::vector<int>> ve;
	//����֮��ı�id
	std::vector<std::map<int , int>> vve;
	//�߳�
	std::vector<double> el;
	//�ߵ���������
	std::vector<std::vector<int>> ev;
	
	//����
	std::vector<std::vector<std::vector<int>>> Neighbour;
	std::vector<std::priority_queue<node>> dijkstra_cache;
	std::vector<std::vector<int>> dijkstra_isvisited;
	std::vector<std::vector<int>> V_VP;
	std::vector<std::vector<double>> V_D;
	std::vector<double> Vx, Vy, Vz;
	
	int n_vertices;
	int n_edges;

	//������
	//��¼�Ѿ������˶�������
	int capacity;

	void updataCapacity();
};