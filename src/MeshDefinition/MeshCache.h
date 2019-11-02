#pragma once
#include "MeshDefinition.h"
#include <iostream>
#include <queue>

class MeshCache
{
public:
	MeshCache(Mesh& mesh);
	~MeshCache();
	//存点的one-ring 
	std::vector<std::vector<int>> vv;
	//存点的邻边
	std::vector<std::vector<int>> ve;
	//两点之间的边id
	std::vector<std::map<int , int>> vve;
	//边长
	std::vector<double> el;
	//边的两个顶点
	std::vector<std::vector<int>> ev;
	
	//邻域
	std::vector<std::vector<std::vector<int>>> Neighbour;
	std::vector<std::priority_queue<node>> dijkstra_cache;
	std::vector<std::vector<int>> dijkstra_isvisited;
	std::vector<std::vector<int>> V_VP;
	std::vector<std::vector<double>> V_D;
	std::vector<double> Vx, Vy, Vz;
	
	int n_vertices;
	int n_edges;
	double avg_el;

	//辅助用
	//记录已经缓存了多少数据
	int capacity;

	void updataCapacity();
};