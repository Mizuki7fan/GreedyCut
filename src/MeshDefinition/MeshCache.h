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
	//面的面积
	std::vector<double> fa;
	std::vector<std::vector<int>> fv;
	
	//邻域
	std::vector<std::priority_queue<node>> dijkstra_cache;
	std::vector<std::vector<int>> dijkstra_isvisited;
	std::vector<std::vector<int>> V_VP;
	std::vector<std::vector<double>> V_D;
	std::vector<std::vector<int>> Neighbour;//存放顶点邻域信息
	std::vector<int> Max_Neighbour;
	std::vector<double> Vx, Vy, Vz;
	
	int n_vertices;
	int n_edges;
	int n_faces;
	double avg_el;
	double avg_fa;

	//辅助用
	//记录已经缓存了多少数据
	int capacity;

	void updataCapacity();
};