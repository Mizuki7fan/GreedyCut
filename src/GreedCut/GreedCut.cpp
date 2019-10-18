#include "GreedCut.h"

GreedCut::GreedCut(Mesh& mesh, MeshCache& MC, std::vector<int>& landmark,double trimming_rate)
	:mesh(mesh),MC(MC),landmark(landmark),trimming_rate(trimming_rate)
{
	std::cout<<"Trimming Rate =" <<trimming_rate<<std::endl;
	for (int i = 0; i < landmark.size(); i++)
	{
		Algorithm::Dijkstra_all(MC, landmark[i]);
	}
	back_landmark = landmark;



}

void GreedCut::GenGraph()
{
	for (int i = 0; i < landmark.size() - 1; i++)
		for (int j = i + 1; j < landmark.size(); j++)
		{
			PresavedPath.push(PathInfo(landmark[i], landmark[j], MC.V_D[landmark[i]][landmark[j]]));
		}
	tree_length = this->CalcTreeLength(PresavedPath, std::vector<int>());
}

double GreedCut::CalcTreeLength(std::priority_queue<PathInfo> GenTreePath, std::vector<int> k)
{
	std::vector<int> t_landmark;
	t_landmark.insert(t_landmark.end(), landmark.begin(), landmark.end());
	t_landmark.insert(t_landmark.end(), k.begin(), k.end());
	int nv = static_cast<int>(t_landmark.size());
	std::vector<int> spanning_tree_current(nv);
	for (int i = 0; i < nv; ++i)
	{
		spanning_tree_current[i] = i;
	}
	int index = 0;
	int j = 0;
	double tree_length = 0;
	while (index < nv - 1)
	{
		PathInfo tmp = GenTreePath.top();
		GenTreePath.pop();
		int m_id = -1, n_id = -1;
		for (int u = 0; u < t_landmark.size(); u++)
		{
			if (t_landmark[u] == tmp.s_p)
				m_id = u;
			if (t_landmark[u] == tmp.e_p)
				n_id = u;
		}
		m_id = spanning_tree_current[m_id];
		n_id = spanning_tree_current[n_id];
		if (m_id < n_id)
			std::swap(m_id, n_id);
		if (m_id != n_id)
		{
			tree_length += tmp.length;
			index++;
			for (int i = 0; i < nv; ++i)
			{
				if (spanning_tree_current[i] == n_id)
					spanning_tree_current[i] = m_id;
			}
		}
	}
	return tree_length;
}

void GreedCut::OP()
{
	int add_count = 0;
	greed1_ok = false;
	while (!greed1_ok)
	{
		if (add_count > 100)
			break;
		Greed1_Op();
		add_count++;
	}
}


void GreedCut::Greed1_Op()
{
	std::vector<double> vertex_weight(MC.n_vertices);
	int para_count = omp_get_num_procs();
#pragma omp parallel for num_threads(para_count)
	for (int i = 0; i < MC.n_vertices; i++)
	{
		if (std::find(landmark.begin(), landmark.end(), i) != landmark.end())
			vertex_weight[i] = this->tree_length;
		else
		{
			vertex_weight[i] = CalcExtraLandmark(i);
		}
	}
	int min_id = -1;
//	double min_dis = tree_length;
	double min_dis=DBL_MAX;
	
	for (int i = 0; i < MC.n_vertices; i++)
	{
		if (vertex_weight[i]<min_dis)
		{
			min_dis=vertex_weight[i];
			min_id=i;
		}
	}
	if (tree_length-min_dis<trimming_rate+0.0000001)
		min_id=-1;
	//else
		//std::cout << "add point " << min_id << "dischange :" << min_dis << std::endl;
	if (min_id == -1)
		greed1_ok = true;
	else
	{
		landmark.push_back(min_id);
		greed1_ok = false;
		Algorithm::Dijkstra_all(MC, min_id);
		for (int i = 0; i < landmark.size() - 1; i++)
		{
			PresavedPath.push(PathInfo(min_id, landmark[i], MC.V_D[landmark[i]][min_id]));
		}
		tree_length = CalcTreeLength(PresavedPath, std::vector<int>());
	}
}

double GreedCut::CalcExtraLandmark(int p)
{
	std::vector<int> e_landmark(1, p);
	std::priority_queue<PathInfo> GenTreePath(PresavedPath);
	for (int i = 0; i < landmark.size(); i++)
	{
		GenTreePath.push(PathInfo(landmark[i], p, MC.V_D[landmark[i]][p]));
	}
	double length = CalcTreeLength(GenTreePath, e_landmark);
	return length;
}

