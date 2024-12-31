#include "Auxiliary.h"

Option::Option(std::string modelPath, std::string optPath) {
  this->modelPath = modelPath;
  this->modelName = modelPath.substr(modelPath.find_last_of('/') + 1,
                                     modelPath.length()); // for windows
  std::ifstream input(optPath);
  std::string line, key, value;
  int position;
  while (std::getline(input, line)) {
    if (line[0] == '#')
      continue;
    position = line.find('=');
    if (position == line.npos)
      continue;
    key = line.substr(0, position);
    value = line.substr(position + 1, line.length());
    if (key == "PointSampling_method")
      PS_method = value;
    else if (key == "PointFinding_method")
      PF_method = value;
    else if (key == "BanArea_Method")
      BanArea_Method = value;
    else if (key == "BanArea_Dn")
      BanArea_Dn = std::stoi(value);
    else if (key == "BanArea_Alpha")
      BanArea_alpha = std::stod(value);
    else if (key == "BanArea_Metric")
      BanArea_Metric = value;
    else if (key == "BanArea_ShrinkRate")
      BanArea_ShrinkRate = std::stod(value);
    else if (key == "Influence_Threshold")
      Influence_Threshold = std::stoi(value);
    else if (key == "Distortion_Threshold")
      Distortion_Threshold = std::stod(value);
    else if (key == "Trimming_Rate")
      Trimming_Rate = std::stod(value);
    else if (key == "Max_AddCount")
      Max_AddCount = std::stoi(value);
  }
}

Algorithm::~Algorithm() {}

void Algorithm::Dijkstra_all(MeshCache &MC, int s_p) {
  std::vector<double> &distance = MC.Vd[s_p];
  std::vector<int> &is_visited = MC.DijkstraIsVisited[s_p];
  std::priority_queue<node> &que = MC.DijkstraCache[s_p];
  std::vector<int> &v_p = MC.VVp[s_p];
  if (distance.size() == 0) {
    is_visited.resize(MC.NVertices, 0);
    distance.resize(MC.NVertices, DBL_MAX);
    distance[s_p] = 0;
    que.push(node(s_p, 0));
    v_p.resize(MC.NVertices, -1);
    v_p[s_p] = s_p;
  }
  while (!que.empty()) {
    node tmp = que.top();
    que.pop();
    if (is_visited[tmp.id])
      continue;
    is_visited[tmp.id] = 1;
    for (int u = 0; u < MC.VV[tmp.id].size(); u++) {
      int vid = MC.VV[tmp.id][u];
      int eid = MC.VE[tmp.id][u];
      if (distance[tmp.id] + MC.EL[eid] < distance[vid]) {
        distance[vid] = distance[tmp.id] + MC.EL[eid];
        que.push(node(vid, distance[vid]));
        v_p[vid] = tmp.id;
      }
    }
  }
}

void Algorithm::Dijkstra_group(MeshCache &MC, std::vector<int> &lmk) {
  std::sort(lmk.begin(), lmk.end());
#pragma omp parallel for
  for (int i = 0; i < lmk.size() - 1; i++) {
    std::vector<int> is_lmk(MC.NVertices, 0);
    for (int j = i; j < lmk.size(); j++)
      is_lmk[lmk[j]] = 1;
    int count = lmk.size() - i;
    int s_p = lmk[i];
    std::vector<int> &is_visited = MC.DijkstraIsVisited[s_p];
    std::vector<double> &distance = MC.Vd[s_p];
    std::priority_queue<node> &que = MC.DijkstraCache[s_p];
    std::vector<int> &v_p = MC.VVp[s_p];
    if (is_visited.size() == 0) {
      is_visited.resize(MC.NVertices, 0);
      distance.resize(MC.NVertices, DBL_MAX);
      distance[s_p] = 0;
      que.push(node(s_p, 0));
      v_p.resize(MC.NVertices, -1);
      v_p[s_p] = s_p;
    }
    for (int i = 0; i < lmk.size(); i++) {
      if (is_visited[lmk[i]] != 0)
        count--;
    }
    if (count < 0)
      continue;
    while (count != 0) {
      if (que.size() == 0)
        printf("debug\n");
      node tmp = que.top();
      que.pop();

      if (is_visited[tmp.id] == 1) {
        continue;
      }
      if (is_lmk[tmp.id] == 1) {
        count--;
        is_lmk[tmp.id] = 0;
      }
      is_visited[tmp.id] = 1;
      for (int u = 0; u < MC.VV[tmp.id].size(); u++) {
        int v_id = MC.VV[tmp.id][u];
        int e_id = MC.VE[tmp.id][u];
        if (distance[tmp.id] + MC.EL[e_id] < distance[v_id]) {
          distance[v_id] = distance[tmp.id] + MC.EL[e_id];
          que.push(node(v_id, distance[v_id]));
          v_p[v_id] = tmp.id;
        }
      }
    }
  }
}

void Algorithm::Kruskal(MeshCache &MCache, std::vector<int> &lmk,
                        std::vector<int> &cutvertex,
                        std::vector<int> &cutedge) {
  std::sort(lmk.begin(), lmk.end());
  int nv = lmk.size();
  std::vector<int> spanning_tree_current(nv);
  std::priority_queue<PathInfo> EdgeInfo;
  for (int i = 0; i < lmk.size(); i++)
    for (int j = i + 1; j < lmk.size(); j++) {
      EdgeInfo.push(PathInfo(lmk[i], lmk[j], MCache.Vd[lmk[i]][lmk[j]]));
    }
  std::vector<int> edge;
  std::set<int> vertex;
  for (int i = 0; i < nv; ++i) {
    spanning_tree_current[i] = i;
  }
  int index = 0;
  int j = 0;
  while (index < nv - 1) {
    PathInfo tmp = EdgeInfo.top();
    EdgeInfo.pop();
    int m_id = -1, n_id = -1;
    for (int u = 0; u < lmk.size(); u++) {
      if (lmk[u] == tmp.s_p)
        m_id = u;
      if (lmk[u] == tmp.e_p)
        n_id = u;
    }
    m_id = spanning_tree_current[m_id];
    n_id = spanning_tree_current[n_id];
    if (m_id < n_id)
      std::swap(m_id, n_id);
    if (m_id != n_id) {
      std::vector<int> path;
      Algorithm::FindPath(MCache.VVp[tmp.s_p], tmp.e_p, path);
      vertex.insert(path[0]);
      for (int i = 0; i < path.size() - 1; i++) {
        vertex.insert(path[i + 1]);
        edge.push_back(MCache.VVE[path[i]][path[i + 1]]);
      }
      index++;
      for (int i = 0; i < nv; ++i) {
        if (spanning_tree_current[i] == n_id)
          spanning_tree_current[i] = m_id;
      }
    }
  }
  std::set<int> e;
  for (auto a : edge)
    e.insert(a);
  for (auto a : e) {
    cutedge.push_back(a);
  }
  for (auto a : vertex)
    cutvertex.push_back(a);
}

void Algorithm::Kruskal(std::vector<int> &lmk,
                        std::priority_queue<PathInfo> que,
                        std::vector<PathInfo> &Res) {
  int nv = lmk.size();
  std::vector<int> spanning_tree_current(nv);
  std::set<int> reV, reE;

  std::vector<int> edge;
  std::set<int> vertex;
  for (int i = 0; i < nv; ++i) {
    spanning_tree_current[i] = i;
  }
  int index = 0;
  int j = 0;
  while (index < nv - 1) {
    PathInfo tmp = que.top();
    que.pop();
    int m_id = -1, n_id = -1;
    for (int u = 0; u < nv; u++) {
      if (lmk[u] == tmp.s_p)
        m_id = u;
      if (lmk[u] == tmp.e_p)
        n_id = u;
    }
    m_id = spanning_tree_current[m_id];
    n_id = spanning_tree_current[n_id];
    if (m_id < n_id)
      std::swap(m_id, n_id);
    if (m_id != n_id) {
      Res.push_back(tmp);
      index++;
      for (int i = 0; i < nv; ++i) {
        if (spanning_tree_current[i] == n_id)
          spanning_tree_current[i] = m_id;
      }
    }
  }
}

void Algorithm::FindPath(std::vector<int> &v_p, int e_p,
                         std::vector<int> &path) {
  path.clear();
  path.push_back(e_p);
  do {
    e_p = v_p[e_p];
    path.push_back(e_p);
  } while (e_p != v_p[e_p]);
}

void Algorithm::UpdateNeighbourhood(MeshCache &MC, int k, int v) {
  if (k <= MC.Max_Neighbour[v])
    return;
  if (k == 1) {
    if (MC.Neighbour[v].size() == 0) {
      MC.Neighbour[v].resize(MC.NVertices, INT_MAX);
      MC.Neighbour[v][v] = 0;
      for (auto a : MC.VV[v])
        MC.Neighbour[v][a] = 1;
      MC.Max_Neighbour[v] = 1;
    }
  } else {
    UpdateNeighbourhood(MC, k - 1, v);
    for (int u = 0; u < MC.NVertices; u++) {
      if (MC.Neighbour[v][u] == k - 1) {
        for (int s = 0; s < MC.VV[u].size(); s++) {
          if (MC.Neighbour[v][MC.VV[u][s]] > k)
            MC.Neighbour[v][MC.VV[u][s]] = k;
        }
      }
    }
    MC.Max_Neighbour[v] = k;
  }
}

void Algorithm::Dijkstra_with_restrict(MeshCache &MCache, int s_p,
                                       std::vector<double> &weight,
                                       std::vector<int> &v_p,
                                       std::vector<double> &d) {
  double total_length = MCache.AVG_EL * MCache.NEdges;
  std::vector<int> is_visited(MCache.NVertices, 0);
  d[s_p] = 0;
  is_visited[s_p] = 0;
  v_p[s_p] = s_p;
  std::priority_queue<node> que;
  que.push(node(s_p, 0));
  while (!que.empty()) {
    node tmp = que.top();
    que.pop();
    if (is_visited[tmp.id] == 1)
      continue;
    is_visited[tmp.id] = 1;
    for (int u = 0; u < MCache.VV[tmp.id].size(); u++) {
      int vid = MCache.VV[tmp.id][u];
      int eid = MCache.VE[tmp.id][u];
      if (d[tmp.id] + weight[eid] < d[vid]) {
        d[vid] = d[tmp.id] + weight[eid];
        que.push(node(vid, d[vid]));
        v_p[vid] = tmp.id;
      }
    }
  }
}

void Algorithm::Dijkstra_with_nearest2(MeshCache &MC, int s_p,
                                       std::vector<int> &is_target,
                                       std::vector<int> &path) {
  int e_p = -1;
  std::vector<double> distance = MC.Vd[s_p];
  std::vector<int> is_visited = MC.DijkstraIsVisited[s_p];
  std::priority_queue<node> que = MC.DijkstraCache[s_p];
  std::vector<int> v_p = MC.VVp[s_p];
  if (distance.size() != 0) {
    double min_dis = DBL_MAX;
    int min_id = -1;
    for (int i = 0; i < is_target.size(); i++) {
      if (is_target[i] == 1) {
        if (is_visited[i] == 0)
          continue;
        if (distance[i] < min_dis) {
          min_dis = distance[i];
          min_id = i;
        }
      }
    }
    if (min_id != -1) {
      e_p = min_id;
    }
  }
  if (e_p == -1) {

    if (is_visited.size() == 0) {
      is_visited.resize(MC.NVertices, 0);
      distance.resize(MC.NVertices, DBL_MAX);
      distance[s_p] = 0;
      que.push(node(s_p, 0));
      v_p.resize(MC.NVertices, -1);
      v_p[s_p] = s_p;
    }
    while (!que.empty()) {
      node tmp = que.top();
      que.pop();
      if (is_visited[tmp.id] == 1)
        continue;
      if (is_target[tmp.id] == 1) {
        e_p = tmp.id;
        break;
      }
      is_visited[tmp.id] = 1;
      for (int u = 0; u < MC.VV[tmp.id].size(); u++) {
        int v_id = MC.VV[tmp.id][u];
        int e_id = MC.VE[tmp.id][u];
        if (distance[tmp.id] + MC.EL[e_id] < distance[v_id]) {
          distance[v_id] = distance[tmp.id] + MC.EL[e_id];
          que.push(node(v_id, distance[v_id]));
          v_p[v_id] = tmp.id;
        }
      }
    }
  }
  FindPath(v_p, e_p, path);
}
