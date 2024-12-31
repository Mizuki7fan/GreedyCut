#pragma once
#include "MeshCache.h"
#include <fstream>
#include <string>

struct Option {
  Option(std::string modelPath, std::string optPath);

  std::string modelName;
  std::string modelPath;
  std::string PS_method = "Dijkstra";
  std::string PF_method = "Dijkstra";
  std::string BanArea_Method = "NonConnect";
  int BanArea_Dn = 10;
  double BanArea_alpha = 0.1;
  std::string BanArea_Metric = "Dijkstra";
  double BanArea_ShrinkRate = 0.9;
  int Influence_Threshold = 20;
  double Distortion_Threshold = 0.01;
  double Trimming_Rate = 0.01;
  int Max_AddCount = 30;
};

class Algorithm {
  Algorithm();
  ~Algorithm();

public:
  static void Dijkstra_all(MeshCache &MCache, int k);
  static void Dijkstra_group(MeshCache &MCache, std::vector<int> &lmk);
  static void Kruskal(MeshCache &MCache, std::vector<int> &lmk,
                      std::vector<int> &cutvertex, std::vector<int> &cutedge);
  static void Kruskal(std::vector<int> &lmk, std::priority_queue<PathInfo> que,
                      std::vector<PathInfo> &Result);
  static void FindPath(std::vector<int> &v_p, int e_p,
                       std::vector<int> &path_v);
  static void UpdateNeighbourhood(MeshCache &MCache, int k, int v);
  static void Dijkstra_with_restrict(MeshCache &MCache, int s_p,
                                     std::vector<double> &weight,
                                     std::vector<int> &, std::vector<double> &);
  static void Dijkstra_with_nearest2(MeshCache &MCache, int s_p,
                                     std::vector<int> &is_target,
                                     std::vector<int> &path);
};
