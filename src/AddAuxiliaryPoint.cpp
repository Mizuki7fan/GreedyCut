#include "AddAuxiliaryPoint.h"

AAP::AAP(Mesh &mesh, MeshCache &MC, std::vector<int> &landmark)
    : mesh(mesh), MC(MC), landmark(landmark) {}

void AAP::Set(double trimming_rate, int max_add_count) {
  this->TrimmingRate = trimming_rate;
  this->MaxAddCount = max_add_count;
}

void AAP::Run() {
  GenGraph();
  int add_count = 0;
  greed1_ok = false;
  while (!greed1_ok) {
    if (add_count == MaxAddCount)
      break;
    Greed1_Op();
    add_count++;
  }
}

AAP::~AAP() {}

void AAP::GenGraph() {
  for (int i = 0; i < landmark.size(); i++)
    Algorithm::Dijkstra_all(MC, landmark[i]);
  for (int i = 0; i < landmark.size(); i++)
    for (int j = i + 1; j < landmark.size(); j++) {
      PresavedPath.push(
          PathInfo(landmark[i], landmark[j], MC.Vd[landmark[i]][landmark[j]]));
    }
  std::vector<PathInfo> Mst;
  std::priority_queue<PathInfo> Back(PresavedPath);
  Algorithm::Kruskal(landmark, Back, Mst);
  double TLength = 0;
  for (auto a : Mst) {
    TLength += a.length;
  }
  this->TreeLength = TLength;
}

void AAP::Greed1_Op() {
  std::vector<double> vertex_weight(MC.NVertices);
#pragma omp parallel for
  for (int i = 0; i < MC.NVertices; i++) {
    if (std::find(landmark.begin(), landmark.end(), i) != landmark.end())
      vertex_weight[i] = this->TreeLength;
    else {
      double TLength = 0;
      std::priority_queue<PathInfo> Path(PresavedPath);
      std::vector<int> e_landmarks(landmark);
      e_landmarks.push_back(i);
      for (int j = 0; j < landmark.size(); j++) {
        Path.push(PathInfo(landmark[j], i, MC.Vd[landmark[j]][i]));
      }
      std::vector<PathInfo> Mst;
      Algorithm::Kruskal(e_landmarks, Path, Mst);
      for (auto a : Mst)
        TLength += a.length;
      vertex_weight[i] = TLength;
    }
  }
  int min_id = -1;
  double min_dis = DBL_MAX;
  for (int i = 0; i < MC.NVertices; i++) {
    if (vertex_weight[i] < min_dis) {
      min_dis = vertex_weight[i];
      min_id = i;
    }
  }
  if (TreeLength - min_dis < TrimmingRate * MC.Lbb + 1e-10)
    min_id = -1;
  if (min_id == -1)
    greed1_ok = true;
  else {
    landmark.push_back(min_id);
    greed1_ok = false;
    Algorithm::Dijkstra_all(MC, min_id);
    for (int i = 0; i < landmark.size() - 1; i++) {
      PresavedPath.push(
          PathInfo(min_id, landmark[i], MC.Vd[landmark[i]][min_id]));
    }
    std::vector<PathInfo> Mst;
    Algorithm::Kruskal(landmark, PresavedPath, Mst);
    double TLength = 0;
    for (auto a : Mst) {
      TLength += a.length;
    }
    this->TreeLength = TLength;
  }
}
