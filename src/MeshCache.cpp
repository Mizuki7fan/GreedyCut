#include "MeshCache.h"

MeshCache::MeshCache(Mesh &mesh) {
  NVertices = mesh.n_vertices();
  NEdges = mesh.n_edges();
  NFaces = mesh.n_faces();

  Vd.resize(NVertices);
  DijkstraIsVisited.resize(NVertices);
  DijkstraCache.resize(NVertices);
  VVp.resize(NVertices);

  VV.resize(NVertices, std::vector<int>());
  VE.resize(NVertices, std::vector<int>());
  VVE.resize(NVertices);

  Vx.resize(NVertices);
  Vy.resize(NVertices);
  Vz.resize(NVertices);
  for (const auto &v : mesh.vertices()) {
    Mesh::Point p = mesh.point(v);
    Vx[v.idx()] = p.data()[0];
    Vy[v.idx()] = p.data()[1];
    Vz[v.idx()] = p.data()[2];

    for (const auto &viheh : mesh.vih_range(v)) {
      Mesh::VertexHandle v1 = mesh.from_vertex_handle(viheh);
      Mesh::EdgeHandle e = mesh.edge_handle(viheh);
      VV[v.idx()].push_back(v1.idx());
      VE[v.idx()].push_back(e.idx());
      VVE[v.idx()][v1.idx()] = e.idx();
    }
  }

  EL.resize(NEdges);
  EV.resize(NEdges, std::vector<int>(2));
  AVG_EL = 0;
  for (const auto &e : mesh.edges()) {
    Mesh::HalfedgeHandle he = mesh.halfedge_handle(e, 0);
    EL[e.idx()] = mesh.calc_edge_length(e);
    EV[e.idx()][0] = mesh.to_vertex_handle(he).idx();
    EV[e.idx()][1] = mesh.from_vertex_handle(he).idx();
    AVG_EL += EL[e.idx()];
  }
  AVG_EL /= NEdges;

  Neighbour.resize(NVertices);
  Max_Neighbour.resize(NVertices, 0);

  FA.resize(NFaces);
  FV.resize(NFaces, std::vector<int>(3));
  AVG_FA = 0;

  for (int i = 0; i < mesh.n_faces(); i++) {
    const auto &f = mesh.face_handle(i);
    auto heh = mesh.halfedge_handle(f);
    const auto &p0 = mesh.point(mesh.from_vertex_handle(heh));
    const auto &p1 = mesh.point(mesh.to_vertex_handle(heh));
    const auto &p2 =
        mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)));
    FA[i] = 0.5 * ((p1 - p0) % (p2 - p0)).norm();
    AVG_FA += FA[i];
    Mesh::FaceVertexIter fvh = mesh.fv_begin(f);
    FV[i][0] = fvh->idx();
    fvh++;
    FV[i][1] = fvh->idx();
    fvh++;
    FV[i][2] = fvh->idx();
  }

  AVG_FA /= NFaces;

  Mesh::Point ptMin, ptMax;
  ptMin[0] = ptMin[1] = ptMin[2] = DBL_MAX;
  ptMax[0] = ptMax[1] = ptMax[2] = -DBL_MAX;
  for (const auto &vh : mesh.vertices()) {
    ptMin.minimize(mesh.point(vh));
    ptMax.maximize(mesh.point(vh));
  }
  Lbb = (ptMax - ptMin).norm();
}

MeshCache::~MeshCache() {}

double MeshTools::Area(const Mesh &mesh) {
  std::vector<double> area(mesh.n_faces(), 0.0);
#pragma omp parallel for
  for (int i = 0; i < mesh.n_faces(); i++) {
    const auto &f = mesh.face_handle(i);
    auto heh = mesh.halfedge_handle(f);
    const auto &p0 = mesh.point(mesh.from_vertex_handle(heh));
    const auto &p1 = mesh.point(mesh.to_vertex_handle(heh));
    const auto &p2 =
        mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)));
    area[i] = 0.5 * ((p1 - p0) % (p2 - p0)).norm();
  }
  double a = 0.0;
  for (auto b : area)
    a += b;
  return a;
}
