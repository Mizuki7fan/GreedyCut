#include "MeshDefinition.h"
double MeshTools::Area(const Mesh& mesh)
{


	std::vector<double> area(mesh.n_faces(),0.0);
#pragma omp parallel for num_threads(16)
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		const auto& f = mesh.face_handle(i);
		auto heh = mesh.halfedge_handle(f);
		const auto& p0 = mesh.point(mesh.from_vertex_handle(heh));
		const auto& p1 = mesh.point(mesh.to_vertex_handle(heh));
		const auto& p2 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)));
		area[i]= 0.5 * ((p1 - p0) % (p2 - p0)).norm();
	}
	double a = 0.0;
	for (auto b : area)
		a += b;
	return a;
}