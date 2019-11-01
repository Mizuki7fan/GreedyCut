#include "PointFinding.h"

PointFinding::PointFinding(const Mesh&m, const Mesh&pm):OriMesh(m),ParaedMesh(pm)
{
}

PointFinding::~PointFinding()
{
}

void PointFinding::Set(std::string metric)
{
	this->Metric = metric;
}

void PointFinding::Find(std::vector<std::pair<int, double>>& result)
{//主体函数，用于找点
	PrepareComputeDistortion();
	ComputeFaceDistortion();
	ComputeVertexDistortion();
	if (Metric == "RealDis")
		FindByRealDis(result);










}

void PointFinding::PrepareComputeDistortion()
{
	facearea.resize(OriMesh.n_faces());
	fpx1.resize(OriMesh.n_faces());
	fpx2.resize(OriMesh.n_faces());
	fpy2.resize(OriMesh.n_faces());
	for (const auto& fh : OriMesh.faces())
	{
		auto heh = OriMesh.halfedge_handle(fh);
		const auto& p0 = OriMesh.point(OriMesh.from_vertex_handle(heh));
		const auto& p1 = OriMesh.point(OriMesh.to_vertex_handle(heh));
		const auto& p2 = OriMesh.point(OriMesh.to_vertex_handle(OriMesh.next_halfedge_handle(heh)));
		auto np = (p1 - p0) % (p2 - p0);
		double area = np.norm() * 0.5;
		np.normalize();
		fpx1[fh.idx()] = (p1 - p0).norm();
		auto pe1 = (p1 - p0).normalized();
		auto pe2 = (np % pe1).normalized();;
		fpx2[fh.idx()] = (p2 - p0) | pe1;
		fpy2[fh.idx()] = (p2 - p0) | pe2;
		facearea[fh.idx()] = area;
	}

}

void PointFinding::ComputeFaceDistortion()
{
	double alpha = 0.5;
	if (OriMesh.n_faces() != ParaedMesh.n_faces())
	{
		return;
	}
	facedistortion.resize(OriMesh.n_faces());
	for (const auto& fh : ParaedMesh.faces())
	{
		auto fid = fh.idx();
		auto heh = ParaedMesh.halfedge_handle(fh);
		const auto& p0 = ParaedMesh.point(ParaedMesh.from_vertex_handle(heh));
		const auto& p1 = ParaedMesh.point(ParaedMesh.to_vertex_handle(heh));
		const auto& p2 = ParaedMesh.point(ParaedMesh.to_vertex_handle(ParaedMesh.next_halfedge_handle(heh)));
		auto np = (p1 - p0) % (p2 - p0);
		double area = np.norm() * 0.5;
		np.normalize();
		auto qx1 = (p1 - p0).norm();
		auto pe1 = (p1 - p0).normalized();
		auto pe2 = (np % pe1).normalized();;
		auto qx2 = (p2 - p0) | pe1;
		auto qy2 = (p2 - p0) | pe2;
		double a = qx1 / fpx1[fid];
		double b = qy2 / fpy2[fid];
		double det = a * b;
		if (det > 0.0)
		{
			double c = ((-qx1 * fpx2[fid] + qx2 * fpx1[fid]) / fpx1[fid] / fpy2[fid]);
			facedistortion[fid] = alpha * (a * a + b * b + c * c) / det * 0.5 + (1 - alpha) * (det + 1.0 / det) * 0.5;
			//facedistortion[fid] = (det + 1.0 / det) * 0.5;
			//facedistortion[fid] = 1.0 / det;
		}
		else
		{
			facedistortion[fid] = std::numeric_limits<double>::infinity();
		}
	}
}

void PointFinding::ComputeVertexDistortion()
{
	//根据面片扭曲计算顶点扭曲
	vertex_distortion.resize(OriMesh.n_vertices());
	for (int i = 0; i < vertex_distortion.size(); i++)
	{
		Mesh::VertexHandle v = OriMesh.vertex_handle(i);
		double v_dis = 0;
		int v_valence = 0;
		for (const auto& vfh : OriMesh.vf_range(v))
		{
			v_dis += facedistortion[vfh.idx()];
			++v_valence;
		}
		vertex_distortion[i] = v_dis / v_valence;
	}
}

void PointFinding::FindByRealDis(std::vector<std::pair<int, double>>&)
{
	//输入



}

