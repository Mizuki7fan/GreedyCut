#include "KPNewton.h"
#include <unordered_map>
#define M_2_SQRTPI 1.12837916709551257390   // 2/sqrt(pi)

KPNewton::KPNewton(Mesh& m)
	:mesh(m),
	isbv(m.n_vertices(), 0)
{
	for (const auto& heh : mesh.halfedges())
	{
		if (mesh.is_boundary(heh))
		{
			auto nextheh = heh;
			do
			{
				isbv[mesh.from_vertex_handle(nextheh).idx()] = 1;
				nextheh = mesh.prev_halfedge_handle(nextheh);
			} while (nextheh != heh);
			break;
		}
	}

}

KPNewton::~KPNewton()
{
}

void KPNewton::Tutte()
{
	auto boundaryvhs = GetBoundary();
	double delta_angle = 2 * M_PI / boundaryvhs.size();
	double area_1_factor = 0.5 * M_2_SQRTPI;
	auto nv = mesh.n_vertices();
	std::vector<double> posx(nv), posy(nv);
	for (size_t i = 0; i < boundaryvhs.size(); ++i)
	{
		posx[boundaryvhs[i].idx()] = area_1_factor * cos(i * delta_angle);
		posy[boundaryvhs[i].idx()] = area_1_factor * sin(i * delta_angle);
	}
	auto& pardiso_it = solver.ia;
	auto& pardiso_jt = solver.ja;
	auto& pardiso_t = solver.a;
	std::vector<double> pardiso_tu;
	std::vector<double> pardiso_tv;

	pardiso_it.reserve(nv + 1);
	pardiso_jt.reserve(6 * nv);
	pardiso_t.reserve(6 * nv);
	pardiso_tu.resize(nv, 0.0);
	pardiso_tv.resize(nv, 0.0);
	//这里其实不太好并行，无法预先知道一个pardiso_jt[i]该选择哪个分支
	for (const auto& vh : mesh.vertices())
	{
		pardiso_it.push_back(static_cast<int>(pardiso_jt.size()));
		auto vid = vh.idx();
		if (mesh.is_boundary(vh))
		{
			pardiso_jt.push_back(vid);
			pardiso_t.push_back(1);
			pardiso_tu[vh.idx()] = posx[vid];
			pardiso_tv[vh.idx()] = posy[vid];
		}
		else
		{
			pardiso_jt.push_back(vid);
			pardiso_t.push_back(mesh.valence(vh));
			std::vector<int> row_id;
			row_id.reserve(mesh.valence(vh));
			for (const auto& vvh : mesh.vv_range(vh))
			{
				int vvid = vvh.idx();
				if (mesh.is_boundary(vvh))
				{
					pardiso_tu[vid] += posx[vvid];
					pardiso_tv[vid] += posy[vvid];
				}
				else
				{
					if (vvid > vid)
					{
						row_id.push_back(vvid);
					}
				}
			}
			std::sort(row_id.begin(), row_id.end(), std::less<int>());
			for (size_t j = 0; j < row_id.size(); j++)
			{
				pardiso_jt.push_back(row_id[j]);
				pardiso_t.push_back(-1);
			}
		}
	}
	pardiso_it.push_back(static_cast<int>(pardiso_jt.size()));
	solver.num = static_cast<int>(nv);
	solver.nnz = pardiso_jt.size();

	solver.pardiso_init();
	solver.factorize();
	solver.rhs = pardiso_tu;
	solver.pardiso_solver();
	posx = solver.result;
	solver.rhs = pardiso_tv;
	solver.pardiso_solver();
	posy = solver.result;
	result.clear();
	result.reserve(nv);
	for (size_t i = 0; i < nv; i++)
	{
		result.push_back({ posx[i], posy[i] });
	}
}

void KPNewton::PrepareDataFree(void)
{
	auto&& ia = solver.ia;
	auto&& ja = solver.ja;
	int nv = static_cast<int>(mesh.n_vertices());
	ia.clear();
	ja.clear();
	ia.reserve(2 * nv + 1);
	ja.reserve(16 * nv);
	std::vector<std::vector<int>> tmp_ja(nv);

#pragma omp parallel for num_threads(16)
	for (int vid = 0; vid < nv; vid++)
	{
		std::vector<int> rowid;

		rowid.push_back(vid + 1);
		rowid.push_back(vid + nv + 1);
		auto vh = mesh.vertex_handle(vid);
		for (const auto& vvh : mesh.vv_range(vh))
		{
			int vvid = vvh.idx();
			if (vvid > vid)
			{
				rowid.push_back(vvid + 1);
			}
			rowid.push_back(vvid + nv + 1);
		}
		std::sort(rowid.begin(), rowid.end(), std::less<int>());
		tmp_ja[vid].insert(tmp_ja[vid].end(), rowid.begin(), rowid.end());
	}
	ia.reserve(tmp_ja.size());
	int cnt = 1;
	for (int i = 0; i < tmp_ja.size(); i++)
	{
		ia.push_back(cnt);
		cnt += tmp_ja[i].size();
		ja.insert(ja.end(), tmp_ja[i].begin(), tmp_ja[i].end());

	}

	tmp_ja.clear();
	tmp_ja.resize(nv);
#pragma omp parallel for num_threads(16)
	for (int vid = 0; vid < mesh.n_vertices(); vid++)
	{
		std::vector<int> rowid;
		rowid.push_back(vid + nv + 1);
		auto vh = mesh.vertex_handle(vid);
		for (const auto& vvh : mesh.vv_range(vh))
		{
			int vvid = vvh.idx();
			if (vvid > vid)
			{
				rowid.push_back(vvid + nv + 1);
			}
		}
		std::sort(rowid.begin(), rowid.end(), std::less<int>());
		tmp_ja[vid].insert(tmp_ja[vid].end(), rowid.begin(), rowid.end());
	}
	ia.reserve(ia.size() + tmp_ja.size());
	for (int i = 0; i < tmp_ja.size(); i++)
	{
		ia.push_back(cnt);
		cnt += tmp_ja[i].size();
		ja.insert(ja.end(), tmp_ja[i].begin(), tmp_ja[i].end());
	}
	ia.push_back(static_cast<int>(ja.size()) + 1);
	solver.a.resize(ja.size());
	solver.num = static_cast<int>(2 * mesh.n_vertices());
	solver.nnz = ja.size();

	std::vector<std::unordered_map<int, int>> spid(2 * nv);
	size_t j = 0;
#pragma omp parallel for num_threads(16)
	for (int i = 0; i < ia.size() - 1; i++)
	{
		for (int j = ia[i]; j < ia[i + 1]; j++)
		{
			spid[i][ja[j - 1] - 1] = j - 1;
		}
	}
	tri.resize(mesh.n_faces());
#pragma omp parallel for num_threads(16)
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		auto fh = mesh.face_handle(i);
		for (const auto& fvh : mesh.fv_range(fh))
		{
			tri[i].push_back(fvh.idx());
		}
	}
	assembleorder.clear();
	std::vector<std::vector<int>> tmp_assembleorder(tri.size());
#pragma omp parallel for num_threads(16)
	for (int fid = 0; fid < tri.size(); fid++)
	{
		const auto& vhs = tri[fid];
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (vhs[i] <= vhs[j])
				{
					tmp_assembleorder[fid].push_back(spid[vhs[i]][vhs[j]]);
					tmp_assembleorder[fid].push_back(spid[vhs[i] + nv][vhs[j] + nv]);
				}
				tmp_assembleorder[fid].push_back(spid[vhs[i]][vhs[j] + nv]);
			}
		}
	}
	for (auto a : tmp_assembleorder)
	{
		assembleorder.insert(assembleorder.end(), a.begin(), a.end());
	}
	for (auto& a : ia)
	{
		a--;
	}
	for (auto& a : ja)
	{
		a--;
	}
	solver.pardiso_init();
}

void KPNewton::RunFree(const EnergyType& etype)
{
	if (EnergyIsNan)
		return;
	//mode决定是否要加面积权
	auto maxiter = 100;
	double Iscale = 1e-8;
	double ls_StepFac = 0.5;
	double ls_GradFac = 0.02;
	double ls_MinStep = 1e-7;
	switch (etype)
	{
	default:
		break;
	case MIPS:
		computeall = ComputeMIPS;
		computeenergy = ComputeEnergyMIPS;
		break;
	case AMIPS:
		computeall = ComputeAMIPS;
		computeenergy = ComputeEnergyAMIPS;
		break;
	}
	if (result.empty())
	{
		std::cerr << "Error: No initial position." << std::endl;
		return;
	}
	// prepare inputs
	auto nV = mesh.n_vertices();
	auto nF = mesh.n_faces();

	std::vector<std::vector<double>> faceElens = MeshFaceEdgeLen2s();
	std::vector<std::vector<double>> faceAngles = MeshAnglesFromFaceEdgeLen2(faceElens);

	std::vector<std::vector< std::complex<double>>> localframe(nF, std::vector<std::complex<double>>(3));
	std::vector<std::vector< std::complex<double>>> D(nF, std::vector<std::complex<double>>(3));
	std::vector<std::vector< std::complex<double>>> DC(nF, std::vector<std::complex<double>>(3));
	std::vector<double> area(nF);
	double totalarea = 0;
#pragma omp parallel for num_threads(16)
	for (int i = 0; i < nF; i++)
	{
		auto&& frame = localframe[i];
		frame[1] = std::sqrt(faceElens[i][2]);
		frame[2] = std::polar(std::sqrt(faceElens[i][1]), faceAngles[i][0]);
		//是否使用面积权
		area[i] = frame[1].real() * frame[2].imag() / 2;
		//area[i] = 1;
		DC[i][0] = std::complex<double>(0.0, -0.25) * (frame[1] - frame[2]) / area[i];
		DC[i][1] = std::complex<double>(0.0, -0.25) * (frame[2] - frame[0]) / area[i];
		DC[i][2] = std::complex<double>(0.0, -0.25) * (frame[0] - frame[1]) / area[i];
		D[i][0] = std::conj(DC[i][0]);
		D[i][1] = std::conj(DC[i][1]);
		D[i][2] = std::conj(DC[i][2]);
	}
	for (int i = 0; i < area.size(); i++)
		totalarea += area[i];
	for (int i = 0; i < area.size(); i++)
		area[i] /= totalarea;

	for (int iter = 0; iter < maxiter; iter++)
	{
		auto&& a = solver.a;
		a.clear();
		a.resize(solver.ja.size(), 0);
		double e = 0;
		std::vector<double> g(2 * nV, 0);
		std::vector<double> b(2 * nV, 0);
		int nele = 0;

		for (size_t fid = 0; fid < tri.size(); fid++)
		{
			std::complex<double> fz, fzb;
			const auto& vhs = tri[fid];
			for (size_t i = 0; i < 3; i++)
			{
				fz += D[fid][i] * result[vhs[i]];
				fzb += DC[fid][i] * result[vhs[i]];
			}
			double x = std::norm(fz);
			double y = std::norm(fzb);
			double etemp, alpha1, alpha2, beta1, beta2, beta3;
			computeall(etemp, alpha1, alpha2, beta1, beta2, beta3, x, y);
			e += etemp * area[fid];
			for (size_t i = 0; i < 3; i++)
			{
				auto gi = 2.0 * (alpha1 * DC[fid][i] * fz + alpha2 * D[fid][i] * fzb) * area[fid];
				g[vhs[i]] += gi.real();
				g[vhs[i] + nV] += gi.imag();
			}
			{// fix alpha and beta
				double temp1 = alpha1 + 2 * beta1 * x;
				double temp2 = alpha2 + 2 * beta2 * y;
				double s1 = temp1 + temp2;
				double s2 = temp1 - temp2;
				double lambda3 = s1 + std::sqrt(s2 * s2 + 16 * beta3 * beta3 * x * y);
				double lambda4 = s1 - std::sqrt(s2 * s2 + 16 * beta3 * beta3 * x * y);
				double t1 = (lambda3 - 2 * alpha1 - 4 * beta1 * x) / (4 * beta3 * y);
				double t2 = (lambda4 - 2 * alpha1 - 4 * beta1 * x) / (4 * beta3 * y);
				double e3norm2 = x + y * t1 * t1;
				double e4norm2 = x + y * t2 * t2;
				lambda3 = (lambda3 > 0 ? lambda3 : 0) / e3norm2;
				lambda4 = (lambda4 > 0 ? lambda4 : 0) / e4norm2;
				if (y > 1e-50 && fabs(beta1 * beta2 * beta3) > 1e-50)
				{
					alpha1 = alpha1 > 0 ? alpha1 : 0;
					alpha2 = alpha2 > 0 ? alpha2 : 0;
					beta1 = (lambda3 + lambda4 - alpha1 * 2.0 / x) / 4.0;
					beta2 = (lambda3 * t1 * t1 + lambda4 * t2 * t2 - alpha2 * 2.0 / y) / 4.0;
					beta3 = (lambda3 * t1 + lambda4 * t2) / 4.0;
				}
			}

			double a1 = area[fid] * (2.0 * alpha1 + 4.0 * beta1 * fz.real() * fz.real());
			double a2 = area[fid] * 4.0 * beta1 * fz.real() * fz.imag();
			double a4 = area[fid] * (2.0 * alpha1 + 4.0 * beta1 * fz.imag() * fz.imag());
			double b1 = area[fid] * 4.0 * beta3 * fzb.real() * fz.real();
			double b2 = area[fid] * 4.0 * beta3 * fzb.real() * fz.imag();
			double b3 = area[fid] * 4.0 * beta3 * fzb.imag() * fz.real();
			double b4 = area[fid] * 4.0 * beta3 * fzb.imag() * fz.imag();
			double d1 = area[fid] * (2.0 * alpha2 + 4.0 * beta2 * fzb.real() * fzb.real());
			double d2 = area[fid] * 4.0 * beta2 * fzb.real() * fzb.imag();
			double d4 = area[fid] * (2.0 * alpha2 + 4.0 * beta2 * fzb.imag() * fzb.imag());
			double s1 = a1 + 2.0 * b1 + d1;
			double s2 = a2 + b2 - b3 - d2;
			double s3 = a4 - 2.0 * b4 + d4;
			double s4 = a2 + b2 + b3 + d2;
			double s5 = d1 - a1;
			double s6 = a4 - d4;
			double s7 = -a2 + b2 + b3 - d2;
			double s8 = a4 + 2.0 * b4 + d4;
			double s9 = -a2 + b2 - b3 + d2;
			double s0 = a1 - 2.0 * b1 + d1;
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					double RR = D[fid][i].real() * D[fid][j].real();
					double RI = D[fid][i].real() * D[fid][j].imag();
					double IR = D[fid][i].imag() * D[fid][j].real();
					double II = D[fid][i].imag() * D[fid][j].imag();

					if (vhs[i] <= vhs[j])
					{
						a[assembleorder[nele++]] += s1 * RR + s2 * (RI + IR) + s3 * II;
						a[assembleorder[nele++]] += s8 * RR + s9 * (RI + IR) + s0 * II;
					}
					a[assembleorder[nele++]] += s4 * RR + s5 * RI + s6 * IR + s7 * II;
				}
			}
		}
		solver.factorize();
		double normg = 0.0;
		std::vector<double> tmp_normg(g.size());

#pragma omp parallel for num_threads(16)
		for (int i = 0; i < g.size(); i++)
		{
			b[i] = -g[i];
			tmp_normg[i] = g[i] * g[i];
		}
		for (auto a : tmp_normg)
			normg += a;
		normg /= std::sqrt(normg);
		std::vector<double> x;
		solver.rhs = b;
		solver.pardiso_solver();
		x = solver.result;
		std::vector<std::complex<double>> sDC(nV);
		std::vector<std::complex<double>> newpos(nV);
		double stepDirlen = 0.0;
		double energyGrad = 0.0;
#pragma omp parallel for num_threads(16)
		for (int i = 0; i < nV; i++)
		{
			sDC[i] = { x[i], x[i + nV] };
			newpos[i] = result[i] + sDC[i];
		}
		for (size_t i = 0; i < nV; i++)
		{
			stepDirlen += x[i] * x[i] + x[i + nV] * x[i + nV];
			energyGrad -= b[i] * x[i] + b[i + nV] * x[i + nV];
		}
		energyGrad *= ls_GradFac;
		stepDirlen = std::sqrt(stepDirlen);
		double lsa = ComputeTMax(result, sDC) * 0.9;
		lsa = lsa > 1 ? 1 : lsa;
#pragma omp parallel for num_threads(16)
		for (int i = 0; i < nV; i++)
		{
			newpos[i] = result[i] + sDC[i] * lsa;
		}
		double newEnergy = ComputeEnergy(newpos, D, DC, area);
		//能量值为正无穷
		if (isnan(newEnergy))
		{
			EnergyIsNan = true;
			break;

		}
			

		while (stepDirlen * lsa > ls_MinStep && e + lsa * energyGrad < newEnergy)
		{
			lsa = ls_StepFac * lsa;
#pragma omp parallel for num_threads(16)
			for (int i = 0; i < nV; i++)
			{
				newpos[i] = result[i] + sDC[i] * lsa;
			}
			newEnergy = ComputeEnergy(newpos, D, DC, area);
		}
		if (newEnergy > e || std::isnan(newEnergy) || std::isinf(newEnergy))
		{
		}
		else
		{
#pragma omp parallel for num_threads(16)
			for (int i = 0; i < nV; i++)
			{
				result[i] = newpos[i];
			}
		}
		if (normg < 1e-4 || (e - newEnergy) / e < 1e-4 || std::isnan(newEnergy) || std::isinf(newEnergy))
		{
			break;
		}

	}

}

void KPNewton::ComputeMIPS(double& energy, double& alpha1, double& alpha2, double& beta1, double& beta2, double& beta3, const double& x, const double& y)
{
	energy = (x + y) / (x - y);
	alpha1 = -2 * y / (x - y) / (x - y);
	alpha2 = 2 * x / (x - y) / (x - y);
	beta1 = -2 * alpha1 / (x - y);
	beta2 = 2 * alpha2 / (x - y);
	beta3 = (alpha1 - alpha2) / (x - y);
}

double KPNewton::ComputeEnergyMIPS(const double& x, const double& y)
{
	return (x + y) / (x - y);
}

void KPNewton::ComputeAMIPS(double& energy, double& alpha1, double& alpha2, double& beta1, double& beta2, double& beta3, const double& x, const double& y)
{
	energy = std::exp((x + y) / (x - y));
	alpha1 = -2 * y / (x - y) / (x - y) * energy;
	alpha2 = 2 * x / (x - y) / (x - y) * energy;
	beta1 = -2 * x * alpha1 / (x - y) / (x - y);
	beta2 = (4 * x - 2 * y) * alpha2 / (x - y) / (x - y);
	beta3 = ((2 * x - y) * alpha1 - x * alpha2) / (x - y) / (x - y);
}

double KPNewton::ComputeEnergyAMIPS(const double& x, const double& y)
{
	return std::exp((x + y) / (x - y));
}

double KPNewton::ComputeTMax(const std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& d) const
{
	double temp_t = std::numeric_limits<double>::infinity();
	std::vector<double> a(tri.size()), b(tri.size()), c(tri.size()), b1(tri.size()), b2(tri.size()), tt(tri.size()), tt1(tri.size()), tt2(tri.size());
	int V_N = (int)mesh.n_vertices();
#pragma omp parallel for num_threads(16)
	for (int i = 0; i < tri.size(); i++)
	{
		const auto f = tri[i];
		const auto& f0 = f[0];
		const auto& f1 = f[1];
		const auto& f2 = f[2];
		const auto& x0 = x[f0].real();
		const auto& x1 = x[f1].real();
		const auto& x2 = x[f2].real();
		const auto& x3 = x[f0].imag();
		const auto& x4 = x[f1].imag();
		const auto& x5 = x[f2].imag();
		const auto& d0 = d[f0].real();
		const auto& d1 = d[f1].real();
		const auto& d2 = d[f2].real();
		const auto& d3 = d[f0].imag();
		const auto& d4 = d[f1].imag();
		const auto& d5 = d[f2].imag();
		a[i] = (d1 - d0) * (d5 - d3) - (d4 - d3) * (d2 - d0);
		b1[i] = (d1 - d0) * (x5 - x3) + (x1 - x0) * (d5 - d3);
		b2[i] = (x4 - x3) * (d2 - d0) + (x2 - x0) * (d4 - d3);
		b[i] = b1[i] - b2[i];
		c[i] = (x1 - x0) * (x5 - x3) - (x4 - x3) * (x2 - x0);
		tt[i] = std::numeric_limits<double>::infinity();
		//tt = 10000;
		if (b[i] * b[i] - 4 * a[i] * c[i] >= 0)
		{
			tt1[i] = 1 / (2 * a[i]) * (-b[i] + sqrt(b[i] * b[i] - 4 * a[i] * c[i]));
			tt2[i] = 1 / (2 * a[i]) * (-b[i] - sqrt(b[i] * b[i] - 4 * a[i] * c[i]));
			if (tt1[i] > 0 && tt2[i] > 0)
			{
				tt[i] = std::min(tt1[i], tt2[i]);
			}
			if (tt1[i] > 0 && tt2[i] < 0)
			{
				tt[i] = tt1[i];
			}
			if (tt1[i] < 0 && tt2[i] > 0)
			{
				tt[i] = tt2[i];
			}
		}
	}

	for (int i = 0; i < tt.size(); i++)
	{
		if (temp_t > tt[i])
		{
			temp_t = tt[i];
		}
	}
	return temp_t;
}

double KPNewton::ComputeEnergy(const std::vector<std::complex<double>>& pos,
	const std::vector<std::vector<std::complex<double>>>& D,
	const std::vector<std::vector<std::complex<double>>>& DC,
	const std::vector<double>& area) const
{
	double energy = 0;
	std::vector<double> tmp_energy(tri.size(), 0);
	int flag = 0;

#pragma omp parallel for num_threads(16)
	for (int fid = 0; fid < tri.size(); fid++)
	{
		if (flag == 0)
		{
			std::complex<double> fz = 0;
			std::complex<double> fzb = 0;
			for (size_t i = 0; i < 3; i++)
			{
				fz += D[fid][i] * pos[tri[fid][i]];
				fzb += DC[fid][i] * pos[tri[fid][i]];
			}
			double x = std::norm(fz);
			double y = std::norm(fzb);
			if (x < y)
			{
				flag = 1;
			}
			tmp_energy[fid] = computeenergy(x, y) * area[fid];
		}
	}
	if (flag == 1)
		return std::numeric_limits<double>::infinity();
	for (auto a : tmp_energy)
	{
		energy += a;
	}

	return energy;
}

void KPNewton::UpdateMesh(void)
{
#pragma omp parallel for num_threads(16)
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		const auto& vh = mesh.vertex_handle(i);
		mesh.set_point(vh, { result[vh.idx()].real(), result[vh.idx()].imag(), 0.0 });
	}
}

void KPNewton::ResultRect(const Mesh& orimesh)
{
	double scale = std::sqrt(MeshTools::Area(orimesh) / MeshTools::Area(mesh));
#pragma omp parallel for num_threads(16)
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		const auto& vh = mesh.vertex_handle(i);
		mesh.point(vh) *= scale;
	}
	OpenMesh::IO::write_mesh(mesh, "C:\\meigyoku\\exp5\\1.obj");
}

std::vector<Mesh::VertexHandle> KPNewton::GetBoundary(void) const
{
	std::vector<Mesh::VertexHandle> vhs;
	for (const auto& heh : mesh.halfedges())
	{
		if (mesh.is_boundary(heh))
		{
			auto nextheh = heh;
			do
			{
				vhs.push_back(mesh.from_vertex_handle(nextheh));
				nextheh = mesh.prev_halfedge_handle(nextheh);
			} while (nextheh != heh);
			break;
		}
	}
	return vhs;
}

std::vector<std::vector<double>> KPNewton::MeshFaceEdgeLen2s(void) const
{
	std::vector<std::vector<double>> len(mesh.n_faces(), std::vector<double>(3));
#pragma omp parallel for num_threads(16)
	for (int i = 0; i < len.size(); i++)
	{
		auto fh = mesh.face_handle(i);
		auto fvit = mesh.fv_iter(fh);
		const auto& p0 = mesh.point(*fvit);
		fvit++;
		const auto& p1 = mesh.point(*fvit);
		fvit++;
		const auto& p2 = mesh.point(*fvit);
		len[i] = { (p2 - p1).sqrnorm(), (p0 - p2).sqrnorm(), (p1 - p0).sqrnorm() };
	}
	return len;
}

std::vector<std::vector<double>> KPNewton::MeshAnglesFromFaceEdgeLen2(const std::vector<std::vector<double>>& len2) const
{
	std::vector<std::vector<double>> ang;
	ang.resize(len2.size(), std::vector<double>(3));
#pragma omp parallel for num_threads(16)
	for (int i = 0; i < len2.size(); i++)
	{
		ang[i] =
		{
		std::acos((len2[i][1] + len2[i][2] - len2[i][0]) / std::sqrt(len2[i][1] * len2[i][2]) / 2.0),
		std::acos((len2[i][2] + len2[i][0] - len2[i][1]) / std::sqrt(len2[i][2] * len2[i][0]) / 2.0),
		std::acos((len2[i][0] + len2[i][1] - len2[i][2]) / std::sqrt(len2[i][0] * len2[i][1]) / 2.0)
		};
	}
	return ang;
}