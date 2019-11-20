#include "GAP.h"

GAP::GAP(Mesh& mesh, MeshCache& MCache) :
	Closedmesh(mesh), MCache(MCache)
{

}

void GAP::Set(const std::vector<std::pair<int, double>>& FP,
	std::string PriorityMetric,
	int FixThreshold,
	int ForbiddenRadius,
	double FilteringRate,
	int ParrCount)
{
	FeaturePoints = FP;
	this->VertexPriorityMetric = PriorityMetric;
	this->FixThreshold = FixThreshold;
	this->GAPForBiddenRadius = ForbiddenRadius;
	this->filtering_rate = FilteringRate;
	this->parr_count = ParrCount;
}

void GAP::SetPardiso(double conv_rate, int MaxIter, int bound_distortion_K)
{
	pardiso = NULL;
	convgence_con_rate = conv_rate;
	MAX_ITER_NUM = MaxIter;
	this->bound_distortion_K = bound_distortion_K;
}

void GAP::gradually_addp_pipeline()
{
	run_bpe();
	bool is_moveon = true;
	int i = 0;
	while (is_moveon)
	{
		i++;
		is_moveon = add1p();
	}
	//OpenMesh::IO::write_mesh(mesh, "result.obj");
	//for (auto it1 = mesh.vertices_begin(); it1 != mesh.vertices_end(); ++it1)
	//{
		//int id = it1->idx();
		//OpenMesh::Vec3d pos(originmesh_area_sqrt * position_of_mesh(id), originmesh_area_sqrt * position_of_mesh(id + V_N), 0);
		//mesh.set_point(*it1, pos);
	//}
	//OpenMesh::IO::write_mesh(mesh, "result.obj");
	//std::cout << "132" << std::endl;
}

void GAP::run_bpe()
{
	bpe_init();
	BPE();
	energy_prev_seam = energy_area;
}

void GAP::bpe_init()
{
	double area_sum = 0.0;
	for (auto f_h = mesh.faces_begin(); f_h != mesh.faces_end(); f_h++)
	{
		auto he_h = mesh.halfedge_handle(*f_h);
		area_sum += mesh.calc_sector_area(he_h);
	}
	originmesh_area_sqrt = std::sqrt(area_sum);
	double area_same_factor = 1.0 / originmesh_area_sqrt;
	VV_ids.resize(V_N);
	for (auto it1 = mesh.vertices_begin(); it1 != mesh.vertices_end(); it1++)
	{
		mesh.set_point(*it1, area_same_factor * mesh.point(*it1));
	}

	for (int i = 0; i < F_N; ++i)
	{
		auto face = mesh.face_handle(i);
		int dd = 0;
		for (Mesh::FaceVertexIter it2 = mesh.fv_begin(face); it2 != mesh.fv_end(face); ++it2)
		{
			auto vertex_ = *it2;
			switch (dd)
			{
			case 0: F0[i] = vertex_.idx(); break;
			case 1: F1[i] = vertex_.idx(); break;
			case 2: F2[i] = vertex_.idx(); break;
			default:
				break;
			}
			dd++;
		}
	}
	init_src_matirx();
	Pre_calculate();
	TutteOp();
	init_area();
}

void GAP::init_src_matirx()
{
	source_p00.resize(F_N);
	source_p01.resize(F_N);
	source_p10.resize(F_N);
	source_p11.resize(F_N);

	for (int i = 0; i < F_N; ++i)
	{
		double p00, p01, p10, p11;
		local_coordinate_inverse(i, p00, p01, p10, p11);
		source_p00[i] = p00;
		source_p01[i] = p01;
		source_p10[i] = p10;
		source_p11[i] = p11;
	}
	update_p00 = source_p00;
	update_p01 = source_p01;
	update_p10 = source_p10;
	update_p11 = source_p11;
}

void GAP::local_coordinate_inverse(int i, double& p00, double& p01, double& p10, double& p11)
{
	int f0 = F0[i];
	int f1 = F1[i];
	int f2 = F2[i];

	OpenMesh::Vec3d x_ = (mesh.point(mesh.vertex_handle(f1)) - mesh.point(mesh.vertex_handle(f0)));
	double x1_0 = x_.length();
	OpenMesh::Vec3d l_ = mesh.point(mesh.vertex_handle(f2)) - mesh.point(mesh.vertex_handle(f0));
	OpenMesh::Vec3d y_ = mesh.normal(mesh.face_handle(i)) % (1 / x1_0 * x_);
	double x2_0 = 1 / x1_0 * l_ | x_;
	double y2_0 = l_ | y_;

	p00 = 1 / x1_0;
	p01 = -x2_0 / (x1_0 * y2_0);
	p10 = 0;
	p11 = 1 / y2_0;
}

void GAP::Pre_calculate()
{//准备Pardiso的数据
	pardiso_i.clear(); pardiso_i.reserve(2 * V_N + 1);
	pardiso_ia.clear(); pardiso_ia.reserve(2 * V_N + 1);
	pardiso_ja.clear(); pardiso_ja.reserve(8 * V_N);
	typedef Eigen::Triplet<int> T;
	std::vector<T> tripletlist;
	for (int i = 0; i < 2 * V_N; i++)
	{
		pardiso_ia.push_back(pardiso_ja.size());
		if (i < V_N)
		{
			auto vertex = mesh.vertex_handle(i);
			vector<int> row_id;

			row_id.push_back(i);
			row_id.push_back(i + V_N);

			for (auto it = mesh.vv_begin(vertex); it != mesh.vv_end(vertex); ++it)
			{
				int id_neighbor = it->idx();
				row_id.push_back(id_neighbor);
				row_id.push_back(id_neighbor + V_N);
			}
			std::sort(row_id.begin(), row_id.end(), less<int>());
			vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);

			int dd = 0;

			for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
			{
				pardiso_ja.push_back(row_id[k]);
				pardiso_i.push_back(i);
				tripletlist.push_back(T(i, row_id[k], dd));
				++dd;
			}
		}
		else
		{
			auto vertex = mesh.vertex_handle(i - V_N);

			vector<int> row_id;

			row_id.push_back(i);

			for (auto it = mesh.vv_begin(vertex); it != mesh.vv_end(vertex); ++it)
			{
				int id_neighbor = it->idx() + V_N;
				row_id.push_back(id_neighbor);
			}
			std::sort(row_id.begin(), row_id.end(), less<int>());
			vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);

			int dd = 0;

			for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
			{
				pardiso_ja.push_back(row_id[k]);
				pardiso_i.push_back(i);
				tripletlist.push_back(T(i, row_id[k], dd));
				++dd;
			}
		}
	}
	Eigen::SparseMatrix<int> find_id_in_rows;
	find_id_in_rows.resize(2 * V_N, 2 * V_N);
	find_id_in_rows.setFromTriplets(tripletlist.begin(), tripletlist.end());

	pardiso_ia.push_back(pardiso_ja.size());

	id_h00.resize(F_N); id_h01.resize(F_N); id_h02.resize(F_N); id_h03.resize(F_N); id_h04.resize(F_N); id_h05.resize(F_N);
	id_h11.resize(F_N); id_h12.resize(F_N); id_h13.resize(F_N); id_h14.resize(F_N); id_h15.resize(F_N);
	id_h22.resize(F_N); id_h23.resize(F_N); id_h24.resize(F_N); id_h25.resize(F_N);
	id_h33.resize(F_N); id_h34.resize(F_N); id_h35.resize(F_N);
	id_h44.resize(F_N); id_h45.resize(F_N);
	id_h55.resize(F_N);

	for (int i = 0; i < F_N; i++)
	{
		int f0 = F0[i]; int f1 = F1[i]; int f2 = F2[i]; int f3 = F0[i] + V_N; int f4 = F1[i] + V_N; int f5 = F2[i] + V_N;

		int min01 = std::min(f0, f1); int max01 = f0 + f1 - min01;
		int min02 = std::min(f0, f2); int max02 = f0 + f2 - min02;
		int min12 = std::min(f1, f2); int max12 = f1 + f2 - min12;

		id_h00[i] = pardiso_ia[f0]; id_h01[i] = pardiso_ia[min01] + find_id_in_rows.coeff(min01, max01); id_h02[i] = pardiso_ia[min02] + find_id_in_rows.coeff(min02, max02);
		id_h03[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f3); id_h04[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f4); id_h05[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f5);

		id_h11[i] = pardiso_ia[f1]; id_h12[i] = pardiso_ia[min12] + find_id_in_rows.coeff(min12, max12);
		id_h13[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f3); id_h14[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f4); id_h15[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f5);

		id_h22[i] = pardiso_ia[f2];
		id_h23[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f3); id_h24[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f4); id_h25[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f5);

		id_h33[i] = pardiso_ia[f3]; id_h34[i] = pardiso_ia[min01 + V_N] + find_id_in_rows.coeff(min01 + V_N, max01 + V_N); id_h35[i] = pardiso_ia[min02 + V_N] + find_id_in_rows.coeff(min02 + V_N, max02 + V_N);

		id_h44[i] = pardiso_ia[f4]; id_h45[i] = pardiso_ia[min12 + V_N] + find_id_in_rows.coeff(min12 + V_N, max12 + V_N);

		id_h55[i] = pardiso_ia[f5];
	}
}

void GAP::TutteOp()
{
	int boundary_num = 0;
	auto it1 = mesh.halfedges_begin();
	while (!mesh.is_boundary(*it1))
		it1++;
	auto he_start = *it1;
	auto he_it = he_start;
	do
	{
		he_it = mesh.next_halfedge_handle(he_it);
		boundary_num++;
	} while (he_it != he_start);

	double delta_angle = 2 * M_PI / boundary_num;
	double area_1_factor = sqrt(1.0 / M_PI);

	for (int i = 0; i < boundary_num; ++i)
	{
		auto v_h = mesh.to_vertex_handle(he_start);
		position_of_mesh(v_h.idx()) = area_1_factor * cos(i * delta_angle);
		position_of_mesh(v_h.idx() + V_N) = area_1_factor * sin(-i * delta_angle);
		he_start = mesh.next_halfedge_handle(he_start);
	}

	std::vector<int> pardiso_it;
	std::vector<int> pardiso_jt;
	std::vector<double> pardiso_t;
	std::vector<double> pardiso_tu;
	std::vector<double> pardiso_tv;

	pardiso_it.reserve(V_N + 1);
	pardiso_jt.reserve(6 * V_N);
	pardiso_t.reserve(6 * V_N);
	pardiso_tu.resize(V_N, 0.0);
	pardiso_tv.resize(V_N, 0.0);
	for (size_t i = 0; i < V_N; i++)
	{
		pardiso_it.push_back(pardiso_jt.size());

		auto v_h = mesh.vertex_handle(i);
		if (mesh.is_boundary(v_h))
		{
			pardiso_jt.push_back(i);
			pardiso_t.push_back(1);

			pardiso_tu[i] = position_of_mesh(i);
			pardiso_tv[i] = position_of_mesh(i + V_N);

		}
		else
		{
			pardiso_jt.push_back(i);
			pardiso_t.push_back(mesh.valence(v_h));
			std::vector<int> row_id;
			row_id.reserve(mesh.valence(v_h));
			double bu = 0.0; double bv = 0.0;
			for (auto it2 = mesh.vv_begin(v_h); it2 != mesh.vv_end(v_h); ++it2)
			{
				int vv_id = it2->idx();
				if (mesh.is_boundary(*it2))
				{
					bu += position_of_mesh(vv_id);
					bv += position_of_mesh(vv_id + V_N);
				}
				else
				{
					if (vv_id > i)
					{
						row_id.push_back(vv_id);
					}
				}
			}
			sort(row_id.begin(), row_id.end(), std::less<int>());
			for (size_t j = 0; j < row_id.size(); j++)
			{
				pardiso_jt.push_back(row_id[j]);
				pardiso_t.push_back(-1);
			}
			pardiso_tu[i] = bu;
			pardiso_tv[i] = bv;
		}
	}
	pardiso_it.push_back(pardiso_jt.size());

	if (pardiso != NULL)
	{
		delete pardiso;
		pardiso = NULL;
	}
	pardiso = new PardisoSolver();
	pardiso->ia = pardiso_it;
	pardiso->ja = pardiso_jt;
	pardiso->nnz = pardiso_jt.size();
	pardiso->num = V_N;

	pardiso->pardiso_init();

	pardiso->a = pardiso_t;

	pardiso->rhs = pardiso_tu;
	pardiso->factorize();

	pardiso->pardiso_solver();

	for (size_t i = 0; i < V_N; i++)
	{
		position_of_mesh(i) = (pardiso->result)[i];
	}

	pardiso->rhs = pardiso_tv;
	pardiso->pardiso_solver();
	for (size_t i = 0; i < V_N; i++)
	{
		position_of_mesh(i + V_N) = (pardiso->result)[i];
	}

	delete pardiso;
	pardiso = NULL;
}

void GAP::init_area()
{
	area.resize(F_N);
	area_uniform.resize(F_N);
	area_src.resize(F_N);
	for (int i = 0; i < F_N; ++i)
	{
		auto face = mesh.face_handle(i);
		//area_src[i] = mesh.calc_sector_area(mesh.halfedge_handle(face));
		area_src[i] = 1.0 / F_N;
	}
	std::fill(area_uniform.begin(), area_uniform.end(), 1.0 / F_N);
	area = area_src;
}

void GAP::BPE()
{
	if (pardiso != NULL)
	{
		delete pardiso;
		pardiso = NULL;
	}
	pardiso = new PardisoSolver();
	pardiso->ia = pardiso_ia;
	pardiso->ja = pardiso_ja;
	pardiso->a.resize(pardiso_ja.size());
	pardiso->nnz = pardiso_ja.size();
	pardiso->num = 2 * V_N;
	pardiso->pardiso_init();
	std::vector<double> energy_area_process;
	energy_area_process.reserve(MAX_ITER_NUM);
	Energysource();

	energy_area_process.push_back(energy_area);
	double energy_pre = 0;
	double energy_cur = energy_uniform;

	int iter_num_cur = 0;
	Intp_T_Min = 0;
	changetocm_flag = 0;

	int slim_iter_num = 0;
	int sum_iter_num = 0;

	double conv_percent = 1;

	g_norm = 1.0;

	long time_beg, time_end;
	time_beg = clock();

	while (iter_num_cur < MAX_ITER_NUM)
	{
		iter_num_cur++;
		energy_pre = energy_cur;
		Update_source_same_t();
		if (changetocm_flag < 0.99 && conv_percent>0.1 && Intp_T_Min < 0.999)
		{
			SLIM();
			energy_area_process.push_back(energy_area);
			slim_iter_num++;
			sum_iter_num++;

			energy_cur = energy_uniform;
			conv_percent = abs(energy_cur - energy_pre) / energy_pre;
			//calc_gradient_norm(position_of_mesh);
			if (conv_percent <= convgence_con_rate /*|| g_norm <= convgence_con_rate*/)
			{
				break;
			}
		}
		else
		{
			break;
		}
	}
	int cm_iter_num = 0;
	while (iter_num_cur < MAX_ITER_NUM)
	{
		iter_num_cur++;
		energy_pre = energy_cur;
		Update_source_same_t();
		if (conv_percent > 0.01 && Intp_T_Min < 0.999)
		{
			CM();
			energy_area_process.push_back(energy_area);
			cm_iter_num++;
			sum_iter_num++;

			energy_cur = energy_uniform;
			conv_percent = abs(energy_cur - energy_pre) / energy_pre;
			//calc_gradient_norm(position_of_mesh);
			if (conv_percent <= convgence_con_rate /*|| g_norm <= convgence_con_rate*/)
			{
				break;
			}
		}
		else
		{
			recover_to_src();
			energy_cur = energy_area;
			while (iter_num_cur < MAX_ITER_NUM)
			{
				iter_num_cur++;
				energy_pre = energy_cur;

				CM();
				energy_area_process.push_back(energy_area);
				sum_iter_num++;

				energy_cur = energy_area;
				conv_percent = abs(energy_cur - energy_pre) / energy_pre;
				//calc_gradient_norm(position_of_mesh);
				if (conv_percent <= convgence_con_rate /*|| g_norm <= convgence_con_rate*/)
				{
					break;
				}
			}
			break;
		}
	}
	time_end = clock();
	time_consumption = (time_end - time_beg) / 1000.0;

	delete pardiso;
	pardiso = NULL;
}

void GAP::Energysource()
{
	double end_e_one_temp = 0, end_e_area = 0;
	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_1, E_2;

	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	const double* pos = position_of_mesh.data();
	for (int i = 0; i < F_N; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + V_N];

		x1 = pos[f1];
		y1 = pos[f1 + V_N];

		x2 = pos[f2];
		y2 = pos[f2 + V_N];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		det = j00 * j11 - j01 * j10;
		if (det <= 0)
			cout << "det " << det << endl;
		E_1 = (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);
		E_2 = 1.0 / (det * det) * E_1;

		end_e_one_temp += E_1;
		end_e_one_temp += E_2;
		end_e_area += ((E_1 + E_2) * area_src[i]);
	}
	energy_uniform = end_e_one_temp / F_N;

	energy_area = end_e_area;
}

void GAP::Update_source_same_t()
{
	double t_min = 1;
	int geqK = 0;

	vector<double> all_s0; all_s0.resize(F_N);
	vector<double> all_s1; all_s1.resize(F_N);

	vector<double> all_w00; all_w00.resize(F_N);
	vector<double> all_w01; all_w01.resize(F_N);
	vector<double> all_w10; all_w10.resize(F_N);
	vector<double> all_w11; all_w11.resize(F_N);

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det;
	double E_d;
	double tt;
	double new_sig0, new_sig1;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	double* position = position_of_mesh.data();

	for (int i = 0; i < F_N; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = position[f0];
		y0 = position[f0 + V_N];

		x1 = position[f1];
		y1 = position[f1 + V_N];

		x2 = position[f2];
		y2 = position[f2 + V_N];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;


		det = j00 * j11 - j01 * j10;
		E_d = (1 + 1 / (det * det)) * (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);

		double alpha_0 = j00 + j11; double alpha_1 = j10 - j01;
		double beta_0 = j00 - j11; double beta_1 = j10 + j01;

		double alpha_norm = 0.5 * sqrt(alpha_0 * alpha_0 + alpha_1 * alpha_1);
		double beta_norm = 0.5 * sqrt(beta_0 * beta_0 + beta_1 * beta_1);

		double sig0 = alpha_norm + beta_norm;
		double sig1 = alpha_norm - beta_norm;
		all_s0[i] = sig0;
		all_s1[i] = sig1;

		double temp = 1 / (sig1 * sig1 - sig0 * sig0);

		all_w00[i] = temp * (j00 * j00 + j10 * j10 - 0.5 * (sig0 * sig0 + sig1 * sig1));
		all_w01[i] = temp * (j00 * j01 + j10 * j11);
		all_w10[i] = temp * (j01 * j00 + j11 * j10);
		all_w11[i] = temp * (j01 * j01 + j11 * j11 - 0.5 * (sig0 * sig0 + sig1 * sig1));

		if (E_d <= bound_distortion_K)
		{
			geqK++;
		}
		else
		{
			tt = newton_equation(sig0, sig1, bound_distortion_K);
			if (tt < t_min)
			{
				t_min = tt;
			}
		}
	}

	changetocm_flag = (double)geqK / F_N;

	for (int i = 0; i < F_N; ++i)
	{
		double sig0 = all_s0[i];
		double sig1 = all_s1[i];

		new_sig0 = pow(sig0, t_min - 1);
		new_sig1 = pow(sig1, t_min - 1);

		double delta_new = new_sig1 - new_sig0;
		double plus_new = 0.5 * (new_sig1 + new_sig0);

		double w00 = delta_new * all_w00[i] + plus_new;
		double w01 = delta_new * all_w01[i];
		double w10 = delta_new * all_w10[i];
		double w11 = delta_new * all_w11[i] + plus_new;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		update_p00[i] = p00 * w00 + p01 * w10;
		update_p01[i] = p00 * w01 + p01 * w11;
		update_p10[i] = p10 * w00 + p11 * w10;
		update_p11[i] = p10 * w01 + p11 * w11;
	}

	Intp_T_Min = t_min;
}

double GAP::newton_equation(const double& a, const double& b, const double& K)
{
	double tt = 1;
	double E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K;
	while (abs(E_d) > 1e-5)
	{
		tt = tt - 1 / (2 * log(a) * pow(a, 2 * tt) + 2 * log(b) * pow(b, 2 * tt) + 2 * log(1 / a) * pow(1 / a, 2 * tt) + 2 * log(1 / b) * pow(1 / b, 2 * tt)) * (pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K);
		E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K;
	}
	return tt;
}

void GAP::SLIM()
{
	double area_now;
	int f0, f1, f2;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;
	double x0, y0, x1, y1, x2, y2;

	double alpha_norm, beta_norm;

	double alpha_0, alpha_1, beta_0, beta_1;

	double sig0, sig1;

	double det, tr;
	double r0, r1, r2, r3;
	double d00, d01, d02,
		d10, d11, d12;

	double new_sig0, new_sig1;
	double temp;
	double w00, w01, w10, w11;
	double p1, p2, p3, w1, w2, w3;

	double h00, h01, h02, h03, h04, h05,
		h11, h12, h13, h14, h15,
		h22, h23, h24, h25,
		h33, h34, h35,
		h44, h45,
		h55;
	double* position = position_of_mesh.data();

	int nnz = pardiso_ja.size();
	pardiso_a.clear(); pardiso_b.clear();
	pardiso_a.resize(nnz, 0.0);
	pardiso_b.resize(2 * V_N, 0.0);

	for (int i = 0; i < F_N; ++i)
	{
		area_now = area[i];
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];


		x0 = position[f0];
		y0 = position[f0 + V_N];

		x1 = position[f1];
		y1 = position[f1 + V_N];

		x2 = position[f2];
		y2 = position[f2 + V_N];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = update_p00[i]; p01 = update_p01[i]; p10 = update_p10[i]; p11 = update_p11[i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11; beta_1 = j10 + j01;

		alpha_norm = 0.5 * sqrt(alpha_0 * alpha_0 + alpha_1 * alpha_1);
		beta_norm = 0.5 * sqrt(beta_0 * beta_0 + beta_1 * beta_1);

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		new_sig0 = sqrt(1 + 1 / sig0 + 1 / (sig0 * sig0) + 1 / (sig0 * sig0 * sig0)); new_sig1 = sqrt(1 + 1 / sig1 + 1 / (sig1 * sig1) + 1 / (sig1 * sig1 * sig1));

		temp = (new_sig1 - new_sig0) / (sig1 * sig1 - sig0 * sig0);

		w00 = temp * (j00 * j00 + j01 * j01 - 0.5 * (sig0 * sig0 + sig1 * sig1)) + 0.5 * (new_sig0 + new_sig1);
		w01 = temp * (j00 * j10 + j01 * j11);
		w10 = temp * (j10 * j00 + j11 * j01);
		w11 = temp * (j10 * j10 + j11 * j11 - 0.5 * (sig0 * sig0 + sig1 * sig1)) + 0.5 * (new_sig0 + new_sig1);

		p1 = p00 * p00 + p01 * p01; p2 = p00 * p10 + p01 * p11; p3 = p10 * p10 + p11 * p11;
		w1 = w00 * w00 + w10 * w10; w2 = w00 * w01 + w10 * w11; w3 = w01 * w01 + w11 * w11;

		area_now *= 2;

		h00 = area_now * (p1 + p2 + p2 + p3) * w1; h01 = -area_now * (p1 + p2) * w1; h02 = -area_now * (p2 + p3) * w1; h03 = area_now * (p1 + p2 + p2 + p3) * w2; h04 = -area_now * (p1 + p2) * w2; h05 = -area_now * (p2 + p3) * w2;
		h11 = area_now * p1 * w1;                  h12 = area_now * p2 * w1;    	 h13 = -area_now * (p1 + p2) * w2; h14 = area_now * p1 * w2;                  h15 = area_now * p2 * w2;
		h22 = area_now * p3 * w1;                  h23 = -area_now * (p2 + p3) * w2; h24 = area_now * p2 * w2;         h25 = area_now * p3 * w2;
		h33 = area_now * (p1 + p2 + p2 + p3) * w3; h34 = -area_now * (p1 + p2) * w3; h35 = -area_now * (p2 + p3) * w3;
		h44 = area_now * p1 * w3;                  h45 = area_now * p2 * w3;
		h55 = area_now * p3 * w3;


		det = j00 * j11 - j01 * j10;
		tr = (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);

		d00 = -p00 - p10; d01 = p00; d02 = p10;
		d10 = -p01 - p11; d11 = p01; d12 = p11;

		r0 = area_now * ((1 + 1 / (det * det)) * j00 - tr * j11 / (det * det * det));
		r1 = area_now * ((1 + 1 / (det * det)) * j01 + tr * j10 / (det * det * det));
		r2 = area_now * ((1 + 1 / (det * det)) * j10 + tr * j01 / (det * det * det));
		r3 = area_now * ((1 + 1 / (det * det)) * j11 - tr * j00 / (det * det * det));


		pardiso_b[f0] -= r0 * d00 + r1 * d10;
		pardiso_b[f1] -= r0 * d01 + r1 * d11;
		pardiso_b[f2] -= r0 * d02 + r1 * d12;
		pardiso_b[f0 + V_N] -= r2 * d00 + r3 * d10;
		pardiso_b[f1 + V_N] -= r2 * d01 + r3 * d11;
		pardiso_b[f2 + V_N] -= r2 * d02 + r3 * d12;

		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;

	}

	pardiso->a = pardiso_a;
	pardiso->rhs = pardiso_b;

	pardiso->factorize();
	pardiso->pardiso_solver();

	std::vector<double> result_d = pardiso->result;

	Eigen::VectorXd negative_grad(2 * V_N), d(2 * V_N);
	for (int i = 0; i < 2 * V_N; i++)
	{
		negative_grad(i) = pardiso_b[i];
		d(i) = result_d[i];
	}

	double temp_t;
	max_step(position_of_mesh, d, temp_t);

	double alpha = min(1.0, 0.8 * temp_t);
	backtracking_line_search(position_of_mesh, d, negative_grad, alpha);
	position_of_mesh += alpha * d;

	Energysource();
}

void GAP::max_step(const Eigen::VectorXd& xx, const Eigen::VectorXd& dd, double& step)
{
	double temp_t = numeric_limits<double>::infinity();
	int f0, f1, f2;
	double a, b, c, b1, b2, tt, tt1, tt2;
	double x0, x1, x2, x3, x4, x5, d0, d1, d2, d3, d4, d5;
	const double* x = xx.data();
	const double* d = dd.data();
	for (int i = 0; i < F_N; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = x[f0]; x1 = x[f1]; x2 = x[f2]; x3 = x[f0 + V_N]; x4 = x[f1 + V_N]; x5 = x[f2 + V_N];
		d0 = d[f0]; d1 = d[f1]; d2 = d[f2]; d3 = d[f0 + V_N]; d4 = d[f1 + V_N]; d5 = d[f2 + V_N];

		a = (d1 - d0) * (d5 - d3) - (d4 - d3) * (d2 - d0);
		b1 = (d1 - d0) * (x5 - x3) + (x1 - x0) * (d5 - d3);
		b2 = (x4 - x3) * (d2 - d0) + (x2 - x0) * (d4 - d3);
		b = b1 - b2;
		c = (x1 - x0) * (x5 - x3) - (x4 - x3) * (x2 - x0);
		tt = get_smallest_pos_quad_zero(a, b, c);
		if (temp_t > tt)
		{
			temp_t = tt;
		}

	}
	step = temp_t;
}

double GAP::get_smallest_pos_quad_zero(double a, double b, double c)
{
	double t1, t2;
	if (std::abs(a) <= 1.0e-10)
	{
		a *= 1e6;
		b *= 1e6;
		c *= 1e6;
	}

	if (std::abs(a) > 1.0e-10)
	{
		double delta_in = pow(b, 2) - 4 * a * c;
		if (delta_in <= 0)
		{
			return INFINITY;
		}

		double delta = sqrt(delta_in); // delta >= 0
		if (b >= 0) // avoid subtracting two similar numbers
		{
			double bd = -b - delta;
			t1 = 2 * c / bd;
			t2 = bd / (2 * a);
		}
		else
		{
			double bd = -b + delta;
			t1 = bd / (2 * a);
			t2 = (2 * c) / bd;
		}

		assert(std::isfinite(t1));
		assert(std::isfinite(t2));

		if (a < 0) std::swap(t1, t2); // make t1 > t2
		// return the smaller positive root if it exists, otherwise return infinity
		if (t1 > 0)
		{
			return t2 > 0 ? t2 : t1;
		}
		else
		{
			return INFINITY;
		}
	}
	else
	{
		if (b == 0) return INFINITY; // just to avoid divide-by-zero
		t1 = -c / b;
		return t1 > 0 ? t1 : INFINITY;
	}
}

void GAP::backtracking_line_search(const Eigen::VectorXd& x, const Eigen::VectorXd& d, const Eigen::VectorXd& negetive_grad, double& alpha)
{
	double h = 0.5;
	double tt = -(negetive_grad.transpose() * d)(0, 0);
	double c = 0.2;
	double ex;
	Energy(x, ex);
	double e;
	Eigen::VectorXd x_new = x + alpha * d;
	Energy(x_new, e);
	while (e > ex + alpha * c * tt)
	{
		alpha = h * alpha;
		x_new = x + alpha * d;
		Energy(x_new, e);
	}
}

void GAP::Energy(const Eigen::VectorXd& position, double& energyupdate)
{
	double energy = 0;
	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_d;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;
	const double* pos = position.data();
	for (int i = 0; i < F_N; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + V_N];

		x1 = pos[f1];
		y1 = pos[f1 + V_N];

		x2 = pos[f2];
		y2 = pos[f2 + V_N];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = update_p00[i]; p01 = update_p01[i]; p10 = update_p10[i]; p11 = update_p11[i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;


		det = j00 * j11 - j01 * j10;
		E_d = (1 + 1 / (det * det)) * (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);

		energy += area[i] * E_d;
	}
	energyupdate = energy;
}

void GAP::CM()
{
	double area_now;
	int f0, f1, f2;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	double x0, y0, x1, y1, x2, y2;

	double hi_0, hi_1;

	double alpha_0, alpha_1, beta_0, beta_1;

	double s1, s2, sig0, sig1;

	double alpha_norm, beta_norm;
	double h_u, h_v, walpha, wbeta;

	double a1x0, a1x1, a1x2, a1x3, a1x4, a1x5,
		a2x0, a2x1, a2x2, a2x3, a2x4, a2x5;

	double aa, bb;
	double uu, vv, uv;
	double u, v;

	double h00, h01, h02, h03, h04, h05,
		h11, h12, h13, h14, h15,
		h22, h23, h24, h25,
		h33, h34, h35,
		h44, h45,
		h55;
	double* position = position_of_mesh.data();
	int nnz = pardiso_ja.size();
	pardiso_a.clear(); pardiso_b.clear();
	pardiso_a.resize(nnz, 0.0);
	pardiso_b.resize(2 * V_N, 0.0);

	for (int i = 0; i < F_N; ++i)
	{
		area_now = area[i];
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = position[f0];
		y0 = position[f0 + V_N];

		x1 = position[f1];
		y1 = position[f1 + V_N];

		x2 = position[f2];
		y2 = position[f2 + V_N];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = update_p00[i]; p01 = update_p01[i]; p10 = update_p10[i]; p11 = update_p11[i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11;  beta_1 = j10 + j01;

		alpha_norm = 0.5 * sqrt(alpha_0 * alpha_0 + alpha_1 * alpha_1);
		beta_norm = 0.5 * sqrt(beta_0 * beta_0 + beta_1 * beta_1);

		if (beta_norm < 1e-30)
			beta_norm = 1e-15;
		s1 = (p00) * (p00 + p10) + (p01) * (p01 + p11);
		s2 = (p10) * (p00 + p10) + (p11) * (p01 + p11);

		double h1 = p00 * p00 + p01 * p01;
		double h2 = p00 * p10 + p01 * p11;
		double h3 = p10 * p10 + p11 * p11;
		double h4 = p00 * p11 - p01 * p10;

		a1x0 = alpha_0 * (-p00 - p10) + alpha_1 * (p01 + p11);  a1x1 = alpha_0 * p00 - alpha_1 * p01; a1x2 = alpha_0 * p10 - alpha_1 * p11;
		a1x3 = alpha_0 * (-p01 - p11) + alpha_1 * (-p00 - p10); a1x4 = alpha_0 * p01 + alpha_1 * p00; a1x5 = alpha_0 * p11 + alpha_1 * p10;

		a2x0 = beta_0 * (-p00 - p10) + beta_1 * (-p01 - p11);   a2x1 = beta_0 * p00 + beta_1 * p01;   a2x2 = beta_0 * p10 + beta_1 * p11;
		a2x3 = beta_0 * (p01 + p11) + beta_1 * (-p00 - p10);    a2x4 = -beta_0 * p01 + beta_1 * p00;  a2x5 = -beta_0 * p11 + beta_1 * p10;

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		hi_0 = 2 + 6 * 1 / (sig0 * sig0 * sig0 * sig0); hi_1 = 2 + 6 * 1 / (sig1 * sig1 * sig1 * sig1);

		aa = 0.25 / alpha_norm; bb = 0.25 / beta_norm;

		uu = aa * aa * (area_now * hi_0 + area_now * hi_1);
		vv = bb * bb * (area_now * hi_0 + area_now * hi_1);
		uv = aa * bb * (area_now * hi_0 - area_now * hi_1);

		h_u = area_now * (2 * sig0 - 2 * 1 / (sig0 * sig0 * sig0));
		h_v = area_now * (2 * sig1 - 2 * 1 / (sig1 * sig1 * sig1));

		walpha = h_u + h_v;
		wbeta = h_u - h_v;

		double hwa1 = (walpha * 0.25 / alpha_norm); double hwa2 = -(walpha * 0.25 * 0.25 / (alpha_norm * alpha_norm * alpha_norm));
		double hwb1 = (wbeta * 0.25 / beta_norm); double hwb2 = -(wbeta * 0.25 * 0.25 / (beta_norm * beta_norm * beta_norm));


		h00 = uu * a1x0 * a1x0 + vv * a2x0 * a2x0 + uv * a1x0 * a2x0 + uv * a2x0 * a1x0; h01 = uu * a1x0 * a1x1 + vv * a2x0 * a2x1 + uv * a1x0 * a2x1 + uv * a2x0 * a1x1; h02 = uu * a1x0 * a1x2 + vv * a2x0 * a2x2 + uv * a1x0 * a2x2 + uv * a2x0 * a1x2; h03 = uu * a1x0 * a1x3 + vv * a2x0 * a2x3 + uv * a1x0 * a2x3 + uv * a2x0 * a1x3; h04 = uu * a1x0 * a1x4 + vv * a2x0 * a2x4 + uv * a1x0 * a2x4 + uv * a2x0 * a1x4; h05 = uu * a1x0 * a1x5 + vv * a2x0 * a2x5 + uv * a1x0 * a2x5 + uv * a2x0 * a1x5;

		h11 = uu * a1x1 * a1x1 + vv * a2x1 * a2x1 + uv * a1x1 * a2x1 + uv * a2x1 * a1x1; h12 = uu * a1x1 * a1x2 + vv * a2x1 * a2x2 + uv * a1x1 * a2x2 + uv * a2x1 * a1x2; h13 = uu * a1x1 * a1x3 + vv * a2x1 * a2x3 + uv * a1x1 * a2x3 + uv * a2x1 * a1x3; h14 = uu * a1x1 * a1x4 + vv * a2x1 * a2x4 + uv * a1x1 * a2x4 + uv * a2x1 * a1x4; h15 = uu * a1x1 * a1x5 + vv * a2x1 * a2x5 + uv * a1x1 * a2x5 + uv * a2x1 * a1x5;

		h22 = uu * a1x2 * a1x2 + vv * a2x2 * a2x2 + uv * a1x2 * a2x2 + uv * a2x2 * a1x2; h23 = uu * a1x2 * a1x3 + vv * a2x2 * a2x3 + uv * a1x2 * a2x3 + uv * a2x2 * a1x3; h24 = uu * a1x2 * a1x4 + vv * a2x2 * a2x4 + uv * a1x2 * a2x4 + uv * a2x2 * a1x4; h25 = uu * a1x2 * a1x5 + vv * a2x2 * a2x5 + uv * a1x2 * a2x5 + uv * a2x2 * a1x5;

		h33 = uu * a1x3 * a1x3 + vv * a2x3 * a2x3 + uv * a1x3 * a2x3 + uv * a2x3 * a1x3; h34 = uu * a1x3 * a1x4 + vv * a2x3 * a2x4 + uv * a1x3 * a2x4 + uv * a2x3 * a1x4; h35 = uu * a1x3 * a1x5 + vv * a2x3 * a2x5 + uv * a1x3 * a2x5 + uv * a2x3 * a1x5;

		h44 = uu * a1x4 * a1x4 + vv * a2x4 * a2x4 + uv * a1x4 * a2x4 + uv * a2x4 * a1x4; h45 = uu * a1x4 * a1x5 + vv * a2x4 * a2x5 + uv * a1x4 * a2x5 + uv * a2x4 * a1x5;

		h55 = uu * a1x5 * a1x5 + vv * a2x5 * a2x5 + uv * a1x5 * a2x5 + uv * a2x5 * a1x5;

		if (walpha >= 0)
		{
			h00 += hwa1 * (s1 + s2) + hwa2 * a1x0 * a1x0; h01 += hwa1 * (-s1) + hwa2 * a1x0 * a1x1; h02 += hwa1 * (-s2) + hwa2 * a1x0 * a1x2; h03 += hwa2 * a1x0 * a1x3; h04 += hwa1 * (h4)+hwa2 * a1x0 * a1x4; h05 += hwa1 * (-h4) + hwa2 * a1x0 * a1x5;
			h11 += hwa1 * (h1)+hwa2 * a1x1 * a1x1;        h12 += hwa1 * (h2)+hwa2 * a1x1 * a1x2;    h13 += hwa1 * (-h4) + hwa2 * a1x1 * a1x3; h14 += hwa2 * a1x1 * a1x4; h15 += hwa1 * (h4)+hwa2 * a1x1 * a1x5;
			h22 += hwa1 * (h3)+hwa2 * a1x2 * a1x2;        h23 += hwa1 * (h4)+hwa2 * a1x2 * a1x3;    h24 += hwa1 * (-h4) + hwa2 * a1x2 * a1x4; h25 += hwa2 * a1x2 * a1x5;
			h33 += hwa1 * (s1 + s2) + hwa2 * a1x3 * a1x3; h34 += hwa1 * (-s1) + hwa2 * a1x3 * a1x4; h35 += hwa1 * (-s2) + hwa2 * a1x3 * a1x5;
			h44 += hwa1 * (h1)+hwa2 * a1x4 * a1x4;        h45 += hwa1 * (h2)+hwa2 * a1x4 * a1x5;
			h55 += hwa1 * (h3)+hwa2 * a1x5 * a1x5;

		}
		h00 += hwb1 * (s1 + s2) + hwb2 * a2x0 * a2x0; h01 += hwb1 * (-s1) + hwb2 * a2x0 * a2x1; h02 += hwb1 * (-s2) + hwb2 * a2x0 * a2x2; h03 += hwb2 * a2x0 * a2x3; h04 += hwb1 * (-h4) + hwb2 * a2x0 * a2x4; h05 += hwb1 * (h4)+hwb2 * a2x0 * a2x5;
		h11 += hwb1 * (h1)+hwb2 * a2x1 * a2x1;        h12 += hwb1 * (h2)+hwb2 * a2x1 * a2x2;    h13 += hwb1 * (h4)+hwb2 * a2x1 * a2x3;    h14 += hwb2 * a2x1 * a2x4; h15 += hwb1 * (-h4) + hwb2 * a2x1 * a2x5;
		h22 += hwb1 * (h3)+hwb2 * a2x2 * a2x2;        h23 += hwb1 * (-h4) + hwb2 * a2x2 * a2x3; h24 += hwb1 * (h4)+hwb2 * a2x2 * a2x4;    h25 += hwb2 * a2x2 * a2x5;
		h33 += hwb1 * (s1 + s2) + hwb2 * a2x3 * a2x3; h34 += hwb1 * (-s1) + hwb2 * a2x3 * a2x4; h35 += hwb1 * (-s2) + hwb2 * a2x3 * a2x5;
		h44 += hwb1 * (h1)+hwb2 * a2x4 * a2x4;        h45 += hwb1 * (h2)+hwb2 * a2x4 * a2x5;
		h55 += hwb1 * (h3)+hwb2 * a2x5 * a2x5;

		u = aa * walpha; v = bb * wbeta;

		pardiso_b[f0] -= (u * a1x0 + v * a2x0);
		pardiso_b[f1] -= (u * a1x1 + v * a2x1);
		pardiso_b[f2] -= (u * a1x2 + v * a2x2);
		pardiso_b[f0 + V_N] -= (u * a1x3 + v * a2x3);
		pardiso_b[f1 + V_N] -= (u * a1x4 + v * a2x4);
		pardiso_b[f2 + V_N] -= (u * a1x5 + v * a2x5);

		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;
	}

	pardiso->a = pardiso_a;
	pardiso->rhs = pardiso_b;

	pardiso->factorize();
	pardiso->pardiso_solver();

	vector<double> result_d = pardiso->result;

	Eigen::VectorXd negative_grad(2 * V_N), d(2 * V_N);
	for (int i = 0; i < 2 * V_N; i++)
	{
		negative_grad(i) = pardiso_b[i];
		d(i) = result_d[i];
	}

	pardiso->free_numerical_factorization_memory();

	double temp_t;
	max_step(position_of_mesh, d, temp_t);

	double alpha = 0.95 * temp_t;

	backtracking_line_search(position_of_mesh, d, negative_grad, alpha);

	double e1;
	double s;
	position_of_mesh += alpha * d;
	Energysource();
}

void GAP::recover_to_src()
{
	area = area_src;
	update_p00 = source_p00;
	update_p01 = source_p01;
	update_p10 = source_p10;
	update_p11 = source_p11;
}

bool GAP::add1p(bool startfromtutte)
{
	if (CanditatePoints.size() == 0)
		return false;
	parr_count = this->parr_count;
	int step;
	if (parr_count > CanditatePoints.size())
	{
		step = CanditatePoints.size();
	}
	else
	{
		step = parr_count;
	}
	p_v_seam.clear();
	p_e_seam.clear();
	p_idx2meshvid.clear();
	p_he2idx.clear();
	p_F0.clear();
	p_F1.clear();
	p_F2.clear();
	p_VV_ids.clear();
	p_position_of_mesh.clear();
	p_V_N.clear();
	p_F_N.clear();
	p_pardiso_i.clear();
	p_pardiso_ia.clear();
	p_pardiso_ja.clear();
	p_pardiso_a.clear();
	p_pardiso_b.clear();
	p_id_h00.clear(); p_id_h01.clear(); p_id_h02.clear(); p_id_h03.clear(); p_id_h04.clear(); p_id_h05.clear();
	p_id_h11.clear(); p_id_h12.clear(); p_id_h13.clear(); p_id_h14.clear(); p_id_h15.clear();
	p_id_h22.clear(); p_id_h23.clear(); p_id_h24.clear(); p_id_h25.clear();
	p_id_h33.clear(); p_id_h34.clear(); p_id_h35.clear();
	p_id_h44.clear(); p_id_h45.clear();
	p_id_h55.clear();

	p_source_p00.clear(); p_source_p01.clear(); p_source_p10.clear(); p_source_p11.clear();
	p_update_p00.clear(); p_update_p01.clear(); p_update_p10.clear(); p_update_p11.clear();
	p_Intp_T_Min.clear();
	p_changetocm_flag.clear();
	p_g_norm.clear();
	p_energy_uniform.clear();
	p_energy_area.clear();
	p_energy_prev_seam.clear();
	p_pardiso.clear();

	p_v_seam.resize(step, v_seam);
	p_e_seam.resize(step, e_seam);
	p_idx2meshvid.resize(step, idx2meshvid);
	p_he2idx.resize(step, he2idx);
	p_F0.resize(step, F0);
	p_F1.resize(step, F1);
	p_F2.resize(step, F2);
	p_VV_ids.resize(step, VV_ids);
	p_position_of_mesh.resize(step, position_of_mesh);
	p_V_N.resize(step, V_N);
	p_F_N.resize(step, F_N);
	p_pardiso_i.resize(step, pardiso_i);
	p_pardiso_ia.resize(step, pardiso_ia);
	p_pardiso_ja.resize(step, pardiso_ja);
	p_pardiso_a.resize(step, pardiso_a);
	p_pardiso_b.resize(step, pardiso_b);

	p_id_h00.resize(step, id_h00); p_id_h01.resize(step, id_h01); p_id_h02.resize(step, id_h02); p_id_h03.resize(step, id_h03); p_id_h04.resize(step, id_h04); p_id_h05.resize(step, id_h05);
	p_id_h11.resize(step, id_h11); p_id_h12.resize(step, id_h12); p_id_h13.resize(step, id_h13); p_id_h14.resize(step, id_h14); p_id_h15.resize(step, id_h15);
	p_id_h22.resize(step, id_h22); p_id_h23.resize(step, id_h23); p_id_h24.resize(step, id_h24); p_id_h25.resize(step, id_h25);
	p_id_h33.resize(step, id_h33); p_id_h34.resize(step, id_h34); p_id_h35.resize(step, id_h35);
	p_id_h44.resize(step, id_h44); p_id_h45.resize(step, id_h45);
	p_id_h55.resize(step, id_h55);

	p_source_p00.resize(step, source_p00); p_source_p01.resize(step, source_p01); p_source_p10.resize(step, source_p10); p_source_p11.resize(step, source_p11);
	p_update_p00.resize(step, update_p00); p_update_p01.resize(step, update_p01); p_update_p10.resize(step, update_p10); p_update_p11.resize(step, update_p11);
	p_Intp_T_Min.resize(step, Intp_T_Min);
	p_changetocm_flag.resize(step, changetocm_flag);
	p_g_norm.resize(step, g_norm);
	p_energy_uniform.resize(step, energy_uniform);
	p_energy_area.resize(step, energy_area);
	p_energy_prev_seam.resize(step, energy_area);
	p_pardiso.resize(step, pardiso);

	std::vector<std::vector<int>> p_path_1p2seam_v(step);
	std::vector<std::vector<int>> p_path_1p2seam_e(step);

	std::vector<double> des_value(step, -1);

#pragma omp parallel for num_threads(step)
	for (int ii = 0; ii < step; ii++)
	{
		if (p_v_seam[ii][CanditatePoints[ii]] == 1)
		{
			des_value[ii] = 0;
			continue;
		}
		p_path_1p2seam_v[ii].clear();
		p_path_1p2seam_e[ii].clear();

		Algorithm::Dijkstra_with_nearest2(MCache, CanditatePoints[ii], p_v_seam[ii], p_path_1p2seam_v[ii]);
		int tmp_V_N = p_V_N[ii] + p_path_1p2seam_v[ii].size() - 1;
		Eigen::VectorXd pos_bk_tmp;
		pos_bk_tmp.resize(2 * tmp_V_N);
		pos_bk_tmp.block(0, 0, p_V_N[ii], 1) = p_position_of_mesh[ii].topRows(p_V_N[ii]);
		pos_bk_tmp.block(tmp_V_N, 0, p_V_N[ii], 1) = p_position_of_mesh[ii].bottomRows(p_V_N[ii]);
		position_of_mesh;
		p_idx2meshvid[ii].resize(tmp_V_N);
		int v_n_tmp = p_V_N[ii];
		OpenMesh::HalfedgeHandle h_iter;
		for (const auto& heh : Closedmesh.voh_range(Closedmesh.vertex_handle(p_path_1p2seam_v[ii].back())))
		{
			if (Closedmesh.to_vertex_handle(heh).idx() == p_path_1p2seam_v[ii][p_path_1p2seam_v[ii].size() - 2])
			{
				h_iter = heh;
				break;
			}
		}
		int i_ = p_path_1p2seam_v[ii].size() - 2;
		int to_vid;

		do
		{
			p_path_1p2seam_e[ii].push_back(h_iter.idx() >> 1);
			auto fh = Closedmesh.face_handle(h_iter);
			pos_bk_tmp(v_n_tmp) = p_position_of_mesh[ii](p_he2idx[ii][h_iter.idx()]);
			pos_bk_tmp(v_n_tmp + tmp_V_N) = p_position_of_mesh[ii](p_he2idx[ii][h_iter.idx()] + p_V_N[ii]);
			p_he2idx[ii][h_iter.idx()] = v_n_tmp;
			p_idx2meshvid[ii][v_n_tmp] = Closedmesh.to_vertex_handle(h_iter).idx();
			auto itfh = Closedmesh.fh_begin(fh);
			p_F0[ii][fh.idx()] = p_he2idx[ii][itfh->idx()];
			itfh++;
			p_F1[ii][fh.idx()] = p_he2idx[ii][itfh->idx()];
			itfh++;
			p_F2[ii][fh.idx()] = p_he2idx[ii][itfh->idx()];
			h_iter = Closedmesh.next_halfedge_handle(h_iter);
			to_vid = Closedmesh.to_vertex_handle(h_iter).idx();

			if (i_ == 0)
			{
				while (!p_e_seam[ii][Closedmesh.edge_handle(h_iter).idx()])
				{
					h_iter = Closedmesh.opposite_halfedge_handle(h_iter);
					p_he2idx[ii][h_iter.idx()] = v_n_tmp;

					fh = Closedmesh.face_handle(h_iter);
					itfh = Closedmesh.fh_begin(fh);
					p_F0[ii][fh.idx()] = p_he2idx[ii][itfh->idx()];
					itfh++;
					p_F1[ii][fh.idx()] = p_he2idx[ii][itfh->idx()];
					itfh++;
					p_F2[ii][fh.idx()] = p_he2idx[ii][itfh->idx()];
					h_iter = Closedmesh.next_halfedge_handle(h_iter);
					to_vid = Closedmesh.to_vertex_handle(h_iter).idx();
				}
				break;
			}

			while (to_vid != p_path_1p2seam_v[ii][i_ - 1])
			{
				h_iter = Closedmesh.opposite_halfedge_handle(h_iter);
				p_he2idx[ii][h_iter.idx()] = v_n_tmp;
				fh = Closedmesh.face_handle(h_iter);
				itfh = Closedmesh.fh_begin(fh);
				p_F0[ii][fh.idx()] = p_he2idx[ii][itfh->idx()];
				itfh++;
				p_F1[ii][fh.idx()] = p_he2idx[ii][itfh->idx()];
				itfh++;
				p_F2[ii][fh.idx()] = p_he2idx[ii][itfh->idx()];

				h_iter = Closedmesh.next_halfedge_handle(h_iter);
				to_vid = Closedmesh.to_vertex_handle(h_iter).idx();
			}

			v_n_tmp++;
			i_--;
		} while (true);
		p_position_of_mesh[ii] = pos_bk_tmp;
		p_V_N[ii] = tmp_V_N;

		std::vector<int> e_seam_tmp;
		e_seam_tmp = p_e_seam[ii];
		for (auto& var : p_path_1p2seam_e[ii])
		{
			e_seam_tmp[var] = 1;
		}

		{
			Mesh::HalfedgeHandle h_begin;
			for (int i = 0; i < e_seam_tmp.size(); i++)
			{
				if (e_seam_tmp[i])
				{
					h_begin = Closedmesh.halfedge_handle(Closedmesh.edge_handle(i), 0);
					break;
				}
			}
			auto h_iter = h_begin;
			std::vector<int> boundary_vs;
			std::vector<int> hiter_id;
			do
			{
				hiter_id.push_back(h_iter.idx());
				boundary_vs.push_back(p_he2idx[ii][h_iter.idx()]);
				h_iter = Closedmesh.next_halfedge_handle(h_iter);
				int tmp = Closedmesh.edge_handle(h_iter).idx();
				while (!e_seam_tmp[Closedmesh.edge_handle(h_iter).idx()])
				{
					h_iter = Closedmesh.opposite_halfedge_handle(h_iter);
					h_iter = Closedmesh.next_halfedge_handle(h_iter);
				}
			} while (h_iter != h_begin);

			p_VV_ids[ii].clear();
			p_VV_ids[ii].resize(tmp_V_N);
			for (size_t i = 0; i < F_N; i++)
			{
				p_VV_ids[ii][p_F0[ii][i]].push_back(p_F2[ii][i]);
				p_VV_ids[ii][p_F1[ii][i]].push_back(p_F0[ii][i]);
				p_VV_ids[ii][p_F2[ii][i]].push_back(p_F1[ii][i]);
			}
			for (size_t i = 0; i < boundary_vs.size() - 1; i++)
			{
				p_VV_ids[ii][boundary_vs[i]].push_back(boundary_vs[i + 1]);
			}
			p_VV_ids[ii][boundary_vs.back()].push_back(boundary_vs.front());
		}
		rePre_calculate(ii);
		BPE(ii);
		double energy_decrease = p_energy_prev_seam[ii] - p_energy_area[ii];
		des_value[ii] = energy_decrease;
	}

	int max_id = -1;
	double max_dec = 0;
	std::vector<int> value;
	for (int a = 0; a < step; a++)
		value.push_back(CanditatePoints[a]);
	for (int m = 0; m < step; m++)
	{
		if (des_value[m] < filtering_rate)
		{
			auto it = std::find(CanditatePoints.begin(), CanditatePoints.end(), value[m]);
			CanditatePoints.erase(it);
		}
		else
		{
			if (des_value[m] > max_dec)
			{
				max_id = m;
				max_dec = DBL_MAX;
			}
		}
	}
	if (max_id != -1)
	{
		std::cout << "addpoint   " << value[max_id] << std::endl;
		result.push_back(value[max_id]);
		auto it = std::find(CanditatePoints.begin(), CanditatePoints.end(), value[max_id]);
		CanditatePoints.erase(it);

		v_seam = p_v_seam[max_id];
		e_seam = p_e_seam[max_id];
		for (size_t i = 0; i < p_path_1p2seam_v[max_id].size(); i++)
		{
			v_seam[p_path_1p2seam_v[max_id][i]] = 1;
		}
		for (size_t i = 0; i < p_path_1p2seam_e[max_id].size(); i++)
		{
			e_seam[p_path_1p2seam_e[max_id][i]] = 1;
		}

		idx2meshvid = p_idx2meshvid[max_id];
		he2idx = p_he2idx[max_id];
		F0 = p_F0[max_id];
		F1 = p_F1[max_id];
		F2 = p_F2[max_id];
		VV_ids = p_VV_ids[max_id];
		position_of_mesh = p_position_of_mesh[max_id];
		V_N = p_V_N[max_id];
		F_N = p_F_N[max_id];

		energy_uniform = p_energy_uniform[max_id];
		energy_area = p_energy_area[max_id];
		energy_prev_seam = p_energy_prev_seam[max_id];

	}
	return true;
}

void GAP::rePre_calculate(int N)
{
	p_pardiso_i[N].clear(); p_pardiso_i[N].reserve(2 * p_V_N[N] + 1);
	p_pardiso_ia[N].clear(); p_pardiso_ia[N].reserve(2 * p_V_N[N] + 1);
	p_pardiso_ja[N].clear(); p_pardiso_ja[N].reserve(8 * p_V_N[N]);

	typedef Eigen::Triplet<int> T;
	std::vector<T> tripletlist;
	for (int i = 0; i < 2 * p_V_N[N]; i++)
	{
		p_pardiso_ia[N].push_back(p_pardiso_ja[N].size());
		if (i < p_V_N[N])
		{
			std::vector<int> row_id;
			row_id.push_back(i);
			row_id.push_back(i + p_V_N[N]);
			for (auto& vv : p_VV_ids[N][i])
			{
				row_id.push_back(vv);
				row_id.push_back(vv + p_V_N[N]);
			}
			std::sort(row_id.begin(), row_id.end(), less<int>());
			std::vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);
			int dd = 0;
			for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
			{
				p_pardiso_ja[N].push_back(row_id[k]);
				p_pardiso_i[N].push_back(i);
				tripletlist.push_back(T(i, row_id[k], dd));
				++dd;
			}
		}
		else
		{
			std::vector<int> row_id;
			row_id.push_back(i);
			for (auto& vv : p_VV_ids[N][i - p_V_N[N]])
			{
				row_id.push_back(vv + p_V_N[N]);
			}
			std::sort(row_id.begin(), row_id.end(), less<int>());
			std::vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);
			int dd = 0;
			for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
			{
				p_pardiso_ja[N].push_back(row_id[k]);
				p_pardiso_i[N].push_back(i);
				tripletlist.push_back(T(i, row_id[k], dd));
				++dd;
			}
		}
	}
	Eigen::SparseMatrix<int> find_id_in_rows;
	find_id_in_rows.resize(2 * p_V_N[N], 2 * p_V_N[N]);
	find_id_in_rows.setFromTriplets(tripletlist.begin(), tripletlist.end());

	p_pardiso_ia[N].push_back(p_pardiso_ja[N].size());

	p_id_h00[N].resize(p_F_N[N]); p_id_h01[N].resize(p_F_N[N]); p_id_h02[N].resize(p_F_N[N]); p_id_h03[N].resize(p_F_N[N]); p_id_h04[N].resize(p_F_N[N]); p_id_h05[N].resize(p_F_N[N]);
	p_id_h11[N].resize(p_F_N[N]); p_id_h12[N].resize(p_F_N[N]); p_id_h13[N].resize(p_F_N[N]); p_id_h14[N].resize(p_F_N[N]); p_id_h15[N].resize(p_F_N[N]);
	p_id_h22[N].resize(p_F_N[N]); p_id_h23[N].resize(p_F_N[N]); p_id_h24[N].resize(p_F_N[N]); p_id_h25[N].resize(p_F_N[N]);
	p_id_h33[N].resize(p_F_N[N]); p_id_h34[N].resize(p_F_N[N]); p_id_h35[N].resize(p_F_N[N]);
	p_id_h44[N].resize(p_F_N[N]); p_id_h45[N].resize(p_F_N[N]);
	p_id_h55[N].resize(p_F_N[N]);

	for (int i = 0; i < p_F_N[N]; i++)
	{
		int f0 = p_F0[N][i]; int f1 = p_F1[N][i]; int f2 = p_F2[N][i]; int f3 = p_F0[N][i] + p_V_N[N]; int f4 = p_F1[N][i] + p_V_N[N]; int f5 = p_F2[N][i] + p_V_N[N];

		int min01 = min(f0, f1); int max01 = f0 + f1 - min01;
		int min02 = min(f0, f2); int max02 = f0 + f2 - min02;
		int min12 = min(f1, f2); int max12 = f1 + f2 - min12;

		p_id_h00[N][i] = p_pardiso_ia[N][f0]; p_id_h01[N][i] = p_pardiso_ia[N][min01] + find_id_in_rows.coeff(min01, max01); p_id_h02[N][i] = p_pardiso_ia[N][min02] + find_id_in_rows.coeff(min02, max02);
		p_id_h03[N][i] = p_pardiso_ia[N][f0] + find_id_in_rows.coeff(f0, f3); p_id_h04[N][i] = p_pardiso_ia[N][f0] + find_id_in_rows.coeff(f0, f4); p_id_h05[N][i] = p_pardiso_ia[N][f0] + find_id_in_rows.coeff(f0, f5);

		p_id_h11[N][i] = p_pardiso_ia[N][f1]; p_id_h12[N][i] = p_pardiso_ia[N][min12] + find_id_in_rows.coeff(min12, max12);
		p_id_h13[N][i] = p_pardiso_ia[N][f1] + find_id_in_rows.coeff(f1, f3); p_id_h14[N][i] = p_pardiso_ia[N][f1] + find_id_in_rows.coeff(f1, f4); p_id_h15[N][i] = p_pardiso_ia[N][f1] + find_id_in_rows.coeff(f1, f5);

		p_id_h22[N][i] = p_pardiso_ia[N][f2];
		p_id_h23[N][i] = p_pardiso_ia[N][f2] + find_id_in_rows.coeff(f2, f3); p_id_h24[N][i] = p_pardiso_ia[N][f2] + find_id_in_rows.coeff(f2, f4); p_id_h25[N][i] = p_pardiso_ia[N][f2] + find_id_in_rows.coeff(f2, f5);

		p_id_h33[N][i] = p_pardiso_ia[N][f3]; p_id_h34[N][i] = p_pardiso_ia[N][min01 + p_V_N[N]] + find_id_in_rows.coeff(min01 + p_V_N[N], max01 + p_V_N[N]); p_id_h35[N][i] = p_pardiso_ia[N][min02 + p_V_N[N]] + find_id_in_rows.coeff(min02 + p_V_N[N], max02 + p_V_N[N]);

		p_id_h44[N][i] = p_pardiso_ia[N][f4]; p_id_h45[N][i] = p_pardiso_ia[N][min12 + p_V_N[N]] + find_id_in_rows.coeff(min12 + p_V_N[N], max12 + p_V_N[N]);

		p_id_h55[N][i] = p_pardiso_ia[N][f5];
	}
}

void GAP::BPE(int N)
{
	if (p_pardiso[N] != NULL)
	{
		delete p_pardiso[N];
		p_pardiso[N] = NULL;
	}
	p_pardiso[N] = new PardisoSolver();
	p_pardiso[N]->ia = p_pardiso_ia[N];
	p_pardiso[N]->ja = p_pardiso_ja[N];
	p_pardiso[N]->a.resize(p_pardiso_ja[N].size());
	p_pardiso[N]->nnz = p_pardiso_ja[N].size();
	p_pardiso[N]->num = 2 * p_V_N[N];

	p_pardiso[N]->pardiso_init();

	std::vector<double> energy_area_process;
	energy_area_process.reserve(MAX_ITER_NUM);
	Energysource(N);

	energy_area_process.push_back(p_energy_area[N]);
	double energy_pre = 0;
	double energy_cur = p_energy_uniform[N];

	int iter_num_cur = 0;
	p_Intp_T_Min[N] = 0;
	p_changetocm_flag[N] = 0;

	int slim_iter_num = 0;
	int sum_iter_num = 0;

	double conv_percent = 1;

	p_g_norm[N] = 1.0;

	long time_beg, time_end;
	time_beg = clock();

	while (iter_num_cur < MAX_ITER_NUM)
	{
		iter_num_cur++;
		energy_pre = energy_cur;
		Update_source_same_t(N);
		if (p_changetocm_flag[N] < 0.99 && conv_percent>0.1 && p_Intp_T_Min[N] < 0.999)
		{
			SLIM(N);
			energy_area_process.push_back(p_energy_area[N]);
			slim_iter_num++;
			sum_iter_num++;

			energy_cur = energy_uniform;
			conv_percent = abs(energy_cur - energy_pre) / energy_pre;
			//calc_gradient_norm(position_of_mesh);
			if (conv_percent <= convgence_con_rate /*|| g_norm <= convgence_con_rate*/)
			{
				break;
			}
		}
		else
		{
			break;
		}
	}

	int cm_iter_num = 0;
	while (iter_num_cur < MAX_ITER_NUM)
	{
		iter_num_cur++;
		energy_pre = energy_cur;
		Update_source_same_t(N);
		if (conv_percent > 0.01 && p_Intp_T_Min[N] < 0.999)
		{
			CM(N);
			energy_area_process.push_back(p_energy_area[N]);
			cm_iter_num++;
			sum_iter_num++;

			energy_cur = energy_uniform;
			conv_percent = abs(energy_cur - energy_pre) / energy_pre;

			if (conv_percent <= convgence_con_rate /*|| g_norm <= convgence_con_rate*/)
			{
				break;
			}
		}
		else
		{
			recover_to_src(N);
			energy_cur = p_energy_area[N];
			while (iter_num_cur < MAX_ITER_NUM)
			{
				iter_num_cur++;
				energy_pre = energy_cur;

				CM(N);
				energy_area_process.push_back(p_energy_area[N]);
				sum_iter_num++;

				energy_cur = p_energy_area[N];
				conv_percent = abs(energy_cur - energy_pre) / energy_pre;
				//calc_gradient_norm(position_of_mesh);
				if (conv_percent <= convgence_con_rate /*|| g_norm <= convgence_con_rate*/)
				{
					break;
				}
			}
			break;
		}

	}

	time_end = clock();
	time_consumption = (time_end - time_beg) / 1000.0;

	delete p_pardiso[N];
	p_pardiso[N] = NULL;
}

void GAP::Energysource(int N)
{
	double end_e_one_temp = 0, end_e_area = 0;

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_1, E_2;

	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	const double* pos = p_position_of_mesh[N].data();
	for (int i = 0; i < p_F_N[N]; ++i)
	{
		f0 = p_F0[N][i];
		f1 = p_F1[N][i];
		f2 = p_F2[N][i];

		x0 = pos[f0];
		y0 = pos[f0 + p_V_N[N]];

		x1 = pos[f1];
		y1 = pos[f1 + p_V_N[N]];

		x2 = pos[f2];
		y2 = pos[f2 + p_V_N[N]];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = p_source_p00[N][i]; p01 = p_source_p01[N][i]; p10 = p_source_p10[N][i]; p11 = p_source_p11[N][i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		det = j00 * j11 - j01 * j10;
		if (det <= 0)
			cout << "det " << det << endl;
		E_1 = (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);
		E_2 = 1.0 / (det * det) * E_1;

		end_e_one_temp += E_1;
		end_e_one_temp += E_2;
		end_e_area += ((E_1 + E_2) * area_src[i]);
	}
	p_energy_uniform[N] = end_e_one_temp / p_F_N[N];

	p_energy_area[N] = end_e_area;
}

void GAP::Update_source_same_t(int N)
{
	double t_min = 1;
	int geqK = 0;

	vector<double> all_s0; all_s0.resize(p_F_N[N]);
	vector<double> all_s1; all_s1.resize(p_F_N[N]);

	vector<double> all_w00; all_w00.resize(p_F_N[N]);
	vector<double> all_w01; all_w01.resize(p_F_N[N]);
	vector<double> all_w10; all_w10.resize(p_F_N[N]);
	vector<double> all_w11; all_w11.resize(p_F_N[N]);


	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det;
	double E_d;
	double tt;
	double new_sig0, new_sig1;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	double* position = p_position_of_mesh[N].data();

	for (int i = 0; i < p_F_N[N]; ++i)
	{
		f0 = p_F0[N][i];
		f1 = p_F1[N][i];
		f2 = p_F2[N][i];

		x0 = position[f0];
		y0 = position[f0 + p_V_N[N]];

		x1 = position[f1];
		y1 = position[f1 + p_V_N[N]];

		x2 = position[f2];
		y2 = position[f2 + p_V_N[N]];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = p_source_p00[N][i]; p01 = p_source_p01[N][i]; p10 = p_source_p10[N][i]; p11 = p_source_p11[N][i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;


		det = j00 * j11 - j01 * j10;
		E_d = (1 + 1 / (det * det)) * (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);

		double alpha_0 = j00 + j11; double alpha_1 = j10 - j01;
		double beta_0 = j00 - j11; double beta_1 = j10 + j01;

		double alpha_norm = 0.5 * sqrt(alpha_0 * alpha_0 + alpha_1 * alpha_1);
		double beta_norm = 0.5 * sqrt(beta_0 * beta_0 + beta_1 * beta_1);

		double sig0 = alpha_norm + beta_norm;
		double sig1 = alpha_norm - beta_norm;
		all_s0[i] = sig0;
		all_s1[i] = sig1;

		double temp = 1 / (sig1 * sig1 - sig0 * sig0);

		all_w00[i] = temp * (j00 * j00 + j10 * j10 - 0.5 * (sig0 * sig0 + sig1 * sig1));
		all_w01[i] = temp * (j00 * j01 + j10 * j11);
		all_w10[i] = temp * (j01 * j00 + j11 * j10);
		all_w11[i] = temp * (j01 * j01 + j11 * j11 - 0.5 * (sig0 * sig0 + sig1 * sig1));

		if (E_d <= bound_distortion_K)
		{
			geqK++;
		}
		else
		{
			tt = newton_equation(sig0, sig1, bound_distortion_K);
			if (tt < t_min)
			{
				t_min = tt;
			}
		}
	}

	p_changetocm_flag[N] = (double)geqK / p_F_N[N];

	for (int i = 0; i < p_F_N[N]; ++i)
	{
		double sig0 = all_s0[i];
		double sig1 = all_s1[i];

		new_sig0 = pow(sig0, t_min - 1);
		new_sig1 = pow(sig1, t_min - 1);

		double delta_new = new_sig1 - new_sig0;
		double plus_new = 0.5 * (new_sig1 + new_sig0);

		double w00 = delta_new * all_w00[i] + plus_new;
		double w01 = delta_new * all_w01[i];
		double w10 = delta_new * all_w10[i];
		double w11 = delta_new * all_w11[i] + plus_new;

		p00 = p_source_p00[N][i]; p01 = p_source_p01[N][i]; p10 = p_source_p10[N][i]; p11 = p_source_p11[N][i];

		p_update_p00[N][i] = p00 * w00 + p01 * w10;
		p_update_p01[N][i] = p00 * w01 + p01 * w11;
		p_update_p10[N][i] = p10 * w00 + p11 * w10;
		p_update_p11[N][i] = p10 * w01 + p11 * w11;
	}

	p_Intp_T_Min[N] = t_min;
}

void GAP::SLIM(int N)
{
	double area_now;
	int f0, f1, f2;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	double x0, y0, x1, y1, x2, y2;

	double alpha_norm, beta_norm;

	double alpha_0, alpha_1, beta_0, beta_1;

	double sig0, sig1;

	double det, tr;
	double r0, r1, r2, r3;
	double d00, d01, d02,
		d10, d11, d12;

	double new_sig0, new_sig1;
	double temp;
	double w00, w01, w10, w11;
	double p1, p2, p3, w1, w2, w3;

	double h00, h01, h02, h03, h04, h05,
		h11, h12, h13, h14, h15,
		h22, h23, h24, h25,
		h33, h34, h35,
		h44, h45,
		h55;
	double* position = p_position_of_mesh[N].data();

	int nnz = p_pardiso_ja[N].size();
	p_pardiso_a[N].clear(); p_pardiso_b[N].clear();
	p_pardiso_a[N].resize(nnz, 0.0);
	p_pardiso_b[N].resize(2 * p_V_N[N], 0.0);

	for (int i = 0; i < p_F_N[N]; ++i)
	{
		area_now = area[i];
		f0 = p_F0[N][i];
		f1 = p_F1[N][i];
		f2 = p_F2[N][i];

		x0 = position[f0];
		y0 = position[f0 + p_V_N[N]];

		x1 = position[f1];
		y1 = position[f1 + p_V_N[N]];

		x2 = position[f2];
		y2 = position[f2 + p_V_N[N]];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = p_update_p00[N][i]; p01 = p_update_p01[N][i]; p10 = p_update_p10[N][i]; p11 = p_update_p11[N][i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11; beta_1 = j10 + j01;

		alpha_norm = 0.5 * sqrt(alpha_0 * alpha_0 + alpha_1 * alpha_1);
		beta_norm = 0.5 * sqrt(beta_0 * beta_0 + beta_1 * beta_1);

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		new_sig0 = sqrt(1 + 1 / sig0 + 1 / (sig0 * sig0) + 1 / (sig0 * sig0 * sig0)); new_sig1 = sqrt(1 + 1 / sig1 + 1 / (sig1 * sig1) + 1 / (sig1 * sig1 * sig1));

		temp = (new_sig1 - new_sig0) / (sig1 * sig1 - sig0 * sig0);

		w00 = temp * (j00 * j00 + j01 * j01 - 0.5 * (sig0 * sig0 + sig1 * sig1)) + 0.5 * (new_sig0 + new_sig1);
		w01 = temp * (j00 * j10 + j01 * j11);
		w10 = temp * (j10 * j00 + j11 * j01);
		w11 = temp * (j10 * j10 + j11 * j11 - 0.5 * (sig0 * sig0 + sig1 * sig1)) + 0.5 * (new_sig0 + new_sig1);

		p1 = p00 * p00 + p01 * p01; p2 = p00 * p10 + p01 * p11; p3 = p10 * p10 + p11 * p11;
		w1 = w00 * w00 + w10 * w10; w2 = w00 * w01 + w10 * w11; w3 = w01 * w01 + w11 * w11;

		area_now *= 2;

		h00 = area_now * (p1 + p2 + p2 + p3) * w1; h01 = -area_now * (p1 + p2) * w1; h02 = -area_now * (p2 + p3) * w1; h03 = area_now * (p1 + p2 + p2 + p3) * w2; h04 = -area_now * (p1 + p2) * w2; h05 = -area_now * (p2 + p3) * w2;
		h11 = area_now * p1 * w1;                  h12 = area_now * p2 * w1;    	 h13 = -area_now * (p1 + p2) * w2; h14 = area_now * p1 * w2;                  h15 = area_now * p2 * w2;
		h22 = area_now * p3 * w1;                  h23 = -area_now * (p2 + p3) * w2; h24 = area_now * p2 * w2;         h25 = area_now * p3 * w2;
		h33 = area_now * (p1 + p2 + p2 + p3) * w3; h34 = -area_now * (p1 + p2) * w3; h35 = -area_now * (p2 + p3) * w3;
		h44 = area_now * p1 * w3;                  h45 = area_now * p2 * w3;
		h55 = area_now * p3 * w3;


		det = j00 * j11 - j01 * j10;
		tr = (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);

		d00 = -p00 - p10; d01 = p00; d02 = p10;
		d10 = -p01 - p11; d11 = p01; d12 = p11;

		r0 = area_now * ((1 + 1 / (det * det)) * j00 - tr * j11 / (det * det * det));
		r1 = area_now * ((1 + 1 / (det * det)) * j01 + tr * j10 / (det * det * det));
		r2 = area_now * ((1 + 1 / (det * det)) * j10 + tr * j01 / (det * det * det));
		r3 = area_now * ((1 + 1 / (det * det)) * j11 - tr * j00 / (det * det * det));


		p_pardiso_b[N][f0] -= r0 * d00 + r1 * d10;
		p_pardiso_b[N][f1] -= r0 * d01 + r1 * d11;
		p_pardiso_b[N][f2] -= r0 * d02 + r1 * d12;
		p_pardiso_b[N][f0 + p_V_N[N]] -= r2 * d00 + r3 * d10;
		p_pardiso_b[N][f1 + p_V_N[N]] -= r2 * d01 + r3 * d11;
		p_pardiso_b[N][f2 + p_V_N[N]] -= r2 * d02 + r3 * d12;

		p_pardiso_a[N][p_id_h00[N][i]] += h00; p_pardiso_a[N][p_id_h01[N][i]] += h01; p_pardiso_a[N][p_id_h02[N][i]] += h02; p_pardiso_a[N][p_id_h03[N][i]] += h03; p_pardiso_a[N][p_id_h04[N][i]] += h04; p_pardiso_a[N][p_id_h05[N][i]] += h05;
		p_pardiso_a[N][p_id_h11[N][i]] += h11; p_pardiso_a[N][p_id_h12[N][i]] += h12; p_pardiso_a[N][p_id_h13[N][i]] += h13; p_pardiso_a[N][p_id_h14[N][i]] += h14; p_pardiso_a[N][p_id_h15[N][i]] += h15;
		p_pardiso_a[N][p_id_h22[N][i]] += h22; p_pardiso_a[N][p_id_h23[N][i]] += h23; p_pardiso_a[N][p_id_h24[N][i]] += h24; p_pardiso_a[N][p_id_h25[N][i]] += h25;
		p_pardiso_a[N][p_id_h33[N][i]] += h33; p_pardiso_a[N][p_id_h34[N][i]] += h34; p_pardiso_a[N][p_id_h35[N][i]] += h35;
		p_pardiso_a[N][p_id_h44[N][i]] += h44; p_pardiso_a[N][p_id_h45[N][i]] += h45;
		p_pardiso_a[N][p_id_h55[N][i]] += h55;

	}

	p_pardiso[N]->a = p_pardiso_a[N];
	p_pardiso[N]->rhs = p_pardiso_b[N];

	p_pardiso[N]->factorize();
	p_pardiso[N]->pardiso_solver();

	std::vector<double> result_d = p_pardiso[N]->result;

	Eigen::VectorXd negative_grad(2 * p_V_N[N]), d(2 * p_V_N[N]);
	for (int i = 0; i < 2 * p_V_N[N]; i++)
	{
		negative_grad(i) = p_pardiso_b[N][i];
		d(i) = result_d[i];
	}

	double temp_t;
	max_step(p_position_of_mesh[N], d, temp_t, N);

	double alpha = min(1.0, 0.8 * temp_t);
	backtracking_line_search(p_position_of_mesh[N], d, negative_grad, alpha, N);
	p_position_of_mesh[N] += alpha * d;

	Energysource(N);
}

void GAP::max_step(const Eigen::VectorXd& xx, const Eigen::VectorXd& dd, double& step, int N)
{
	double temp_t = numeric_limits<double>::infinity();
	int f0, f1, f2;
	double a, b, c, b1, b2, tt, tt1, tt2;
	double x0, x1, x2, x3, x4, x5, d0, d1, d2, d3, d4, d5;
	const double* x = xx.data();
	const double* d = dd.data();
	for (int i = 0; i < F_N; ++i)
	{
		f0 = p_F0[N][i];
		f1 = p_F1[N][i];
		f2 = p_F2[N][i];

		x0 = x[f0]; x1 = x[f1]; x2 = x[f2]; x3 = x[f0 + p_V_N[N]]; x4 = x[f1 + p_V_N[N]]; x5 = x[f2 + p_V_N[N]];
		d0 = d[f0]; d1 = d[f1]; d2 = d[f2]; d3 = d[f0 + p_V_N[N]]; d4 = d[f1 + p_V_N[N]]; d5 = d[f2 + p_V_N[N]];

		a = (d1 - d0) * (d5 - d3) - (d4 - d3) * (d2 - d0);
		b1 = (d1 - d0) * (x5 - x3) + (x1 - x0) * (d5 - d3);
		b2 = (x4 - x3) * (d2 - d0) + (x2 - x0) * (d4 - d3);
		b = b1 - b2;
		c = (x1 - x0) * (x5 - x3) - (x4 - x3) * (x2 - x0);
		tt = get_smallest_pos_quad_zero(a, b, c);
		if (temp_t > tt)
		{
			temp_t = tt;
		}

	}
	step = temp_t;
}

void GAP::backtracking_line_search(const Eigen::VectorXd& x, const Eigen::VectorXd& d, const Eigen::VectorXd& negetive_grad, double& alpha, int N)
{
	double h = 0.5;
	double tt = -(negetive_grad.transpose() * d)(0, 0);
	double c = 0.2;
	double ex;
	Energy(x, ex, N);
	double e;
	Eigen::VectorXd x_new = x + alpha * d;
	Energy(x_new, e, N);
	while (e > ex + alpha * c * tt)
	{
		alpha = h * alpha;
		x_new = x + alpha * d;
		Energy(x_new, e, N);
	}
}

void GAP::Energy(const Eigen::VectorXd& position, double& energyupdate, int N)
{
	double energy = 0;
	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_d;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;
	const double* pos = position.data();
	for (int i = 0; i < p_F_N[N]; ++i)
	{
		f0 = p_F0[N][i];
		f1 = p_F1[N][i];
		f2 = p_F2[N][i];

		x0 = pos[f0];
		y0 = pos[f0 + p_V_N[N]];

		x1 = pos[f1];
		y1 = pos[f1 + p_V_N[N]];

		x2 = pos[f2];
		y2 = pos[f2 + p_V_N[N]];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = p_update_p00[N][i]; p01 = p_update_p01[N][i]; p10 = p_update_p10[N][i]; p11 = p_update_p11[N][i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;


		det = j00 * j11 - j01 * j10;
		E_d = (1 + 1 / (det * det)) * (j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);

		energy += area[i] * E_d;
	}
	energyupdate = energy;
}

void GAP::CM(int N)
{
	double area_now;
	int f0, f1, f2;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	double x0, y0, x1, y1, x2, y2;

	double hi_0, hi_1;

	double alpha_0, alpha_1, beta_0, beta_1;

	double s1, s2, sig0, sig1;

	double alpha_norm, beta_norm;
	double h_u, h_v, walpha, wbeta;

	double a1x0, a1x1, a1x2, a1x3, a1x4, a1x5,
		a2x0, a2x1, a2x2, a2x3, a2x4, a2x5;

	double aa, bb;
	double uu, vv, uv;
	double u, v;

	double h00, h01, h02, h03, h04, h05,
		h11, h12, h13, h14, h15,
		h22, h23, h24, h25,
		h33, h34, h35,
		h44, h45,
		h55;
	double* position = p_position_of_mesh[N].data();
	int nnz = p_pardiso_ja[N].size();
	p_pardiso_a[N].clear(); p_pardiso_b[N].clear();
	p_pardiso_a[N].resize(nnz, 0.0);
	p_pardiso_b[N].resize(2 * p_V_N[N], 0.0);

	for (int i = 0; i < p_F_N[N]; ++i)
	{
		area_now = area[i];
		f0 = p_F0[N][i];
		f1 = p_F1[N][i];
		f2 = p_F2[N][i];

		x0 = position[f0];
		y0 = position[f0 + p_V_N[N]];

		x1 = position[f1];
		y1 = position[f1 + p_V_N[N]];

		x2 = position[f2];
		y2 = position[f2 + p_V_N[N]];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = p_update_p00[N][i]; p01 = p_update_p01[N][i]; p10 = p_update_p10[N][i]; p11 = p_update_p11[N][i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11;  beta_1 = j10 + j01;

		alpha_norm = 0.5 * sqrt(alpha_0 * alpha_0 + alpha_1 * alpha_1);
		beta_norm = 0.5 * sqrt(beta_0 * beta_0 + beta_1 * beta_1);

		if (beta_norm < 1e-30)
			beta_norm = 1e-15;
		s1 = (p00) * (p00 + p10) + (p01) * (p01 + p11);
		s2 = (p10) * (p00 + p10) + (p11) * (p01 + p11);

		double h1 = p00 * p00 + p01 * p01;
		double h2 = p00 * p10 + p01 * p11;
		double h3 = p10 * p10 + p11 * p11;
		double h4 = p00 * p11 - p01 * p10;

		a1x0 = alpha_0 * (-p00 - p10) + alpha_1 * (p01 + p11);  a1x1 = alpha_0 * p00 - alpha_1 * p01; a1x2 = alpha_0 * p10 - alpha_1 * p11;
		a1x3 = alpha_0 * (-p01 - p11) + alpha_1 * (-p00 - p10); a1x4 = alpha_0 * p01 + alpha_1 * p00; a1x5 = alpha_0 * p11 + alpha_1 * p10;

		a2x0 = beta_0 * (-p00 - p10) + beta_1 * (-p01 - p11);   a2x1 = beta_0 * p00 + beta_1 * p01;   a2x2 = beta_0 * p10 + beta_1 * p11;
		a2x3 = beta_0 * (p01 + p11) + beta_1 * (-p00 - p10);    a2x4 = -beta_0 * p01 + beta_1 * p00;  a2x5 = -beta_0 * p11 + beta_1 * p10;

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		hi_0 = 2 + 6 * 1 / (sig0 * sig0 * sig0 * sig0); hi_1 = 2 + 6 * 1 / (sig1 * sig1 * sig1 * sig1);

		aa = 0.25 / alpha_norm; bb = 0.25 / beta_norm;

		uu = aa * aa * (area_now * hi_0 + area_now * hi_1);
		vv = bb * bb * (area_now * hi_0 + area_now * hi_1);
		uv = aa * bb * (area_now * hi_0 - area_now * hi_1);

		h_u = area_now * (2 * sig0 - 2 * 1 / (sig0 * sig0 * sig0));
		h_v = area_now * (2 * sig1 - 2 * 1 / (sig1 * sig1 * sig1));

		walpha = h_u + h_v;
		wbeta = h_u - h_v;

		double hwa1 = (walpha * 0.25 / alpha_norm); double hwa2 = -(walpha * 0.25 * 0.25 / (alpha_norm * alpha_norm * alpha_norm));
		double hwb1 = (wbeta * 0.25 / beta_norm); double hwb2 = -(wbeta * 0.25 * 0.25 / (beta_norm * beta_norm * beta_norm));


		h00 = uu * a1x0 * a1x0 + vv * a2x0 * a2x0 + uv * a1x0 * a2x0 + uv * a2x0 * a1x0; h01 = uu * a1x0 * a1x1 + vv * a2x0 * a2x1 + uv * a1x0 * a2x1 + uv * a2x0 * a1x1; h02 = uu * a1x0 * a1x2 + vv * a2x0 * a2x2 + uv * a1x0 * a2x2 + uv * a2x0 * a1x2; h03 = uu * a1x0 * a1x3 + vv * a2x0 * a2x3 + uv * a1x0 * a2x3 + uv * a2x0 * a1x3; h04 = uu * a1x0 * a1x4 + vv * a2x0 * a2x4 + uv * a1x0 * a2x4 + uv * a2x0 * a1x4; h05 = uu * a1x0 * a1x5 + vv * a2x0 * a2x5 + uv * a1x0 * a2x5 + uv * a2x0 * a1x5;

		h11 = uu * a1x1 * a1x1 + vv * a2x1 * a2x1 + uv * a1x1 * a2x1 + uv * a2x1 * a1x1; h12 = uu * a1x1 * a1x2 + vv * a2x1 * a2x2 + uv * a1x1 * a2x2 + uv * a2x1 * a1x2; h13 = uu * a1x1 * a1x3 + vv * a2x1 * a2x3 + uv * a1x1 * a2x3 + uv * a2x1 * a1x3; h14 = uu * a1x1 * a1x4 + vv * a2x1 * a2x4 + uv * a1x1 * a2x4 + uv * a2x1 * a1x4; h15 = uu * a1x1 * a1x5 + vv * a2x1 * a2x5 + uv * a1x1 * a2x5 + uv * a2x1 * a1x5;

		h22 = uu * a1x2 * a1x2 + vv * a2x2 * a2x2 + uv * a1x2 * a2x2 + uv * a2x2 * a1x2; h23 = uu * a1x2 * a1x3 + vv * a2x2 * a2x3 + uv * a1x2 * a2x3 + uv * a2x2 * a1x3; h24 = uu * a1x2 * a1x4 + vv * a2x2 * a2x4 + uv * a1x2 * a2x4 + uv * a2x2 * a1x4; h25 = uu * a1x2 * a1x5 + vv * a2x2 * a2x5 + uv * a1x2 * a2x5 + uv * a2x2 * a1x5;

		h33 = uu * a1x3 * a1x3 + vv * a2x3 * a2x3 + uv * a1x3 * a2x3 + uv * a2x3 * a1x3; h34 = uu * a1x3 * a1x4 + vv * a2x3 * a2x4 + uv * a1x3 * a2x4 + uv * a2x3 * a1x4; h35 = uu * a1x3 * a1x5 + vv * a2x3 * a2x5 + uv * a1x3 * a2x5 + uv * a2x3 * a1x5;

		h44 = uu * a1x4 * a1x4 + vv * a2x4 * a2x4 + uv * a1x4 * a2x4 + uv * a2x4 * a1x4; h45 = uu * a1x4 * a1x5 + vv * a2x4 * a2x5 + uv * a1x4 * a2x5 + uv * a2x4 * a1x5;

		h55 = uu * a1x5 * a1x5 + vv * a2x5 * a2x5 + uv * a1x5 * a2x5 + uv * a2x5 * a1x5;

		if (walpha >= 0)
		{
			h00 += hwa1 * (s1 + s2) + hwa2 * a1x0 * a1x0; h01 += hwa1 * (-s1) + hwa2 * a1x0 * a1x1; h02 += hwa1 * (-s2) + hwa2 * a1x0 * a1x2; h03 += hwa2 * a1x0 * a1x3; h04 += hwa1 * (h4)+hwa2 * a1x0 * a1x4; h05 += hwa1 * (-h4) + hwa2 * a1x0 * a1x5;
			h11 += hwa1 * (h1)+hwa2 * a1x1 * a1x1;        h12 += hwa1 * (h2)+hwa2 * a1x1 * a1x2;    h13 += hwa1 * (-h4) + hwa2 * a1x1 * a1x3; h14 += hwa2 * a1x1 * a1x4; h15 += hwa1 * (h4)+hwa2 * a1x1 * a1x5;
			h22 += hwa1 * (h3)+hwa2 * a1x2 * a1x2;        h23 += hwa1 * (h4)+hwa2 * a1x2 * a1x3;    h24 += hwa1 * (-h4) + hwa2 * a1x2 * a1x4; h25 += hwa2 * a1x2 * a1x5;
			h33 += hwa1 * (s1 + s2) + hwa2 * a1x3 * a1x3; h34 += hwa1 * (-s1) + hwa2 * a1x3 * a1x4; h35 += hwa1 * (-s2) + hwa2 * a1x3 * a1x5;
			h44 += hwa1 * (h1)+hwa2 * a1x4 * a1x4;        h45 += hwa1 * (h2)+hwa2 * a1x4 * a1x5;
			h55 += hwa1 * (h3)+hwa2 * a1x5 * a1x5;

		}
		h00 += hwb1 * (s1 + s2) + hwb2 * a2x0 * a2x0; h01 += hwb1 * (-s1) + hwb2 * a2x0 * a2x1; h02 += hwb1 * (-s2) + hwb2 * a2x0 * a2x2; h03 += hwb2 * a2x0 * a2x3; h04 += hwb1 * (-h4) + hwb2 * a2x0 * a2x4; h05 += hwb1 * (h4)+hwb2 * a2x0 * a2x5;
		h11 += hwb1 * (h1)+hwb2 * a2x1 * a2x1;        h12 += hwb1 * (h2)+hwb2 * a2x1 * a2x2;    h13 += hwb1 * (h4)+hwb2 * a2x1 * a2x3;    h14 += hwb2 * a2x1 * a2x4; h15 += hwb1 * (-h4) + hwb2 * a2x1 * a2x5;
		h22 += hwb1 * (h3)+hwb2 * a2x2 * a2x2;        h23 += hwb1 * (-h4) + hwb2 * a2x2 * a2x3; h24 += hwb1 * (h4)+hwb2 * a2x2 * a2x4;    h25 += hwb2 * a2x2 * a2x5;
		h33 += hwb1 * (s1 + s2) + hwb2 * a2x3 * a2x3; h34 += hwb1 * (-s1) + hwb2 * a2x3 * a2x4; h35 += hwb1 * (-s2) + hwb2 * a2x3 * a2x5;
		h44 += hwb1 * (h1)+hwb2 * a2x4 * a2x4;        h45 += hwb1 * (h2)+hwb2 * a2x4 * a2x5;
		h55 += hwb1 * (h3)+hwb2 * a2x5 * a2x5;

		u = aa * walpha; v = bb * wbeta;

		p_pardiso_b[N][f0] -= (u * a1x0 + v * a2x0);
		p_pardiso_b[N][f1] -= (u * a1x1 + v * a2x1);
		p_pardiso_b[N][f2] -= (u * a1x2 + v * a2x2);
		p_pardiso_b[N][f0 + p_V_N[N]] -= (u * a1x3 + v * a2x3);
		p_pardiso_b[N][f1 + p_V_N[N]] -= (u * a1x4 + v * a2x4);
		p_pardiso_b[N][f2 + p_V_N[N]] -= (u * a1x5 + v * a2x5);

		p_pardiso_a[N][p_id_h00[N][i]] += h00; p_pardiso_a[N][p_id_h01[N][i]] += h01; p_pardiso_a[N][p_id_h02[N][i]] += h02; p_pardiso_a[N][p_id_h03[N][i]] += h03; p_pardiso_a[N][p_id_h04[N][i]] += h04; p_pardiso_a[N][p_id_h05[N][i]] += h05;
		p_pardiso_a[N][p_id_h11[N][i]] += h11; p_pardiso_a[N][p_id_h12[N][i]] += h12; p_pardiso_a[N][p_id_h13[N][i]] += h13; p_pardiso_a[N][p_id_h14[N][i]] += h14; p_pardiso_a[N][p_id_h15[N][i]] += h15;
		p_pardiso_a[N][p_id_h22[N][i]] += h22; p_pardiso_a[N][p_id_h23[N][i]] += h23; p_pardiso_a[N][p_id_h24[N][i]] += h24; p_pardiso_a[N][p_id_h25[N][i]] += h25;
		p_pardiso_a[N][p_id_h33[N][i]] += h33; p_pardiso_a[N][p_id_h34[N][i]] += h34; p_pardiso_a[N][p_id_h35[N][i]] += h35;
		p_pardiso_a[N][p_id_h44[N][i]] += h44; p_pardiso_a[N][p_id_h45[N][i]] += h45;
		p_pardiso_a[N][p_id_h55[N][i]] += h55;
	}

	p_pardiso[N]->a = p_pardiso_a[N];
	p_pardiso[N]->rhs = p_pardiso_b[N];

	p_pardiso[N]->factorize();
	p_pardiso[N]->pardiso_solver();

	vector<double> result_d = p_pardiso[N]->result;

	Eigen::VectorXd negative_grad(2 * p_V_N[N]), d(2 * p_V_N[N]);
	for (int i = 0; i < 2 * p_V_N[N]; i++)
	{
		negative_grad(i) = p_pardiso_b[N][i];
		d(i) = result_d[i];
	}

	p_pardiso[N]->free_numerical_factorization_memory();

	double temp_t;
	max_step(p_position_of_mesh[N], d, temp_t, N);

	double alpha = 0.95 * temp_t;

	backtracking_line_search(p_position_of_mesh[N], d, negative_grad, alpha, N);

	double e1;
	double s;
	p_position_of_mesh[N] += alpha * d;
	Energysource(N);
}

void GAP::recover_to_src(int N)
{
	area = area_src;
	p_update_p00[N] = p_source_p00[N];
	p_update_p01[N] = p_source_p01[N];
	p_update_p10[N] = p_source_p10[N];
	p_update_p11[N] = p_source_p11[N];
}


void GAP::GenFirstCut()
{
	double unit=0;
	if (VertexPriorityMetric == "Neighbourhood")
		unit = 1;
	else if (VertexPriorityMetric == "RealDis")
		unit = MCache.avg_el;

	ClassifyFeaturePoints(FixThreshold*unit);
	std::ofstream fix("fix.txt");
	for (auto a : FixPoints)
		fix << a << std::endl;
	fix.close();
	std::ofstream cand("cand.txt");
	for (auto a : CanditatePoints)
		cand << a << std::endl;
	cand.close();

	std::vector<std::vector<int>> ForbiddenArea;
	std::vector<double> RadiusInfo(CanditatePoints.size(), GAPForBiddenRadius*unit);
	for (int i = 0; i < CanditatePoints.size(); i++)
	{
		if (VertexPriorityMetric == "Neighbourhood")
		{
			Algorithm::UpdateNeighbourhood(MCache, RadiusInfo[i], CanditatePoints[i]);
			for (int u = 0; u < FixPoints.size(); u++)
			{
				if (MCache.Neighbour[CanditatePoints[i]][FixPoints[u]] < RadiusInfo[i])
				{
					RadiusInfo[i] = MCache.Neighbour[CanditatePoints[i]][FixPoints[u]];
				}
			}
		}
		else if (VertexPriorityMetric == "RealDis")
		{
			Algorithm::Dijkstra_all(MCache, CanditatePoints[i]);
			for (int u = 0; u < FixPoints.size(); u++)
			{
				if (MCache.V_D[CanditatePoints[i]][FixPoints[u]] < RadiusInfo[i])
				{
					RadiusInfo[i] = MCache.V_D[CanditatePoints[i]][FixPoints[u]];
				}
			}
		}
	}
	int m = 0;
	while (1)
	{
		ForbiddenArea.clear();
		ForbiddenArea.resize(MCache.n_vertices, std::vector<int>());
		for (int i = 0; i < CanditatePoints.size(); i++)
		{
			int v = CanditatePoints[i];
			for (int u = 0; u < MCache.n_vertices; u++)
			{
				if ((VertexPriorityMetric == "Neighbourhood") && (MCache.Neighbour[v][u] < RadiusInfo[i]))
				{
					ForbiddenArea[u].push_back(i);
				}
				else 	if ((VertexPriorityMetric == "RealDis") && (MCache.V_D[CanditatePoints[i]][u]< RadiusInfo[i]))
				{
					ForbiddenArea[u].push_back(i);
				}
			}
		}
		bool close = true;
		std::vector<int> flag(CanditatePoints.size(), 0);
		for (int i = 0; i < ForbiddenArea.size(); i++)
		{
			if (ForbiddenArea[i].size() >= 2)
			{
				for (auto a : ForbiddenArea[i])
				{
					close = false;
					flag[a] = 1;
				}
			}
		}
		if (close)
			break;
		for (int i = 0; i < flag.size(); i++)
		{
			RadiusInfo[i] -= flag[i]*unit;
		}
		std::ofstream out(std::to_string(m) + ".txt");
		out << "VERTICES" << std::endl;
		for (int i = 0; i < ForbiddenArea.size(); i++)
		{
			if (ForbiddenArea[i].size() != 0)
				out << i << std::endl;
		}
		out.close();
		m++;
		}
	for (int i = 0; i < RadiusInfo.size(); i++)
		RadiusInfo[i] -= unit;
	ForbiddenArea.clear();
	ForbiddenArea.resize(MCache.n_vertices, std::vector<int>());
	for (int i = 0; i < CanditatePoints.size(); i++)
	{
		int v = CanditatePoints[i];
		for (int u = 0; u < MCache.n_vertices; u++)
		{
			if ((VertexPriorityMetric == "Neighbourhood") && (MCache.Neighbour[v][u] < RadiusInfo[i]))
			{
				ForbiddenArea[u].push_back(i);
			}
			else 	if ((VertexPriorityMetric == "RealDis") && (MCache.V_D[CanditatePoints[i]][u])< RadiusInfo[i])
			{
				ForbiddenArea[u].push_back(i);
			}
		}
	}
	std::vector<double> weight(MCache.n_edges, -1);
	for (int i = 0; i < MCache.n_edges; i++)
		weight[i] = MCache.el[i];
	for (int i = 0; i < ForbiddenArea.size(); i++)
	{
		if (ForbiddenArea[i].size() != 0)
		{
			for (int u = 0; u < MCache.ve[i].size(); u++)
			{
				weight[MCache.ve[i][u]] = DBL_MAX;
			}
		}
	}
	std::priority_queue<PathInfo> que;
	for (int i = 0; i < FixPoints.size(); i++)
	{//给定所有路径获得其生成树
		std::vector<int> v_p(MCache.n_vertices, -1);
		std::vector<double> d(MCache.n_vertices, DBL_MAX);
		Algorithm::Dijkstra_with_restrict(MCache, FixPoints[i], weight, v_p, d);
		for (int j = i + 1; j < FixPoints.size(); j++)
		{
			std::vector<int> path;
			int s_p = FixPoints[i];
			int e_p = FixPoints[j];
			double length = 0;
			while (e_p != s_p)
			{
				path.push_back(e_p);
				length += MCache.el[MCache.vve[e_p][v_p[e_p]]];
				e_p = v_p[e_p];
			}
			path.push_back(e_p);
			PathInfo tmp(FixPoints[i], FixPoints[j], length);
			tmp.path = path;
			que.push(tmp);
		}
	}
	std::vector<PathInfo> mst;
	//使用Kruskal算法求出最小生成树
	Algorithm::Kruskal(FixPoints, que, mst);
	std::set<int> e, v;
	for (auto a : mst)
	{
		for (int i = 0; i < a.path.size() - 1; i++)
		{
			e.insert(MCache.vve[a.path[i]][a.path[i + 1]]);
			v.insert(a.path[i]);
			v.insert(a.path[i + 1]);
		}
	}
	std::vector<int> cv, ce;
	for (auto a : v)
		cv.push_back(a);
	for (auto a : e)
		ce.push_back(a);
	/*
	std::ofstream m2("mst.txt");
	m2 << "VERTICES" << std::endl;
	for (auto a : v)
		m2 << a << std::endl;
	m2 << "EDGES" << std::endl;
	for (auto a : e)
		m2 << a << std::endl;
	m2.close();
	*/
	meshcut = std::make_unique<MeshCut>(const_cast<Mesh&>(Closedmesh),MCache);
	meshcut->SetCut(cv, ce);
	v_seam.resize(Closedmesh.n_vertices());
	e_seam.resize(Closedmesh.n_edges());
	for (int i = 0; i < cv.size(); i++)
	{
		v_seam[cv[i]] = 1;
	}
	for (int i = 0; i < ce.size(); i++)
	{
		e_seam[ce[i]] = 1;
	}
	meshcut->MakeSeam();
	//meshcut->get_correspondence();
	mesh = meshcut->GetCutedMesh();
	meshcut->get_correspondence(he2idx, idx2meshvid);
	OpenMesh::IO::write_mesh(mesh,"tmp.obj");
	V_N = mesh.n_vertices();
	F_N =mesh.n_faces();
	area.resize(F_N);
	F0.resize(F_N);
	F1.resize(F_N);
	F2.resize(F_N);
	position_of_mesh.resize(2 * V_N);
	negative_grad_norm.resize(2 * V_N);
	mesh.request_face_normals();
	mesh.update_normals();
	//设置Pardiso的相关参数
}

void GAP::ClassifyFeaturePoints(double threshold)
{
	FixPoints.clear();
	CanditatePoints.clear();
	//按照优先级进行排序,匿名函数
	std::sort(FeaturePoints.begin(), FeaturePoints.end(), [&](std::pair<int, double>& x, std::pair<int, double>& y) {return x.second > y.second; });
	for (auto a : FeaturePoints)
	{
		if (a.second > threshold)
			FixPoints.push_back(a.first);
		else
			CanditatePoints.push_back(a.first);
	}
}