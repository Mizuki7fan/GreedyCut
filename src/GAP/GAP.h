#pragma once
#include "../MeshDefinition/MeshDefinition.h"
#include "../MeshDefinition/MeshCache.h"
#include "../Auxiliary/Algorithm.h"
#include "../MeshCut/MeshCut.h"
#include"Eigen/Dense"
#include"Eigen/Sparse"
#include"Eigen/SparseLU"
#include "Eigen/SVD"
#include "..\PardisoSolver\PardisoSolver.h"

class GAP
{
public:
	GAP(Mesh& mesh,MeshCache& MCache);
	void Set(const std::vector<std::pair<int, double>>& FP,
		std::string PriorityMetric, 
		int FixThreshold, 
		int ForbiddenRadius,
		double FilteringRate);
	void SetPardiso(double conv_rate,int MaxIter,int bound_distortion_K);

	void gradually_addp_pipeline();
	//函数：run_bpe()
	void run_bpe();
	//run_bpe()-
	void bpe_init();
	//run_bpe()-bpe_init()-
	void init_src_matirx();
	//run_bpe()-bpe_init()-init_src_matrix()-
	void local_coordinate_inverse(int i, double& p00, double& p01, double& p10, double& p11);
	//run_bpe()-bpe_init()-
	void Pre_calculate();
	//run_bpe()-bpe_init()-
	void TutteOp();
	//run_bpe()-bpe_init()-
	void init_area();
	//run_bpe()-
	void BPE();
	//run_bpe()-BPE()-
	void Energysource();
	//run_bpe()-BPE()-
	void Update_source_same_t();
	//run_bpe()-BPE()-Update_source_same_t()-
	double newton_equation(const double& a, const double& b, const double& K);
	//run_bpe()-BPE()-
	void SLIM();
	//run_bpe()-BPE()-SLIM()-
	void max_step(const Eigen::VectorXd& xx, const Eigen::VectorXd& dd, double& step);
	//run_bpe()-BPE()-SLIM()-max_step()-
	double get_smallest_pos_quad_zero(double a, double b, double c);
	//run_bpe()-BPE()-SLIM()-
	void backtracking_line_search(const Eigen::VectorXd& x, const Eigen::VectorXd& d, const Eigen::VectorXd& negetive_grad, double& alpha);
	//run_bpe()-BPE()-SLIM()-backtracking_line_search()-
	void Energy(const Eigen::VectorXd& position, double& energyupdate);
	//run_bpe()-BPE()-
	void CM();
	//run_bpe()-BPE()-
	void recover_to_src();




	//计算将所有固定点链接而得到的初始切割
	void GenFirstCut();
	void ClassifyFeaturePoints(double threshold);
	void InitCut();

private:
	Mesh& mesh;
	MeshCache MCache;
	//点的优先级的度量
	std::string VertexPriorityMetric = "Neighbourhood";
	int FixThreshold = 20;
	int GAPForBiddenRadius = 5;
	//MeshCut类需要构造函数,这里只能先创建其指针
	std::unique_ptr<MeshCut> meshcut;

	//备选的特征点
	std::vector<std::pair<int, double>> FeaturePoints;
	std::vector<int> FixPoints;
	std::vector<int> CanditatePoints;

	std::vector<double> source_p00;
	std::vector<double> source_p01;
	std::vector<double> source_p10;
	std::vector<double> source_p11;

	std::vector<double> update_p00;
	std::vector<double> update_p01;
	std::vector<double> update_p10;
	std::vector<double> update_p11;


private:
	int F_N, V_N;
	double g_norm;
	std::vector<double> area;
	std::vector<double> area_uniform;
	std::vector<double> area_src;

	std::vector<int> F0, F1, F2;
	std::vector<std::vector<int>> VV_ids;

	Eigen::VectorXd position_of_mesh;
	Eigen::VectorXd negative_grad_norm;
	PardisoSolver* pardiso;
	std::vector<int> pardiso_i, pardiso_ia, pardiso_ja;
	std::vector<double> pardiso_a, pardiso_b;

	double time_consumption;

	double Intp_T_Min;
	double changetocm_flag;
	double energy_prev_seam;
	double energy_uniform;
	double energy_area;

	double convgence_con_rate;
	int MAX_ITER_NUM;
	double bound_distortion_K;
	//过滤阈值
	double filtering_rate = 0.01;

	std::vector<int> id_h00; std::vector<int> id_h01; std::vector<int> id_h02; std::vector<int> id_h03; std::vector<int> id_h04; std::vector<int> id_h05;
	std::vector<int> id_h11; std::vector<int> id_h12; std::vector<int> id_h13; std::vector<int> id_h14; std::vector<int> id_h15;
	std::vector<int> id_h22; std::vector<int> id_h23; std::vector<int> id_h24; std::vector<int> id_h25;
	std::vector<int> id_h33; std::vector<int> id_h34; std::vector<int> id_h35;
	std::vector<int> id_h44; std::vector<int> id_h45;
	std::vector<int> id_h55;

	


};