#pragma once
#include "MeshCache.h"
#include "Solver/Solver.h"
#include <complex>
#include <iostream>

class KPNewton
{
public:
	KPNewton(Mesh& m);
	~KPNewton();
	enum EnergyType
	{
		MIPS, AMIPS
	};

	void Tutte();
	// Prepare data for free boundary parameterization
	void PrepareDataFree(void);
	// Optimize parameterization using a type of energy
	void RunFree(const EnergyType& etype);
	bool EnergyNan() { return EnergyIsNan; };

	// Compute parameters using MIPS
	static void ComputeMIPS(double& energy, double& alpha1, double& alpha2,
		double& beta1, double& beta2, double& beta3,
		const double& x, const double& y);

	// Compute MIPS energy
	static double ComputeEnergyMIPS(const double& x, const double& y);

	// Compute parameters using AMIPS
	static void ComputeAMIPS(double& energy, double& alpha1, double& alpha2,
		double& beta1, double& beta2, double& beta3,
		const double& x, const double& y);

	// Compute AMIPS energy
	static double ComputeEnergyAMIPS(const double& x, const double& y);
	// Compute maximum step size in line search
	double ComputeTMax(const std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& d) const;

	double ComputeEnergy(const std::vector<std::complex<double>>& pos,
		const std::vector<std::vector<std::complex<double>>>& D,
		const std::vector<std::vector<std::complex<double>>>& DC,
		const std::vector<double>& area) const;
	// Update mesh vertices
	void UpdateMesh(void);

	void ResultRect(const Mesh& mesh);

private:
	std::vector<Mesh::VertexHandle> GetBoundary(void) const;
	// Compute edge length of each face
	std::vector<std::vector<double>> MeshFaceEdgeLen2s(void) const;
	// Compute angles in each triangles
	std::vector<std::vector<double>> MeshAnglesFromFaceEdgeLen2(const std::vector<std::vector<double>>& len2) const;

	std::vector< std::complex<double>> result;
	std::vector<std::vector<int>> tri;
	std::vector<int> assembleorder;

	Mesh& mesh;
	std::vector<int> isbv;

	decltype(&ComputeMIPS) computeall; // parameter functional
	decltype(&ComputeEnergyMIPS) computeenergy; // energy functional
	bool EnergyIsNan = false;

public:
	Solver* solver;

};