#include <io.h>
#include <time.h>
#include <string>
#include "Auxiliary.h"
#include "MeshCache.h"
#include "PointSampling.h"
#include "MeshCut.h"
#include "KPNewton.h"
#include "PointFinding.h"
#include "GAP.h"
#include "AddAuxiliaryPoint.h"

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		printf("error parameter count!\n");
		return -1;
	}
	time_t st = clock();
	Option opt(argv[1], argv[2]);

	Mesh mesh;
	OpenMesh::IO::read_mesh(mesh, argv[1]);
	MeshCache MC(mesh);

	if (MC.NVertices == 0)
	{
		printf("n_vertices is 0");
		return -1;
	}
	printf("Find Cut 1\n");
	PointSampling PS(MC);
	PS.Set(opt.PS_method);
	std::vector<int> SamplePoints;
	PS.ComputeSamples(SamplePoints);
	
	MeshCut Mcut(mesh,MC);
	Mcut.Set(SamplePoints);
	Mcut.Connect();
	Mcut.MakeSeam();
	Mesh CutedMesh = Mcut.GetCutedMesh();
	printf("KPNewton1\n");

	KPNewton kpn(CutedMesh);
	kpn.Tutte();
	kpn.PrepareDataFree();
	kpn.RunFree(KPNewton::MIPS);
	kpn.RunFree(KPNewton::AMIPS);
	kpn.UpdateMesh();
	kpn.ResultRect(mesh);
	if (kpn.EnergyNan())
	{
		printf("fail in KPNewton!\n");
		return -1;
	}
	printf("PF1\n");
	PointFinding PF(mesh, CutedMesh, MC);
	//PF.Set(opt.PF_vertex_priority_metric);
	PF.FindLocalMaximizer();

	printf("Cut2\n");
	MeshCut Mcut2(mesh, MC);
	Mcut2.SetBanCondition(PF.GetLocalMaximizer(), Mcut.GetCutvertex(),opt.BanArea_Method);
	Mcut2.CalcBanArea(opt.BanArea_Dn, opt.BanArea_Metric, opt.BanArea_alpha, opt.BanArea_ShrinkRate);
	std::vector<int> AllowedArea = Mcut2.GetMaxConnectedRegion();
	std::vector<int> cut2Edges, cut2Vertices;
	PS.ComputeSamplesFromSelectedRegion(opt.PS_method, AllowedArea, cut2Vertices, cut2Edges);
	Mcut2.SetCut(cut2Vertices, cut2Edges);
	Mcut2.MakeSeam();
	Mesh CutedMesh2 = Mcut2.GetCutedMesh();

	printf("KPNewton2\n");
	KPNewton kpn2(CutedMesh2);
	kpn2.Tutte();
	kpn2.PrepareDataFree();
	kpn2.RunFree(KPNewton::MIPS);
	kpn2.RunFree(KPNewton::AMIPS);
	kpn2.UpdateMesh();
	kpn2.ResultRect(mesh);
	if (kpn2.EnergyNan())
	{
		std::cerr << "fail in KPNewton2!" << std::endl;
		return -1;
	}

	printf("PF2\n");
	PointFinding PF2(mesh, CutedMesh2, MC);
	PF2.Set(opt.PS_method);
	std::vector<std::pair<int, double>> VertexPriority;
	PF2.FindLocalMaximizer();
	PF2.Find(VertexPriority);

	printf("GAP\n");
	GAP gap(mesh, MC,VertexPriority);
	gap.Set(opt.Influence_Threshold,opt.Distortion_Threshold);
	gap.SetSolver(1e-6, 500, 250);
	gap.Run();
	std::vector<int> gapResult = gap.getResult();

	printf("AAP\n");
	AAP aap(mesh, MC, gapResult);
	aap.Set(opt.Trimming_Rate, opt.Max_AddCount);
	aap.Run();

	std::vector<int> Result = aap.GetResult();
	for (auto a : Result)
		printf("%i\n", a);
	printf("Fin\n");

}