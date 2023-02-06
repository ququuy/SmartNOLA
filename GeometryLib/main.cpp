#include "IO.h"
#include "GeometryLib.h"
#include "cmake_definition.h"



int main() {
	// ----- Region Growing & Coloring Result 
	GEO::PN3_Range points = IO::read_PLY((std::string(SOLUTION_ROOT_PATH) + "/data/test.ply").c_str());
	std::vector<GEO::PN3_Range> points_of_planes = GEO::detect_planes_growing(points);
	GEO::PNC3_Range pnc_planes = GEO::coloring_PN3(points_of_planes);
	IO::write_PNC3((std::string(SOLUTION_ROOT_PATH) + "/data/output/test_planes.ply").c_str(), pnc_planes);

	// ----- Plane Clustering & Coloring Result
	ALGO::PlaneData_Range plane_datas;
	for (const auto& pn3_range : points_of_planes) {
		GEO::Point3_Range p3_range = GEO::PN3Range_to_Point3Range(pn3_range);
		plane_datas.push_back(ALGO::project_points(p3_range));
	}
	std::vector<ALGO::PlaneData_Range> planes_clustered;
	ALGO::plane_clustering(plane_datas, planes_clustered);

	GEO::PNC3_Range pnc_clustered_planes = ALGO::extract_cluster_result(planes_clustered);
	IO::write_PNC3((std::string(SOLUTION_ROOT_PATH) + "/data/output/test_planes_clustered.ply").c_str(), pnc_clustered_planes);


	// ----- Translation Clustering
	ALGO::trans_cluster_test(planes_clustered);
	return 0;
}