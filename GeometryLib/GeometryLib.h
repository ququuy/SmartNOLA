// GeometryLib.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>
#include <cmath>
#include "cmake_definition.h" // temp
#define _USE_MATH_DEFINES

// ---------- BASIC ------------- //
namespace MA {
	inline float RandomReal() {
		return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	}
}

// ----------- GLM ----------- //
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>


// ----------  CGAL ------------- //
#include <CGAL/boost/graph/IO/PLY.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Timer.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/write_points.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>
#include "CGAL/Aff_transformation_3.h"
#include "CGAL/aff_transformation_tags.h"
#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/vcm_estimate_edges.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/intersections.h>
/* Alpha Shape */
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>

/* Shape Regularization */
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_contours.h>


namespace GEO {
	typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
	typedef CGAL::Parallel_if_available_tag Concurrency_tag;

	typedef Kernel::FT FT;

	typedef Kernel::Point_3 Point_3;
	typedef Kernel::Point_2 Point_2;
	typedef Kernel::Vector_3 Vector_3;
	typedef Kernel::Vector_2 Vector_2;
	typedef Kernel::Line_2                   Line_2;
	typedef Kernel::Line_3                   Line_3;
	typedef Kernel::Segment_2				 Segment_2;
	typedef Kernel::Segment_3				 Segment_3;
	typedef Kernel::Plane_3                  Plane_3;
	typedef Kernel::Direction_2				 Direction_2;

	typedef std::array<unsigned char, 4> Color;

	typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;

	typedef std::tuple<Point_3, Vector_3> PN_3;
	typedef CGAL::Nth_of_tuple_property_map<0, PN_3> PN3_Point_map;
	typedef CGAL::Nth_of_tuple_property_map<1, PN_3> PN3_Normal_map;
	typedef std::vector<PN_3> PN3_Range;


	typedef std::tuple<Point_2, Vector_2> PN_2;
	typedef CGAL::Nth_of_tuple_property_map<0, PN_2> PN2_Point_map;
	typedef CGAL::Nth_of_tuple_property_map<1, PN_2> PN2_Normal_map;
	typedef std::vector<PN_2> PN2_Range;

	typedef std::tuple<Point_3, Vector_3, Color> PNC_3;
	typedef CGAL::Nth_of_tuple_property_map<0, PNC_3> PNC3_Point_map;
	typedef CGAL::Nth_of_tuple_property_map<1, PNC_3> PNC3_Normal_map;
	typedef CGAL::Nth_of_tuple_property_map<2, PNC_3> PNC3_Color_map;
	typedef std::vector<PNC_3> PNC3_Range;
	typedef std::vector<Point_3> Point3_Range;
	typedef std::vector<Point_2> Point2_Range;
	typedef std::vector<Segment_3> Segment3_Range;
	typedef std::vector<Segment_2> Segment2_Range;
	typedef std::vector<Direction_2> Direction2_Range;
	typedef Point2_Range Contour;

	typedef CGAL::Convex_hull_traits_adapter_2<Kernel, CGAL::Pointer_property_map<Point_2>::type > Convex_hull_traits_2;


	// ---- plane detection
	using Sphere_Neighbor_query_3d = CGAL::Shape_detection::Point_set::Sphere_neighbor_query<Kernel, PN3_Range, PN3_Point_map>;
	using K_Neighbor_query_3d = CGAL::Shape_detection::Point_set::K_neighbor_query<Kernel, PN3_Range, PN3_Point_map>;
	using Region_type_3d = CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region<Kernel, PN3_Range, PN3_Point_map, PN3_Normal_map>;
	//using Region_growing_3d = CGAL::Shape_detection::Region_growing<PN3_Range, Sphere_Neighbor_query_3d, Region_type_3d>;
	using Region_growing_3d = CGAL::Shape_detection::Region_growing<PN3_Range, K_Neighbor_query_3d, Region_type_3d>;

	namespace PlaneDetection {
		struct RegionGrowingParams {
			int k;//TODO 
		};
	}


	// ---- line detection (2d)
	using Sphere_Neighbor_query_2d = CGAL::Shape_detection::Point_set::Sphere_neighbor_query<Kernel, PN2_Range, PN2_Point_map>;
	using K_Neighbor_query_2d = CGAL::Shape_detection::Point_set::K_neighbor_query<Kernel, PN2_Range, PN2_Point_map>;
	using Region_type_2d = CGAL::Shape_detection::Point_set::Least_squares_line_fit_region<Kernel, PN2_Range, PN2_Point_map, PN2_Normal_map>;
	//using Region_growing_2d = CGAL::Shape_detection::Region_growing<PN2_Range, Sphere_Neighbor_query_2d, Region_type_2d>;
	using Region_growing_2d = CGAL::Shape_detection::Region_growing<PN2_Range, K_Neighbor_query_2d, Region_type_2d>;




	// --- alpha shape`
	using AlphaShape_Vb = CGAL::Alpha_shape_vertex_base_2<Kernel>;//                   Vb;
	using AlphaShape_Fb = CGAL::Alpha_shape_face_base_2<Kernel>;//                     Fb;
	using AlphaShape_Tds = CGAL::Triangulation_data_structure_2<AlphaShape_Vb, AlphaShape_Fb>;//          Tds;
	using Triangulation_2 = CGAL::Delaunay_triangulation_2<Kernel, AlphaShape_Tds>;//                Triangulation_2;
	using AlphaShape_2 = CGAL::Alpha_shape_2<Triangulation_2>;//                 Alpha_shape_2;

	// --- shape regularization
	using Contour_Directions =
		CGAL::Shape_regularization::Contours::Multiple_directions_2<Kernel, Contour>;
	using Contour_Directions_Custom =
		CGAL::Shape_regularization::Contours::User_defined_directions_2<Kernel, Contour>;

}

namespace GEO {
	//
	inline glm::vec3 p3_to_glm(const Point_3& p) { return glm::vec3(p.x(), p.y(), p.z()); }
	inline glm::vec2 p2_to_glm(const Point_2& p) { return glm::vec2(p.x(), p.y()); }
	inline Point3_Range PN3Range_to_Point3Range(const PN3_Range& pns) {
		Point3_Range result;
		for (const auto& pn : pns) {
			result.push_back(std::get<0>(pn));
		}
		return result;
	}
	inline FT cross_2(const glm::vec2& v0, const glm::vec2& v1) { return v0.x * v1.y - v0.y * v1.x; }
	inline Color rand_color() {
		const Color color = {
				static_cast<unsigned char>(rand() % 256),
				static_cast<unsigned char>(rand() % 256),
				static_cast<unsigned char>(rand() % 256), 255 };
		return color;
	}
	

	void normal_estimate(PN3_Range& pointsv);
	Point2_Range normalize_points(const Point2_Range& range);
	Point2_Range normalize_points_T(const Point2_Range& range);
	Point3_Range normalize_points(const Point3_Range& range);
	PNC3_Range coloring_PN3(const std::vector<PN3_Range>& ranges);
	PNC3_Range coloring_PN3(const PN3_Range& pn3_range);
	PNC3_Range coloring_Point3(const std::vector<Point3_Range>& ranges);
	std::tuple<size_t, FT> nearest(Point_2 p, const Point2_Range& range);
	std::tuple<size_t, FT> nearest(Point_3 p, const Point3_Range& range);
	FT chamfer_distance(const Point2_Range& X, const Point2_Range& Y);
	Point_2 center_point(const Point2_Range& range);
	Point_3 center_point(const Point3_Range& range);
	FT convex_hull_area(const Point2_Range& range);


	// ----------- plane detection
	class Config_RegionGrowing {
	public:
		FT          search_sphere_radius;// = FT(1);
		std::size_t k;// = 12;
		FT          max_distance_to_shape;// = FT(0.5);
		FT          max_accepted_angle;// = FT(20);
		std::size_t min_region_size;// = 50;
		Config_RegionGrowing(
			FT _s,
			std::size_t _k,
			FT _md,
			FT _ma,
			std::size_t _mr
		);

	};// rg_default;
	//Config_RegionGrowing rg_default;

	//std::vector<PN3_Range> detect_planes_growing(PN3_Range& points);
	std::vector<PN3_Range> detect_planes_growing(PN3_Range& points, const Config_RegionGrowing& conf);
	std::vector<Kernel::Plane_3> extract_planes(const std::vector<PN3_Range>& ranges);
	//void project_points(const Kernel::Plane_3 plane, const Point3_Range& p3_range, Point3_Range& p3_in_plane, Point2_Range& p2_in_plane, Point2_Range& p2_in_convex);


	// ---------- line detection
	std::vector<PN2_Range> detect_lines_growing(PN2_Range& points, const Config_RegionGrowing& conf);
	std::vector<Kernel::Line_2> extract_lines(const std::vector<Point2_Range>& ranges);
	Segment2_Range extract_segments(const std::vector<Point2_Range>& ranges);


	// ------------ AlphaShape

	Segment2_Range compute_alphashape(const std::vector<Point_2>& p2);
	Mesh compute_alphashape_mesh(const std::vector<Point_2>& p2);
	void compute_alphashape_mesh(const std::vector<Point_2>& p2,
		const Plane_3 plane,
		std::vector<glm::vec3>& poses,
		std::vector<unsigned int>& indices
		);
	Mesh contour_to_mesh(const Segment2_Range& contour);
	Mesh poly_segs_2_mesh(const Segment2_Range& segments);
	void mesh_to_3d(const Plane_3& plane, Mesh& mesh);
	void extract_mesh(const Mesh& mesh, 
		std::vector<glm::vec3>& poses,
		std::vector<unsigned int>& indices
		);

}

namespace CGAL {
	template< class F >
	struct Output_rep< GEO::Color, F > {
		const GEO::Color& c;
		static const bool is_specialized = true;
		Output_rep(const GEO::Color& c) : c(c)
		{ }
		std::ostream& operator() (std::ostream& out) const
		{
			if (IO::is_ascii(out))
				out << int(c[0]) << " " << int(c[1]) << " " << int(c[2]) << " " << int(c[3]);
			else
				out.write(reinterpret_cast<const char*>(&c), sizeof(c));
			return out;
		}
	};
}



namespace ALGO {
	using namespace GEO;
	
	// ----------- Plane Clustering
	//Config_RegionGrowing rg_facade; 
	//Config_RegionGrowing rg_default; 
	//void init_configs();

	struct PlaneData {
		std::shared_ptr<Point3_Range> points_3;
		std::shared_ptr<Point2_Range> points_2;
		std::shared_ptr<Point2_Range> points_2_convex;
		Kernel::Plane_3 plane;

		FT confidence;
	};
	typedef std::vector<PlaneData> PlaneData_Range;

	PlaneData project_points(const Point3_Range& p3_range);
	PlaneData get_plane_data(const Point3_Range& p3_range);
	glm::vec3 center_distance(const PlaneData& plane_1, const PlaneData& plane_2);
	FT point_plane_distance(const glm::vec3 point, const PlaneData& plane_data);
	glm::vec3 point_plane_distance_v(const glm::vec3 point, const PlaneData& plane_data);
	//PlaneData project_points(const PN3_Range& pn3_range);
	//PlaneData project_points(const Kernel::Plane_3 plane, const Point3_Range& p3_range);
	FT plane_distance_square(const Point2_Range& p_convex1, const Point2_Range& p_convex2);
	FT plane_distance_chamfer(const Point2_Range& p_1, const Point2_Range& p_2);
	void calc_plane_confidence(PlaneData& plane_data);
	void plane_clustering(const PlaneData_Range& planes_input, std::vector<PlaneData_Range>& planes_clustered);
	void plane_clustering_fuzzy(const PlaneData_Range& planes_input, std::vector<PlaneData_Range>& planes_clustered);
	std::vector<PlaneData_Range> translate_clustering(std::vector<PlaneData_Range>& planes_clustered);
	PlaneData_Range simular_translations(const PlaneData& plane_1, const PlaneData& plane_2, const std::vector<PlaneData_Range>& planes_clustered);
	void trans_cluster_test(std::vector<PlaneData_Range>& planes_clustered);

	PNC3_Range extract_cluster_result(std::vector<PlaneData_Range>& planes_clustered);
	std::vector<Point3_Range> extract_plane_points(const PlaneData_Range& range);

	// ------- Plane Processing
	
	class PlaneProxy {
	public:
		std::shared_ptr<PlaneData> plane_data;
		glm::vec3 center;

		PlaneProxy(std::shared_ptr<PlaneData> _plane_data);
		void bbox_3d(glm::vec3& a, glm::vec3& b, glm::vec3& c, glm::vec3& d);
		glm::vec3 calc_translate(const glm::vec3& p);
		static void translation_clustering(const std::vector<std::shared_ptr<PlaneProxy>>& plane_proxies,
			std::vector<glm::vec3>& clustered_centers,
			std::vector<std::vector<glm::vec3>>& clustered_translations
			);

	protected:
	};


	void translation_clustering(const std::vector<glm::vec3>& poses,
		std::vector<glm::vec3>& clustered_centers,
		std::vector<std::vector<glm::vec3>>& clustered_translations
	);

	// --------------- Mesh Generating

	void regularize_alpha_contour(const Segment2_Range& segs_in, Segment2_Range& segs_out);
	// --- regularization methods
	PN2_Range alpha2pn(const Segment2_Range& segs_in);
	void sort_points2d(Contour& points);
	void merge_lines(std::vector<Line_2>& lines);
	void merge_lines_quad(std::vector<Line_2>& lines);
	void merge_segments(Segment2_Range& segments);
	Point2_Range lines_polygon(const std::vector<Line_2>& lines);
	Direction2_Range search_dom_dirs_2div(const Segment2_Range& segs);
	Direction2_Range search_dom_dirs_ransac(const Segment2_Range& segs);
	void reg_segs_cgal(const Segment2_Range& segs_in, Segment2_Range& segs_out);
	void reg_segs_line_fitting(const Segment2_Range& segs_in, Segment2_Range& segs_out);


	class TemplateProxy {
	public:
		std::vector<std::shared_ptr<PlaneProxy>> planes;

		//TemplateProxy() = default;
		// Q: const params
		TemplateProxy(const std::vector<std::shared_ptr<PlaneProxy>>& planes_);

	};

}



namespace FZ {
	typedef std::vector<float> Feature;
	class Fuzzy {
	public:
		Fuzzy(const std::vector<Feature>& _features);
		void set_config(float _miu);
		void do_cluster(std::vector<size_t> &cluster_indices, size_t& c_num);

		virtual ~Fuzzy();
	protected:
		void transitive_closure(); 
		//TODO
		void normalize_simularity(bool positive = true);

		float distance_L2(const Feature& f1, const Feature& f2);
		size_t find(size_t x);


		// --- Configs
		float miu = 0.7;

		// ------

		std::vector<Feature> features;
		std::vector<size_t> cluster;
		std::vector<size_t> f;
		size_t N, M;
		float** s_matrix;
		bool** b_matrix;
		bool** c_matrix;

	};
}
