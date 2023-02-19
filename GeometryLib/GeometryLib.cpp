// GeometryLib.cpp : Defines the entry point for the application.
//

#include "GeometryLib.h"
#include "IO.h"

using namespace GEO;

ALGO::Config_RegionGrowing::Config_RegionGrowing(
	FT _s,
	std::size_t _k,
	FT _md,
	FT _ma,
	std::size_t _mr
):
search_sphere_radius(_s),
k(_k),
max_distance_to_shape(_md),
max_accepted_angle(_ma),
min_region_size(_mr)
{};

//void ALGO::init_configs() {
//	rg_default = Config_RegionGrowing(
//		float(1),
//		12,
//		float(0.5),
//		float(20),
//		50
//	);
//
//	rg_facade = Config_RegionGrowing(
//		float(1),
//		12,
//		float(1.0),
//		float(30),
//		500
//	);
//
//
//}


void GEO::normal_estimate(PN3_Range& pointsv) {

	std::list<PN_3> points(pointsv.begin(), pointsv.end());

	// Estimates normals direction.
	// Note: pca_estimate_normals() requiresa range of points
	// as well as property maps to access each point's position and normal.
	const int nb_neighbors = 16; // K-nearest neighbors = 3 rings
	{
		// First compute a spacing using the K parameter
		double spacing
			= CGAL::compute_average_spacing<Concurrency_tag>
			(points, nb_neighbors,
				CGAL::parameters::point_map(PN3_Point_map()));
		//double spacing = 0.1;
		// Then, estimate normals with a fixed radius
		CGAL::pca_estimate_normals<Concurrency_tag>
			(points,
				0, // when using a neighborhood radius, K=0 means no limit on the number of neighbors returns
				CGAL::parameters::point_map(PN3_Point_map())
				.normal_map(PN3_Normal_map())
				.neighbor_radius(2. * spacing)); // use 2*spacing as neighborhood radius
	}

	// Orients normals.
	// Note: mst_orient_normals() requires a range of points
	// as well as property maps to access each point's position and normal.
	std::list<PN_3>::iterator unoriented_points_begin =
		CGAL::mst_orient_normals(points, nb_neighbors,
			CGAL::parameters::point_map(PN3_Point_map())
			.normal_map(PN3_Normal_map()));
	// Optional: delete points with an unoriented normal
	// if you plan to call a reconstruction algorithm that expects oriented normals.
	points.erase(unoriented_points_begin, points.end());

	pointsv = PN3_Range(points.begin(), points.end());
}

std::tuple<size_t, FT> GEO::nearest(Point_2 p, const Point2_Range& range) {
	FT distance = FLT_MAX;//std::numeric_limits<float>::max();
	size_t index = -1;
	for (size_t i = 0; i < range.size(); ++i) {
		glm::vec2 pvec0 = p2_to_glm(p);
		glm::vec2 pvec1 = p2_to_glm(range[i]);
		FT distance_ = glm::length((pvec0 - pvec1));
		if (distance > distance_) {
			distance = distance_;
			index = i;
		}
	}
	return std::tuple<size_t, FT>(index, distance);
}

std::tuple<size_t, FT> GEO::nearest(Point_3 p, const Point3_Range& range) { return std::tuple<size_t, FT>(1, 1); } // TODO

//std::vector<PN3_Range> GEO::detect_planes_growing(PN3_Range& points) {
std::vector<PN3_Range> GEO::detect_planes_growing(PN3_Range& points, const Config_RegionGrowing& conf) {
	//auto t1 = coloring_PN3(points);
	//IO::write_PNC3((std::string(SOLUTION_ROOT_PATH) + "/data/output/before.ply").c_str(), t1);
	normal_estimate(points);
	//auto t2 = coloring_PN3(points);
	//IO::write_PNC3((std::string(SOLUTION_ROOT_PATH) + "/data/output/after.ply").c_str(), t2);


	// Default parameter values for the data file point_set_2.xyz.
	const FT          search_sphere_radius = conf.search_sphere_radius;
	const std::size_t k = conf.k;
	const FT          max_distance_to_plane = conf.max_distance_to_shape;
	const FT          max_accepted_angle = conf.max_accepted_angle;
	const std::size_t min_region_size = conf.min_region_size;
	// Create instances of the classes Neighbor_query and Region_type.
	//Sphere_Neighbor_query_3d neighbor_query(
	//	points,
	//	search_sphere_radius);
	K_Neighbor_query_3d neighbor_query(
		points,
		k,
		PN3_Point_map());
	Region_type_3d region_type(
		points,
		max_distance_to_plane, max_accepted_angle, min_region_size);
	// Create an instance of the region growing class.
	Region_growing_3d region_growing(
		points, neighbor_query, region_type);
	// Run the algorithm.
	std::vector<std::vector<size_t>> regions;
	region_growing.detect(std::back_inserter(regions));
	// Print the number of found regions.
	std::cout << "* " << regions.size() <<
		" regions have been found"
		<< std::endl;
	srand(static_cast<unsigned int>(time(nullptr)));
	// Iterate through all regions.
	std::vector<PN3_Range> result;
	for (const auto& region : regions) {
		// Generate a random color.
		const Color color = {
				static_cast<unsigned char>(rand() % 256),
				static_cast<unsigned char>(rand() % 256),
				static_cast<unsigned char>(rand() % 256) };
		std::vector<PN_3> pns_in_plane;
		std::vector<Point_2> points2d;
		result.push_back(PN3_Range());

		for (const auto index : region) {
			const auto& key = *(points.begin() + index);
			//const Point_2& point = get(Point_map(), key);
			//pwc.push_back(std::make_pair(Point_3(point.x(), point.y(), 0), color));
			Point_3 point = std::get<0>(key);
			Vector_3 normal = std::get<1>(key);
			result.back().push_back(key);
			//points_in_plane[tot_plane_count].push_back(point);
			//pns_in_plane.push_back(key);
		}


		//CGAL::Shape_detection::

	//	Plane plane;
	//	linear_least_squares_fitting_3(points_in_plane[tot_plane_count].begin(), points_in_plane[tot_plane_count].end(), plane, CGAL::Dimension_tag<0>());
	//	
	//	Vector p_normal = plane.plane_normal();
	//	FT p_d = plane.d();
	//	if (p_normal[2] < 1e-8) p_normal = -p_normal, p_d = -p_d;
	//
	//    for (const auto index : region) {
	//		Points_PN& p = *(points.begin() + (index));
	//		Point point = std::get<0>(p);
	//
	//		//points_in_plane[tot_plane_count].push_back(point);
	//
	//		Point2 point2d = plane.to_2d(plane.projection(point));
	//		points2d.push_back(point2d);
	//		pns_in_plane.push_back(p);
	//
	//    }

	}

	return result;
}



// --- Normals needed
std::vector<PN2_Range> GEO::detect_lines_growing(PN2_Range& points, const Config_RegionGrowing& conf) {

	// Default parameter values for the data file point_set_2.xyz.
	const FT          search_sphere_radius = conf.search_sphere_radius;
	const std::size_t k = conf.k;
	const FT          max_distance_to_line = conf.max_distance_to_shape;
	const FT          max_accepted_angle = conf.max_accepted_angle;
	const std::size_t min_region_size = conf.min_region_size;
	// Create instances of the classes Neighbor_query and Region_type.
	//Sphere_Neighbor_query_2d neighbor_query(
	//	points,
	//	search_sphere_radius);
	K_Neighbor_query_2d neighbor_query(
		points,
		k,
		PN2_Point_map());
	Region_type_2d region_type(
		points,
		max_distance_to_line, max_accepted_angle, min_region_size);
	// Create an instance of the region growing class.
	Region_growing_2d region_growing(
		points, neighbor_query, region_type);
	// Run the algorithm.
	std::vector<std::vector<size_t>> regions;
	region_growing.detect(std::back_inserter(regions));
	// Print the number of found regions.
	std::cout << "* " << regions.size() <<
		" regions have been found"
		<< std::endl;
	srand(static_cast<unsigned int>(time(nullptr)));
	// Iterate through all regions.
	std::vector<PN2_Range> result;
	for (const auto& region : regions) {
		// Generate a random color.
		//const Color color = {
		//		static_cast<unsigned char>(rand() % 256),
		//		static_cast<unsigned char>(rand() % 256),
		//		static_cast<unsigned char>(rand() % 256) };
		std::vector<PN_2> pns_in_plane;
		std::vector<Point_2> points2d;
		result.push_back(PN2_Range());

		for (const auto index : region) {
			const auto& key = *(points.begin() + index);
			Point_2 point = std::get<0>(key);
			Vector_2 normal = std::get<1>(key);
			result.back().push_back(key);
		}

	}

	return result;
}


std::vector<Kernel::Line_2> GEO::extract_lines(const std::vector<Point2_Range>& ranges) {
	std::vector<Kernel::Line_2> lines(ranges.size());
	for (size_t i = 0; i < ranges.size(); ++i) {
		const auto& p2_range = ranges[i];
		CGAL::linear_least_squares_fitting_2(p2_range.begin(), p2_range.end(), lines[i], CGAL::Dimension_tag<0>());
	}
	return lines;
}


Segment2_Range GEO::extract_segments(const std::vector<Point2_Range>& ranges) {
	Segment2_Range segs;
	for (size_t i = 0; i < ranges.size(); ++i) {
		Line_2 line;
		const auto& points_of_line = ranges[i];
		FT fitting_quality;
		fitting_quality = linear_least_squares_fitting_2(points_of_line.begin(), points_of_line.end(), line, CGAL::Dimension_tag<0>());

		Point_2 seg_a, seg_b;
		FT mn = FLT_MAX, mx = -FLT_MAX;
		FT point_count = points_of_line.size();
		for (const auto point : points_of_line) {
			Point_2 point_in_line = line.projection(point);
			int axis = line.is_vertical();
			if (mn > point_in_line[axis]) {
				mn = point_in_line[axis];
				seg_a = point_in_line;
			}
			if (mx < point_in_line[axis]) {
				mx = point_in_line[axis];
				seg_b = point_in_line;
			}
		}
		FT seg_length = (seg_a - seg_b) * (seg_a - seg_b);
		FT density = point_count / seg_length;

		segs.push_back(
			Segment_2(seg_a, seg_b)
		);
	}
	return segs;
}


PNC3_Range GEO::coloring_PN3(const std::vector<PN3_Range>& ranges) {
	PNC3_Range pncs;
	for (const auto &pn3_range : ranges) {
		const Color color = {
				static_cast<unsigned char>(rand() % 256),
				static_cast<unsigned char>(rand() % 256),
				static_cast<unsigned char>(rand() % 256), 255 };
		for (const auto& pn3 : pn3_range) {
			Point_3 point = std::get<0>(pn3);
			Vector_3 normal = std::get<1>(pn3);
			pncs.emplace_back(
				point,
				normal,
				color
			);
		}
	}
	return pncs;
}


PNC3_Range GEO::coloring_PN3(const PN3_Range& pn3_range) {
	PNC3_Range pncs;
	for (const auto& pn3 : pn3_range) {
		Point_3 point = std::get<0>(pn3);
		Vector_3 normal = std::get<1>(pn3);
		normal /= normal.squared_length();
		Color color = {
			normal.x() * 128+128,
			normal.y() * 128+128,
			normal.z() * 128+128, 255
		};
		pncs.emplace_back(point, normal, color);
	}
	return pncs;
}

PNC3_Range GEO::coloring_Point3(const std::vector<Point3_Range>& ranges) {
	PNC3_Range pncs;
	for (const auto &p3_range : ranges) {
		const Color color = {
				static_cast<unsigned char>(rand() % 256),
				static_cast<unsigned char>(rand() % 256),
				static_cast<unsigned char>(rand() % 256), 255 };
		for (const auto& p3 : p3_range) {
			Point_3 point = p3;
			Vector_3 normal; // not set normal temporally
			pncs.emplace_back(
				point,
				normal,
				color
			);
		}
	}
	return pncs;
}

//std::vector<Kernel::Plane_3> GEO::extract_planes(const std::vector<PN3_Range>& ranges) {
//	std::vector<Kernel::Plane_3> planes;
//	for (size_t i = 0; i < ranges.size(); ++i) {
//		const auto& pn3_range = ranges[i];
//		CGAL::linear_least_squares_fitting_3(pn3_range.begin(), pn3_range.end(), planes[i], CGAL::Dimension_tag<0>());
//	}
//	return planes;
//}


//void GEO::project_points(const Kernel::Plane_3 plane, const Point3_Range& p3_range, Point3_Range& p3_in_plane, Point2_Range& p2_in_plane, Point2_Range& p2_in_convex) {
//	// initialize
//	p3_in_plane.clear();
//	p2_in_plane.clear();
//	p2_in_convex.clear();
//
//	for (const auto& p3 : p3_range) {
//		Point_3 p3_ = plane.projection(p3);
//		Point_2 p2_ = plane.to_2d(p3_);
//
//		p3_in_plane.push_back(p3_);
//		p2_in_plane.push_back(p2_);
//	}
//
//	std::vector<std::size_t> indices(p2_in_plane.size()), out;
//	std::iota(indices.begin(), indices.end(), 0);
//	CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
//		Convex_hull_traits_2(CGAL::make_property_map(p2_in_plane)));
//	for (int i = 0; i < out.size(); i++) p2_in_convex.push_back(p2_in_plane[out[i]]);
//}

Point2_Range GEO::normalize_points(const Point2_Range& range) {
	Point2_Range result;
	glm::vec2 center(0,0);
	float scale = -1.0;
	for (const auto& p : range) {
		glm::vec2 pvec = glm::vec2(p.x(), p.y());
		center += pvec;
	}
	center /= range.size();
	for (auto p : range) {
		glm::vec2 pvec = glm::vec2(p.x(), p.y());
		scale = max(scale, abs(pvec-center)[0]);
		scale = max(scale, abs(pvec-center)[1]);
	}

	for (auto& p : range) {
		glm::vec2 pvec = glm::vec2(p.x(), p.y());
		pvec = (pvec - center) / scale;
		result.emplace_back(pvec[0], pvec[1]);
	}

	return result;
}

Point2_Range GEO::normalize_points_T(const Point2_Range& range) {
	Point2_Range result;
	glm::vec2 center(0,0);
	for (const auto& p : range) {
		glm::vec2 pvec = glm::vec2(p.x(), p.y());
		center += pvec;
	}
	center /= range.size();

	for (auto& p : range) {
		glm::vec2 pvec = glm::vec2(p.x(), p.y());
		pvec = (pvec - center);
		result.emplace_back(pvec[0], pvec[1]);
	}
	return result;
}

Point3_Range GEO::normalize_points(const Point3_Range& range) {
	Point3_Range result;
	glm::vec3 center(0,0,0);
	FT scale = -1;
	for (const auto& p : range) {
		glm::vec3 pvec = glm::vec3(p.x(), p.y(), p.z());
		center += pvec;
	}
	center /= range.size();
	for (auto p : range) {
		glm::vec3 pvec = glm::vec3(p.x(), p.y(), p.z());
		scale = max(scale, abs(pvec-center)[0]);
		scale = max(scale, abs(pvec-center)[1]);
		scale = max(scale, abs(pvec-center)[2]);
	}
	glm::mat4 m_translate = glm::translate(glm::mat4(1.0), -center);
    glm::mat4 m_scale = glm::scale(glm::mat4(1.0), glm::vec3(1.0/scale, 1.0/scale, 1.0/scale));
	// Generally, m_scale should be multiplied on the right of m_translate
	// But the transform here is not general transfrom
	// it is [ pvec = (pvec - center) / scale; ]
	// not [ pvec = (pvec / scale) - center; ]
	glm::mat4 m_T = m_scale * m_translate;

	for (auto& p : range) {
		glm::vec3 pvec = glm::vec3(p.x(), p.y(), p.z());
		pvec = m_T * glm::vec4(pvec, 1.0);
		result.emplace_back(pvec[0], pvec[1], pvec[2]);
	}

	return result;
}


FT GEO::chamfer_distance(const Point2_Range& X, const Point2_Range& Y) {
	FT sum0 = 0;
	for (const auto& x : X) {
		auto r = nearest(x, Y);
		sum0 += std::get<1>(r);
	} sum0 /= (FT)X.size();
	FT sum1 = 0;
	for (const auto& y : Y) {
		auto r = nearest(y, X);
		sum1 += std::get<1>(r);
	} sum1 /= (FT)Y.size();
	return sum0 + sum1;
}

GEO::Segment2_Range GEO::compute_alphashape(const std::vector<Point_2>& p2) {

	AlphaShape_2 A(p2.begin(), p2.end(),
		//FT(0.001),
		FT(0.1),
		AlphaShape_2::GENERAL);

	Segment2_Range segments;

	/* get alpha edges */
	AlphaShape_2::Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(),
		end = A.alpha_shape_edges_end();
	for (; it != end; ++it)
		segments.push_back(A.segment(*it));

	return segments;
}



Mesh GEO::compute_alphashape_mesh(const std::vector<Point_2>& p2) {
	AlphaShape_2 A(p2.begin(), p2.end(),
		//FT(0.001),
		FT(0.005),
		AlphaShape_2::GENERAL);

	std::vector<Segment_2> segments;

	/* get alpha edges */
	AlphaShape_2::Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(),
		end = A.alpha_shape_edges_end();
	for (; it != end; ++it)
		segments.push_back(A.segment(*it));

	//A.classify
		/* get alpha faces */
	Mesh mesh;

	//for (Triangulation_2::face_handles::iterator it = ) {

	for (Triangulation_2::Face_handle f : A.Triangulation_2::finite_face_handles()) {

		Point_2 p2;
		Point_3 p3;
		Mesh::Vertex_index idx[3];
		for (int i = 0; i < 3; i++) {

			Triangulation_2::Vertex_handle v0 = f->vertex(i);
			p2 = v0->point();
			p3 = Point_3(p2[0], p2[1], 0);
			idx[i] = mesh.add_vertex(p3);
		}
		if (A.classify(f) == AlphaShape_2::INTERIOR)
			mesh.add_face(idx[0], idx[1], idx[2]);
	}


	//Triangulation_2::Face_iterator it = A.Triangulation_2::faces_begin(),
	//	end = A.Triangulation_2::faces_end();
	//for (; it != end; ++it) {
	//	it->
	//	Mesh::Vertex_index s = mesh.add_vertex(p);
	//	Mesh::Vertex_index u = mesh.add_vertex(q);
	//	Mesh::Vertex_index v = mesh.add_vertex(p);
	//	mesh.add_face(s, u, v);

	//}

	//std::ofstream f("D:\\Codes\\PointsProcessing\\BuildingClustering\\NOLA\\data\\classified\\1_segs\\f_" +
	//	std::to_string(alpha_id) + ".ply");
	////CGAL::IO::set_binary_mode(f); // The PLY file will be written in the binary format
	//CGAL::IO::write_PLY(
	//	f,
	//	mesh
	//);
	//}


	//write_segments(segments, "D:\\Codes\\PointsProcessing\\BuildingClustering\\NOLA\\data\\classified\\1_segs\\alpha" +
	//	std::to_string(alpha_id++) + ".ply");
	//regularize_contours(segments);

	return mesh;
}


void GEO::compute_alphashape_mesh(const std::vector<Point_2>& p2,
	const Plane_3 plane,
	std::vector<glm::vec3>& poses,
	std::vector<unsigned int>& indices
) {
	//poses.clear();
	//indices.clear();
	AlphaShape_2 A(p2.begin(), p2.end(),
		FT(0.05),
		//FT(0.005),
		AlphaShape_2::GENERAL);

	std::vector<Segment_2> segments;

	/* get alpha edges */
	AlphaShape_2::Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(),
		end = A.alpha_shape_edges_end();
	for (; it != end; ++it)
		segments.push_back(A.segment(*it));

	//A.classify
		/* get alpha faces */
	Mesh mesh;

	//for (Triangulation_2::face_handles::iterator it = ) {

	// TODO: simplify
	for (Triangulation_2::Face_handle f : A.Triangulation_2::finite_face_handles()) {

		Point_2 p2;
		Point_3 p3;
		Mesh::Vertex_index idx_[3];
		unsigned int idx[3];
		for (int i = 0; i < 3; i++) {

			Triangulation_2::Vertex_handle v0 = f->vertex(i);
			p2 = v0->point();
			p3 = plane.to_3d(p2);// Point_3(p2[0], p2[1], 0);
			poses.push_back(p3_to_glm(p3));
			idx[i] = poses.size()-1;

			idx_[i] = mesh.add_vertex(Point_3(p2[0], p2[1], 0));
		}
		if (A.classify(f) == AlphaShape_2::INTERIOR) {
			mesh.add_face(idx_[0], idx_[1], idx_[2]);
			for (auto i : idx) 
				indices.push_back(i);
		}
	}


	std::ofstream f((std::string(SOLUTION_ROOT_PATH) + "/data/output/alpha.ply").c_str());
	//CGAL::IO::set_binary_mode(f); // The PLY file will be written in the binary format
	CGAL::IO::write_PLY(
		f,
		mesh
	);


	//write_segments(segments, "D:\\Codes\\PointsProcessing\\BuildingClustering\\NOLA\\data\\classified\\1_segs\\alpha" +
	//	std::to_string(alpha_id++) + ".ply");
	//regularize_contours(segments);

	//return mesh;

}


GEO::Mesh GEO::contour_to_mesh(const Segment2_Range& contour) {
	GEO::Mesh mesh;
	// TODO

	return mesh;
}


GEO::Mesh GEO::poly_segs_2_mesh(const Segment2_Range& segments) {
	Mesh mesh;
	std::vector<Mesh::Vertex_index> vid;
	for (const auto& seg : segments) {
		Point_2 p2 = seg.source();
		vid.push_back(
			mesh.add_vertex(Point_3(p2.x(), p2.y(), 0))
		);
	}
	for (size_t i = 0, j = 1, k = 2; k < segments.size(); ++j, ++k) {
		mesh.add_face(vid[i], vid[j], vid[k]);
	}

	return mesh;
}


void GEO::mesh_to_3d(const Plane_3& plane, Mesh& mesh) {
	// TODO
	//Mesh new_mesh;
	for (auto& vid : mesh.vertices()) {
		Point_3& p = mesh.point(vid);
		p = plane.to_3d(Point_2(p.x(), p.y()));
	}

}


void GEO::extract_mesh(const Mesh& mesh,
	std::vector<glm::vec3>& poses,
	std::vector<unsigned int>& indices) {

	unsigned int idx_offset = poses.size();

	for (auto& vid : mesh.vertices()) {
		const Point_3& p = mesh.point(vid);
		poses.push_back(p3_to_glm(p));
	}

	for (const auto& f : mesh.faces()) {
		for (const auto& halfedge : mesh.halfedges_around_face(mesh.halfedge(f))) {
			auto vertex = mesh.target(halfedge);
			indices.push_back(vertex.idx() + idx_offset);
		}
	}

}




//ALGO::PlaneData ALGO::project_points(const Kernel::Plane_3 plane, const Point3_Range& p3_range) {
//	PlaneData result;
//	result.plane = plane;
//	result.points_2 = std::make_shared<Point2_Range>();
//	result.points_2_convex = std::make_shared<Point2_Range>();
//	result.points_3 = std::make_shared<Point3_Range>();
//
//	auto& p3_in_plane = *result.points_3;
//	auto& p2_in_plane = *result.points_2;
//	auto& p2_in_convex = *result.points_2_convex;
//
//	for (const auto& p3 : p3_range) {
//		Point_3 p3_ = plane.projection(p3);
//		Point_2 p2_ = plane.to_2d(p3_);
//
//		p3_in_plane.push_back(p3_);
//		p2_in_plane.push_back(p2_);
//	}
//
//	std::vector<std::size_t> indices(p2_in_plane.size()), out;
//	std::iota(indices.begin(), indices.end(), 0);
//	CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
//		Convex_hull_traits_2(CGAL::make_property_map(p2_in_plane)));
//	for (int i = 0; i < out.size(); i++) p2_in_convex.push_back(p2_in_plane[out[i]]);
//
//	return result;
//}


ALGO::PlaneData ALGO::get_plane_data(const Point3_Range& p3_range) {

	PlaneData result;
	result.points_2 = std::make_shared<Point2_Range>();
	result.points_2_convex = std::make_shared<Point2_Range>();
	result.points_3 = std::make_shared<Point3_Range>();
	CGAL::linear_least_squares_fitting_3(p3_range.begin(), p3_range.end(), result.plane, CGAL::Dimension_tag<0>());

	const auto& plane = result.plane;
	auto& p3_in_world = *result.points_3;
	auto& p2_in_plane = *result.points_2;
	auto& p2_in_convex = *result.points_2_convex;

	for (const auto& p3 : p3_range) {
		Point_3 p3_ = plane.projection(p3);
		Point_2 p2_ = plane.to_2d(p3_);

		//p3_in_plane.push_back(p3_);
		p3_in_world.push_back(p3);
		p2_in_plane.push_back(p2_);
	}

	std::vector<std::size_t> indices(p2_in_plane.size()), out;
	std::iota(indices.begin(), indices.end(), 0);
	CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
		Convex_hull_traits_2(CGAL::make_property_map(p2_in_plane)));
	for (int i = 0; i < out.size(); i++) p2_in_convex.push_back(p2_in_plane[out[i]]);

	return result;
}


ALGO::PlaneData ALGO::project_points(const Point3_Range& p3_range) {
	PlaneData result;
	result.points_2 = std::make_shared<Point2_Range>();
	result.points_2_convex = std::make_shared<Point2_Range>();
	result.points_3 = std::make_shared<Point3_Range>();
	CGAL::linear_least_squares_fitting_3(p3_range.begin(), p3_range.end(), result.plane, CGAL::Dimension_tag<0>());

	const auto& plane = result.plane;
	auto& p3_in_plane = *result.points_3;
	auto& p2_in_plane = *result.points_2;
	auto& p2_in_convex = *result.points_2_convex;

	for (const auto& p3 : p3_range) {
		Point_3 p3_ = plane.projection(p3);
		Point_2 p2_ = plane.to_2d(p3_);

		p3_in_plane.push_back(p3_);
		p2_in_plane.push_back(p2_);
	}

	std::vector<std::size_t> indices(p2_in_plane.size()), out;
	std::iota(indices.begin(), indices.end(), 0);
	CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
		Convex_hull_traits_2(CGAL::make_property_map(p2_in_plane)));
	for (int i = 0; i < out.size(); i++) p2_in_convex.push_back(p2_in_plane[out[i]]);

	return result;

}

//ALGO::PlaneData ALGO::project_points(const PN3_Range& pn3_range) {
//	PlaneData result;
//	result.points_2 = std::make_shared<Point2_Range>();
//	result.points_2_convex = std::make_shared<Point2_Range>();
//	result.points_3 = std::make_shared<Point3_Range>();
//	CGAL::linear_least_squares_fitting_3(pn3_range.begin(), pn3_range.end(), result.plane, CGAL::Dimension_tag<0>());
//
//	const auto& plane = result.plane;
//	auto& p3_in_plane = *result.points_3;
//	auto& p2_in_plane = *result.points_2;
//	auto& p2_in_convex = *result.points_2_convex;
//
//	for (const auto& pn3 : pn3_range) {
//		Point_3 p3 = std::get<0>(pn3);
//		Point_3 p3_ = plane.projection(p3);
//		Point_2 p2_ = plane.to_2d(p3_);
//
//		p3_in_plane.push_back(p3_);
//		p2_in_plane.push_back(p2_);
//	}
//
//	std::vector<std::size_t> indices(p2_in_plane.size()), out;
//	std::iota(indices.begin(), indices.end(), 0);
//	CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
//		Convex_hull_traits_2(CGAL::make_property_map(p2_in_plane)));
//	for (int i = 0; i < out.size(); i++) p2_in_convex.push_back(p2_in_plane[out[i]]);
//
//	return result;
//}


FT ALGO::plane_distance_square(const Point2_Range& p_convex1, const Point2_Range& p_convex2) {
	Point2_Range points[2];
	points[0] = normalize_points_T(p_convex1);
	points[1] = normalize_points_T(p_convex2);

	glm::vec2 left_up[2], left_down[2], right_up[2], right_down[2]; // I think, init as (0,0) should be ok?
	for (size_t i = 0; i < 2; ++i) {
		left_up[i] = { 0, 0 };
		left_down[i] = { 0, 0 };
		right_up[i] = { 0, 0 };
		right_down[i] = { 0, 0 };
	}
	for (size_t i = 0; i < 2; ++i) {
		for (size_t j = 0; j < points[i].size(); ++j) {
			left_up[i].x = min(left_up[i].x, points[i][j].x());
			left_up[i].y = min(left_up[i].y, points[i][j].y());

			left_down[i].x = min(left_down[i].x, points[i][j].x());
			left_down[i].y = max(left_down[i].y, points[i][j].y());

			right_up[i].x = max(right_up[i].x, points[i][j].x());
			right_up[i].y = min(right_up[i].y, points[i][j].y());

			right_down[i].x = max(left_down[i].x, points[i][j].x());
			right_down[i].y = max(left_down[i].y, points[i][j].y());
		}
	}

	FT distance = 0;
	distance += glm::length(left_up[0] - left_up[1]) +
		glm::length(left_down[0] - left_down[1]) +
		glm::length(right_up[0] - right_up[1]) +
		glm::length(right_down[0] - right_down[1]);
	return distance;
}


FT ALGO::plane_distance_chamfer(const Point2_Range& p_1, const Point2_Range& p_2) {
	Point2_Range points[2];
	points[0] = normalize_points_T(p_1);
	points[1] = normalize_points_T(p_2);

	return chamfer_distance(points[0], points[1]);
}


void ALGO::plane_clustering(const PlaneData_Range& planes_input, std::vector<PlaneData_Range>& planes_clustered) {
	// Args Config:
	// for wall.ply
	// FT _threshold_orient = 0.3;
	// FT _threshold_shape = 0.1;
	// for wall_54.ply
	FT _threshold_orient = 0.3;
	FT _threshold_shape = 0.05;

	// Variables Initialize
	planes_clustered.clear();
	std::vector<PlaneData> cluster_center;

	// TODO:
	// Sort planes_input by confidence firstly

	// Clustering : Orientation Filter -> Distance Filter
	for (const auto& plane_data : planes_input) {
		const auto& plane = plane_data.plane;
		glm::vec3 normal = glm::normalize(glm::vec3(plane.a(), plane.b(), plane.c()));
		size_t i = 0;
		size_t index = cluster_center.size();
		FT mindist = 999999999999;
		for (; i < cluster_center.size(); ++i) {
			const auto& plane_data_ = cluster_center[i];
			const auto& plane_ = plane_data_.plane;
			glm::vec3 normal_ = glm::normalize(glm::vec3(plane_.a(), plane_.b(), plane_.c()));

			FT distance_orient = glm::length(glm::cross(normal, normal_));
			if (distance_orient > _threshold_orient) continue;

			//FT distance_shape = plane_distance_square(*plane_data.points_2_convex, *plane_data_.points_2_convex);
			FT distance_shape = plane_distance_chamfer(*plane_data.points_2, *plane_data_.points_2);
			if (distance_shape > _threshold_shape) continue;

			if (distance_shape < mindist) {
				index = i;
				mindist = distance_shape;
			}

			printf("Yes ds %f\n", distance_shape);

			break;

		}
		if (index == cluster_center.size()) { // add a new cluster
			//printf("No %d %d\n", i, cluster_center.size());

			cluster_center.push_back(plane_data);
			planes_clustered.push_back(PlaneData_Range());
		}
		planes_clustered[index].push_back(plane_data);
	}


	// -----   remove outlier clusters
	//for (int i = planes_clustered.size() - 1; ~i; --i) {
	//	if (planes_clustered[i].size() <= 2) {
	//		swap(planes_clustered[i], planes_clustered.back());
	//		planes_clustered.pop_back();
	//	}
	//}

	std::cout << "Plane Clusters Number : " << planes_clustered.size() << std::endl;
}


void ALGO::plane_clustering_fuzzy(const PlaneData_Range& planes_input, std::vector<PlaneData_Range>& planes_clustered) {
	// Args Config:
	// for wall.ply
	// FT _threshold_orient = 0.3;
	// FT _threshold_shape = 0.1;
	// for wall_54.ply
	//FT _miu = 0.995;
	FT _miu = 0.996;
	//FT _threshold_shape = 0.05;

	// Variables Initialize
	planes_clustered.clear();
	std::vector<PlaneData> cluster_center;

	// TODO:
	// Sort planes_input by confidence firstly

	// Clustering : Orientation Filter -> Distance Filter
	std::vector<FZ::Feature> features;
	for (const auto& plane_data : planes_input) {
		features.push_back(FZ::Feature());
		const auto& convex_p2 = normalize_points_T(*plane_data.points_2_convex);
		float mx[2] = {-99999,-99999}, mn[2] = {99999, 99999};
		for (const auto& p2 : convex_p2) {
			mx[0] = max(mx[0], p2.x());
			mx[1] = max(mx[1], p2.y());
			mn[0] = min(mn[0], p2.x());
			mn[1] = min(mn[1], p2.y());
		}
		features.back().push_back(mx[0] - mn[0]);
		features.back().push_back(mx[1] - mn[1]);

	}
	std::vector<size_t> cluster_indices;
	size_t cluster_number;
	FZ::Fuzzy fuzzy(features);
	fuzzy.set_config((float)_miu);
	fuzzy.do_cluster(cluster_indices, cluster_number);
	std::cout << "Plane Clusters Number : " << cluster_number << std::endl;


	planes_clustered.assign(cluster_number, PlaneData_Range());
	for (size_t i = 0; i < planes_input.size(); ++i) {
		size_t cid = cluster_indices[i+1]-1;
		planes_clustered[cid].push_back(planes_input[i]);
	}


	// -----   remove outlier clusters
	//for (int i = planes_clustered.size() - 1; ~i; --i) {
	//	if (planes_clustered[i].size() <= 2) {
	//		swap(planes_clustered[i], planes_clustered.back());
	//		planes_clustered.pop_back();
	//	}
	//}


}



Point_2 GEO::center_point(const Point2_Range& range) {
	glm::vec2 c(0, 0);
	for (const auto& p : range) {
		c += glm::vec2(p.x(), p.y());
	} c /= range.size();
	return Point_2(c.x, c.y);
}
Point_3 GEO::center_point(const Point3_Range& range) {
	glm::vec3 c(0, 0,0);
	for (const auto& p : range) {
		c += glm::vec3(p.x(), p.y(), p.z());
	} c /= range.size();
	return Point_3(c.x, c.y, c.z);

}


glm::vec3 ALGO::center_distance(const PlaneData& plane_1, const PlaneData& plane_2) {
	Point_3 c_1 = center_point(*plane_1.points_3);
	Point_3 c_2 = center_point(*plane_2.points_3);
	return glm::vec3(c_2.x() - c_1.x(), c_2.y() - c_1.y(), c_2.z() - c_1.z());
}


FT ALGO::point_plane_distance(const glm::vec3 point, const PlaneData& plane_data) {
	FT min_dist = 9999999;
	for (const auto& p2 : *plane_data.points_2_convex) {
		auto plane = plane_data.plane;
		Point_3 p3 = plane.to_3d(p2);
		glm::vec3 p3vec = p3_to_glm(p3);
		min_dist = min(min_dist, glm::length(point - p3vec));
	}
	return min_dist;
}


// -----------------------
// for T that T · plane_1 = plane_2
// search for simular T', T' · plane_x = plane_y, (plane_x, plane_y is from the same cluster)
// push_back plane_x, plane_y to result
ALGO::PlaneData_Range ALGO::simular_translations(const PlaneData& plane_1, const PlaneData& plane_2, const std::vector<PlaneData_Range>& planes_clustered) {
	// Args Config:
	FT _threshold_distance = 0.3;

	PlaneData_Range result;
	glm::vec3 T = center_distance(plane_1, plane_2);
	result.push_back(plane_1);
	result.push_back(plane_2);

	for (const auto& cluster : planes_clustered) {
		for (size_t i = 0; i < cluster.size(); ++i) {
			for (size_t j = 0; j < cluster.size(); ++j) {
				if (i == j) continue; // T(i, j) \neq T(j, i)
				const auto& plane_x = cluster[i];
				const auto& plane_y = cluster[j];
				glm::vec3 T_ = center_distance(plane_x, plane_y);
				if (glm::length(T - T_) < _threshold_distance) {
					result.push_back(plane_x);
					result.push_back(plane_y);
				}
			}
		}
	}

	return result;
}

void ALGO::trans_cluster_test(std::vector<PlaneData_Range>& planes_clustered) {
	auto input_point3 = extract_plane_points(PlaneData_Range(planes_clustered[0].begin(), planes_clustered[0].begin() + 2));
	auto input_pnc3 = coloring_Point3(input_point3);
	IO::write_PNC3((std::string(SOLUTION_ROOT_PATH) + "/data/output/trans_cls_input.ply").c_str(), input_pnc3);

	auto output_plane_range = simular_translations(planes_clustered[0][0], planes_clustered[0][1], planes_clustered);
	auto output_point3 = extract_plane_points(output_plane_range);
	auto output_pnc3 = coloring_Point3(output_point3);
	IO::write_PNC3((std::string(SOLUTION_ROOT_PATH) + "/data/output/trans_cls_output.ply").c_str(), input_pnc3);
}

std::vector<ALGO::PlaneData_Range> ALGO::translate_clustering(std::vector<PlaneData_Range>& planes_clustered) {
	std::vector<PlaneData_Range> result;
	std::vector<std::vector<glm::vec3>> translates;
	for (size_t i = 0; i < planes_clustered.size(); ++i) {
		translates.push_back(std::vector<glm::vec3>());
		const auto& range = planes_clustered[i];
		for (size_t j = 0; j < range.size(); ++j) {
			Point_3 cj = center_point(*range[j].points_3);
			for (size_t k = 0; k < range.size(); ++k) {
				Point_3 ck = center_point(*range[k].points_3);
				translates.back().push_back(glm::vec3(ck.x() - cj.x(), ck.y() - cj.y(), ck.z() - cj.z()));
			}
		}
	}

	for (size_t i = 0; i < translates.size(); ++i) {
		result.push_back(PlaneData_Range());
		for (size_t j = i+1; j < translates.size(); ++j) {
			FT mindist = 999999999;
			size_t iii = -1, jjj = -1;
			for (size_t ii = 0; ii < translates[i].size(); ++ii) {
				size_t ia = ii / planes_clustered[i].size();
				size_t ib = ii % planes_clustered[i].size();
				if (ia == ib) continue;
				for (size_t jj = 0; jj < translates[j].size(); ++jj) {
					size_t ja = jj / planes_clustered[j].size();
					size_t jb = jj % planes_clustered[j].size();
					if (ja == jb) continue;
					if (glm::length(translates[i][ii] - translates[j][jj]) < mindist) {
						mindist = glm::length(translates[i][ii] - translates[j][jj]);
						iii = ii;
						jjj = jj;
					}
				}
			}
			if (iii == (size_t) - 1 || jjj == (size_t) - 1) continue;
			size_t ia = iii / planes_clustered[i].size();
			size_t ib = iii % planes_clustered[i].size();
			size_t ja = jjj / planes_clustered[j].size();
			size_t jb = jjj % planes_clustered[j].size();
			result.back().push_back(planes_clustered[i][ia]);
			result.back().push_back(planes_clustered[i][ib]);
			result.back().push_back(planes_clustered[j][ja]);
			result.back().push_back(planes_clustered[j][jb]);
		}
	}
	return result;
}




PNC3_Range ALGO::extract_cluster_result(std::vector<PlaneData_Range>& planes_clustered) {
	std::vector<Point3_Range> result_point3;
	for (const auto& cluster : planes_clustered) {
		result_point3.push_back(Point3_Range());
		for (const auto& plane_data: cluster) {
			for (const auto& p : *plane_data.points_3) {
				result_point3.back().push_back(p);
			}
		}
	}
	PNC3_Range colored_result = coloring_Point3(result_point3);
	return colored_result;
}


std::vector<Point3_Range> ALGO::extract_plane_points(const PlaneData_Range& range) {
	std::vector<Point3_Range> result_point3;
	for (const auto& plane_data : range) {
		result_point3.push_back(*plane_data.points_3);
	}
	return result_point3;
}


ALGO::PlaneProxy::PlaneProxy(std::shared_ptr<PlaneData> _plane_data) :
	plane_data(_plane_data)
{
	center = p3_to_glm(center_point(*(_plane_data->points_3)));
}



void ALGO::PlaneProxy::bbox_3d(glm::vec3& a, glm::vec3& b, glm::vec3& c, glm::vec3& d) {
	auto plane = plane_data->plane;
	auto points2d = *plane_data->points_2;
	glm::vec2 bbox[2];
	bbox[0].x = 999999;
	bbox[0].y = 999999;
	bbox[1].x = -999999;
	bbox[1].y = -999999;
	for (const auto& p2d : points2d) {
		bbox[0].x = min(p2d.x(), bbox[0].x);
		bbox[0].y = min(p2d.y(), bbox[0].y);

		bbox[1].x = max(p2d.x(), bbox[1].x);
		bbox[1].y = max(p2d.y(), bbox[1].y);
	}
	Point_3 a_ = plane.to_3d(
		Point_2(bbox[0].x, bbox[0].y)
	);
	Point_3 b_ = plane.to_3d(
		Point_2(bbox[1].x, bbox[0].y)
	);
	Point_3 c_ = plane.to_3d(
		Point_2(bbox[1].x, bbox[1].y)
	);
	Point_3 d_ = plane.to_3d(
		Point_2(bbox[0].x, bbox[1].y)
	);
	a = p3_to_glm(a_);
	b = p3_to_glm(b_);
	c = p3_to_glm(c_);
	d = p3_to_glm(d_);

	return;
}

glm::vec3 ALGO::PlaneProxy::calc_translate(const glm::vec3& p) {
	return center - p;
}


void ALGO::PlaneProxy::translation_clustering(const std::vector<std::shared_ptr<PlaneProxy>>& plane_proxies,
	std::vector<glm::vec3>& clustered_centers,
	std::vector<std::vector<glm::vec3>>& clustered_translations) {
	// Args Config:
	FT _threshold_distance = 0.2;

	std::vector<glm::vec3> translations;
	for (size_t i = 0; i < plane_proxies.size(); ++i) {
		for (size_t j = 0; j < plane_proxies.size(); ++j) {
			if (i == j) continue;
			translations.push_back(plane_proxies[i]->center - plane_proxies[j]->center);
		}
	}

	for (const auto& t0 : translations) {
		size_t c_id = 0;
		for (; c_id < clustered_centers.size(); ++c_id) {
			if (glm::length(t0 - clustered_centers[c_id]) < _threshold_distance) break;
		}
		if (c_id == clustered_centers.size()) clustered_centers.push_back(t0);
		clustered_translations[c_id].push_back(t0);
	}
	// TODO:
	// calcu center for every cluster in the end
}


void ALGO::translation_clustering(const std::vector<glm::vec3>& poses,
	std::vector<glm::vec3>& clustered_centers,
	std::vector<std::vector<glm::vec3>>& clustered_translations
) {
	// Args Config:
	FT _threshold_distance = 0.2;

	std::vector<glm::vec3> translations;
	for (size_t i = 0; i < poses.size(); ++i) {
		for (size_t j = 0; j < poses.size(); ++j) {
			if (i == j) continue;
			translations.push_back(poses[i] - poses[j]);
		}
	}

	for (const auto& t0 : translations) {
		size_t c_id = 0;
		for (; c_id < clustered_centers.size(); ++c_id) {
			if (glm::length(t0 - clustered_centers[c_id]) < _threshold_distance) break;
		}
		if (c_id == clustered_centers.size()) {
			clustered_centers.push_back(t0);
			clustered_translations.push_back(std::vector<glm::vec3>());
		}
		clustered_translations[c_id].push_back(t0);
	}

	std::cout << poses.size() << std::endl;
	std::cout << translations.size() << std::endl;
	std::cout << clustered_centers.size() << std::endl;
	// TODO:
	// calcu center for every cluster in the end

}



void ALGO::regularize_alpha_contour(const Segment2_Range& segs_in, Segment2_Range& segs_out) {
	IO::write_Seg2(IO::FAST_PATH("alpha_before.ply"), segs_in);

	//reg_segs_cgal(segs_in, segs_out);
	reg_segs_line_fitting(segs_in, segs_out);

	IO::write_Seg2(IO::FAST_PATH("alpha_after.ply"), segs_out);
}



void ALGO::sort_points2d(Contour& points) {
	Point_2 center = center_point(points);

	std::vector<std::pair<FT, Point_2>> angle_with_points;
	for (const auto& p : points) {
		angle_with_points.emplace_back(
			atan2(p.y() - center.y(), p.x() - center.x()),
			p
		);
	}
	std::sort(angle_with_points.begin(), angle_with_points.end());
	points.clear();
	for (const auto& awp : angle_with_points) {
		points.push_back(awp.second);
	}
}


PN2_Range ALGO::alpha2pn(const Segment2_Range& segs_in) {
	//auto segs = segs_in;
	PN2_Range result;
	Contour contour;
	for (const auto& seg : segs_in) contour.push_back(seg.source());
	sort_points2d(contour);

	size_t n = contour.size();
	for (size_t i = 0; i < n; ++i) {
		auto t = contour[(i + 1)%n];
		auto s = contour[i];
		Vector_2 dir(s, t);
		dir /= dir.squared_length(); // normalize
		dir = Vector_2(dir.y(), -dir.x()); // rotate 90 degree
		result.emplace_back(
			s,
			dir
		);
	}

	return result;
}


Direction2_Range ALGO::search_dom_dirs_2div(const Segment2_Range& segs) {
	FT div_eps = 1e-9;
	FT angle = M_PI / 2;
	FT angle_low  = angle - M_PI / 6.0;
	FT angle_high = angle + M_PI / 6.0;

	while ((angle_high - angle_low) > div_eps) {
		FT loss = 0;
		angle = (angle_high + angle_low) * .5;
		glm::vec2 dom_v[2] = {
			glm::vec2(cos(angle), sin(angle)),
			glm::vec2(sin(angle), -cos(angle))
		};

		//auto 

		for (const auto& seg : segs) {
			glm::vec2 s = p2_to_glm(seg.source());
			glm::vec2 t = p2_to_glm(seg.target());
			glm::vec2 v(t - s); v = glm::normalize(v);
			size_t id = 0;
			if (abs(glm::dot(dom_v[1], v)) >
				abs(glm::dot(dom_v[0], v)) - div_eps) {
				id = 1; // be more close to dom_v[1] than dom_v[0]
			}
			// Not necessary because sin(x) is symmetric about pi/2
			// if (glm::dot(dom_v[id], v) < div_eps) v = -v;
			FT contri = cross_2(dom_v[id], v) * glm::length(t-s);
			loss += contri;
		}

		if (loss > -div_eps) angle_high = angle;
		else angle_low = angle;
		printf(
			"loss %f, angle %f\n", loss, angle
		);
	}

	angle = (angle_high + angle_low) * .5;
	Direction2_Range ans = {
		Direction_2(cos(angle), sin(angle)),
		Direction_2(sin(angle), -cos(angle))
	};
	return ans;
}


void ALGO::reg_segs_cgal(const Segment2_Range& segs_in, Segment2_Range& segs_out) {
	// Args:
	const FT min_length_2 = FT(0.5);
	const FT  max_angle_2 = FT(10);
	const FT max_offset_2 = FT(0.5);

	Contour contour;
	for (const auto& seg : segs_in) contour.push_back(seg.source());
	sort_points2d(contour);

	//std::reverse(contour.begin(), contour.end());
	std::vector<GEO::PNC_3> points_before;
	for (size_t i = 0; i < contour.size(); ++i) {
		auto& p = contour[i];
		points_before.emplace_back(
			GEO::Point_3(p.x(), p.y(), 0),
			GEO::Vector_3(0, 0, 0),
			GEO::Color{ (unsigned char)((float)i / (float)contour.size() * 255) , 0, 0, 255}
		);
	}
	IO::write_PNC3(IO::FAST_PATH("alpha_contour_before.ply"), points_before);



	const bool is_closed = true;
	//Contour_Directions directions(
	//	contour, is_closed, CGAL::parameters::
	//	minimum_length(min_length_2).maximum_angle(max_angle_2));
	
	//std::vector<Direction_2> dirs = {
	//	Direction_2(Vector_2(0, 1)),
	//	Direction_2(Vector_2(1, 0))
	//};
	std::vector<Direction_2> dirs = search_dom_dirs_2div(segs_in);
	Contour_Directions_Custom directions(
		contour, is_closed, dirs);

	Contour contour_regularized;
	CGAL::Shape_regularization::Contours::regularize_closed_contour(
		contour, directions, std::back_inserter(contour_regularized),
		CGAL::parameters::maximum_offset(max_offset_2));

	std::cout << "* number of directions = " <<
		directions.number_of_directions() << std::endl;

	segs_out.clear();
	size_t n = contour_regularized.size();
	for (size_t i = 0; i < n; ++i) {
		segs_out.emplace_back(contour_regularized[i],
			contour_regularized[(i + 1) % n]);
	}

	std::vector<GEO::PNC_3> points;
	for (const auto& p : contour_regularized) {
		points.emplace_back(
			GEO::Point_3(p.x(), p.y(), 0),
			GEO::Vector_3(0, 0, 0),
			GEO::Color{ 255, 0, 0, 255 }
		);
	}
	IO::write_PNC3(IO::FAST_PATH("alpha_contour_after.ply"), points);
}



void ALGO::reg_segs_line_fitting(const Segment2_Range& segs_in, Segment2_Range& segs_out) {
	//Config_RegionGrowing rg_edges(0.5, 3, 0.2, 30, 5);
	Config_RegionGrowing rg_edges(0.5, 3, 0.5, 60, 5);


	auto pn2_range = alpha2pn(segs_in);
	auto line_points_ranges = detect_lines_growing(pn2_range, rg_edges);

	// --- for check lines
	PNC3_Range lines_pn3;
	for (const auto& range : line_points_ranges) {
		Color color = rand_color();
		for (const auto& pn2 : range) {
			auto p = std::get<0>(pn2);
			auto n = std::get<1>(pn2);
			lines_pn3.emplace_back(
				Point_3(p.x(), p.y(), 0),
				Vector_3(n.x(), n.y(), 0),
				color
			);
		}
	}
	IO::write_PNC3(IO::FAST_PATH("lines.ply"), lines_pn3);


	std::vector<Point2_Range> line_points_ranges0;
	Point2_Range poly_points;
	for (const auto& range : line_points_ranges) {
		line_points_ranges0.push_back(Point2_Range());
		for (const auto& pn2 : range) {
			line_points_ranges0.back().push_back(std::get<0>(pn2));
		}
	}

	if (line_points_ranges.size() == 2) {
		auto segs_inter = extract_segments(line_points_ranges0);
		poly_points.push_back(segs_inter[0].source());
		poly_points.push_back(segs_inter[0].target());
		poly_points.push_back(segs_inter[1].source());
		poly_points.push_back(segs_inter[1].target());
	}
	else {
		auto lines = extract_lines(line_points_ranges0);
		merge_lines_quad(lines);
		//merge_lines(lines);
		poly_points = lines_polygon(lines);
		sort_points2d(poly_points);
	}

	for (size_t i = 0; i < poly_points.size(); ++i) {
		segs_out.emplace_back(
			poly_points[i],
			poly_points[(i + 1) % poly_points.size()]
		);
	}

	//std::vector<Point2_Range> line_points_ranges0;
	//for (const auto& range : line_points_ranges) {
	//	line_points_ranges0.push_back(Point2_Range());
	//	for (const auto& pn2 : range) {
	//		line_points_ranges0.back().push_back(std::get<0>(pn2));
	//	}
	//}
	//auto segments = extract_segments(line_points_ranges0);
	//merge_segments(segments);

	//segs_out = segs_in;

}


void ALGO::merge_lines(std::vector<Line_2>& lines) {
	// Args:
	FT threshold_orient = 0.1;
	FT threshold_offset = 1.0;

	std::vector<std::vector<Line_2>> result;
	for (const auto& line : lines) {
		glm::vec2 p(line.a(), line.b());
		FT c = line.c();
		bool fg = false;
		for (auto& cluster : result) {
			glm::vec2 p0(0, 0);
			FT c0 = 0;
			for (const auto& line0 : cluster) {
				glm::vec2 pj(line0.a(), line0.b());
				FT cj = line0.c();
				p0 += pj;
				c0 += cj;
			}
			p0 /= cluster.size();
			c0 /= cluster.size();
			FT distance_orient = glm::length(p0 - p);
			FT distance_offset = c0 - c;
			if (distance_orient < threshold_orient + 1e-8 &&
				distance_offset < threshold_offset + 1e-8
				) {
				fg = true;
				cluster.push_back(line);
			}
		}
		if (!fg) {
			result.push_back(std::vector<Line_2>());
			result.back().push_back(line);
		}
	}
	lines.clear();
	for (const auto& cluster : result) {
		glm::vec3 p0(0, 0, 0);
		for (const auto& line0 : cluster) {
			glm::vec3 pj(line0.a(), line0.b(), line0.c());
			p0 += pj;
		}
		p0 /= cluster.size();
		lines.emplace_back(p0.x, p0.y, p0.z);
	}
}


void ALGO::merge_lines_quad(std::vector<Line_2>& lines) {
	// Args:
	FT threshold_orient = 0.15;
	FT threshold_offset = 1.0;

	std::vector<std::vector<Line_2>> result;
	for (const auto& line : lines) {
		glm::vec2 p(line.a(), line.b());
		bool fg = false;
		for (auto& cluster : result) {
			glm::vec2 p0(0, 0);
			for (const auto& line0 : cluster) {
				glm::vec2 pj(line0.a(), line0.b());
				FT cj = line0.c();
				p0 += pj;
			}
			p0 /= cluster.size();
			FT distance_orient = glm::length(p0 - p);
			if (distance_orient < threshold_orient + 1e-8) {
				fg = true;
				cluster.push_back(line);
			}
		}
		if (!fg) {
			result.push_back(std::vector<Line_2>());
			result.back().push_back(line);
		}
	}
	assert(result.size() == 2);

	lines.clear();
	for (const auto& cluster : result) {
		FT mx_c = -FLT_MAX;
		FT mn_c =  FLT_MAX;
		for (const auto& line : cluster) {
			mx_c = max(mx_c, line.c());
			mn_c = min(mn_c, line.c());
		}
		FT mid_c = (mx_c + mn_c) * .5;
		std::vector<glm::vec3> ls[2];
		for (const auto& line : cluster) {
			FT c = line.c();
			if (c < mid_c) {
				ls[0].emplace_back(line.a(), line.b(), line.c());
			}
			else {
				ls[1].emplace_back(line.a(), line.b(), line.c());
			}
		}

		for (size_t i = 0; i < 2; ++i) {
			glm::vec3 line_vector(0, 0, 0);
			for (const auto& l : ls[i]) {
				line_vector += l;
			} line_vector /= ls[i].size();
			lines.emplace_back(line_vector.x, line_vector.y, line_vector.z);
		}
	}

}



void ALGO::merge_segments(Segment2_Range& segments) {
	// Args:
	const FT min_length_2 = FT(0.5);
	const FT  max_angle_2 = FT(40);
	const FT max_offset_2 = FT(0.5);


	Contour contour;
	for (const auto& seg : segments) {
		contour.push_back(seg.source());
		contour.push_back(seg.target());
	}

	std::vector<GEO::PNC_3> points_before;
	for (size_t i = 0; i < contour.size(); ++i) {
		auto& p = contour[i];
		points_before.emplace_back(
			GEO::Point_3(p.x(), p.y(), 0),
			GEO::Vector_3(0, 0, 0),
			GEO::Color{ (unsigned char)((float)i / (float)contour.size() * 255) , 0, 0, 255}
		);
	}
	IO::write_PNC3(IO::FAST_PATH("segs_contour_before.ply"), points_before);



	const bool is_closed = true;
	Contour_Directions directions(
		contour, is_closed, CGAL::parameters::
		minimum_length(min_length_2).maximum_angle(max_angle_2));
	

	Contour contour_regularized;
	CGAL::Shape_regularization::Contours::regularize_closed_contour(
		contour, directions, std::back_inserter(contour_regularized),
		CGAL::parameters::maximum_offset(max_offset_2));

	std::cout << "* number of directions = " <<
		directions.number_of_directions() << std::endl;

	Segment2_Range segs_out;
	segs_out.clear();
	size_t n = contour_regularized.size();
	for (size_t i = 0; i < n; ++i) {
		segs_out.emplace_back(contour_regularized[i],
			contour_regularized[(i + 1) % n]);
	}

	std::vector<GEO::PNC_3> points;
	for (const auto& p : contour_regularized) {
		points.emplace_back(
			GEO::Point_3(p.x(), p.y(), 0),
			GEO::Vector_3(0, 0, 0),
			GEO::Color{ 255, 0, 0, 255 }
		);
	}
	IO::write_PNC3(IO::FAST_PATH("segs_contour_after.ply"), points);


}

Point2_Range ALGO::lines_polygon(const std::vector<Line_2>& lines) {
	// Args:
	FT threshold = 0.4;
	assert(lines.size() == 4);
	Point2_Range  quadrilateral;
	//for (size_t i = 0; i < 4; ++i) {
	for (size_t i = 0; i < 2; ++i) {
		Line_2 l0 = lines[i];
		//for (size_t j = i + 1; j < 3; ++j) {
		for (size_t j = 2; j < 4; ++j) {
			Line_2 l1 = lines[j];
			glm::vec2 p0(l0.a(), l0.b());
			glm::vec2 p1(l1.a(), l1.b());
			FT distance = glm::length(p0 - p1);
			//if (distance < threshold) continue;

			const auto result = intersection(l0, l1);
			if (result) {
				if (const Point_2* point = boost::get<Point_2>(&*result)) {
					quadrilateral.push_back(*point);
				}
			}
		}
	}
	assert(quadrilateral.size() == 4);
	sort_points2d(quadrilateral);
	return quadrilateral;
}






FZ::Fuzzy::Fuzzy(const std::vector<Feature>& _features) :
	features(_features)
{
	N = features.size();
	M = features[0].size();
	// ----------- construct simularity matrix
	s_matrix = new float* [N];
	for (size_t i = 0; i < N; ++i) s_matrix[i] = new float[N];
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			s_matrix[i][j] = distance_L2(features[i], features[j]);
		}
	}
	normalize_simularity(false);

	// ----
	b_matrix = new bool* [N];
	c_matrix = new bool* [N];
	for (size_t i = 0; i < N; ++i) {
		b_matrix[i] = new bool[N];
		c_matrix[i] = new bool[N];
	}
}

float FZ::Fuzzy::distance_L2(const Feature& f1, const Feature& f2) {
	float ans = 0;
	for (size_t i = 0; i < M; ++i) ans += (f1[i] - f2[i]) * (f1[i] - f2[i]);
	return sqrt(ans);
}

FZ::Fuzzy::~Fuzzy() {
	for (size_t i = 0; i < N; ++i) {
		delete[] s_matrix[i];
		delete[] b_matrix[i];
		delete[] c_matrix[i];
	}
	delete[] s_matrix;
	delete[] b_matrix;
	delete[] c_matrix;
}


void FZ::Fuzzy::transitive_closure() {
	// ------ construct c_matrix
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			for (size_t k = 0; k < N; ++k) {
				c_matrix[i][j] = b_matrix[i][k] && b_matrix[k][j];
			}
		}
	}
}

void FZ::Fuzzy::set_config(float _miu) { miu = _miu; }


size_t FZ::Fuzzy::find(size_t x) {
	f[x] = f[x] == x ? f[x] : find(f[x]);
	return f[x];
}


void FZ::Fuzzy::do_cluster(std::vector<size_t>& cluster_indices, size_t& c_num) {
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			b_matrix[i][j] = s_matrix[i][j] >= miu;// 0.996;// miu - 1e-9;
			//b_matrix[i][j] = s_matrix[i][j] >= 0.95;// miu - 1e-9;
		}
	}


	f.assign(N+1, 0);
	cluster.assign(N+1, 0);
	for (int i = 1; i <= N; ++i) f[i] = i;
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			if (b_matrix[i][j]) {
				f[j+1] = find(i+1);
			}
		}
	}
	c_num = 0;
	for (size_t i = 1; i <= N; ++i) {
		size_t fi = find(i);
		printf("ci%d\n", cluster[fi]);
		if (cluster[fi]) cluster[i] = cluster[fi];
		else cluster[i] = ++c_num;
	}
	cluster_indices = cluster;
}


void FZ::Fuzzy::normalize_simularity(bool positive) {
	float mx = -1;
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			mx = max(mx, s_matrix[i][j]);
		}
	}
	float factor = 1.0 / mx;
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			if (positive) {
				s_matrix[i][j] *= factor;
			}
			else {
				s_matrix[i][j] = 1.0 - (factor * s_matrix[i][j]);
			}
		}
	}
}

