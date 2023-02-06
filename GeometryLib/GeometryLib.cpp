﻿// GeometryLib.cpp : Defines the entry point for the application.
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
max_distance_to_plane(_md),
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
	const FT          max_distance_to_plane = conf.max_distance_to_plane;
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
		Color color = {
			normal.x() * 255,
			normal.y() * 255,
			normal.z() * 255, 255
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
