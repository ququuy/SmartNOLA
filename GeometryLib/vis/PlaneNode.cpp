#include "PlaneNode.h"
#include <iostream>

TriangleNode::TriangleNode(const glm::vec3& p0, const glm::vec3& p1, const glm::vec3& p2)
{
	p[0] = p0;
	p[1] = p1;
	p[2] = p2;
}


void TriangleNode::print() {
	std::cout << p[0].x << " " << p[0].y << " " << p[0].z << std::endl;
	std::cout << p[1].x << " " << p[1].y << " " << p[1].z << std::endl;
	std::cout << p[2].x << " " << p[2].y << " " << p[2].z << std::endl;
}


bool TriangleNode::ray_hit(glm::vec3 ray_origin, glm::vec3 ray_direction, float& ray_t, float tmin, float tmax) {
	float a = p[0][0] - p[1][0];
	float b = p[0][1] - p[1][1];
	float c = p[0][2] - p[1][2];
	float d = p[0][0] - p[2][0];
	float e = p[0][1] - p[2][1];
	float f = p[0][2] - p[2][2];
	float g = ray_direction[0];
	float h = ray_direction[1];
	float i = ray_direction[2];
	float j = p[0][0] - ray_origin[0];
	float k = p[0][1] - ray_origin[1];
	float l = p[0][2] - ray_origin[2];
	float ei_minus_hf = e * i - h * f;
	float gf_minus_di = g * f - d * i;
	float dh_minus_eg = d * h - e * g;
	float ak_minus_jb = a * k - j * b;
	float jc_minus_al = j * c - a * l;
	float bl_minus_kc = b * l - k * c;
	float M = a * ei_minus_hf + b * gf_minus_di + c * dh_minus_eg;

	// compute t
	float ray_t_ = -(f * ak_minus_jb + e * jc_minus_al + d * bl_minus_kc) / M;
	if (ray_t_ < tmin || ray_t_ > tmax)
		return false;

	// compute gamma
	float gamma = (i * ak_minus_jb + h * jc_minus_al + g * bl_minus_kc) / M;
	if (gamma < 0 || gamma > 1)
		return false;

	// compute beta
	float beta = (j * ei_minus_hf + k * gf_minus_di + l * dh_minus_eg) / M;
	if (beta < 0 || beta > 1 - gamma)
		return false;

	ray_t = ray_t_;
	return true;
}



PlaneNode::PlaneNode(
	std::shared_ptr<ALGO::PlaneProxy> plane_proxy_,
	std::shared_ptr<PointCloud> pcd_,
	unsigned int cluster_id_
) :
plane_proxy(plane_proxy_),
pcd(pcd_),
cluster_id(cluster_id_)
{
	stat = normal;
	setup_face();
}

PlaneNode::PlaneNode(
	std::shared_ptr<ALGO::PlaneData> plane_data_,
	std::shared_ptr<Shader> shader_
) {
	glm::vec4 color = glm::vec4{
			MA::RandomReal(),
			MA::RandomReal(),
			MA::RandomReal(), 1.0 };
	raw_color = color;
	// construct PointCloud
	std::vector<glm::vec3> pvecs;
	for (const auto& p : *(plane_data_->points_3)) {
		pvecs.push_back(GEO::p3_to_glm(p));
	}
	pcd = std::make_shared<PointCloud>(pvecs, shader_);

	// construct PlaneProxy
	plane_proxy = std::make_shared<ALGO::PlaneProxy>(plane_data_);

	stat = normal;
	setup_face();
	cluster_id = -1;
}


void PlaneNode::Draw() {
	if (stat == STAT::normal) {
		pcd->SetColor(raw_color);
	}
	else if (stat == STAT::selected) {
		pcd->SetColor(glm::vec4(1.0, 1.0, 0.0, 1.0));
	}
	else if (stat == STAT::tied) { //hide
		pcd->SetColor(glm::vec4(0.0, 1.0, 0.0, 1.0));
		return;
	}
	//pcd->SetColor(glm::vec4(1.0, 0.0, 0.0, 1.0));

	pcd->Draw();
	return;
}


bool PlaneNode::ray_hit(glm::vec3 ray_origin, glm::vec3 ray_direction, float& ray_t, float tmin, float tmax) {
	float result_t;
	if (faces[0].ray_hit(ray_origin, ray_direction, result_t, tmin, tmax)) ray_t = result_t;
	else if (faces[1].ray_hit(ray_origin, ray_direction, result_t, tmin, tmax)) ray_t = result_t;
	else return false;
	return true;
}



void PlaneNode::setup_face() {
	glm::vec3 p_square[4];
	plane_proxy->bbox_3d(p_square[0], p_square[1], p_square[2], p_square[3]);

	faces[0] = TriangleNode(p_square[0], p_square[1], p_square[2]);
	faces[1] = TriangleNode(p_square[0], p_square[2], p_square[3]);

	area = 0.0;
	for (auto& face : faces) {
		area += glm::length(glm::cross(face.p[1] - face.p[0], face.p[2] - face.p[0])) * .5;
	}
}




void TemplateNode::add_plane(std::shared_ptr<PlaneNode> p_plane_) {
	//center = (center * (float)p_planes.size() + p_plane_->plane_proxy->center) / (float)(p_planes.size() + 1);
	p_planes.push_back(p_plane_);
}


void TemplateNode::setup(std::shared_ptr<Shader> shader_) {
	main_plane(); // set main plane to index 0
	// ---- calcu center ( avg of all planes )
	// center = glm::vec3(0, 0, 0);
	// for (const auto& pn : p_planes) {
	// 	center += pn->plane_proxy->center;
	// }
	// center /= (float)p_planes.size();

	// ---- calcu center ( center of plane[0] )
	center = p_planes[0]->plane_proxy->center;

	// ------ alpha shape
	// TODO: visualize test

	std::vector<glm::vec3> poses;
	GEO::Point3_Range points3d;

	for (auto& pn : p_planes) {
		auto& pp = pn->plane_proxy;
		auto& pd = pp->plane_data;
		auto& plane = (pd->plane);
		auto& p2 = *(pd->points_2);
		//GEO::compute_alphashape_mesh(p2, plane, poses, indices);
		//GEO::GEN_MESH(p2, plane, poses, indices);
		
		//GEO::Mesh mesh = GEO::compute_alphashape_mesh(p2);
		//GEO::mesh_to_3d(plane, mesh);
		//GEO::extract_mesh(mesh, poses, indices);

		GEO::Segment2_Range segs = GEO::compute_alphashape(p2);
		GEO::Segment2_Range segs_regularized;
		ALGO::regularize_alpha_contour(segs, segs_regularized);

		GEO::Mesh mesh = GEO::poly_segs_2_mesh(segs_regularized);
		GEO::mesh_to_3d(plane, mesh);
		GEO::extract_mesh(mesh, poses, indices);
		if (mesh.faces().size() >= 1)
			GEO::sample_points(mesh, points3d);

		pn->stat = PlaneNode::tied;
	}
	for (const auto& pos : poses) {
		vertices.push_back({ pos });
	}


	// ----- rectangle
	//std::vector<unsigned int> indices;
	//std::vector<Vertex> vertices;

	//for (auto& pn : p_planes) {
	//	//auto& pn = p_planes[0];
	//	auto& faces = pn->faces;
	//	for (int i = 0; i < 2; ++i) {
	//		for (int j = 0; j < 3; ++j) {
	//			auto& p = faces[i].p[j];
	//			indices.push_back(vertices.size());
	//			vertices.push_back({ p });
	//		}
	//	}
	//	pn->stat = PlaneNode::tied;
	//}


	// ---- normalize mesh
	for (auto& v : vertices) {
		v.Position -= center;
	}
	for (const auto& p3d : points3d) {
		glm::vec3 p = GEO::p3_to_glm(p3d);
		points.push_back(p - center);
	}

	rg = std::make_shared<RenderGeometry>(vertices, indices, shader_);

}


void TemplateNode::Draw() {
	if (rg != nullptr) {
		glm::mat4 m_model = glm::translate(glm::mat4(1.0), center);
		rg->SetModelMatrix(m_model);
		//rg->SetModelMatrix(glm::mat4(1.0));
		rg->Draw();
	}
}

void TemplateNode::Draw(glm::vec3 pos) {
	if (rg) {
		glm::mat4 m_model = glm::translate(glm::mat4(1.0), pos);
		rg->SetModelMatrix(m_model);
		//rg->SetModelMatrix(glm::mat4(1.0));
		rg->Draw();
	}

}

void TemplateNode::reset() {
	p_planes.clear();
	indices.clear();
	vertices.clear();
	points.clear();
	center = glm::vec3(0, 0, 0);
	rg = nullptr;
}



void TemplateNode::main_plane() {
	size_t m_i = 0;
	float max_area = -FLT_MAX;
	for (size_t i = 0; i < p_planes.size(); ++i) {
		float area = GEO::convex_hull_area(*(p_planes[i]->plane_proxy->plane_data->points_2_convex));
		if (area > max_area) {
			max_area = area;
			m_i = i;
		}
	}
	swap(p_planes[m_i], p_planes[0]);
}



float TemplateNode::distance_p2tri(const glm::vec3& point,
	const glm::vec3& tp0, const glm::vec3& tp1, const glm::vec3& tp2) {
	auto p = GEO::glm_to_p3(point); 
	auto p0 = GEO::glm_to_p3(tp0);
	auto p1 = GEO::glm_to_p3(tp1);
	auto p2 = GEO::glm_to_p3(tp2);
	GEO::Triangle_3 triangle(p0, p1, p2);

	return CGAL::squared_distance(p, triangle);
}


float TemplateNode::distance_p2m(const glm::vec3& p) {
	float distance = FLT_MAX;
	for (size_t i = 0; i < indices.size(); i += 3) {
		float di = distance_p2tri(p,
			vertices[indices[i]].Position,
			vertices[indices[i+1]].Position,
			vertices[indices[i+2]].Position
			);
		distance = std::min(di, distance);
	}
	return distance;
}





FacadeNode::FacadeNode(
	std::shared_ptr<ALGO::PlaneData> plane_data_,
	std::shared_ptr<Shader> shader_
) : 
PlaneNode(plane_data_, shader_)
{
	//pcd->SetPointScale(0.03);
}


//void FacadeNode::extract_planes(std::vector<std::shared_ptr<PlaneNode>> &planes) const {
//	GEO::Config_RegionGrowing rg_default = GEO::Config_RegionGrowing(
//		float(1),
//		12,
//		//5,
//		float(0.1),//float(0.5),
//		float(10),//float(20),
//		//float(0.001),//float(0.5),
//		//float(5),//float(20),
//		//50
//		20
//	);
//	auto shader = pcd->shader;
//	auto plane_data = plane_proxy->plane_data;
//
//	GEO::Point3_Range points_3 = *plane_data->points_3;
//	GEO::PN3_Range pns;
//	for (const auto& p3 : points_3) {
//		GEO::Point_3 point = p3;
//		GEO::Vector_3 normal;
//		pns.emplace_back(p3, normal);
//	}
//
//	// normal estimateh
//	//auto t1 = GEO::coloring_PN3(pns);
//	//IO::write_PNC3((std::string(SOLUTION_ROOT_PATH) + "/data/output/before.ply").c_str(), t1);
//	GEO::normal_estimate(pns);
//	//auto t2 = GEO::coloring_PN3(pns);
//	//IO::write_PNC3((std::string(SOLUTION_ROOT_PATH) + "/data/output/after.ply").c_str(), t2);
//	auto points_of_planes = GEO::detect_planes_growing(pns, rg_default);
//
//	for (const auto& pn3_range : points_of_planes) {
//		GEO::Point3_Range p3_range = GEO::PN3Range_to_Point3Range(pn3_range);
//		std::shared_ptr<ALGO::PlaneData> plane_data =
//			std::make_shared<ALGO::PlaneData>(ALGO::project_points(p3_range));
//
//		planes.push_back(std::make_shared<PlaneNode>(plane_data, shader));
//	}
//
//}


void FacadeNode::extract_planes() {
	// Small b28
	//GEO::Config_RegionGrowing rg_default = GEO::Config_RegionGrowing(
	//	float(1),
	//	8,
	//	float(0.01),//float(0.5),
	//	float(50),//float(20),
	//	30
	//);
	//b56
	//GEO::Config_RegionGrowing rg_default = GEO::Config_RegionGrowing(
	//	float(1),
	//	8,
	//	float(0.08),//float(0.5),
	//	float(30),//float(20),
	//	30
	//);

	//dc1
	GEO::Config_RegionGrowing rg_default = GEO::Config_RegionGrowing(
		float(1),
		12,
		float(0.15),//float(0.2),
		//float(0.15),//float(0.2),
		float(40),//float(20),
		50
	);
	// Middle f..
	//GEO::Config_RegionGrowing rg_default = GEO::Config_RegionGrowing(
	//	float(1),
	//	12,
	//	float(0.2),//float(0.5),
	//	float(40),//float(20),
	//	50
	//);
	// large (C54
	//GEO::Config_RegionGrowing rg_default = GEO::Config_RegionGrowing(
	//	float(1),
	//	12,
	//	float(0.8),//float(0.5),
	//	float(30),//float(20),
	//	//float(40),//float(20),
	//	50
	//);
	auto shader = pcd->shader;
	auto plane_data = plane_proxy->plane_data;

	GEO::Point3_Range points_3 = *plane_data->points_3;
	GEO::PN3_Range pns;
	for (const auto& p3 : points_3) {
		GEO::Point_3 point = p3;
		GEO::Vector_3 normal;
		pns.emplace_back(p3, normal);
	}
	GEO::normal_estimate(pns);
	auto points_of_planes = GEO::detect_planes_growing(pns, rg_default);

	for (const auto& pn3_range : points_of_planes) {
		GEO::Point3_Range p3_range = GEO::PN3Range_to_Point3Range(pn3_range);
		std::shared_ptr<ALGO::PlaneData> plane_data =
			//std::make_shared<ALGO::PlaneData>(ALGO::project_points(p3_range));
			std::make_shared<ALGO::PlaneData>(ALGO::get_plane_data(p3_range));

		plane_nodes.push_back(std::make_shared<PlaneNode>(plane_data, shader));
	}

}



SelectNode::SelectNode(std::shared_ptr<UIRectangle> rect_) :
rect(rect_)
{}


void SelectNode::Draw() {
	rect->update();
	rect->Draw();
}


void SelectNode::reset() {
	tied = false;
	rect->A = glm::vec2(0, 0);
	rect->B = glm::vec2(0, 0);
}


bool SelectNode::inside(glm::vec2 point) {
	return point.x >= std::min(rect->A.x, rect->B.x) &&
		point.x <= std::max(rect->A.x, rect->B.x) &&
		point.y >= std::min(rect->A.y, rect->B.y) &&
		point.y <= std::max(rect->A.y, rect->B.y);
}


void SelectNode::drag(glm::vec2 pos) {
	if (tied) rect->B = pos;
	else {
		rect->A = pos; tied = true;
	}
}
