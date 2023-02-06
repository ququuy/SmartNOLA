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
	else if (stat == STAT::tied) {
		pcd->SetColor(glm::vec4(0.0, 1.0, 0.0, 1.0));
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
}




FacadeNode::FacadeNode(
	std::shared_ptr<ALGO::PlaneData> plane_data_,
	std::shared_ptr<Shader> shader_
) : 
PlaneNode(plane_data_, shader_)
{
	//pcd->SetPointScale(0.03);
}


void FacadeNode::extract_planes(std::vector<std::shared_ptr<PlaneNode>> &planes) const {
	GEO::Config_RegionGrowing rg_default = GEO::Config_RegionGrowing(
		float(1),
		12,
		float(0.5),//float(0.5),
		float(20),//float(20),
		50
	);
	auto shader = pcd->shader;
	auto plane_data = plane_proxy->plane_data;

	GEO::Point3_Range points_3 = *plane_data->points_3;
	GEO::PN3_Range pns;
	for (const auto& p3 : points_3) {
		GEO::Point_3 point = p3;
		GEO::Vector_3 normal;
		pns.emplace_back(p3, normal);
	}
	auto points_of_planes = GEO::detect_planes_growing(pns, rg_default);

	for (const auto& pn3_range : points_of_planes) {
		GEO::Point3_Range p3_range = GEO::PN3Range_to_Point3Range(pn3_range);
		std::shared_ptr<ALGO::PlaneData> plane_data =
			std::make_shared<ALGO::PlaneData>(ALGO::project_points(p3_range));

		planes.push_back(std::make_shared<PlaneNode>(plane_data, shader));
	}

}


void FacadeNode::extract_planes() {

	GEO::Config_RegionGrowing rg_default = GEO::Config_RegionGrowing(
		float(1),
		12,
		float(0.2),//float(0.5),
		float(20),//float(20),
		50
	);
	auto shader = pcd->shader;
	auto plane_data = plane_proxy->plane_data;

	GEO::Point3_Range points_3 = *plane_data->points_3;
	GEO::PN3_Range pns;
	for (const auto& p3 : points_3) {
		GEO::Point_3 point = p3;
		GEO::Vector_3 normal;
		pns.emplace_back(p3, normal);
	}
	auto points_of_planes = GEO::detect_planes_growing(pns, rg_default);

	for (const auto& pn3_range : points_of_planes) {
		GEO::Point3_Range p3_range = GEO::PN3Range_to_Point3Range(pn3_range);
		std::shared_ptr<ALGO::PlaneData> plane_data =
			std::make_shared<ALGO::PlaneData>(ALGO::project_points(p3_range));

		plane_nodes.push_back(std::make_shared<PlaneNode>(plane_data, shader));
	}

}
