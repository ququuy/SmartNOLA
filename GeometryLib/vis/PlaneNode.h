#pragma once

#include "../GeometryLib.h"
#include "RenderObject.h"


class TriangleNode {
public:
	glm::vec3 p[3];
	TriangleNode(const glm::vec3& p0, const glm::vec3& p1, const glm::vec3& p2);
	TriangleNode() = default;
	bool ray_hit(glm::vec3 ray_origin, glm::vec3 ray_direction, float& ray_t, float tmin=0.01, float tmax=9999);
	
	void print();

};


class PlaneNode {
public:
	enum STAT {
		normal,
		selected,
		tied,//hide
	};
	std::shared_ptr<ALGO::PlaneProxy> plane_proxy;
	std::shared_ptr<PointCloud> pcd;
	STAT stat;
	glm::vec4 raw_color;
	TriangleNode faces[2]; // use 2 triangles to represent a square ( bbox of projected plane points )
	unsigned int cluster_id;
	float area;

	PlaneNode(
		std::shared_ptr<ALGO::PlaneProxy> plane_proxy_,
		std::shared_ptr<PointCloud> pcd_,
		unsigned int cluster_id_
	);
	PlaneNode(
		std::shared_ptr<ALGO::PlaneData> plane_data_,
		std::shared_ptr<Shader> shader_
		);

	void Draw();
	bool ray_hit(glm::vec3 ray_origin, glm::vec3 ray_direction, float& ray_t, float tmin=0.01, float tmax=9999);

protected:
	void setup_face();
};

typedef std::vector<PlaneNode> PlaneNode_Range;


class TemplateNode {
	// not used yet
public:
	TemplateNode() = default;
	std::vector<std::shared_ptr<PlaneNode>> p_planes;

	glm::vec3 center;
	std::shared_ptr<RenderGeometry> rg = nullptr;
	std::vector<unsigned int> indices;
	std::vector<Vertex> vertices;
	std::vector<glm::vec3> points;

	void add_plane(std::shared_ptr<PlaneNode> p_plane_);
	void setup(std::shared_ptr<Shader> shader_);
	void Draw();
	void Draw(glm::vec3 pos);
	void reset();
	void main_plane();
	float distance_p2m(const glm::vec3& p);

protected:
	float distance_p2tri(const glm::vec3& point,
		const glm::vec3& tp0, const glm::vec3& tp1, const glm::vec3& tp2);
};

class InstanceNode {
};




class FacadeNode : public PlaneNode {
public:
	FacadeNode(
		std::shared_ptr<ALGO::PlaneData> plane_data_,
		std::shared_ptr<Shader> shader_
		);

	// TODO
	void extract_planes(std::vector<std::shared_ptr<PlaneNode>> &planes) const;
	void extract_planes();


	std::vector<std::shared_ptr<PlaneNode>> plane_nodes;
protected:

};


class SelectNode {
public:

	std::shared_ptr<UIRectangle> rect;
	bool tied = false;

	SelectNode(std::shared_ptr<UIRectangle> rect_);

	void Draw();
	void reset();
	bool inside(glm::vec2 point);
	void drag(glm::vec2 pos);

protected:

};

