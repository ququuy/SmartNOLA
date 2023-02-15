#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "PlaneNode.h"
#include "camera.h"
#include "../GeometryLib.h"

enum GLOB_STATUS { // action
	viewing,
	selecting,
	dragging,
	generating_1,
	generating_2
};

enum DISP_STATUS {
	facades,
	planes
};


class StatusManager {
public:

	void _update();
	void _initialize();
	void _draw();
	void _setup(std::shared_ptr<Camera> camera_);


	// facadeNode -> planeNode
	void collect_facade_planes();

	void select(glm::vec3 ray_origin, glm::vec3 ray_dir);
	void finish_selecting();
	void search_same_row();

	// generators:
	void generate();
	void generate_t_cluster();
	void cluster();
	

	void drag();
	void search();
	void set_window_info(unsigned int SCR_WIDTH_, unsigned int SCR_HEIGHT_);
	void clear();
	std::shared_ptr<Camera> get_camera();

	GLOB_STATUS glob_stat = viewing;
	DISP_STATUS disp_stat = facades;

protected:

	std::vector<size_t> m_cluster; // merge-set of cluster_id
	std::vector<size_t> cluster_size;

	// global variables
	std::shared_ptr<Sphere> test_sphere;
	std::vector<std::shared_ptr<FacadeNode>> static_facade_nodes;
	std::vector<std::shared_ptr<FacadeNode>> facade_nodes;
	std::vector<std::shared_ptr<PlaneNode>> plane_nodes;
	std::vector<std::shared_ptr<PlaneNode>> static_plane_nodes;
	std::shared_ptr<Shader> shader_default;
	std::shared_ptr<Camera> camera;

	unsigned int SCR_WIDTH, SCR_HEIGHT;


	bool select_success;
	//std::shared_ptr<PlaneNode> template_selected;
	std::shared_ptr<TemplateNode> template_selected;
	std::shared_ptr<PlaneNode> template_drag;
	glm::vec3 ray_hit_point;
	glm::vec3 template_point;
	glm::vec3 drag_point;
	glm::vec3 gap_vector;

	std::vector<glm::vec3> generated_poses;

	glm::vec3 bbox[2];

	void draw_poses();
	void gen_boundingbox();

	size_t find_mc(size_t c);

};