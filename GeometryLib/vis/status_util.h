#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "PlaneNode.h"
#include "camera.h"
#include "../GeometryLib.h"

enum GLOB_STATUS { // action
	viewing,
	selecting_p,
	instancing,
	selecting_i,
	dragging,
	result,
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
	void select_r(glm::vec2 scr_pos);
	void finish_selecting();
	void search_same_row();

	// generators:
	void generate();
	void generate_t_cluster();
	void cluster();


	std::pair<glm::vec3, float> nearby_instance(glm::vec3 position, size_t cid);
	void back_to_facades_view();
	bool copy_area(glm::vec3 start_pos);
	void copy_copy();
	void stash_instances();
	void final_step();
	

	void drag();
	void search();
	void set_window_info(unsigned int SCR_WIDTH_, unsigned int SCR_HEIGHT_);
	void clear();
	std::shared_ptr<Camera> get_camera();

	GLOB_STATUS glob_stat = viewing;
	DISP_STATUS disp_stat = facades;
	glm::vec2 mouse_pos;

protected:

	std::vector<size_t> m_cluster; // merge-set of cluster_id
	std::vector<size_t> cluster_size;

	// global variables
	std::shared_ptr<Sphere> test_sphere;
	//std::shared_ptr<UIRectangle> test_rect;
	std::shared_ptr<SelectNode> select_node;
	std::vector<std::shared_ptr<FacadeNode>> static_facade_nodes;
	std::vector<std::shared_ptr<FacadeNode>> tied_facade_nodes;
	std::vector<std::shared_ptr<FacadeNode>> facade_nodes;
	std::vector<std::shared_ptr<PlaneNode>> plane_nodes;
	std::vector<std::shared_ptr<PlaneNode>> static_plane_nodes;
	std::shared_ptr<Shader> shader_default;
	std::shared_ptr<Shader> shader_ui;
	std::shared_ptr<Camera> camera;
	std::vector<glm::vec3> result_points;
	std::shared_ptr<PointCloud> result_pcd;

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
	std::vector<glm::vec3> area_poses; // relative poses
	glm::vec3 ref_pos, delta_pos; // relative poses


	// final result:
	std::vector<std::shared_ptr<TemplateNode>> stashed_templates;
	std::vector<std::vector<glm::vec3>> stashed_poses;


	glm::vec3 bbox[2];

	void draw_poses();
	void gen_boundingbox();
	void ndc2worldray(const Camera& camera, const glm::vec2& scr_pos,
		const unsigned int scr_width,
		const unsigned int scr_height,
		glm::vec3& w_origin, glm::vec3& w_direction);
	void select_planes();
	void select_instances();

	size_t find_mc(size_t c);

};