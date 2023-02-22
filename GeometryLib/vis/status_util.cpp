#include "status_util.h"
#include "../IO.h"
#include <iostream>


void StatusManager::collect_facade_planes() {
	plane_nodes.clear();
	static_plane_nodes.clear();
	for (const auto& facade_node : static_facade_nodes) {
		//facade_node->extract_planes(static_plane_nodes);
		facade_node->extract_planes();
	}
}

void StatusManager::cluster() {
	ALGO::PlaneData_Range plane_datas;
	for (const auto& plane_node : plane_nodes) plane_datas.push_back(*(plane_node->plane_proxy->plane_data));
	plane_nodes.clear();
	std::vector<ALGO::PlaneData_Range> planes_clustered;
	//ALGO::plane_clustering(plane_datas, planes_clustered);
	ALGO::plane_clustering_fuzzy(plane_datas, planes_clustered);

	GEO::PNC3_Range pnc_clustered_planes = ALGO::extract_cluster_result(planes_clustered);
	IO::write_PNC3((std::string(SOLUTION_ROOT_PATH) + "/data/output/wall_planes_clustered.ply").c_str(), pnc_clustered_planes);

	unsigned int cluster_id = 1;
	cluster_size.push_back(0);
    for (const auto& cluster : planes_clustered) {
        glm::vec4 color = glm::vec4{
				MA::RandomReal(),
				MA::RandomReal(),
				MA::RandomReal(), 1.0};
        for (const auto& plane : cluster) {
            // construct PointCloud
			std::vector<glm::vec3> pvecs;
            for (const auto& p : *plane.points_3) {
                pvecs.push_back(GEO::p3_to_glm(p));
            }
            auto pcd = std::make_shared<PointCloud>(pvecs, shader_default);

            // construct PlaneProxy
			auto plane_data = std::make_shared<ALGO::PlaneData>(plane);
            auto plane_proxy = std::make_shared<ALGO::PlaneProxy>( plane_data );
			auto plane_node = std::make_shared<PlaneNode>(plane_proxy, pcd, cluster_id);

			//if (cluster.size() > 2) {
			if (plane_node->area <= 8.0) {
				plane_nodes.push_back(plane_node);
				plane_nodes.back()->raw_color = color;
			}
			else {
				static_plane_nodes.push_back(plane_node);
				static_plane_nodes.back()->raw_color = color;
			}
        }
        cluster_id++;
		cluster_size.push_back(cluster.size());
    }

	// ---- init merge set
	// for (size_t i = 0; i < cluster_id; ++i) {
	// 	m_cluster.push_back(i);
	// }
	m_cluster.resize(cluster_id);
	std::iota(m_cluster.begin(), m_cluster.end(), 0);

}


size_t StatusManager::find_mc(size_t c) {
	if (m_cluster[c] != c) m_cluster[c] = find_mc(m_cluster[c]);
	return m_cluster[c];
}


void StatusManager::select_r(glm::vec2 scr_pos) {
	if (disp_stat == planes) {
		if (glob_stat == viewing) {
			//select_node->drag(ndc);
			glob_stat = selecting;
		}
		else if (glob_stat == selecting) {
			select_planes();
			glob_stat = viewing;
			select_node->reset();
		}

	}
	else if (disp_stat == facades) {
		glm::vec3 ray_origin, ray_direction;
		ndc2worldray(*camera, scr_pos, SCR_WIDTH, SCR_HEIGHT, ray_origin, ray_direction);

		std::shared_ptr<FacadeNode> selected_facade_node = nullptr;
		float ray_t = -1;
		float tmax = 1000;
		float tmin = 0.01;
		size_t index = 0, selected_index = -1;
		for (auto& facade_node : static_facade_nodes) {
			float temp_t;
			if (facade_node->ray_hit(ray_origin, ray_direction, temp_t, tmin, tmax)) {
				tmax = temp_t;
				ray_t = temp_t;
				selected_facade_node = facade_node;
				selected_index = index;
			}
			index++;
		}
		if (selected_facade_node) { // select true
			selected_facade_node->stat = PlaneNode::selected;
			facade_nodes.push_back(selected_facade_node);
			for (const auto& plane_node : selected_facade_node->plane_nodes) {
				plane_nodes.push_back(plane_node);
			}
			swap(static_facade_nodes[selected_index], static_facade_nodes.back());
			static_facade_nodes.pop_back();
		}
	}
}




void StatusManager::select(glm::vec3 ray_origin, glm::vec3 ray_dir) {
	std::cout << "calling select" << std::endl;
	if (disp_stat == planes) {
		std::shared_ptr<PlaneNode> selected_plane_node = nullptr;
		float ray_t = -1;
		float tmax = 1000;
		float tmin = 0.01;
		for (auto& plane_node : plane_nodes) {
			float temp_t;
			if (plane_node->ray_hit(ray_origin, ray_dir, temp_t, tmin, tmax)) {
				tmax = temp_t;
				ray_t = temp_t;
				selected_plane_node = plane_node;
			}
		}
		if (selected_plane_node) { // select true
			//std::cout << "select true" << std::endl;
			ray_hit_point = ray_origin + ray_t * ray_dir;
			if (glob_stat == selecting) {
				selected_plane_node->stat = PlaneNode::selected;
				template_point = ray_origin + ray_t * ray_dir;
				//template_selected = selected_plane_node;
				template_selected->add_plane(selected_plane_node);
				glob_stat = viewing;
			}
			//else if (glob_stat == dragging) {
				//selected_plane_node->stat = PlaneNode::tied;
				//drag_point = ray_hit_point;
			template_drag = selected_plane_node;
			drag_point = template_drag->plane_proxy->center;
			gap_vector = ray_hit_point - template_point;
			//}

			unsigned int id = selected_plane_node->cluster_id;
			size_t sz = cluster_size[id];
			printf("cluster %d size: %d\n", id, sz);
		}

	}
	else if (disp_stat == facades) {
		std::shared_ptr<FacadeNode> selected_facade_node = nullptr;
		float ray_t = -1;
		float tmax = 1000;
		float tmin = 0.01;
		size_t index = 0, selected_index = -1;
		for (auto& facade_node : static_facade_nodes) {
			float temp_t;
			if (facade_node->ray_hit(ray_origin, ray_dir, temp_t, tmin, tmax)) {
				tmax = temp_t;
				ray_t = temp_t;
				selected_facade_node = facade_node;
				selected_index = index;
			}
			index++;
		}
		if (selected_facade_node) { // select true
			selected_facade_node->stat = PlaneNode::selected;
			facade_nodes.push_back(selected_facade_node);
			for (const auto& plane_node : selected_facade_node->plane_nodes) {
				plane_nodes.push_back(plane_node);
			}
			swap(static_facade_nodes[selected_index], static_facade_nodes.back());
			static_facade_nodes.pop_back();
		}
	}
}


void StatusManager::finish_selecting() {
	template_selected->setup(shader_default);
}



void StatusManager::search_same_row() {
	// for wall54.ply
	//const float line_threshold = 0.03;
	const float line_threshold = 0.5;
	// for wall.ply
	//const float line_threshold = 0.5;

	glm::vec3 p_origin = template_selected->center;
	size_t target_cluster = find_mc(template_selected->p_planes[0]->cluster_id);
	generated_poses.push_back(p_origin);
	for (auto& plane_node : plane_nodes) {

		// ---- filter:
		glm::vec3 translate = plane_node->plane_proxy->calc_translate(p_origin);
		//if (translate.y > line_threshold) continue;
		if (abs(translate.z) > abs(line_threshold)) continue;
		if (abs(translate.x) > abs(line_threshold)) continue;

		size_t s_cluster = find_mc(plane_node->cluster_id);
		if (s_cluster != target_cluster) continue;


		plane_node->stat = PlaneNode::tied;
		generated_poses.push_back(plane_node->plane_proxy->center);
	}
	std::cout << "plane_nodes.size(): " << plane_nodes.size() << std::endl;
	std::cout << "poses.size(): " << generated_poses.size() << std::endl;

}


void StatusManager::generate() {
	size_t n = generated_poses.size();
	if (n == 0) {
		std::cout << " No Poses" << std::endl;
		return;
	}
	glm::vec3 p_origin = generated_poses[0];
	for (size_t i = 1; i < n; ++i) {
		glm::vec3 p_a = -(generated_poses[i] - p_origin) + p_origin;
		generated_poses.push_back(p_a);
	}
}


void StatusManager::generate_t_cluster() {
	std::vector<glm::vec3> clustered_centers;
	std::vector<std::vector<glm::vec3>> clustered_translations;
	size_t n = generated_poses.size();
	ALGO::translation_clustering(generated_poses, clustered_centers, clustered_translations);

	size_t c_id = 0;
	size_t max_size = 0;
	for (size_t i = 0; i < clustered_translations.size(); ++i) {
		//if (max_size < clustered_translations[i].size()) {
		//	max_size = clustered_translations[i].size();
		//	c_id = i;
		//}
		if (clustered_translations.size() < 3) continue;

		glm::vec3 T = clustered_centers[i];
		for (size_t j = 0; j < n; ++j) {
			generated_poses.push_back(generated_poses[j] + T);
		}
	}

}



void StatusManager::drag() {
	size_t n = generated_poses.size();
	//glm::vec3 translate = drag_point - template_selected->plane_proxy->center;
	glm::vec3 translate = drag_point - template_selected->center;
	for (size_t i = 0; i < n; ++i) {
		glm::vec3 newpos = generated_poses[i];
		newpos.x += translate.x;
		generated_poses.push_back(newpos);
	}
}


void StatusManager::search() {
	//glm::vec3 world_point = ray_hit_point + gap_vector;
	//for (unsigned int i = 0; i < 3; ++i) {
	//	drag(world_point);
	//	world_point += gap_vector;
	//}
}


void StatusManager::set_window_info(unsigned int SCR_WIDTH_, unsigned int SCR_HEIGHT_) {
	SCR_WIDTH = SCR_WIDTH_;
	SCR_HEIGHT = SCR_HEIGHT_;
}


void StatusManager::clear() { // reset
	generated_poses.clear();
	for (auto& plane_node : plane_nodes) {
		plane_node->stat = PlaneNode::normal;
	}
	template_selected->reset();
	// TEMP
	select_node->reset();
}


void StatusManager::_draw() {
	// ---- for plane mode
	//for (auto& plane_node : plane_nodes) {
	//	plane_node->Draw();
	//}
	//for (auto& plane_node : static_plane_nodes) {
	//	plane_node->Draw();
	//}
	//draw_poses();

	// ----- for building mode
	if (disp_stat == facades) {
		for (auto& facade_node : facade_nodes) {
			facade_node->Draw();
		}
		for (auto& facade_node : static_facade_nodes) {
			facade_node->Draw();
		}
	}
	else if (disp_stat == planes) {
		for (auto& plane_node : plane_nodes) {
			plane_node->Draw();
		}
		for (auto& facade_node : static_facade_nodes) {
			facade_node->Draw();
		}
		for (auto& plane_node : static_plane_nodes) {
			plane_node->Draw();
		}
		// TEMP
		template_selected->Draw();
		draw_poses();
	}
	if (glob_stat == selecting) select_node->Draw();

}

void StatusManager::draw_poses() {
	for (const auto& pos : generated_poses) {
		//glm::mat4 m_s = glm::scale(glm::mat4(1.0), glm::vec3(2.0,2.0,2.0));
		//glm::mat4 m_model = glm::translate(glm::mat4(1.0), pos);
		//m_model = m_model * m_s;
		//test_sphere->SetModelMatrix(m_model);
		//test_sphere->Draw();
		template_selected->Draw(pos);
	}
}

void StatusManager::_update() {
	if (glob_stat == generating_1) {
		search_same_row();
		glob_stat = viewing;
	}
	else if (glob_stat == selecting) {
		select_node->drag(mouse_pos);
	}
	else if (glob_stat == generating_2) {
		//generate_t_cluster();
		generate();
		glob_stat = viewing;
	}
	else if (glob_stat == dragging) {
		//generate_t_cluster();
		drag();
		glob_stat = viewing;
	}


	glm::mat4 view_matrix = camera->GetViewMatrix();
	glm::mat4 projection_matrix = glm::perspective(glm::radians(camera->Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.01f, 500.0f);
	shader_default->use();
	shader_default->setMat4("view", view_matrix);
	shader_default->setMat4("projection", projection_matrix);
}

void StatusManager::_initialize() {
	template_selected = std::make_shared<TemplateNode>();

	//ALGO::init_configs();
	GEO::Config_RegionGrowing rg_default = GEO::Config_RegionGrowing(
		float(1),
		12,
		float(0.5),
		float(20),
		50
	);

	GEO::Config_RegionGrowing rg_facade = GEO::Config_RegionGrowing(
		float(1),
		12,
		float(0.5),
		float(90),
		500
	);



	shader_default = std::make_shared<Shader>((std::string(SOLUTION_ROOT_PATH) + "/shaders/default_vs.glsl").c_str(),
		(std::string(SOLUTION_ROOT_PATH) + "/shaders/default_fs.glsl").c_str());
	shader_ui = std::make_shared<Shader>((std::string(SOLUTION_ROOT_PATH) + "/shaders/ui_vs.glsl").c_str(),
		(std::string(SOLUTION_ROOT_PATH) + "/shaders/ui_fs.glsl").c_str());


	//   ------- Facades
	// ----- Region Growing & Coloring Result 
	GEO::PN3_Range points = IO::read_PLY((std::string(SOLUTION_ROOT_PATH) + "/data/54_C.ply").c_str());
	//GEO::PN3_Range points = IO::read_PLY((std::string(SOLUTION_ROOT_PATH) + "/data/wall_54.ply").c_str());
	//GEO::PN3_Range points = IO::read_PLY((std::string(SOLUTION_ROOT_PATH) + "/data/wall.ply").c_str());
	std::vector<GEO::PN3_Range> points_of_planes = GEO::detect_planes_growing(points, rg_facade);
	//std::vector<GEO::PN3_Range> points_of_planes = GEO::detect_planes_growing(points, rg_default);
	GEO::PNC3_Range pnc_planes = GEO::coloring_PN3(points_of_planes);
	IO::write_PNC3((std::string(SOLUTION_ROOT_PATH) + "/data/output/wall_planes.ply").c_str(), pnc_planes);
	std::cout << points.size() << std::endl;

	// ----- Plane Clustering & Coloring Result
	std::vector<std::shared_ptr<ALGO::PlaneData>> facade_plane_datas;
	for (const auto& pn3_range : points_of_planes) {
		GEO::Point3_Range p3_range = GEO::PN3Range_to_Point3Range(pn3_range);
		//facade_plane_datas.push_back(std::make_shared<ALGO::PlaneData>(ALGO::project_points(p3_range)));
		facade_plane_datas.push_back(std::make_shared<ALGO::PlaneData>(ALGO::get_plane_data(p3_range)));
	}

	for (const auto& plane_data : facade_plane_datas) {
		static_facade_nodes.push_back(std::make_shared<FacadeNode>(plane_data, shader_default));
	}

	// TEMP
	collect_facade_planes();



	//   ---------------- Planes
	//   // ----- Region Growing & Coloring Result 
	//   GEO::PN3_Range points = IO::read_PLY((std::string(SOLUTION_ROOT_PATH) + "/data/wall_54.ply").c_str());
	//   //GEO::PN3_Range points = IO::read_PLY((std::string(SOLUTION_ROOT_PATH) + "/data/wall.ply").c_str());
	//   std::vector<GEO::PN3_Range> points_of_planes = GEO::detect_planes_growing(points);
	//   GEO::PNC3_Range pnc_planes = GEO::coloring_PN3(points_of_planes);
	//   IO::write_PNC3((std::string(SOLUTION_ROOT_PATH) + "/data/output/wall_planes.ply").c_str(), pnc_planes);
	//   std::cout << points.size() << std::endl;

	//   // ----- Plane Clustering & Coloring Result
	//   ALGO::PlaneData_Range plane_datas;
	//   for (const auto& pn3_range : points_of_planes) {
	//   	GEO::Point3_Range p3_range = GEO::PN3Range_to_Point3Range(pn3_range);
	//   	plane_datas.push_back(ALGO::project_points(p3_range));
	//   }

	//   // ----- collect small planes
	//   for (const auto& plane_data : plane_datas) {
    //       glm::vec4 color = glm::vec4{
	//   			MA::RandomReal(),
	//   			MA::RandomReal(),
	//   			MA::RandomReal(), 1.0};
	//   	// construct PointCloud
	//   	std::vector<glm::vec3> pvecs;
	//   	for (const auto& p : *plane_data.points_3) {
	//   		pvecs.push_back(GEO::p3_to_glm(p));
	//   	}
	//   	auto pcd = std::make_shared<PointCloud>(pvecs, shader);

	//   	// construct PlaneProxy
	//   	// ? auto plane_data = std::make_shared<ALGO::PlaneData>(plane);
	//   	auto plane_proxy = std::make_shared<ALGO::PlaneProxy>(plane_data);
	//   	auto plane_node = std::make_shared<PlaneNode>(plane_proxy, pcd, 0);

	//   	static_plane_nodes.push_back(plane_node);
	//   	static_plane_nodes.back()->raw_color = color;
	//   }



	//// ---- do cluster ( different algos )
	//std::vector<ALGO::PlaneData_Range> planes_clustered;
	////ALGO::plane_clustering(plane_datas, planes_clustered);
	//ALGO::plane_clustering_fuzzy(plane_datas, planes_clustered);

	//GEO::PNC3_Range pnc_clustered_planes = ALGO::extract_cluster_result(planes_clustered);
	//IO::write_PNC3((std::string(SOLUTION_ROOT_PATH) + "/data/output/wall_planes_clustered.ply").c_str(), pnc_clustered_planes);

	//unsigned int cluster_id = 0;
    //for (const auto& cluster : planes_clustered) {
    //    glm::vec4 color = glm::vec4{
	//			MA::RandomReal(),
	//			MA::RandomReal(),
	//			MA::RandomReal(), 1.0};
    //    for (const auto& plane : cluster) {
    //        // construct PointCloud
	//		std::vector<glm::vec3> pvecs;
    //        for (const auto& p : *plane.points_3) {
    //            pvecs.push_back(GEO::p3_to_glm(p));
    //        }
    //        auto pcd = std::make_shared<PointCloud>(pvecs, shader);

    //        // construct PlaneProxy
	//		auto plane_data = std::make_shared<ALGO::PlaneData>(plane);
    //        auto plane_proxy = std::make_shared<ALGO::PlaneProxy>( plane_data );
	//		auto plane_node = std::make_shared<PlaneNode>(plane_proxy, pcd, cluster_id);

	//		if (cluster.size() > 2) {
	//			plane_nodes.push_back(plane_node);
	//			plane_nodes.back()->raw_color = color;
	//		}
	//		else {
	//			static_plane_nodes.push_back(plane_node);
	//			static_plane_nodes.back()->raw_color = color;
	//		}
    //    }
    //    cluster_id++;
    //}

	test_sphere = std::make_shared<Sphere>(shader_default);
	auto rect = std::make_shared<UIRectangle>(shader_ui);
	select_node = std::make_shared<SelectNode>(rect);
	//test_rect->update();

	gen_boundingbox();
}

void StatusManager::_setup(std::shared_ptr<Camera> camera_) {
	camera = camera_;
}

std::shared_ptr<Camera> StatusManager::get_camera() { return camera;  }



void StatusManager::gen_boundingbox() {
	const float MX = 999999.f;
	bbox[0] = glm::vec3(MX, MX, MX);
	bbox[1] = -glm::vec3(MX, MX, MX);
	for (const auto& plane_node : plane_nodes) {
		for (size_t i = 0; i < 2; ++i) {
			for (size_t j = 0; j < 3; ++j) {
				bbox[0].x = min(plane_node->faces[i].p[j].x, bbox[0].x);
				bbox[0].y = min(plane_node->faces[i].p[j].y, bbox[0].y);
				bbox[0].z = min(plane_node->faces[i].p[j].z, bbox[0].z);

				bbox[1].x = max(plane_node->faces[i].p[j].x, bbox[1].x);
				bbox[1].y = max(plane_node->faces[i].p[j].y, bbox[1].y);
				bbox[1].z = max(plane_node->faces[i].p[j].z, bbox[1].z);
			}
		}
	}
}



void StatusManager::ndc2worldray(const Camera& camera, const glm::vec2& scr_pos,
	const unsigned int scr_width,
	const unsigned int scr_height,
	glm::vec3& w_origin, glm::vec3& w_direction) {

		glm::vec2 ndc = scr_pos;
		ndc.x = ndc.x / scr_width - 0.5;
		ndc.y = -(ndc.y / scr_height - 0.5); // [-0.5, 0.5]

		float height = 2.0 * tan(glm::radians(camera.Zoom / 2.0));
		float width = height * scr_width / scr_height;
		w_origin = camera.Position;
		w_direction = glm::vec3(ndc.x * width, ndc.y * height, -1);
		w_direction = glm::inverse(camera.GetViewMatrix()) *
			glm::vec4(w_direction, 1.0);
		w_direction = glm::normalize(w_direction - w_origin);
}


void StatusManager::select_planes() {
	for (auto& plane_node : plane_nodes) {
		glm::vec4 pos = glm::vec4(plane_node->plane_proxy->center, 1.0);
		// MVP transform
		glm::mat4 model_matrix = plane_node->pcd->GetModelMatrix();
		glm::mat4 view_matrix = camera->GetViewMatrix();
		glm::mat4 projection_matrix = glm::perspective(glm::radians(camera->Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.01f, 500.0f);

		glm::vec4 clip = projection_matrix * view_matrix * model_matrix * pos;
		glm::vec2 ndc = glm::vec2(clip.x, clip.y) / clip.w; // TODO: consider depth occlusion

		if (select_node->inside(ndc)) { // select
			plane_node->stat = PlaneNode::selected;
			template_selected->add_plane(plane_node);
		}
	}
}
