#pragma once
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include "Camera.h"
#include "status_util.h"
#include "IO.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"


class WindowManager {
protected:
	float deltaTime;
	float lastFrameTime;

	bool firstMouse;
	float lastX;
	float lastY;

	ImFont* font_tahoma;
	ImFont* font_sans;

	std::shared_ptr<Camera> camera;
	GLFWwindow* window;

	static void destroy(WindowManager*) {}


public:
	std::shared_ptr<StatusManager> m_status;

	//unsigned int SCR_WIDTH = 1024;
	//unsigned int SCR_HEIGHT = 768;
	unsigned int SCR_WIDTH = 1920;
	unsigned int SCR_HEIGHT = 1080;
	glm::vec2 mouse_ndc;

	static std::shared_ptr<WindowManager> instance_m_window;

	static std::shared_ptr<WindowManager> getInstance();

	static void mouse_pos_callback(GLFWwindow* w, double x, double y);

	static void framebuffer_size_callback(GLFWwindow* window, int width, int height);

	static void mouse_click_callback(GLFWwindow* window, int button, int action, int mods);

	static void key_input_callback(GLFWwindow* window, int key, int scancode, int action, int mods);



	//auto a = framebuffer_size_callback;
	//GLFWframebuffersizefun func_framebuffersizefun = framebuffer_size_callback;//std::function<GLFWframebuffersizefun>(framebuffer_size_callback).target<GLFWframebuffersizefun>();

	//std::function<GLFWframebuffersizefun> func_framebuffersizefun = &framebuffer_size_callback;
	//std::function<void(GLFWwindow*,int,int)> func_framebuffersizefun = framebuffer_size_callback;
	//std::function<GLFWcursorposfun> func_cursorpos = &mouse_callback;
	//std::function<GLFWmousebuttonfun> func_mousebutton = &mouse_click_callback;
	//std::function<GLFWframebuffersizefun> func_framebuffersizefun = framebuffer_size_callback;

	void draw_gui();

	void call_cluster();
	void call_change_display_mode();


	bool should_close() {
		return glfwWindowShouldClose(window);
	}

	void swap_buffers() {

		glfwSwapBuffers( window );
	}

	void _setup(std::shared_ptr<StatusManager> m_status_) {
		m_status = m_status_;
		camera = m_status->get_camera();
		firstMouse = true;
		m_status->set_window_info(SCR_WIDTH, SCR_HEIGHT);
	}

		//typedef std::remove_pointer<GLFWframebuffersizefun>::type ftype_framebuffer;
		//typedef std::remove_pointer<GLFWcursorposfun>::type ftype_cursorpos;
		//typedef std::remove_pointer<GLFWmousebuttonfun>::type ftype_mousebutton;

		//std::function<ftype_framebuffer> fun_framebuffer;
		//std::function<ftype_cursorpos> fun_mousepos;
		//std::function<ftype_mousebutton> fun_mouseclick;
	void _initialize() {
	// glfw: initialize and configure
	// ------------------------------
		glfwInit();
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

		// glfw window creation
		// --------------------
		//GLFWwindow* window = glfwCreateWindow(m_window->SCR_WIDTH, m_window->SCR_HEIGHT, "rep det", NULL, NULL);
		window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "BuildingPoints Consolidation", NULL, NULL);
		if (window == NULL)
		{
			std::cout << "Failed to create GLFW window" << std::endl;
			glfwTerminate();
			return;
			//return -1;
		}
		glfwMakeContextCurrent(window);
		//glfwSetFramebufferSizeCallback(window,m_window->func);
		//glfwSetCursorPosCallback(window,*(m_window->func_cursorpos.target<GLFWcursorposfun>()));
		//glfwSetMouseButtonCallback(window,*(m_window->func_mousebutton.target<GLFWmousebuttonfun>()));

		//glfwSetScrollCallback(window, scroll_callback);



		//std::function<void(GLFWwindow*, int, int)> fun_framebuffer1;
		//typedef ftype_framebuffer = std::remove_pointer<GLFWframebuffersizefun>::type ;
		//fun_framebuffer = std::bind(&WindowManager::framebuffer_size_callback,
		//	this,
		//	std::placeholders::_1,
		//	std::placeholders::_2,
		//	std::placeholders::_3
		//);
		//fun_mousepos = std::bind(&WindowManager::mouse_callback,
		//	this,
		//	std::placeholders::_1,
		//	std::placeholders::_2,
		//	std::placeholders::_3
		//);
		//fun_mouseclick = std::bind(&WindowManager::mouse_click_callback,
		//	this,
		//	std::placeholders::_1,
		//	std::placeholders::_2,
		//	std::placeholders::_3,
		//	std::placeholders::_4
		//);
		//auto func = [](GLFWwindow* w, double x, double y) { fun_mousepos(w, x, y); }
		glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
		//glfwSetFramebufferSizeCallback(window, *(fun_framebuffer.target<ftype_framebuffer>()));
		//glfwSetCursorPosCallback(window, *(fun_mousepos.target<ftype_cursorpos>()));
		glfwSetCursorPosCallback(window, mouse_pos_callback);
		//glfwSetMouseButtonCallback(window, *(fun_mouseclick.target<ftype_mousebutton>()));
		glfwSetMouseButtonCallback(window, mouse_click_callback);
		glfwSetKeyCallback(window, key_input_callback);
		std::cout << "glfw window initialization finished." << std::endl;
		// glad: load all OpenGL function pointers
		// ---------------------------------------
		if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
		{
			std::cout << "Failed to initialize GLAD" << std::endl;
			return;
		}


		// init imgui
		// Setup Dear ImGui context
		IMGUI_CHECKVERSION();
		ImGui::CreateContext();
		ImGuiIO& io = ImGui::GetIO(); (void)io;
		//io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
		//io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

		// Setup Dear ImGui style
		//ImGui::StyleColorsDark();
		//ImGui::StyleColorsClassic();
		ImGui::StyleColorsLight();
		font_tahoma = ImGui::GetIO().Fonts->AddFontFromFileTTF(IO::ASSETS_PATH("tahoma.ttf"), 30.0f, nullptr);
		font_sans = ImGui::GetIO().Fonts->AddFontFromFileTTF(IO::ASSETS_PATH("sansserif.otf"), 24.0f, nullptr);

		// Setup Platform/Renderer backends
		ImGui_ImplGlfw_InitForOpenGL(window, true);
		ImGui_ImplOpenGL3_Init("#version 130");
	}

	//WindowManager() = default;

	void updateFrame() {
		processInput(window);
		float currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrameTime;
		lastFrameTime = currentFrame;
	}


	// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
	// ---------------------------------------------------------------------------------------------------------
	void processInput(GLFWwindow* window)
	{
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
			glfwSetWindowShouldClose(window, true);

		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
			camera->ProcessKeyboard(UP, deltaTime);
		if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
			camera->ProcessKeyboard(DOWN, deltaTime);
		if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
			camera->ProcessKeyboard(LEFT, deltaTime);
		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
			camera->ProcessKeyboard(RIGHT, deltaTime);
		if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
			camera->ProcessKeyboard(FORWARD, deltaTime);
		if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
			camera->ProcessKeyboard(BACKWARD, deltaTime);



	}
	void test() {}

	// glfw: whenever the window size changed (by OS or user resize) this callback function executes
	// ---------------------------------------------------------------------------------------------
	void framebuffer_size_impl(GLFWwindow* window, int width, int height) {
		// make sure the viewport matches the new window dimensions; note that width and 
		// height will be significantly larger than specified on retina displays.
		glViewport(0, 0, width, height);
	}


	void mouse_pos_impl(GLFWwindow* window, double xpos, double ypos) {
		m_status->mouse_pos = glm::vec2(xpos/SCR_WIDTH * 2 - 1, -(ypos/SCR_HEIGHT * 2 - 1));
		if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE)
		{
			lastX = xpos;
			lastY = ypos;
			return;
		}
		if (!firstMouse)
		{
			lastX = xpos;
			lastY = ypos;
			firstMouse = true;
		}

		float xoffset = xpos - lastX;
		float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

		lastX = xpos;
		lastY = ypos;

		camera->ProcessMouseMovement(xoffset, yoffset);


	}

	void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {}

	void mouse_click_impl(GLFWwindow* window, int button, int action, int mods) {
		if (
			button == GLFW_MOUSE_BUTTON_RIGHT &&
			//action == GLFW_PRESS)
			action == GLFW_RELEASE)
		{
			//STATUS_SELECTING = true;

			double xpos, ypos;
			glfwGetCursorPos(window, &xpos, &ypos);
			m_status->select_r(glm::vec2(xpos, ypos));
			//mouse_ndc.x = (xpos / SCR_WIDTH) * 2 - 1.0;
			//mouse_ndc.y = (ypos / SCR_HEIGHT) * 2 - 1.0;
			//mouse_ndc.y = -mouse_ndc.y;

			//glm::vec3 w_origin, w_direction;
			//ndc2worldray(*camera, mouse_ndc, w_origin, w_direction);
			//m_status->select(w_origin, w_direction);

			return;
		}

	}

	void key_input_impl(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		if (key == GLFW_KEY_O && action == GLFW_RELEASE);
			
			//m_status->glob_stat = selecting;
		if (key == GLFW_KEY_P && action == GLFW_RELEASE)
			m_status->finish_selecting();
		if (key == GLFW_KEY_R && action == GLFW_RELEASE)
			m_status->clear();
		if (key == GLFW_KEY_I && action == GLFW_RELEASE)
			m_status->search_same_row();
		if (key == GLFW_KEY_U && action == GLFW_RELEASE)
			m_status->copy_copy();
		if (key == GLFW_KEY_F && action == GLFW_RELEASE)
			m_status->final_step();
		if (key == GLFW_KEY_T && action == GLFW_RELEASE)
			m_status->stash_instances();
		//if (key == GLFW_KEY_U && action == GLFW_RELEASE);
		//	m_status->glob_stat = generating_2;
		//if (key == GLFW_KEY_Y && action == GLFW_RELEASE) {
		//	m_status->glob_stat = dragging;
		//}

		if (key == GLFW_KEY_C && action == GLFW_RELEASE) {
			call_cluster();
			//m_status->cluster();
		}

		if (key == GLFW_KEY_B && action == GLFW_RELEASE) {
			call_change_display_mode();
		}


		if (key == GLFW_KEY_Z && action == GLFW_RELEASE) {
			camera->rotateAroundY(0.5);
		}
		if (key == GLFW_KEY_X && action == GLFW_RELEASE) {
			camera->rotateAroundY(-0.5);
		}

	}

	void ndc2worldray(const Camera& camera, glm::vec2 ndc, glm::vec3& w_origin, glm::vec3& w_direction) {
		ndc *= 0.5f; // [-1, 1] to [-.5, .5]

		float height = 2.0 * tan(glm::radians(camera.Zoom / 2.0));
		float width = height * SCR_WIDTH / SCR_HEIGHT;
		w_origin = camera.Position;
		w_direction = glm::vec3(ndc.x * width, ndc.y * height, -1);
		w_direction = glm::inverse(camera.GetViewMatrix()) *
			glm::vec4(w_direction, 1.0);
		w_direction = glm::normalize(w_direction - w_origin);
	}

};


void WindowManager::call_cluster() {
	m_status->cluster();
}
void WindowManager::call_change_display_mode() {
	if (m_status->disp_stat == facades) m_status->disp_stat = planes;
	else if (m_status->disp_stat == planes) {
		m_status->back_to_facades_view();
		m_status->disp_stat = facades;
	}
}

void WindowManager::draw_gui() {
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();
	//ImGui::SetNextWindowPos(ImVec2(10, 600), ImGuiCond_FirstUseEver);
	//ImGui::SetNextWindowSize(ImVec2(600, 600), ImGuiCond_FirstUseEver);
//ImGui::SetNextWindowPos(ImVec2(10, 600));
//ImGui::SetNextWindowSize(ImVec2(600, 600));

// 2. Show a simple window that we create ourselves. We use a Begin/End pair to create a named window.
	{
		ImGui::PushFont(font_tahoma);
		ImGui::SetNextWindowPos(ImVec2(0, 0));
		ImGui::SetNextWindowSize(ImVec2(1920, 65));
		//ImGui::SetWindowFontScale(1.8);
		static float f = 0.0f;
		//static int counter = 0;

		ImGui::Begin("Inspector", NULL, ImGuiWindowFlags_NoTitleBar);                          // Create a window called "Hello, world!" and append into it.
		//ImGui::SetWindowFontScale(2.5);

		//ImGui::Text("Try to select some repeated instances.");               // Display some text (you can use a format strings too)
		ImGui::SameLine();
		ImGui::SetCursorPosX(20);
		if (ImGui::Button("Change Mode", ImVec2(180,50))) {
			call_change_display_mode();
		}
		ImGui::SameLine();
		ImGui::SetCursorPosX(220);
		if (ImGui::Button("Cluster", ImVec2(180, 50))) {
			call_cluster();
		}
		ImGui::SameLine();
		ImGui::SetCursorPosX(420);
		if (ImGui::Button("Calc Shape", ImVec2(180,50))) {
			//call_change_display_mode();
			m_status->finish_selecting();
		}
		ImGui::SameLine();
		ImGui::SetCursorPosX(620);
		if (ImGui::Button("Fit Cluster", ImVec2(180,50))) {
			//call_change_display_mode();
			m_status->search_same_row();
		}
		ImGui::SameLine();
		ImGui::SetCursorPosX(820);
		if (ImGui::Button("Repeat", ImVec2(180,50))) {
			//call_change_display_mode();
			//m_status->search_same_row();
			m_status->copy_copy();
		}
		//ImGui::SameLine();
		//ImGui::Text(" ");            
		//ImGui::SameLine();
		ImGui::SameLine();
		ImGui::SetCursorPosX(1020);
		if (ImGui::Button("Generate", ImVec2(180,50))) {
			//call_change_display_mode();
			//m_status->search_same_row();
			m_status->clear();
		}
		ImGui::SameLine();
		ImGui::SetCursorPosX(1220);
		if (ImGui::Button("Delete", ImVec2(180,50))) {
		}
		ImGui::SameLine();
		ImGui::SetCursorPosX(1420);
		if (ImGui::Button("Reset", ImVec2(180,50))) {
		}

		//ImGui::Checkbox("Demo Window", &show_demo_window);      // Edit bools storing our window open/close state
		//ImGui::Checkbox("Another Window", &show_another_window);

		//ImGui::SliderFloat("Point Scale", &f, 0.0f, 1.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
		//ImGui::ColorEdit3("Clear Color", (float*)&clear_color); // Edit 3 floats representing a color
		//ImGui::SliderFloat("Plane Cluster Radius", &_plane_cluster_radius, 0.001, 0.03);

		ImGui::SameLine();
		ImGui::SetCursorPosX(1750);
		ImGui::Text("(%.1f FPS)", ImGui::GetIO().Framerate);

		ImGui::End();
	}
	{
		ImGui::SetNextWindowPos(ImVec2(0, 970));
		ImGui::SetNextWindowSize(ImVec2(800, 110));
		ImGui::Begin("Description");
		ImGui::PopFont();
		ImGui::PushFont(font_sans);
		ImGui::TextWrapped("Try to select some repeated instances.");
		ImGui::TextWrapped("Right Click to select a structure or grouping.");

		ImGui::End();
		ImGui::PopFont();
	}




	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}



//WindowManager::instance_m_window = nullptr;
std::shared_ptr<WindowManager> WindowManager::instance_m_window(
	new WindowManager(),
	WindowManager::destroy
);

std::shared_ptr<WindowManager> WindowManager::getInstance() {
	if (instance_m_window == nullptr) {
		instance_m_window = std::make_shared<WindowManager>();
		//instance_m_window->_initialize();
	}
	return instance_m_window;
}

void WindowManager::mouse_pos_callback(GLFWwindow* w, double x, double y) {
	getInstance()->mouse_pos_impl(w, x, y);
}

void WindowManager::framebuffer_size_callback(GLFWwindow* window, int width, int height) {
	getInstance()->framebuffer_size_impl(window, width, height);
}

void WindowManager::mouse_click_callback(GLFWwindow* window, int button, int action, int mods) {
	getInstance()->mouse_click_impl(window, button, action, mods);
}

void WindowManager::key_input_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	getInstance()->key_input_impl(window, key, scancode, action, mods);

}

