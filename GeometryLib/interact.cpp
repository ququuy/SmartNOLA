
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "IO.h"
#include "GeometryLib.h"

#include "vis/window_util.h"
#include "vis/camera.h"
#include "vis/PlaneNode.h"
#include "cmake_definition.h"

//Camera camera(glm::vec3(0.0f, 0.0f, 5.0f));

//PlaneNode_Range plane_nodes;

std::shared_ptr<Camera> camera;
std::shared_ptr<StatusManager> m_status;
std::shared_ptr<WindowManager> m_window;

// should be called after gl initialization
void init_() {
    camera = std::make_shared<Camera>(glm::vec3(0.f, 0.f, 45.f));
    m_window->_initialize();
    m_status->_initialize();
    m_status->_setup(camera);
    m_window->_setup(m_status);
}

int main() {
    m_status = std::make_shared<StatusManager>();
    m_window = WindowManager::getInstance();



    srand(time(0));

    // ------------ init data
    init_();

    //Sphere sphere(shader);


	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    while (!m_window->should_close())
    {

        // input
        // -----
        //m_window->processInput(window);

        
        //float currentFrame = glfwGetTime();
        //deltaTime = currentFrame - lastFrameTime;
        //lastFrameTime = currentFrame;
        m_window->updateFrame();

		// Let GLFW process events
		glfwPollEvents();

		// Update state
        m_status->_update();
	
		// Draw scene
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

		(glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT));
		glEnable(GL_DEPTH_TEST);
        // draw
        m_status->_draw();


		// Display results
        m_window->swap_buffers();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    //glDeleteVertexArrays(1, &VAO);
    //glDeleteBuffers(1, &VBO);
    //glDeleteProgram(shaderProgram);

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();


	return 0;
}
