#pragma once
#include <vector>
#include <string>
#include <utility>
#include <learnopengl/shader_s.h>
#include <glm/glm.hpp>

#include "../GeometryLib.h"

struct Vertex {
	glm::vec3 Position;
	glm::vec3 Normal;
	//Vec2f TexCoords;
};

class RenderObject {
public:
    /* Mesh Data */
    const std::vector<Vertex> vertices;
    const std::vector<unsigned int> indices;
    //std::vector<Vec3f> normal;
    unsigned int VAO, VBO, EBO;

    /* Functions */
    //MyMesh(vector<Vertex> _vertices, vector<unsigned int> _indices, vector<Texture> _textures);
    //void readMeshFile(string filename);
    RenderObject(const std::vector<Vertex> &vertices_, const std::vector<unsigned int> &indices_);
    RenderObject() {};
    RenderObject(const RenderObject &renderObject) = default;
    ~RenderObject();
    void Draw(Shader* shader);

    //glm::mat4 m_model_;


private:
    /* Rendering Data */
    //unsigned int VAO, VBO;// , EBO;
    void setupMesh();

};

class RenderGeometry {
public:
	RenderGeometry() = default;
	RenderGeometry(const std::vector<Vertex>& vertices_,
		const std::vector<unsigned int>& indices_,
		std::shared_ptr<Shader> shader_
	);

	void Draw();
	void ApplyTransform(glm::mat4 transform);
	glm::mat4 GetModelMatrix();
	void SetModelMatrix(const glm::mat4& model);
	void SetColor(glm::vec4 color_);

	std::shared_ptr<Shader> shader;

protected: 
	glm::vec4 color = glm::vec4(0.4, 0.7, 0.9, 1.0);
	std::shared_ptr<RenderObject> ro;
	glm::mat4 world_transform = glm::mat4(1.0f);

	virtual void setup();// = 0;
};

class PointCloud : public RenderGeometry
{
public:

	std::vector<glm::vec3> points;
	//std::vector<glm::vec3> normals;
	PointCloud(const std::vector<glm::vec3>& points_,
		std::shared_ptr<Shader> shader_
		);
	PointCloud(const std::vector<glm::vec3>& points_,
		const std::vector<glm::vec3>& normals_,
		std::shared_ptr<Shader> shader_
		);
	void SetPointScale(float s);

protected:
	//float scale = 0.025;
	float scale = 0.04;
	void point2mesh(glm::vec3 point, std::vector<Vertex>& vertices, std::vector<unsigned int>& indices);
	void point2mesh(glm::vec3 point, glm::vec3 normal, std::vector<Vertex>& vertices, std::vector<unsigned int>& indices);
	void setup() override; 
};

//class PlanePointCloud : public PointCloud {
//public:
//	PlanePointCloud()
//
//protected:
//
//};

class Sphere : public RenderGeometry {
public:
	Sphere(std::shared_ptr<Shader> snhader_);
	void SetPosition(const glm::vec3 position_);
	glm::vec3 GetPosition();

protected:
	float radius = 0.3;
	glm::vec3 position;

	void setup() override;
};


class UIRectangle : public RenderGeometry {
public:

	// NDC space [-1, 1]
	glm::vec2 A;
	glm::vec2 B;

	UIRectangle(std::shared_ptr<Shader> shader_);
	void update();

protected:

	void setup() override;

};



