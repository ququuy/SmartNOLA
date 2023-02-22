#include "RenderObject.h"

RenderObject::RenderObject(const std::vector<Vertex> &vertices_, const std::vector<unsigned int>& indices_) :
vertices(vertices_),
indices(indices_)
{
	setupMesh();
}


RenderObject::~RenderObject() {
	glDeleteBuffers(1, & EBO);
	glDeleteBuffers(1, & VBO);
	glDeleteVertexArrays(1, & VAO);
}



void RenderObject::Draw(Shader* shader) {
	//shader.use();
	shader->use();

	// draw mesh
	glBindVertexArray(VAO);
	
	//glPolygonMode( GL_FRONT_AND_BACK, GL_LINE ); 
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL); 
	glDrawElements(GL_TRIANGLES, static_cast<unsigned int>(indices.size()), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

	// always good practice to set everything back to defaults once configured.
	//glActiveTexture(GL_TEXTURE0);
}


void RenderObject::setupMesh() {
	// create buffers/arrays
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	glBindVertexArray(VAO);
	// load data into vertex buffers
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	// A great thing about structs is that their memory layout is sequential for all its items.
	// The effect is that we can simply pass a pointer to the struct and it translates perfectly to a glm::vec3/2 array which
	// again translates to 3/2 floats which translates to a byte array.
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), &vertices[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);

	// set the vertex attribute pointers
	// vertex Positions
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
	// vertex Colors 
	//glEnableVertexAttribArray(1);
	//glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, Normal));
	// vertex Normals
	//glEnableVertexAttribArray(2);
	//glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, TexCoords));

	glBindVertexArray(0);
}

void RenderGeometry::setup() {
	//
}

RenderGeometry::RenderGeometry(const std::vector<Vertex>& vertices_,
	const std::vector<unsigned int>& indices_,
	std::shared_ptr<Shader> shader_
) :
	shader(shader_)
{
	ro = std::make_shared<RenderObject>(vertices_, indices_);
}

void RenderGeometry::Draw() {
	shader->use();
	shader->setMat4("model", world_transform);
	shader->setVec4("bColor", color);

	ro->Draw(shader.get());
}

void RenderGeometry::ApplyTransform(glm::mat4 transform) {
	world_transform = transform * world_transform;
}

glm::mat4 RenderGeometry::GetModelMatrix() {
	return world_transform;
}

void RenderGeometry::SetModelMatrix(const glm::mat4& model) {
	world_transform = model;
}
void RenderGeometry::SetColor(glm::vec4 color_) {
	color = color_;
}



PointCloud::PointCloud(const std::vector<glm::vec3>& points_, 
	std::shared_ptr<Shader> shader_) :
	points(points_)
{
	shader = shader_;
	setup();
}

void PointCloud::point2mesh(glm::vec3 point, std::vector<Vertex> &vertices, std::vector<unsigned int> &indices) {
	unsigned int offset = vertices.size();
	for (int i = -1; i <= 1; i+=2) {
		for (int j = -1; j <= 1; j += 2) {
			for (int k = -1; k <= 1; k += 2) {
				vertices.push_back(
					Vertex{
						glm::vec3(i, j, k) * scale + point,
					}
				);
			}
		}
	}

	{
		indices.push_back(0 + offset);
		indices.push_back(1 + offset);
		indices.push_back(3 + offset);
	
		indices.push_back(0 + offset);
		indices.push_back(3 + offset);
		indices.push_back(2 + offset);
	
		indices.push_back(1 + offset);
		indices.push_back(5 + offset);
		indices.push_back(7 + offset);
	
		indices.push_back(1 + offset);
		indices.push_back(7 + offset);
		indices.push_back(3 + offset);

		indices.push_back(5 + offset);
		indices.push_back(4 + offset);
		indices.push_back(6 + offset);

		indices.push_back(5 + offset);
		indices.push_back(6 + offset);
		indices.push_back(7 + offset);

		indices.push_back(4 + offset);
		indices.push_back(0 + offset);
		indices.push_back(2 + offset);

		indices.push_back(4 + offset);
		indices.push_back(2 + offset);
		indices.push_back(6 + offset);

		indices.push_back(0 + offset);
		indices.push_back(4 + offset);
		indices.push_back(5 + offset);

		indices.push_back(0 + offset);
		indices.push_back(5 + offset);
		indices.push_back(1 + offset);

		indices.push_back(2 + offset);
		indices.push_back(3 + offset);
		indices.push_back(7 + offset);

		indices.push_back(2 + offset);
		indices.push_back(7 + offset);
		indices.push_back(6 + offset);
	}
}

void PointCloud::SetPointScale(float s) {
	scale = s;
}

void PointCloud::setup() {
	std::vector<Vertex> vertices;
	std::vector<unsigned int> indices;
	for (int i = 0; i < points.size(); ++i) {
		point2mesh(points[i], vertices, indices);
	}
	ro = std::make_shared<RenderObject>(vertices, indices);
}




Sphere::Sphere(std::shared_ptr<Shader> shader_) {
	shader = shader_;
	setup();
}

void Sphere::setup() {
	const size_t grid_nx = 10;
	const size_t grid_ny = 10;
	std::vector<Vertex> vertices;
	std::vector<unsigned int> indices;
	unsigned int** positions = new unsigned int* [grid_nx];
	for (size_t i = 0; i < grid_nx; ++i) {
		positions[i] = new unsigned int[grid_ny - 1];
	}
	unsigned int vtop, vbottom;


	float theta = 0.0;
	float phi = 0.0;
	float delta_t = M_PI * 2.0 / grid_nx;
	float delta_p = M_PI / grid_ny;
	glm::vec3 p(0, 0, 1);

	p = p * radius;// +center;
	//vtop = triMesh->add_vertex(TriMesh::Point(p));
	vertices.push_back(Vertex{ p });
	vtop = vertices.size() - 1;

	phi = delta_p;
	for (size_t i = 0; i < grid_nx; ++i, theta += delta_t, phi = delta_p) {
		for (size_t j = 0; j < grid_ny - 1; ++j, phi += delta_p) {
			glm::vec3 p(sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi));
			p = p * radius;// +center;
			vertices.push_back(Vertex{ p });
			positions[i][j] = vertices.size() - 1;// triMesh->add_vertex(TriMesh::Point(p));
		}
	}
	p = glm::vec3(0, 0, -1);
	p = p * radius;// +center;
	//vbottom = triMesh->add_vertex(TriMesh::Point(p));
	vertices.push_back(Vertex{ p });
	vbottom = vertices.size() - 1;


	for (size_t i = 0; i < grid_nx; ++i) {
		size_t next_ = (i + 1) % (grid_nx);
		for (size_t j = 1; j < grid_ny - 1; ++j) {
			indices.push_back(positions[i][j]);
			indices.push_back(positions[next_][j]);
			indices.push_back(positions[next_][j - 1]);

			indices.push_back(positions[i][j]);
			indices.push_back(positions[next_][j - 1]);
			indices.push_back(positions[i][j - 1]);
		}
		indices.push_back(positions[i][0]);
		indices.push_back(positions[next_][0]);
		indices.push_back(vtop);

		indices.push_back(positions[i][grid_ny - 2]);
		indices.push_back(vbottom);
		indices.push_back(positions[next_][grid_ny - 2]);
	}

	ro = std::make_shared<RenderObject>(vertices, indices);

	delete positions;
	return;
}



void Sphere::SetPosition(const glm::vec3 position_) {
	position = position_;
}

glm::vec3 Sphere::GetPosition() {
	return position;
}





UIRectangle::UIRectangle(std::shared_ptr<Shader> shader_) {
	shader = shader_;
	setup();
}

void UIRectangle::update() {

	glm::mat4 mscale = glm::scale(glm::mat4(1.0), 
		glm::vec3(abs(B[0] - A[0])*.5f, abs(B[1] - A[1])*.5f, 1.0));
	glm::mat4 mtranslate = glm::translate(glm::mat4(1.0),
		glm::vec3((A + B) * 0.5f, 0.0)
		);
	world_transform = mtranslate * mscale;

}

void UIRectangle::setup() {
	color.w = 0.5;
	// --- generate mesh
	std::vector<Vertex> vertices;
	std::vector<unsigned int> indices;
	for (int i = -1; i <= 1; i+=2) {
		for (int j = -1; j <= 1; j += 2) {
			vertices.push_back(
				Vertex{
					glm::vec3(i, j, 0)
					}
				);
		}
	}
	indices.push_back(0);
	indices.push_back(1);
	indices.push_back(3);

	indices.push_back(0);
	indices.push_back(3);
	indices.push_back(2);


	// TODO:
	// initialize A, B

	ro = std::make_shared<RenderObject>(vertices, indices);
}
