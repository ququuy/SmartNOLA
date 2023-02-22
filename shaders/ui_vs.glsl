#version 330 core
layout (location = 0) in vec3 aPos;

// declare an interface block; see 'Advanced GLSL' for what these are.
uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;

void main()
{
    gl_Position = model * vec4(aPos, 1.0);
}

