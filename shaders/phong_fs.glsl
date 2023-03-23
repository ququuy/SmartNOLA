#version 330 core
out vec4 FragColor;

in VS_OUT {
    vec3 FragPos;
    vec3 Normal;
} fs_in;

uniform vec3 CameraPos;
uniform vec4 bColor;

void main()
{           
    vec3 L = normalize(-vec3(1, -1, 1));
    vec3 V = normalize(CameraPos - fs_in.FragPos);
    vec3 H = normalize(L + V);
    vec3 N = normalize(fs_in.Normal);

    vec3 diffuse = max(dot(L, N), 0) * vec3(bColor.x, bColor.y, bColor.z);

    vec3 specular = pow(dot(H, N), 128) * vec3(1.0, 1.0, 1.0) * 0.2;

    vec3 color = diffuse + specular;

    FragColor = vec4(color, bColor.w);
}
