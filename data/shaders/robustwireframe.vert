#version 330 core

in vec3 vertexPosition;
in vec3 vertexNormal;
in vec4 vertexColor;

out EyeSpaceVertex {
    vec3 position;
    vec3 normal;
    vec4 color;
} vs_out;

uniform mat4 modelView;
uniform mat3 modelViewNormal;
uniform mat4 mvp;

void main()
{
    vs_out.normal = normalize( modelViewNormal * vertexNormal );
    vs_out.position = vec3( modelView * vec4( vertexPosition, 1.0 ) );
    vs_out.color = vertexColor;
    gl_Position = mvp * vec4( vertexPosition, 1.0 );
}
