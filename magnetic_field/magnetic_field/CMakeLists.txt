# CMake entry point
cmake_minimum_required (VERSION 3.0)
project (magnetic_field)

find_package(OpenGL REQUIRED)

set(ALL_LIBS
	${OPENGL_LIBRARY}
	glfw
	GLEW
)

add_executable(magnetic_field
	magnetic_field.cpp
	../common/shader.cpp
	../common/shader.hpp
	../common/texture.cpp
	../common/texture.hpp
	../common/controls.cpp
	../common/controls.hpp
	Particle.fragmentshader
	Particle.vertexshader
)

target_link_libraries(magnetic_field
	${ALL_LIBS}
)

set(EXECUTABLE_OUTPUT_PATH "/home/sushranth/Downloads/Advanced_Computer_Graphics/magnetic_field/magnetic_field") 


