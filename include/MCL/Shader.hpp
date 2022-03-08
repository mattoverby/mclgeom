// Copyright Matt Overby 2022.
// Distributed under the MIT License.
// adapted from (https://r3dux.org)

// Example use:
//	Shader myshader;
//	... OpenGL context creation
//	myshader.init_from_files("myshader.vert", "myshader.frag");
//	myshader.enable();
//	while (rendering)
//  {
//		glUniform3f(myshader.uniform("eye"), 0.f, 0.f, 0.f);
//		... Draw stuff
//	}
//
// Assumes OpenGL extensions already included
//
#ifndef MCL_SHADER_HPP
#define MCL_SHADER_HPP 1

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <stdexcept>

namespace mcl
{

class Shader
{
public:
	Shader() : program_id(0) {}

	~Shader(){ if(program_id) { glDeleteProgram(program_id); } }

	// Init the shader from files (must create OpenGL context first!)
	void init_from_files(const std::string &vertex_file, const std::string &frag_file);

	// Init the shader from strings (must create OpenGL context first!)
	void init_from_strings(const std::string &vertex_source, const std::string &frag_source);

	// Bind the shader to set uniforms
	void enable();

	// Returns the bound location of a named attribute
	GLuint attribute(const std::string &name);

	// Returns the bound location of a named uniform
	GLuint uniform(const std::string &name);
 
private:
	GLuint program_id;
	GLuint vertex_id;
	GLuint fragment_id;

	std::unordered_map<std::string, GLuint> attributes;
	std::unordered_map<std::string, GLuint> uniforms;

	// Initialize the shader, called by init_from_*
	void init(const std::string &vertex_source, const std::string &frag_source);

	// Compiles the shader, called by init
	GLuint compile(const std::string &shaderSource, GLenum type);

}; // end of shader

//
// Implementation
//

GLuint Shader::compile(const std::string &source, GLenum type)
{
	// Generate a shader id
	// Note: Shader id will be non-zero if successfully created.
	GLuint shaderId = glCreateShader(type);
	if (shaderId == 0){ throw std::runtime_error("\nglCreateShader Error"); }

	// Attach the GLSL source code and compile the shader
	const char *shaderchar = source.c_str();
	glShaderSource(shaderId, 1, &shaderchar, NULL);
	glCompileShader(shaderId);

	// Check the compilation status and throw a runtime_error if shader compilation failed
	GLint shaderStatus;
	glGetShaderiv(shaderId, GL_COMPILE_STATUS, &shaderStatus);
	if (shaderStatus == GL_FALSE)
	{
		// Print compile error
		GLint logSize = 0;
		glGetShaderiv(shaderId, GL_INFO_LOG_LENGTH, &logSize);

		// The logSize includes the NULL character
		GLchar errorLog[logSize];
		glGetShaderInfoLog(shaderId, logSize, &logSize, errorLog);
		printf("\nError compiling shader: %s", errorLog);

		// Exit with failure.
		glDeleteShader(shaderId); // Don't leak the shader.
		throw std::runtime_error("\nglCompileShader Error");
	}

	return shaderId;
}


void Shader::init(const std::string &vertex_source, const std::string &frag_source)
{
	// Create the resource
	program_id = glCreateProgram();
	if (program_id == 0) { throw std::runtime_error("\nglCreateProgram Error"); }

	// Compile the shaders and return their id values
	vertex_id = compile(vertex_source, GL_VERTEX_SHADER);
	fragment_id = compile(frag_source, GL_FRAGMENT_SHADER);

	// Attach and link the shader program
	glAttachShader(program_id, vertex_id);
	glAttachShader(program_id, fragment_id);
	glLinkProgram(program_id);

	// Once the shader program has the shaders attached and linked, the shaders are no longer required.
	// If the linking failed, then we're going to abort anyway so we still detach the shaders.
	glDetachShader(program_id, vertex_id);
	glDetachShader(program_id, fragment_id);

	// Check the program link status and throw a runtime_error if program linkage failed.
	GLint programLinkSuccess = GL_FALSE;
	glGetProgramiv(program_id, GL_LINK_STATUS, &programLinkSuccess);
	if (programLinkSuccess != GL_TRUE) { throw std::runtime_error("\nShader Error: Problem with link"); }

	glUseProgram(0);
}


void Shader::init_from_files(const std::string &vertex_file, const std::string &frag_file)
{
	std::string vert_string, frag_string;

	// Load the vertex shader
	std::ifstream vert_in( vertex_file, std::ios::in | std::ios::binary );
	if( vert_in ){ vert_string = (std::string((std::istreambuf_iterator<char>(vert_in)), std::istreambuf_iterator<char>())); }
	else{ throw std::runtime_error("\nShader Error: failed to load \""+vertex_file+"\"" ); }

	// Load the fragement shader
	std::ifstream frag_in( frag_file, std::ios::in | std::ios::binary );
	if( frag_in ){ frag_string = (std::string((std::istreambuf_iterator<char>(frag_in)), std::istreambuf_iterator<char>())); }
	else{ throw std::runtime_error("\nShader Error: failed to load \""+frag_file+"\"" ); }

	init( vert_string, frag_string );
}


void Shader::init_from_strings(const std::string &vertex_source, const std::string &frag_source)
{
	init(vertex_source, frag_source);
}


void Shader::enable()
{
	if( program_id!=0 ){ glUseProgram(program_id); }
	else{ throw std::runtime_error("\nShader Error: Can't enable, not initialized"); }
}


GLuint Shader::attribute(const std::string &name)
{
	// Add the attribute to the map table if it doesn't already exist
	if( attributes.count(name)==0 ){
		attributes[name] = glGetAttribLocation( program_id, name.c_str() );
		if( (int)attributes[name] == -1 ){ throw std::runtime_error("\nShader Error: bad attribute ("+name+")"); }
	}
	return attributes[name];
}


GLuint Shader::uniform(const std::string &name)
{
	// Add the uniform to the map table if it doesn't already exist
	if( uniforms.count(name)==0 ){ 
		uniforms[name] = glGetUniformLocation( program_id, name.c_str() );
		if( (int)uniforms[name] == -1 ){ throw std::runtime_error("\nShader Error: bad uniform ("+name+")"); }
	}
	return uniforms[name];
}

} // end namespace mcl

#endif
