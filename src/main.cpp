#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

#include <xatlas/xatlas.h>
#include <embree3/rtcore.h>
#include <OpenImageDenoise/oidn.h>
#include <glm/glm.hpp>

#include "system.h"
#include "embree.h"
#include "oidn.h"
#include <lm.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

struct Vertex
{
	float X, Y, Z;
	float NormalX, NormalY, NormalZ;
	float U, V;
};

static void embreeError(void* /*userPtr*/, enum RTCError /*code*/, const char* str)
{
	fatalError(str);
}

static void oidnError(void* /*userPtr*/, OIDNError code, const char* message)
{
	if (code == OIDN_ERROR_CANCELLED)
		return;
	fatalError(message);
}

void saveImage(const char* filename, int width, int height, const std::vector<float>& image)
{
	uint8_t* buffer = (uint8_t*)malloc(sizeof(uint8_t) * width * height * 4);
	for (int y = 0; y < height; ++y)
	{
		for (int x = 0; x < width; ++x)
		{
			int index = y * width + x;
			uint8_t r = image[index * 4 + 0] * 255.f;
			uint8_t g = image[index * 4 + 1] * 255.f;
			uint8_t b = image[index * 4 + 2] * 255.f;
			uint8_t a = image[index * 4 + 3] * 255.f;

			buffer[index * 4 + 0] = r;
			buffer[index * 4 + 1] = g;
			buffer[index * 4 + 2] = b;
			buffer[index * 4 + 3] = a;
		}
	}
	stbi_write_png(filename, width, height, 4, buffer, width * 4);
	free(buffer);
}

std::vector<uint32_t> generateIndicesFromVertices(const std::vector<Vertex>& vertices)
{
	std::vector<uint32_t> indices(vertices.size());
	for (uint32_t i = 0; i < vertices.size(); ++i)
	{
		indices.push_back(i);
	}
	return indices;
}

// https://en.wikipedia.org/wiki/Halton_sequence
static float halton(int index, int base)
{
	float result = 0;
	float f = 1;
	while (index > 0) {
		f /= base;
		result += f * (index % base);
		index = (int)floor(index / (float)base);
	}
	return result;
}

inline float fract(float a)
{
	return a - trunc(a);
}

inline float max(float a, float b)
{
	return (a < b) ? b : a;
}

inline float min(float a, float b)
{
	return (a < b) ? a : b;
}

inline float frand()
{
	uint32_t rnd = rand() & UINT16_MAX;
	return float(rnd) * 1.0f / float(UINT16_MAX);
}

constexpr float kPi2 = 6.2831853071795864769252867665590f;
constexpr int kMaxDepth = 8;

// Precomputed Global Illumination in Frostbite (GDC 2018)
static glm::vec3 randomDirHemisphere(int index, const float* offset, const glm::vec3& normal)
{
	// Generate 2 uniformly-distributed values in range 0 to 1
	float u = halton(index, 3);
	float v = halton(index, 5);
	// Apply per-texel randomization
	u = fract(u + offset[0]);
	v = fract(v + offset[1]);
	// Transform unit square sample to uniform hemisphere direction
	const float cosTheta = u * 2.0f - 1.0f;
	const float sinTheta = sqrtf(1.0f - cosTheta * cosTheta);
	const float vrad = v * kPi2;
	const float sinPhi = sin(vrad);
	const float cosPhi = cos(vrad);
	glm::vec3 dir(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
	if (glm::dot(dir, normal) < 0.0f)
		return glm::vec3(-dir.x, -dir.y, -dir.z);
	return dir;
}

// https://github.com/aras-p/ToyPathTracer
static void traceRays(uint32_t start, uint32_t end, const std::vector<Vertex>& vertices)
{
	const float kNear = 0.01f * modelGetScale();
	// Initialize paths.
	for (uint32_t i = start; i < end; i++) 
	{
		TexelData& texel = s_bake.texels[i];
		const SampleLocation& sample = s_bake.sampleLocations[s_bake.sampleLocationRanks[i]];
		glm::vec3 color (0.0f);
		glm::vec3 throughput (1.0f);
		glm::vec3 rayOrigin = sample.pos;
		glm::vec3 rayDir = randomDirHemisphere(texel.numPathsTraced, texel.randomOffset, sample.normal);
		for (int depth = 0; depth < s_bake.options.maxDepth; depth++) 
		{
			RTCIntersectContext context;
			rtcInitIntersectContext(&context);
			alignas(16) RTCRayHit rh;
			rh.ray.org_x = rayOrigin.x;
			rh.ray.org_y = rayOrigin.y;
			rh.ray.org_z = rayOrigin.z;
			rh.ray.tnear = kNear;
			rh.ray.dir_x = rayDir.x;
			rh.ray.dir_y = rayDir.y;
			rh.ray.dir_z = rayDir.z;
			rh.ray.time = 0.0f;
			rh.ray.tfar = FLT_MAX;
			rh.ray.mask = UINT32_MAX;
			rh.ray.id = 0;
			rh.ray.flags = 0;
			rh.hit.primID = RTC_INVALID_GEOMETRY_ID;
			rh.hit.geomID = RTC_INVALID_GEOMETRY_ID;
			embree::Intersect1(s_bake.embreeScene, &context, &rh);
			s_bake.numRaysTraced++;
			texel.numPathsTraced++;
			if (rh.hit.geomID == RTC_INVALID_GEOMETRY_ID) 
			{
				// Ray missed, use sky color.
				color = color + (throughput * s_bake.options.skyColor);
				break;
			}
			const Vertex& v0 = vertices[rh.hit.primID * 3 + 0];
			const Vertex& v1 = vertices[rh.hit.primID * 3 + 1];
			const Vertex& v2 = vertices[rh.hit.primID * 3 + 2];

			// we got a new ray bounced from the surface; recursively trace it
			glm::vec3 diffuse(0.5f);
			glm::vec3 emission(0.0f);

			//const objzMaterial* mat = s_bake.triMaterials[rh.hit.primID];
			//if (mat) 
			//{
			//	diffuse = glm::vec3(mat->diffuse[0], mat->diffuse[1], mat->diffuse[2]);
			//	emission = glm::vec3(mat->emission[0], mat->emission[1], mat->emission[2]);
			//	float uv[2];
			//	for (int j = 0; j < 2; j++)
			//		uv[j] = v0.texcoord[j] + (v1.texcoord[j] - v0.texcoord[j]) * rh.hit.u + (v2.texcoord[j] - v0.texcoord[j]) * rh.hit.v;
			//	glm::vec3 texelColor;
			//	if (modelSampleMaterialDiffuse(mat, uv, &texelColor))
			//	{
			//		diffuse = diffuse * texelColor;
			//	}
			//	if (modelSampleMaterialEmission(mat, uv, &texelColor))
			//	{
			//		emission = texelColor;
			//	}
			//}

			if (emission.x > 0.0f || emission.y > 0.0f || emission.z > 0.0f)
			{
				color = color + (throughput * emission);
			}
			else
			{
				throughput = throughput * diffuse;
			}

			// Russian Roulette
			// https://computergraphics.stackexchange.com/questions/2316/is-russian-roulette-really-the-answer
			const float p = max(throughput.x, max(throughput.y, throughput.z));
			if (frand() > p)
			{
				break;
			}

			throughput = throughput * (1.0f / p);
			if (depth + 1 < kMaxDepth) 
			{
				// Using barycentrics should be more precise than "origin + dir * rh.ray.tfar".
				glm::vec3 v0pos(v0.X, v0.Y, v0.Z);
				glm::vec3 v1pos(v1.X, v1.Y, v1.Z);
				glm::vec3 v2pos(v2.X, v2.Y, v2.Z);
				rayOrigin = v0pos + (v1pos - v0pos) * rh.hit.u + (v2pos - v0pos) * rh.hit.v;
				const glm::vec3 normal = glm::normalize(glm::vec3(rh.hit.Ng_x, rh.hit.Ng_y, rh.hit.Ng_z));
				rayDir = randomDirHemisphere(texel.numPathsTraced, texel.randomOffset, normal);
			}
		}
		texel.accumColor = texel.accumColor + color;
		texel.numColorSamples++;
	}

	// Copy texel data to lightmap.
	for (uint32_t i = start; i < end; i++)
	{
		TexelData& texel = s_bake.texels[i];
		const SampleLocation& sample = s_bake.sampleLocations[s_bake.sampleLocationRanks[i]];
		float* rgba = &lightmap[(sample.uv[0] + sample.uv[1] * atlas->width) * 4];
		const float invn = 1.0f / (float)texel.numColorSamples;
		rgba[0] = texel.accumColor.x * invn;
		rgba[1] = texel.accumColor.y * invn;
		rgba[2] = texel.accumColor.z * invn;
		rgba[3] = 1.0f;
	}
}


int main(int argc, char** argv)
{
	embree::init();
	oidn::init();

	if (argc < 1)
	{
		fatalError("No input file");
	}

	std::vector<Vertex> vertices;
	// TODO: load mesh

	std::vector<uint32_t> indices = generateIndicesFromVertices(vertices);

	xatlas::Atlas* atlas = xatlas::Create();

	xatlas::MeshDecl meshDecl;
	meshDecl.vertexCount = vertices.size();
	meshDecl.vertexPositionData = &vertices[0];
	meshDecl.vertexPositionStride = sizeof(Vertex);
	meshDecl.vertexNormalData = &vertices[0].NormalX;
	meshDecl.vertexNormalStride = sizeof(Vertex);
	meshDecl.vertexUvData = &vertices[0].U;
	meshDecl.vertexUvStride = sizeof(Vertex);
	meshDecl.indexCount = vertices.size();
	meshDecl.indexData = &indices[0];
	meshDecl.indexFormat = xatlas::IndexFormat::UInt32;

	xatlas::AddMeshError::Enum error = xatlas::AddMesh(atlas, meshDecl, 1);
	if (error != xatlas::AddMeshError::Success)
	{
		fatalError(xatlas::StringForEnum(error));
	}

	xatlas::ChartOptions chartOptions;
	xatlas::ComputeCharts(atlas, chartOptions);

	xatlas::PackOptions packOptions;
	xatlas::PackCharts(atlas, packOptions);

	{
		std::fstream file("output.txt");
		for (uint32_t i = 0; i < atlas->meshCount; i++)
		{
			const xatlas::Mesh& outputMesh = atlas->meshes[i];
			for (uint32_t j = 0; j < outputMesh.indexCount; ++j)
			{
				const xatlas::Vertex& vertex = outputMesh.vertexArray[outputMesh.indexArray[j]];
				float u = vertex.uv[0] / (float)atlas->width;
				float v = vertex.uv[1] / (float)atlas->height;
				file << u << v << std::endl;
			}
		}
	}

	std::vector<float> lightmap(atlas->width * atlas->height * 4);

	RTCDevice embreeDevice = embree::NewDevice(nullptr);
	if (!embreeDevice)
	{
		fatalError("embree device creation failed");
	}

	embree::SetDeviceErrorFunction(embreeDevice, embreeError, nullptr);

	RTCGeometry embreeGeometry = embree::NewGeometry(embreeDevice, RTC_GEOMETRY_TYPE_TRIANGLE);
	embree::SetSharedGeometryBuffer(embreeGeometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, &vertices[0], 0, sizeof(Vertex), vertices.size());
	embree::SetSharedGeometryBuffer(embreeGeometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, &indices[0], 0, sizeof(uint32_t) * 3, indices.size() / 3);
	embree::CommitGeometry(embreeGeometry);

	RTCScene embreeScene = embree::NewScene(embreeDevice);
	embree::AttachGeometry(embreeScene, embreeGeometry);
	embree::CommitScene(embreeScene);

	embree::ReleaseScene(embreeScene);
	embree::ReleaseGeometry(embreeGeometry);
	embree::ReleaseDevice(embreeDevice);

	std::vector<float> denoisedLightmap(atlas->width * atlas->height * 4);

	OIDNDevice device = oidn::NewDevice(OIDN_DEVICE_TYPE_DEFAULT);
	if(!device)
	{
		fatalError("OIDN device creation failed");
	}

	oidn::SetDeviceErrorFunction(device, oidnError, nullptr);
	oidn::SetDevice1b(device, "setAffinity", false);
	oidn::CommitDevice(device);
	OIDNFilter filter = oidn::NewFilter(device, "RTLightmap");
	oidn::SetSharedFilterImage(filter, "color", &lightmap[0], OIDN_FORMAT_FLOAT3, atlas->width, atlas->height, 0, sizeof(float) * 4, 0);
	oidn::SetSharedFilterImage(filter, "output", &denoisedLightmap[0], OIDN_FORMAT_FLOAT3, atlas->width, atlas->height, 0, sizeof(float) * 4, 0);
	oidn::CommitFilter(filter);
	oidn::ExecuteFilter(filter);
	oidn::ReleaseFilter(filter);
	oidn::ReleaseDevice(device);

	for (uint32_t i = 0; i < atlas->width * atlas->height; i++)
	{
		denoisedLightmap[i * 4 + 3] = lightmap[i * 4 + 3];
	}

	xatlas::Destroy(atlas);

	saveImage("lightmap.png", atlas->width, atlas->height, lightmap);

	return 0;
}

#include <xatlas/xatlas.cpp>

#include <windows.h>

void fatalError(const char* error)
{
	std::cerr << error << std::endl;
	exit(EXIT_FAILURE);
}

void* loadLibrary(const char* filename)
{
	return (void*)::LoadLibraryA(filename);
}
void* loadSymbol(void* handle, const char* symbol)
{
	return (void*)::GetProcAddress((HMODULE)handle, symbol);
}

