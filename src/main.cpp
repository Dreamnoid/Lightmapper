#include <vector>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <string>

#include <xatlas/xatlas.h>
#include <embree3/rtcore.h>
#include <OpenImageDenoise/oidn.h>
#include <glm/glm.hpp>

#include "system.h"
#include "embree.h"
#include "oidn.h"
#include "lm.h"
#include "smesh.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

class RngMwc
{
public:
	RngMwc(uint32_t _z = 12345, uint32_t _w = 65435);
	void reset(uint32_t _z = 12345, uint32_t _w = 65435);
	uint32_t gen();

private:
	uint32_t m_z;
	uint32_t m_w;
};

inline RngMwc::RngMwc(uint32_t _z, uint32_t _w)
	: m_z(_z)
	, m_w(_w)
{
}

inline void RngMwc::reset(uint32_t _z, uint32_t _w)
{
	m_z = _z;
	m_w = _w;
}

inline uint32_t RngMwc::gen()
{
	m_z = 36969 * (m_z & 65535) + (m_z >> 16);
	m_w = 18000 * (m_w & 65535) + (m_w >> 16);
	return (m_z << 16) + m_w;
}

template <typename Rng>
inline float frnd(Rng* _rng)
{
	uint32_t rnd = _rng->gen() & UINT16_MAX;
	return float(rnd) * 1.0f / float(UINT16_MAX);
}

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

// https://en.wikipedia.org/wiki/Halton_sequence
static float halton(int index, int base)
{
	float result = 0;
	float f = 1;
	while (index > 0)
	{
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

constexpr float kPi2 = 6.2831853071795864769252867665590f;
constexpr int kMaxDepth = 8;
constexpr int kPassCount = 100;
constexpr float kNear = 0.01f;

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

struct TexelData
{
	glm::vec3 accumColor;
	uint32_t numColorSamples;
	int numPathsTraced; // Seed for randomDirHemisphere.
	float randomOffset[2];
};

struct SampleLocation
{
	glm::vec3 pos;
	glm::vec3 normal;
	uint32_t uv[2];
};

static std::vector<uint32_t> sortSamples(const std::vector<SampleLocation>& sampleLocations, int width, int height)
{
	std::vector<uint32_t> sampleLocationRanks(sampleLocations.size());
	for (uint32_t i = 0; i < sampleLocationRanks.size(); ++i)
	{
		sampleLocationRanks[i] = i;
	}
	auto compareSampleLocations = [sampleLocations, width, height](const uint32_t& i0, const uint32_t& i1)
	{
		const uint32_t offset0 = sampleLocations[i0].uv[0] + sampleLocations[i0].uv[1] * width;
		const uint32_t offset1 = sampleLocations[i1].uv[0] + sampleLocations[i1].uv[1] * width;
		if (offset0 < offset1)
			return -1;
		else if (offset0 > offset1)
			return 1;
		return 0;
	};
	std::sort(sampleLocationRanks.begin(), sampleLocationRanks.end(), compareSampleLocations);
	return sampleLocationRanks;
}

static void traceRays(int width, int height, const std::vector<SMESH::Vertex>& vertices, std::vector<float>& lightmap, RTCScene embreeScene)
{
	int numTrianglesRasterized = 0;
	std::vector<SampleLocation> sampleLocations;

	// Don't allow duplicate samples at the same uv.
	std::vector<bool> sampleExists(width * height);
	for (int i = 0; i < width * height; ++i)
	{
		sampleExists[i] = false;
	}

	for (uint32_t tri = 0; tri < uint32_t(vertices.size() / 3); ++tri)
	{
		lm_context ctx;
		ctx.rasterizer.x = ctx.rasterizer.y = 0;

		glm::vec2 uvMin (FLT_MAX, FLT_MAX);
		glm::vec2 uvMax (-FLT_MAX, -FLT_MAX);

		for (int i = 0; i < 3; ++i) 
		{
			const SMESH::Vertex& vertex = vertices[tri * 3 + i];

			ctx.triangle.p[i] = glm::vec3(vertex.X, vertex.Y, vertex.Z);
			ctx.triangle.uv[i].x = vertex.U2 * width;
			ctx.triangle.uv[i].y = vertex.V2 * height;

			// update bounds on lightmap
			uvMin = glm::min(uvMin, ctx.triangle.uv[i]);
			uvMax = glm::max(uvMax, ctx.triangle.uv[i]);
		}

		// Calculate area of interest (on lightmap) for conservative rasterization.
		glm::vec2 bbMin = glm::floor(uvMin);
		glm::vec2 bbMax = glm::ceil(uvMax);
		ctx.rasterizer.minx = ctx.rasterizer.x = glm::max((int)bbMin.x - 1, 0);
		ctx.rasterizer.miny = ctx.rasterizer.y = glm::max((int)bbMin.y - 1, 0);
		ctx.rasterizer.maxx = glm::min((int)bbMax.x + 1, width);
		ctx.rasterizer.maxy = glm::min((int)bbMax.y + 1, height);

		assert(ctx.rasterizer.minx <= ctx.rasterizer.maxx && ctx.rasterizer.miny <= ctx.rasterizer.maxy);
		if (lm_findFirstConservativeTriangleRasterizerPosition(&ctx)) 
		{
			while(true)
			{
				SampleLocation sample;
				sample.pos = ctx.sample.position;
				sample.normal = ctx.sample.direction;
				sample.uv[0] = ctx.rasterizer.x;
				sample.uv[1] = ctx.rasterizer.y;

				const uint32_t offset = sample.uv[0] + sample.uv[1] * width;
				if (!sampleExists[offset]) 
				{
					sampleExists[offset] = true;
					sampleLocations.push_back(sample);
				}

				if (!lm_findNextConservativeTriangleRasterizerPosition(&ctx))
				{
					break;
				}
			}
		}
	}

	auto sampleLocationRanks = sortSamples(sampleLocations, width, height);

	RngMwc rng;
	rng.reset();

	std::vector<TexelData> texels(sampleLocations.size());
	for (uint32_t i = 0; i < (uint32_t)sampleLocations.size(); i++)
	{
		TexelData& texel = texels[i];
		texel.accumColor = glm::vec3(0.0f);
		texel.numColorSamples = 0;
		texel.numPathsTraced = 0;
		texel.randomOffset[0] = frnd(&rng);
		texel.randomOffset[1] = frnd(&rng);
	}

	int samplesPerTexelCount = 0;
	for (int pass = 0; pass < kPassCount; ++pass)
	{
		if (pass % 10 == 0)
		{
			std::cout << "Pass " << pass << std::endl;
		}
		
		for (uint32_t i = 0; i < sampleLocations.size(); i++)
		{
			TexelData& texel = texels[i];
			const SampleLocation& sample = sampleLocations[sampleLocationRanks[i]];
			glm::vec3 color(0.0f);
			glm::vec3 throughput(1.0f);

			glm::vec3 rayOrigin = sample.pos;
			glm::vec3 rayDir = glm::normalize(randomDirHemisphere(texel.numPathsTraced, texel.randomOffset, sample.normal));
			for (int depth = 0; depth < kMaxDepth; ++depth)
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

				embree::Intersect1(embreeScene, &context, &rh);

				texel.numPathsTraced++;

				if (rh.hit.geomID == RTC_INVALID_GEOMETRY_ID) // Ray missed
				{
					static const glm::vec3 skyColor(1.f);
					color = color + (throughput * skyColor);
					break;
				}

				const SMESH::Vertex& v0 = vertices[rh.hit.primID * 3 + 0];
				const SMESH::Vertex& v1 = vertices[rh.hit.primID * 3 + 1];
				const SMESH::Vertex& v2 = vertices[rh.hit.primID * 3 + 2];

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
					color += (throughput * emission);
				}
				else
				{
					throughput *= diffuse;
				}

				// Russian Roulette
				// https://computergraphics.stackexchange.com/questions/2316/is-russian-roulette-really-the-answer
				const float p = glm::max(throughput.x, glm::max(throughput.y, throughput.z));
				if (frnd(&rng) > p)
				{
					break;
				}

				throughput *= (1.0f / p);
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

			texel.accumColor += color;
			texel.numColorSamples++;
		}

		samplesPerTexelCount++;
	}

	// Copy texel data to lightmap.
	for (uint32_t i = 0; i < sampleLocations.size(); ++i)
	{
		const SampleLocation& sample = sampleLocations[sampleLocationRanks[i]];

		float intensity = 1.f;
		{
			const static glm::vec3 sunPosition = glm::vec3(1000.f, 1000.f, 1000.f);
			glm::vec3 dirToSun = glm::normalize(sunPosition - sample.pos);

			RTCIntersectContext context;
			rtcInitIntersectContext(&context);

			RTCRay ray;
			ray.org_x = sample.pos.x;
			ray.org_y = sample.pos.y;
			ray.org_z = sample.pos.z;
			ray.tnear = kNear;
			ray.dir_x = dirToSun.x;
			ray.dir_y = dirToSun.y;
			ray.dir_z = dirToSun.z;
			ray.tfar = FLT_MAX;
			ray.mask = UINT32_MAX;
			ray.id = 0;
			ray.flags = 0;
			embree::Occluded1(embreeScene, &context, &ray);
			if (ray.tfar != FLT_MAX)
			{
				intensity = 0.5f;
			}
		}

		TexelData& texel = texels[i];
		float* rgba = &lightmap[(sample.uv[0] + sample.uv[1] * width) * 4];
		const float invn = 1.0f / (float)texel.numColorSamples;
		rgba[0] = texel.accumColor.x * invn * intensity;
		rgba[1] = texel.accumColor.y * invn * intensity;
		rgba[2] = texel.accumColor.z * invn * intensity;
		rgba[3] = 1.0f;
	}
}

std::vector<float> denoise(std::vector<float>& lightmap, int width, int height)
{
	std::vector<float> denoisedLightmap(lightmap.size());

	OIDNDevice device = oidn::NewDevice(OIDN_DEVICE_TYPE_DEFAULT);
	if (!device)
	{
		fatalError("OIDN device creation failed");
	}

	oidn::SetDeviceErrorFunction(device, oidnError, nullptr);
	oidn::SetDevice1b(device, "setAffinity", false);
	oidn::CommitDevice(device);
	OIDNFilter filter = oidn::NewFilter(device, "RTLightmap");
	oidn::SetSharedFilterImage(filter, "color", &lightmap[0], OIDN_FORMAT_FLOAT3, width, height, 0, sizeof(float) * 4, 0);
	oidn::SetSharedFilterImage(filter, "output", &denoisedLightmap[0], OIDN_FORMAT_FLOAT3, width, height, 0, sizeof(float) * 4, 0);
	oidn::CommitFilter(filter);
	oidn::ExecuteFilter(filter);
	oidn::ReleaseFilter(filter);
	oidn::ReleaseDevice(device);

	for (int i = 0; i < width * height; i++)
	{
		denoisedLightmap[i * 4 + 3] = lightmap[i * 4 + 3];
	}

	return denoisedLightmap;
}

xatlas::Atlas* generateUVs(std::vector<SMESH::Vertex>& vertices, const std::vector<uint32_t>& indices)
{
	xatlas::Atlas* atlas = xatlas::Create();

	xatlas::MeshDecl meshDecl;
	meshDecl.vertexCount = vertices.size();
	meshDecl.vertexPositionData = &vertices[0];
	meshDecl.vertexPositionStride = sizeof(SMESH::Vertex);
	meshDecl.vertexNormalData = &vertices[0].NormalX;
	meshDecl.vertexNormalStride = sizeof(SMESH::Vertex);
	meshDecl.vertexUvData = &vertices[0].U1;
	meshDecl.vertexUvStride = sizeof(SMESH::Vertex);
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

	for (uint32_t i = 0; i < atlas->meshCount; i++)
	{
		const xatlas::Mesh& outputMesh = atlas->meshes[i];
		for (uint32_t j = 0; j < outputMesh.indexCount; ++j)
		{
			const xatlas::Vertex& vertex = outputMesh.vertexArray[outputMesh.indexArray[j]];
			vertices[j].U2 = vertex.uv[0] / (float)atlas->width;
			vertices[j].V2 = vertex.uv[1] / (float)atlas->height;
		}
	}

	return atlas;
}

std::vector<float> traceLightmap(xatlas::Atlas* atlas, std::vector<SMESH::Vertex>& vertices, const std::vector<uint32_t>& indices)
{
	RTCDevice embreeDevice = embree::NewDevice(nullptr);
	if (!embreeDevice)
	{
		fatalError("embree device creation failed");
	}

	embree::SetDeviceErrorFunction(embreeDevice, embreeError, nullptr);

	RTCGeometry embreeGeometry = embree::NewGeometry(embreeDevice, RTC_GEOMETRY_TYPE_TRIANGLE);
	embree::SetSharedGeometryBuffer(embreeGeometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, &vertices[0], 0, sizeof(SMESH::Vertex), vertices.size());
	embree::SetSharedGeometryBuffer(embreeGeometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, &indices[0], 0, sizeof(uint32_t) * 3, indices.size() / 3);
	embree::CommitGeometry(embreeGeometry);

	RTCScene embreeScene = embree::NewScene(embreeDevice);
	embree::AttachGeometry(embreeScene, embreeGeometry);
	embree::CommitScene(embreeScene);

	std::vector<float> lightmap(atlas->width * atlas->height * 4);
	traceRays(atlas->width, atlas->height, vertices, lightmap, embreeScene);
	return lightmap;
}

int main(int argc, char** argv)
{
	embree::init();
	oidn::init();

	if (argc < 2)
	{
		fatalError("No input file");
	}
	std::string meshFilename(argv[1]);
	std::cout << "Loading mesh " << meshFilename << std::endl;
	SMESH::Mesh mesh = SMESH::load(meshFilename);

	std::vector<SMESH::Vertex> allVertices;
	for (const SMESH::Group& group : mesh)
	{
		for (const SMESH::Vertex& v : group.Vertices)
		{
			allVertices.push_back(v);
		}
	}

	std::vector<uint32_t> indices = SMESH::generateIndicesFromVertices(allVertices);

	std::cout << "Generate UVs" << std::endl;
	std::vector<SMESH::Vertex> allVerticesCopy = allVertices; // generateUVs will slightly move vertices around
	xatlas::Atlas* atlas = generateUVs(allVerticesCopy, indices);

	std::cout << "Trace" << std::endl;
	std::vector<float> lightmap = traceLightmap(atlas, allVerticesCopy, indices);

	std::cout << "Denoise" << std::endl;
	std::vector<float> denoisedLightmap = denoise(lightmap, atlas->width, atlas->height);

	xatlas::Destroy(atlas);
	
	std::cout << "Save" << std::endl;
	std::string lightmapFilename = meshFilename + ".png";
	saveImage(lightmapFilename.c_str(), atlas->width, atlas->height, denoisedLightmap);

	SMESH::copyUV2ToUV(allVerticesCopy, allVertices);

	SMESH::Mesh lightmappedMesh;
	SMESH::Group group;
	group.Name = std::filesystem::path(lightmapFilename).filename().string();
	group.Vertices = allVertices;
	lightmappedMesh.push_back(group);

	SMESH::save(meshFilename + "_lightmap.smesh", lightmappedMesh);

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

