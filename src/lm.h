
/***********************************************************
* A single header file OpenGL lightmapping library         *
* https://github.com/ands/lightmapper                      *
* no warranty implied | use at your own risk               *
* author: Andreas Mantler (ands) | last change: 12.06.2016 *
*                                                          *
* License:                                                 *
* This software is in the public domain.                   *
* Where that dedication is not recognized,                 *
* you are granted a perpetual, irrevocable license to copy *
* and modify this file however you want.                   *
***********************************************************/
#include <glm/glm.hpp>

union FloatUint32
{
	float f;
	uint32_t u;
};

static bool isFinite(float f)
{
	FloatUint32 fu;
	fu.f = f;
	return fu.u != 0x7F800000u && fu.u != 0x7F800001u;
}

static float cross(const glm::vec2& a, const glm::vec2& b)
{
	return a.x * b.y - a.y * b.x;
}

static glm::vec2 toBarycentric(glm::vec2 p1, glm::vec2 p2, glm::vec2 p3, glm::vec2 p)
{
	// http://www.blackpawn.com/texts/pointinpoly/
	// Compute vectors
	glm::vec2 v0 = (p3 - p1);
	glm::vec2 v1 = (p2 - p1);
	glm::vec2 v2 = (p - p1);
	// Compute dot products
	float dot00 = glm::dot(v0, v0);
	float dot01 = glm::dot(v0, v1);
	float dot02 = glm::dot(v0, v2);
	float dot11 = glm::dot(v1, v1);
	float dot12 = glm::dot(v1, v2);
	// Compute barycentric coordinates
	float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
	float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	float v = (dot00 * dot12 - dot01 * dot02) * invDenom;
	return glm::vec2(u, v);
}

static inline int lm_leftOf(glm::vec2 a, glm::vec2 b, glm::vec2 c)
{
	float x = cross((b - a), (c - b));
	return x < 0 ? -1 : x > 0;
}

static bool lm_lineIntersection(glm::vec2 x0, glm::vec2 x1, glm::vec2 y0, glm::vec2 y1, glm::vec2* res)
{
	glm::vec2 dx = (x1- x0);
	glm::vec2 dy = (y1- y0);
	glm::vec2 d = (x0- y0);
	float dyx = cross(dy, dx);
	if (dyx == 0.0f)
		return false;
	dyx = cross(d, dx) / dyx;
	if (dyx <= 0 || dyx >= 1)
		return false;
	res->x = y0.x + dyx * dy.x;
	res->y = y0.y + dyx * dy.y;
	return true;
}

// this modifies the poly array! the poly array must be big enough to hold the result!
// res must be big enough to hold the result!
static int lm_convexClip(glm::vec2 *poly, int nPoly, const glm::vec2*clip, int nClip, glm::vec2 *res)
{
	int nRes = nPoly;
	int dir = lm_leftOf(clip[0], clip[1], clip[2]);
	for (int i = 0, j = nClip - 1; i < nClip && nRes; j = i++)
	{
		if (i != 0)
		{
			for (nPoly = 0; nPoly < nRes; nPoly++)
			{
				poly[nPoly] = res[nPoly];
			}
		}
		nRes = 0;
		glm::vec2 v0 = poly[nPoly - 1];
		int side0 = lm_leftOf(clip[j], clip[i], v0);
		if (side0 != -dir)
		{
			res[nRes++] = v0;
		}
		for (int k = 0; k < nPoly; k++)
		{
			glm::vec2 v1 = poly[k], x;
			int side1 = lm_leftOf(clip[j], clip[i], v1);
			if (side0 + side1 == 0 && side0 && lm_lineIntersection(clip[j], clip[i], v0, v1, &x))
				res[nRes++] = x;
			if (k == nPoly - 1)
				break;
			if (side1 != -dir)
				res[nRes++] = v1;
			v0 = v1;
			side0 = side1;
		}
	}

	return nRes;
}

struct lm_context
{
	struct
	{
		glm::vec3 p[3];
		glm::vec2 uv[3];
	} triangle;

	struct
	{
		int minx, miny;
		int maxx, maxy;
		int x, y;
	} rasterizer;

	struct
	{
		glm::vec3 position;
		glm::vec3 direction;
	} sample;
};

static bool lm_hasConservativeTriangleRasterizerFinished(lm_context *ctx)
{
	return ctx->rasterizer.y >= ctx->rasterizer.maxy;
}

static void lm_moveToNextPotentialConservativeTriangleRasterizerPosition(lm_context *ctx)
{
	if (++ctx->rasterizer.x >= ctx->rasterizer.maxx)
	{
		ctx->rasterizer.x = ctx->rasterizer.minx;
		++ctx->rasterizer.y;
	}
}

static bool lm_trySamplingConservativeTriangleRasterizerPosition(lm_context *ctx)
{
	if (lm_hasConservativeTriangleRasterizerFinished(ctx))
		return false;

	glm::vec2 pixel[16];
	pixel[0] = glm::vec2(ctx->rasterizer.x, ctx->rasterizer.y);
	pixel[1] = glm::vec2(ctx->rasterizer.x + 1, ctx->rasterizer.y);
	pixel[2] = glm::vec2(ctx->rasterizer.x + 1, ctx->rasterizer.y + 1);
	pixel[3] = glm::vec2(ctx->rasterizer.x, ctx->rasterizer.y + 1);

	glm::vec2 res[16];
	int nRes = lm_convexClip(pixel, 4, ctx->triangle.uv, 3, res);
	if (nRes > 0)
	{
		// do centroid sampling
		glm::vec2 centroid = res[0];
		float area = res[nRes - 1].x * res[0].y - res[nRes - 1].y * res[0].x;
		for (int i = 1; i < nRes; i++)
		{
			centroid = (centroid+ res[i]);
			area += res[i - 1].x * res[i].y - res[i - 1].y * res[i].x;
		}
		centroid = (centroid/ (float)nRes);
		area = glm::abs(area / 2.0f);

		if (area > 0.0f)
		{
			// calculate 3D sample position and orientation
			glm::vec2 uv = toBarycentric(
				ctx->triangle.uv[0],
				ctx->triangle.uv[1],
				ctx->triangle.uv[2],
				centroid);

			// sample it only if its's not degenerate
			if (isFinite(uv.x) && isFinite(uv.y))
			{
				const glm::vec3 &p0 = ctx->triangle.p[0];
				const glm::vec3&p1 = ctx->triangle.p[1];
				const glm::vec3&p2 = ctx->triangle.p[2];
				const glm::vec3 v1 = p1 - p0;
				const glm::vec3 v2 = p2 - p0;
				ctx->sample.position = p0 + v2 * uv.x + v1 * uv.y;
				ctx->sample.direction = glm::normalize(glm::cross(v1, v2));

				if (isFinite(ctx->sample.position.x) && isFinite(ctx->sample.position.y) && isFinite(ctx->sample.position.z) &&
					isFinite(ctx->sample.direction.x) && isFinite(ctx->sample.direction.y) && isFinite(ctx->sample.direction.z) &&
					glm::length(ctx->sample.direction) > 0.5f) // don't allow 0.0f. should always be ~1.0f
				{
					return true;
				}
			}
		}
	}
	return false;
}

// returns true if a sampling position was found and
// false if we finished rasterizing the current triangle
static bool lm_findFirstConservativeTriangleRasterizerPosition(lm_context *ctx)
{
	while (!lm_trySamplingConservativeTriangleRasterizerPosition(ctx))
	{
		lm_moveToNextPotentialConservativeTriangleRasterizerPosition(ctx);
		if (lm_hasConservativeTriangleRasterizerFinished(ctx))
			return false;
	}
	return true;
}

static bool lm_findNextConservativeTriangleRasterizerPosition(lm_context *ctx)
{
	lm_moveToNextPotentialConservativeTriangleRasterizerPosition(ctx);
	return lm_findFirstConservativeTriangleRasterizerPosition(ctx);
}