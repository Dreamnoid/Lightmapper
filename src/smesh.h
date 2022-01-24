#pragma once
#include <vector>
#include <fstream>
#include "system.h"

namespace SMESH
{
	struct Vertex
	{
		float X, Y, Z;
		float NormalX, NormalY, NormalZ;
		float U1, V1;
		float U2, V2;
		uint32_t Color;
	};

	struct Group
	{
		std::string Name;
		std::vector<Vertex> Vertices;
	};

	using Mesh = std::vector<Group>;

	int32_t readInt(std::fstream& file)
	{
		int32_t intValue;
		file.read((char*)&intValue, sizeof(int32_t));
		return intValue;
	}

	void writeByte(std::fstream& file, uint8_t value)
	{
		file.write((char*)&value, sizeof(uint8_t));
	}

	int32_t read7BitEncodedInt(std::fstream& file)
	{
		char current;
		int32_t index = 0, result = 0;
		do
		{
			file.read((char*)&current, sizeof(char));
			result |= (current & 127) << index;
			index += 7;
		}   
		while ((current & 128) != 0);
		return result;
	}

	// https://referencesource.microsoft.com/#mscorlib/system/io/binarywriter.cs,2daa1d14ff1877bd,references
	void write7BitEncodedInt(std::fstream& file, int32_t value)
	{
		// Write out an int 7 bits at a time.  The high bit of the byte,
		// when on, tells reader to continue reading more bytes.
		uint32_t v = (uint32_t)value;   // support negative numbers
		while (v >= 0x80)
		{
			writeByte(file, (uint8_t)(v | 0x80));
			v >>= 7;
		}
		writeByte(file, (uint8_t)v);
	}

	void writeString(std::fstream& file, const std::string& str)
	{
		write7BitEncodedInt(file, str.length());
		for (char c : str)
		{
			file.write(&c, sizeof(char));
		}
	}

	std::string readString(std::fstream& file)
	{
		int32_t charsCount = read7BitEncodedInt(file);
		std::string str = "";
		for (int j = 0; j < charsCount; ++j)
		{
			uint8_t character;
			file.read((char*)&character, sizeof(uint8_t));
			str += character;
		}
		return str;
	}

	Mesh load(const std::string& filename)
	{
		std::fstream file(filename, std::ios::in | std::ios::binary);
		if (!file.is_open())
		{
			fatalError("Could not open mesh file");
		}

		int32_t groupsCount = readInt(file);

		std::vector<Group> groups(groupsCount);

		for (int g = 0; g < groupsCount; ++g)
		{
			Group group;
			group.Name = readString(file);

			int32_t verticesCount = readInt(file);

			group.Vertices.reserve(verticesCount);
			for (int i = 0; i < verticesCount; ++i)
			{
				Vertex v;
				file.read((char*)&v, sizeof(Vertex));
				group.Vertices.push_back(v);
			}

			groups.push_back(group);
		}
		file.close();
		return groups;
	}

	void save(const std::string& filename, const Mesh& groups)
	{
		std::fstream file(filename, std::ios::out | std::ios::binary);

		int32_t groupsCount = groups.size();
		file.write((char*)&groupsCount, sizeof(int32_t));

		for(const Group& group : groups)
		{
			writeString(file, group.Name);

			int32_t size = group.Vertices.size();
			file.write((char*)&size, sizeof(int32_t));

			for(const Vertex& v : group.Vertices)
			{
				file.write((char*)&v, sizeof(Vertex));
			}
		}

		file.close();
	}

	std::vector<uint32_t> generateIndicesFromVertices(const std::vector<Vertex>& vertices)
	{
		std::vector<uint32_t> indices;
		indices.reserve(vertices.size());
		for (uint32_t i = 0; i < vertices.size(); ++i)
		{
			indices.push_back(i);
		}
		return indices;
	}

	void copyUV2ToUV(std::vector<Vertex>& from, std::vector<Vertex>& to)
	{
		for (uint32_t i = 0; i < from.size(); ++i)
		{
			Vertex v = from[i];
			to[i].U1 = v.U2;
			to[i].V1 = v.V2;
		}
	}
}