#pragma once

#include "GeometryLib.h"
#include <fstream>

namespace IO {
	GEO::PN3_Range read_PLY(const char* filepath);
	void write_PNC3(const char* filename, const std::vector<GEO::PNC_3>& points);
}