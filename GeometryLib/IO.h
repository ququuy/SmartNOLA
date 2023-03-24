#pragma once

#include "GeometryLib.h"
#include <fstream>

namespace IO {
	inline const char* FAST_PATH(const char* filename) {
		return (std::string(SOLUTION_ROOT_PATH) + "/data/output/" + std::string(filename)).c_str();
	}
	inline const char* ASSETS_PATH(const char* filename) {
		return (std::string(SOLUTION_ROOT_PATH) + "/data/assets/" + std::string(filename)).c_str();
	}
	GEO::PN3_Range read_PLY(const char* filepath);
	void write_PNC3(const char* filename, const std::vector<GEO::PNC_3>& points);
	void write_Point3(const char* filename, const std::vector<GEO::Point_3>& points);

	void write_Seg3(const char* filename, const GEO::Segment3_Range& segs);
	void write_Seg2(const char* filename, const GEO::Segment2_Range& segs);
}