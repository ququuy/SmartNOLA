#include "IO.h"
using namespace IO;
using namespace GEO;

PN3_Range IO::read_PLY(const char* filepath) {
	std::ifstream in(filepath);
	PN3_Range points;
	if (!CGAL::IO::read_PLY_with_properties(in,
		std::back_inserter(points),
		CGAL::IO::make_ply_point_reader(PN3_Point_map())))
	{
		std::cerr << "Error: cannot read file " << filepath << std::endl;
	}
	return points;
}

void IO::write_PNC3(const char* filename, const std::vector<GEO::PNC_3>& points) {
	// std::ofstream f("30-e1.5-c0.8.ply", std::ios::binary);
	std::ofstream f(filename, std::ios::binary);
	//CGAL::IO::set_binary_mode(f); // The PLY file will be written in the binary format
	CGAL::IO::write_PLY_with_properties(f, points,
		CGAL::make_ply_point_writer(PNC3_Point_map()),
		std::make_tuple(PNC3_Normal_map(),
			CGAL::IO::PLY_property<double>("nx"),
			CGAL::IO::PLY_property<double>("ny"),
			CGAL::IO::PLY_property<double>("nz")),
		std::make_tuple(PNC3_Color_map(),
			CGAL::IO::PLY_property<unsigned char>("red"),
			CGAL::IO::PLY_property<unsigned char>("green"),
			CGAL::IO::PLY_property<unsigned char>("blue"),
			CGAL::IO::PLY_property<unsigned char>("alpha")));
}
