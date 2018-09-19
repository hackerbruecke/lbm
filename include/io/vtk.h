#pragma once

#include "domain.h"
#include "boundary.h"
#include "helper.h"

#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <exception>

// Since VTK uses some deprecated headers we want to ignore the related warnings while including
// to avoid compiler output flooding
#pragma GCC diagnostic push
// TODO: Add GCC pragma! The following one is not working :-/
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkStructuredPoints.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkDoubleArray.h>

#pragma GCC diagnostic pop

namespace lbm {
namespace io {

template<typename lattice_model>
void write_vtk_file(const Domain<lattice_model>& domain, const std::string& output_dir,
                    const std::string& output_filename, uint64_t t);

template<typename lattice_model, typename solid_collision_model>
auto read_vtk_point_file(const std::string& filename,
                         FluidCollision<lattice_model>& fluid_collision_model) -> Domain_ptr<lattice_model>;

} //namespace io
} //namespace lbm

#if defined (LEGACY_WRITER)
#include "vtk_legacy.hpp"
#endif

#include "vtk.hpp"
