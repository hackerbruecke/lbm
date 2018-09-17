#pragma once

#include <vtkSmartPointer.h>


namespace lbm
{
namespace io
{
#if !defined (LEGACY_WRITER)

namespace
{

template <typename lattice_model>
void compute_coordinates(const Domain<lattice_model>& domain,
        vtkSmartPointer<vtkPoints>& points)
{
    const auto xl = domain.xlength();
    const auto yl = domain.ylength();
    const auto zl = domain.zlength();
    const auto xs = domain.xspacing();
    const auto ys = domain.yspacing();
    const auto zs = domain.zspacing();
    const auto xo = domain.xorigin();
    const auto yo = domain.yorigin();
    const auto zo = domain.zorigin();

    /* Write lattice coordinates */
    for (auto z = 1u; z < zl + 1; ++z) {
        for (auto y = 1u; y < yl + 1; ++y) {
            for (auto x = 1u; x < xl + 1; ++x) {
                points->InsertNextPoint(xs*x+xo-1, ys*y+yo-1, zs*z+zo-1);
            }
        }
    }
}

}//namespace

template<typename lattice_model>
void write_vtk_file(const Domain<lattice_model>& domain, const std::string& output_dir,
        const std::string& output_filename, uint64_t t)
{
    const auto xl = domain.xlength();
    const auto yl = domain.ylength();
    const auto zl = domain.zlength();
    // Compute point coordinates
    auto points = vtkSmartPointer<vtkPoints>::New();
    compute_coordinates(domain, points);

    // Compute velocity and density vectors
    auto velocities = vtkSmartPointer<vtkDoubleArray>::New();
    auto densities = vtkSmartPointer<vtkDoubleArray>::New();

    velocities->SetNumberOfComponents(lattice_model::D);
    densities->SetNumberOfComponents(1);

    velocities->SetName("Velocity");
    densities->SetName("Density");

    for (auto z = 1u; z < zl + 1; ++z) {
        for (auto y = 1u; y < yl + 1; ++y) {
            for (auto x = 1u; x < xl + 1; ++x) {
                auto current_cell = domain.cell(x, y, z);
                auto density = current_cell.density();
                auto vel = current_cell.velocity(density);

                densities->InsertNextTuple1(density);
                velocities->InsertNextTuple3(vel[0], vel[1], vel[2]);
            }
        }
    }

    // Create a grid and write coordinates and velocity/density
    auto structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
    structuredGrid->SetDimensions(xl, yl, zl);
    structuredGrid->SetPoints(points);
    structuredGrid->GetPointData()->SetVectors(velocities);
    structuredGrid->GetPointData()->SetScalars(densities);
    // Save filename as a combination of passed filename and timestep
    std::stringstream sstr;
    sstr << output_dir << "/" << output_filename << "." << t << ".vts";
    // Write file
    auto writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName(sstr.str().c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}

#endif //!LEGACY_WRITER

template<typename lattice_model, typename solid_collision_model>
auto read_vtk_point_file(const std::string& filename,
        FluidCollision<lattice_model>& fluid_collision_model) -> Domain_ptr<lattice_model>
{
    vtkObject::GlobalWarningDisplayOff();
    auto reader = vtkSmartPointer<vtkDataSetReader>::New();

    reader->SetFileName(filename.c_str());
    reader->Update();

    if (!reader->IsFileValid("structured_points"))
        throw std::logic_error("VTK file \"" + filename + "\" does not exist or does not "
                "seem to be a valid structured grids file!");

    auto output = reader->GetOutput();
    if (output == nullptr) {
        // TODO: Improve exception handling
        throw std::logic_error("Could not read file!");
    }
    output->Register(reader);

    // Now check for point data
    auto dataset = reader->GetOutput();
    auto pd = dataset->GetPointData();

    Domain_ptr<lattice_model> domain;
    if (pd != nullptr) {
        auto p = dynamic_cast<vtkStructuredPoints*>(dataset);

        // Read Dimensions
        assert(p->GetDataDimension() == 3);
        const auto xl = p->GetDimensions()[0];
        const auto yl = p->GetDimensions()[1];
        const auto zl = p->GetDimensions()[2];
        // Read Origins
        double ox, oy, oz;
        p->GetOrigin(ox, oy, oz);
        // Read Spacing
        double sx, sy, sz;
        p->GetSpacing(sx, sy, sz);
        // Create Domain
        domain = make_unique<Domain<lattice_model>>(xl, yl, zl, fluid_collision_model, ox, oy, oz, sx, sy, sz);
        // Get domain points
        const auto points = dynamic_cast<vtkUnsignedCharArray*>(pd->GetArray(0));
        static auto solid_collision = solid_collision_model(*domain);
        // Read only inner cells since we need to apply boundary conditions to the outer cells
        // Of course, we read the cells in a linear manner
        int index = 0;
        for (int z = 1; z < zl + 1; ++z) {
            for (int y = 1; y < yl + 1; ++y) {
                for (int x = 1; x < xl + 1; ++x) {
                    const bool fluid = points->GetValue(index++);
                    if (!fluid)
                        domain->cell(x, y, z).set_collision_handler(&solid_collision);
                }
            }
        }
    }
    else {
        dataset->Delete();
        throw std::logic_error("Point data is null!");
    }
    dataset->Delete();
    return std::move(domain);
}

} //namespace output
} //namespace lbm
