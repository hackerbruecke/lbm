#pragma once

#include <pugixml.hpp>
#include "domain.h"
#include "vtk.h"
#include "boundary.h"
#include "configuration.h"

#include <memory>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

namespace lbm
{
namespace io
{

inline void check_attribute(const pugi::xml_node& node, const std::string& name)
{
    if (!node.attribute(name.c_str()))
        throw std::logic_error("Missing attribute \"" + name +
                "\" for node \"" + node.name() + "\"!");
}

template <typename lattice_model>
auto parse_condition(const pugi::xml_node& boundary, Domain<lattice_model>& domain)
    -> NonFluidCollision<lattice_model>&
{
    using BF = BoundaryKeeper<lattice_model>;
    check_attribute(boundary, "condition");

    const std::string condition = boundary.attribute("condition").value();
    if (condition == "noslip") {
        return BF::template get_collision<NoSlipBoundary<lattice_model>>(domain);
    }
    else if (condition == "movingwall") {
        //wallvelocity[3]
        check_attribute(boundary, "vx");
        check_attribute(boundary, "vy");
        check_attribute(boundary, "vz");
        const lbm::double_array<lattice_model::D> wall_velocity = {
                boundary.attribute("vx").as_double(),
                boundary.attribute("vy").as_double(),
                boundary.attribute("vz").as_double()
        };
        return BF::template get_collision<MovingWallBoundary<lattice_model>>(domain, wall_velocity);
    }
    else if (condition == "freeslip") {
        return BF::template get_collision<FreeSlipBoundary<lattice_model>>(domain);
    }
    else if (condition == "outflow") {
        double rho_ref = 1.0;

        //Ref density optional
        const auto xml_rho_ref = boundary.attribute("rho-ref");
        if (xml_rho_ref)
            rho_ref = xml_rho_ref.as_double();
        return BF::template get_collision<OutflowBoundary<lattice_model>>(domain, rho_ref);
    }
    else if (condition == "inflow") {
        //inflow velocity[3]
        check_attribute(boundary, "vx");
        check_attribute(boundary, "vy");
        check_attribute(boundary, "vz");
        const lbm::double_array<lattice_model::D> inflow_velocity = {
                boundary.attribute("vx").as_double(),
                boundary.attribute("vy").as_double(),
                boundary.attribute("vz").as_double()
        };
        //Ref dens optional
        double rho_ref = 1.0;
        const auto xml_rho_ref = boundary.attribute("rho-ref");
        if (xml_rho_ref)
            rho_ref = xml_rho_ref.as_double();
        return BF::template get_collision<InflowBoundary<lattice_model>>(
                domain, inflow_velocity, rho_ref);
    }
    else if (condition == "pressure") {
        //input density
        check_attribute(boundary, "rho-in");
        double rho_in = boundary.attribute("rho-in").as_double();
        return BF::template get_collision<PressureBoundary<lattice_model>>(domain, rho_in);
    }
    throw std::logic_error(condition + " boundary condition not supported!");
}

template <typename lattice_model>
void parse_boundary(const pugi::xml_node& boundary, Domain<lattice_model>& domain)
{
    check_attribute(boundary, "extent");
    const std::string extent = boundary.attribute("extent").value();
    auto& condition = parse_condition<lattice_model>(boundary, domain);
    const auto xl = domain.xlength();
    const auto yl = domain.ylength();
    const auto zl = domain.zlength();
    if (extent == "z0") {
        domain.setBoundaryCondition(condition, 0, xl+1, 0, yl+1, 0, 0);
    }
    else if (extent == "zmax") {
        domain.setBoundaryCondition(condition, 0, xl+1, 0, yl+1, zl+1, zl+1);
    }
    else if (extent == "x0") {
        domain.setBoundaryCondition(condition, 0, 0, 0, yl+1, 0, zl+1);
    }
    else if (extent == "xmax") {
        domain.setBoundaryCondition(condition, xl+1, xl+1, 0, yl+1, 0, zl+1);
    }
    else if (extent == "y0") {
        domain.setBoundaryCondition(condition, 0, xl+1, 0, 0, 0, zl+1);
    }
    else if (extent == "ymax") {
        domain.setBoundaryCondition(condition, 0, xl+1, yl+1, yl+1, 0, zl+1);
    }
    else {
        std::vector<std::uint64_t> ex;
        boost::char_separator<char> sep(" ");
        boost::tokenizer<boost::char_separator<char>> tokens(extent, sep);
        for (const auto& t : tokens) {
            ex.push_back(boost::lexical_cast<std::uint64_t>(t));
        }
        if (ex.size() < 6)
            throw std::logic_error("Extent \"" + extent + "\" is not complete! Must be six values!");
        domain.setBoundaryCondition(condition, ex[0], ex[1], ex[2], ex[3], ex[4], ex[5]);
    }
}


template <typename lattice_model>
auto parse_scenario_file(const std::string& filename, Config& cfg,
        FluidCollision<lattice_model>& collision) -> std::unique_ptr<Domain<lattice_model>>
{
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());
    if (!result)
        throw std::logic_error("XML file \"" + filename + "\" could not be read properly!");

    auto scenario = doc.child("scenario");
    if (!scenario)
        throw std::logic_error("Scenario node missing!");
    if (!scenario.attribute("name"))
        throw std::logic_error("Scenario name is missing!");
    const auto scenario_name = scenario.attribute("name").value();

    std::cout << "Reading scenario: \"" << scenario_name << "\"..." << std::endl;

    cfg.set_output_filename(scenario_name);
    auto xml_domain = scenario.child("domain");
    if (!xml_domain)
        throw std::logic_error("Domain node is missing!");
    auto vtk_file = xml_domain.attribute("vtk-file");
    auto xml_xl = xml_domain.attribute("xl");
    auto xml_yl = xml_domain.attribute("yl");
    auto xml_zl = xml_domain.attribute("zl");
    double xl, yl, zl;

    std::unique_ptr<Domain<lattice_model>> domain;
    if (vtk_file) {
        domain = lbm::io::read_vtk_point_file<lattice_model,
                lbm::NoSlipBoundary<lattice_model>>(vtk_file.value(), collision);
        xl = domain->xlength();
        yl = domain->ylength();
        zl = domain->zlength();
    }
    else if (xml_xl && xml_yl && xml_zl) {
        assert(xml_xl.as_uint() > 0 && xml_yl.as_uint() > 0 && xml_zl.as_uint() > 0);
        xl = xml_xl.as_uint();
        yl = xml_yl.as_uint();
        zl = xml_zl.as_uint();
        domain = make_unique<lbm::Domain<lattice_model>>(xl, yl, zl, collision);
    }
    else {
        throw std::logic_error(
                "Neither vtk-file nor xl/yl/zl attribute provided to domain node!");
    }
//    std::cout << "X length: " << xl << std::endl;
//    std::cout << "Y length: " << yl << std::endl;
//    std::cout << "Z length: " << zl << std::endl;

    // Iterate through boundary conditions
    for (auto xml_boundary = xml_domain.child("boundary"); xml_boundary;
            xml_boundary = xml_boundary.next_sibling("boundary")) {
        parse_boundary(xml_boundary, *domain);
    }
    return std::move(domain);
}


}//namespace io
}//namespace lbm
