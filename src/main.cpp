/*This file is part of xraymcExample4.

xraymcExample4 is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

xraymcExample4 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with xraymcExample3. If not, see < https://www.gnu.org/licenses/>.

Copyright 2026 Erlend Andersen
*/

#include "xraymc/xraymc.hpp"

#include "hollowsphere.hpp"

// Pretty beam name
using Beam = xraymc::DXBeam<>;

// Construct the world with all relevant items
auto constructWorld()
{

    using Room = xraymc::EnclosedRoom<>;

    // Create world object
    xraymc::World<Room, HollowSphere> world;

    // Pre-allocating avoids internal reallocation that would invalidate the references
    // returned by addItem().
    world.reserveNumberOfItems(2);

    // Make a room of 2 mm lead thickness
    auto& room = world.template addItem<Room>("Room");
    room.setInnerRoomAABB({ -50, -50, -50, 50, 50, 50 }); // Size of inner walls in cm
    room.setWallThickness(0.2); // Wall thickness in cm
    const auto lead = xraymc::Material<>::byZ(82).value(); // Create lead material
    const auto lead_dens = xraymc::AtomHandler::Atom(82).standardDensity; // of standard lead density
    room.setMaterial(lead, lead_dens);

    // add a hollow sphere with center at (0,0,0)
    auto& hsphere = world.template addItem<HollowSphere>("Hollow sphere");

    // Build the world and construct ray tracing structures
    world.build();
    return world;
}

void simulate(auto& world, auto& beam, std::uint32_t nThreads = 0)
{
    xraymc::Transport transport;
    transport.runConsole(world, beam, nThreads);
}

// Generate an image of the world using xraymc path tracing
void visualizeWorld(const auto& world, auto* beam = nullptr, double zoom = 1)
{
    // Set the image resolution
    constexpr int resy = 1024 * 1;
    constexpr int resx = (resy * 3) / 2;

    // Create the visualization object
    xraymc::VisualizeWorld viz(world);
    viz.setDistance(1000); // Setting camera view distance in cm

    // Create an image buffer
    auto buffer = viz.template createBuffer<double>(resx, resy);

    // Add lines illustrating the beam collimation
    if (beam) {
        viz.addLineProp(*beam, 95, 0.2);
    }

    // Set the camera azimuthal angle
    viz.setAzimuthalAngleDeg(60);

    viz.setPolarAngleDeg(45);
    viz.suggestFOV(zoom); // Suggest camera FOV based on world size
    viz.generate(world, buffer); // Raytrace image
    viz.savePNG("hollowsphere.png", buffer); // Save image to png
}

void runSimulation()
{
    std::cout << "Starting simulation\n";

    std::cout << "Building the world...";
    // create world
    auto world = constructWorld();
    std::cout << "Done\n";

    // add a beam
    Beam beam;
    beam.setTubeVoltage(65);
    beam.clearTubeFiltrationMaterials();
    beam.addTubeFiltrationMaterial(13, 3.5);
    beam.addTubeFiltrationMaterial(29, 0.1);

    beam.setCollimationHalfAnglesDeg(2, 2);
    // 128 independent exposures each run on a separate thread; total histories = 128 × 1e6.
    beam.setNumberOfExposures(128);
    beam.setNumberOfParticlesPerExposure(1E6);
    beam.setDAPvalue(1); // Normalization DAP value in mGycm^2

    // Rebuild after rotating the C-arm so the BVH acceleration structures
    // reflect the updated triangle positions.
    world.build();

    // show the simulation
    visualizeWorld(world, &beam);

    // Run the simulation
    std::cout << "Running simulation with " << beam.numberOfParticles() << " particles\n";
    xraymc::Transport transport;
    transport.runConsole(world, beam);
}

int main()
{

    runSimulation();
    return EXIT_SUCCESS;
}