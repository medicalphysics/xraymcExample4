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

// DXBeam models a diagnostic X-ray tube: polychromatic spectrum, filtration, and collimation.
using Beam = xraymc::DXBeam<>;

// constructWorld builds the scene: a lead-lined room containing a concrete hollow sphere.
auto constructWorld()
{
    // EnclosedRoom is a box-shaped shell.
    using Room = xraymc::EnclosedRoom<>;

    // The World template parameters list every item type that can be added; the BVH acceleration
    // structure is built over all of them when world.build() is called.
    xraymc::World<Room, HollowSphere> world;

    // Pre-allocating avoids internal reallocation that would invalidate the references
    // returned by addItem().
    world.reserveNumberOfItems(2);

    // --- Lead-lined room ---
    // Inner walls span a 100 cm cube centred at the origin; 0.2 cm (2 mm) of lead shielding.
    auto& room = world.template addItem<Room>("Room");
    room.setInnerRoomAABB({ -50, -50, -50, 50, 50, 50 });
    room.setWallThickness(0.2);
    const auto lead = xraymc::Material<>::byZ(82).value(); // Z=82 → Pb
    const auto lead_dens = xraymc::AtomHandler::Atom(82).standardDensity;
    room.setMaterial(lead, lead_dens);

    // --- Hollow concrete sphere ---
    // Default: 7 cm outer radius, 0.5 cm wall, centred at the origin.
    // No further configuration needed here; defaults are set in the HollowSphere constructor.
    auto& hsphere = world.template addItem<HollowSphere>("Hollow sphere");

    // Build the BVH acceleration structure over all items before transport or visualization.
    world.build();
    return world;
}

// simulate runs the Monte Carlo transport with a progress readout to the console.
// nThreads == 0 lets the runtime pick the hardware thread count automatically.
void simulate(auto& world, auto& beam, std::uint32_t nThreads = 0)
{
    xraymc::Transport transport;
    transport.runConsole(world, beam, nThreads);
}

// visualizeWorld ray-traces the scene from a fixed camera position and saves a PNG.
// beam is optional: when provided, lines showing the collimation cone are overlaid.
// zoom < 1 zooms out; zoom > 1 zooms in relative to the auto-computed FOV.
void visualizeWorld(const auto& world, auto* beam = nullptr, double zoom = 1)
{
    // 3:2 aspect ratio; increase resy for higher resolution at the cost of render time.
    constexpr int resy = 1024 * 1;
    constexpr int resx = (resy * 3) / 2;

    xraymc::VisualizeWorld viz(world);
    viz.setDistance(1000); // Camera distance from the world origin in cm

    // Buffer holds per-pixel radiance values as doubles before tone-mapping to PNG.
    auto buffer = viz.template createBuffer<double>(resx, resy);

    if (beam) {
        // Draw the beam cone edges out to 95 cm with 0.2 cm line width for clarity.
        viz.addLineProp(*beam, 95, 0.2);
    }

    // Camera angles: 60° azimuth (around Z), 45° polar (elevation) gives an isometric-like view.
    viz.setAzimuthalAngleDeg(60);
    viz.setPolarAngleDeg(45);

    // suggestFOV derives a field-of-view that fits the whole world AABB into the frame.
    viz.suggestFOV(zoom);
    viz.generate(world, buffer);
    viz.savePNG("hollowsphere.png", buffer);
}

void runSimulation()
{
    std::cout << "Starting simulation\n";
    std::cout << "Building the world...";

    // create world
    auto world = constructWorld();
    std::cout << "Done\n";

    // --- Beam configuration ---
    // 65 kVp tube voltage
    Beam beam;
    beam.setTubeVoltage(65);
    beam.clearTubeFiltrationMaterials();
    beam.addTubeFiltrationMaterial(13, 3.5); // 3.5 mm aluminium (Z=13) inherent + added filter
    beam.addTubeFiltrationMaterial(29, 0.1); // 0.1 mm copper (Z=29) additional hardening filter

    // ±2° half-angle in both directions → small field at isocenter (≈ 7 cm radius at 100 cm SID).
    beam.setCollimationHalfAnglesDeg(2, 2);

    beam.setPosition({ 0, 0, -30 }); // X-ray source 30 cm from the origin along the negative Z axis
    beam.setDirectionCosines({ 1, 0, 0 }, { 0, 1, 0 }); // Row = +X, column = +Y; beam central ray points in +Z (towards the sphere)

    // 128 independent exposures each run on a separate thread; total histories = 128 × 1e6.
    beam.setNumberOfExposures(128);
    beam.setNumberOfParticlesPerExposure(1E3);

    // DAP = 1 mGy·cm² normalises output so reported dose values are per unit DAP.
    beam.setDAPvalue(1);

    // Ray-trace a preview image before the Monte Carlo run.
    visualizeWorld(world, &beam);

    // Run the Monte Carlo simulation; runConsole prints progress to stdout.
    std::cout << "Running simulation with " << beam.numberOfParticles() << " particles\n";
    xraymc::Transport transport;
    const auto time = transport.runConsole(world, beam);
    std::cout << "Simulation time: " << time << std::endl;

    // Print some statistics from the sumulation
    // Obtain the type erased hollow sphere object
    auto sphere_cand = world.getItemPointerFromName("Hollow sphere");
    if (const auto* sphere = std::get_if<HollowSphere>(sphere_cand)) {
        const auto& doseScored = sphere->doseScored();
        std::cout << "Dose to hollow sphere wall: " << doseScored.dose();
        std::cout << " (" << doseScored.standardDeviation();
        std::cout << ") mGy/mGycm2 (uncertainty)\n";
        std::cout << "Number of events in sphere wall (photoelectric, coherent and Compton): ";
        std::cout << doseScored.numberOfEvents() << std::endl;
    }
}

int main()
{
    runSimulation();
    return EXIT_SUCCESS;
}