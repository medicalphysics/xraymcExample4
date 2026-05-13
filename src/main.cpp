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



// Set the low energu correction model
// 0: No binding energy correction
// 1: Livermore binding correction in Compton events
// 2: RIA correction of Compton events and emission of characteristic radiation from photoelectric events
constexpr int MODEL = 2;

// Set number of atomic shells to model
constexpr int NSHELLS = 16;

// Specify if we do forced interaction biasing in tetrahedral meshes
constexpr bool FORCEINTERACTIONS = true;

// Making pretty names from templated classes in xraymc
using CarmType = Carm<NSHELLS, MODEL>;
using Mesh = xraymc::TriangulatedMesh<NSHELLS, MODEL>;
using VGrid = xraymc::AAVoxelGrid<NSHELLS, MODEL, 255>;
using Room = xraymc::EnclosedRoom<NSHELLS, MODEL>;
using TetMesh = xraymc::TetrahedalMesh<NSHELLS, MODEL, FORCEINTERACTIONS>;
using Box = xraymc::WorldBox<NSHELLS, MODEL>;
using Surface = xraymc::TriangulatedOpenSurface<NSHELLS, MODEL>;

// Creating a world class from world template
using World = xraymc::World<Mesh, VGrid, Room, TetMesh, CarmType, Box, Surface>;

// Pretty beam name
using Beam = xraymc::DXBeam<false>;

// Making a pretty name for the material class
using Material = xraymc::Material<NSHELLS>;

// Convienence funtion to read the ICRP110 voxel phantoms
static VGrid getICRP110Phantom(bool female = true)
{
    auto phantom_data = female ? ICRP110PhantomReader::readFemalePhantom() : ICRP110PhantomReader::readMalePhantom();
    VGrid phantom;

    // Create the materials from material definitions
    std::vector<Material> materials;
    for (auto& w : phantom_data.mediaComposition()) {
        auto mat_cand = Material::byWeight(w);
        if (mat_cand)
            materials.push_back(mat_cand.value());
        else
            throw std::runtime_error("error");
    }
    // Populate the xraymc voxel object with phantom data
    phantom.setData(phantom_data.dimensions(), phantom_data.densityData(), phantom_data.mediaData(), materials);
    phantom.setSpacing(phantom_data.spacing());
    // Flip and roll the phantom to make it lay on the table
    // ICRP 110 data is stored head-up in the Z direction; roll twice then flip
    // to reorient the phantom so it lies supine along the X axis on the table.
    phantom.rollAxis(2, 0);
    phantom.rollAxis(2, 1);
    phantom.flipAxis(2);
    return phantom;
}

// Convienence funtion to read ICRP 145 phantom
static TetMesh readICRP145Phantom(bool female = true)
{
    // Simply get filenames based on gender
    const std::string name = female ? "MRCP_AF" : "MRCP_AM";
    const std::string elefile = name + ".ele";
    const std::string nodefile = name + ".node";
    const std::string mediafile = name + "_media.dat";
    const std::string organfile = "icrp145organs.csv";

    // Read tetrahedral data from disc
    xraymc::TetrahedalMeshReader reader(nodefile, elefile, mediafile, organfile);

    // ICRP 145 meshes face the -Y direction by default; rotate 180° around Z to face +Y.
    reader.rotate({ 0, 0, 1 }, std::numbers::pi_v<double>);
    // Create tetrahedral object
    return TetMesh { reader.data() };
}

// Construct the world with all relevant items
World constructWorld()
{
    // Create world object
    World world;

    // Pre-allocating avoids internal reallocation that would invalidate the references
    // returned by addItem() below (auto& carm, auto& table, etc.).
    world.reserveNumberOfItems(16);

    // Start by adding the c-arm, CarmType is a custom type that reads a stl file looking like a c-arm
    // the stl file has origin at the source position for zero primary and secondary angles
    auto& carm = world.template addItem<CarmType>({ "carm.stl" }, "C-Arm");
    // Shift the entire C-arm so its isocenter lands at (0,0,0); all other objects
    // are then positioned relative to that isocenter.
    carm.translate(xraymc::vectormath::scale(carm.isoCenter(), -1.0));

    // Adding the patient table, this is a triangulated solid mesh
    xraymc::STLReader stlreader;
    stlreader.setFilePath("table.stl");
    auto& table = world.template addItem<Mesh>({ stlreader() }, "Table");
    // Translate the table (in cm)
    table.translate({ 0, 0, 0 });
    const std::array<double, 6> table_aabb = table.AABB(); // get the bounding box of the table (X0 Y0 Z0 X1 Y1 Z1)
    // Filling the table with carbon
    const auto carbon = Material::byZ(6).value();
    table.setMaterial(carbon);
    table.setMaterialDensity(0.2); // g/cm3

    // add a male ICRP 110 voxel phantom as patient
    auto& patient = world.addItem(getICRP110Phantom(false), "Patient");
    // translate patient to lay on table
    auto patient_aabb = patient.AABB();
    patient.translate({ table_aabb[3] - patient_aabb[3], 0, table_aabb[5] - patient_aabb[2] });

    // Make a room of 2 mm lead thickness
    auto& room = world.template addItem<Room>("Room");
    room.setInnerRoomAABB({ -300, -200, -100, 200, 200, 180 }); // Size of inner walls in cm
    room.setWallThickness(0.2); // Wall thickness in cm
    const auto lead = Material::byZ(82).value(); // Create lead material
    const auto lead_dens = xraymc::AtomHandler::Atom(82).standardDensity; // of standard lead density
    room.setMaterial(lead, lead_dens);

    // Add a concrete floor
    auto& floor = world.template addItem<Box>("Floor");
    floor.setAABB({ -299.99, -199.99, -99.99, 199.99, 199.99, -90 });
    const auto concrete = Material::byNistName("Concrete, Ordinary").value(); // Concrete material
    const auto concrete_dens = xraymc::NISTMaterials::density("Concrete, Ordinary"); // and density
    floor.setMaterial(concrete, concrete_dens);
    auto floor_aabb = floor.AABB();

    // Add a doctor
    auto& doctor = world.template addItem<TetMesh>(readICRP145Phantom(true), "Doctor");
    // Make the doctor stand on the floor and beside the table
    doctor.translate({ -50, table_aabb[1] - doctor.AABB()[4] - 4, floor_aabb[5] - doctor.AABB()[2] });

    // Add a nurse
    auto& nurse = world.template addItem<TetMesh>(readICRP145Phantom(false), "Nurse");
    // Make the doctor stand on the floor and beside the table
    nurse.translate({ doctor.AABB()[0] - (nurse.AABB()[3] - nurse.AABB()[1]) / 2 - 10, table_aabb[1] - nurse.AABB()[4] - 4, floor_aabb[5] - nurse.AABB()[2] });

    // Add ceiling shield
    // Reading mesh data from disc
    stlreader.setFilePath("ceilingShield.stl");
    auto ceilingshield_tri = stlreader();
    // Rotate mesh object data
    std::for_each(std::execution::par_unseq, ceilingshield_tri.begin(), ceilingshield_tri.end(), [](auto& tri) {
        tri.rotate(std::numbers::pi_v<double> / 2, { 1, 0, 0 });
        tri.rotate(std::numbers::pi_v<double> / 2 + xraymc::DEG_TO_RAD() * 60, { 0, 0, 1 });
    });
    // Create a surface mesh object
    auto& ceilingshield = world.template addItem<Surface>({ ceilingshield_tri }, "Ceiling shield");
    // Move the shield to proper position in the world
    ceilingshield.translate({ -15, -35, 30 });
    // Set material and density, in this case just ordinary lead
    ceilingshield.setMaterial(lead, lead_dens);
    ceilingshield.setSurfaceThickness(0.05); // Set the surface thickness on cm

    // Add table shield as a simple box
    auto& table_shield = world.template addItem<Box>("Table shield");
    table_shield.setAABB({ -130, table.AABB()[1] - 1.05, floor.AABB()[5] + 5, 10, table.AABB()[1] - 1.0, table.AABB()[5] + 5.0 });
    table_shield.setMaterial(lead, lead_dens);

    // Add patient shielding blanket as a surface mesh
    stlreader.setFilePath("blanket.stl");
    auto blanket_tri = stlreader();
    auto& blanket = world.template addItem<Surface>({ blanket_tri }, "Blanket");
    auto bc = blanket.center();
    // Position the blanket over the patient
    blanket.translate(xraymc::vectormath::scale(bc, -1.0));
    blanket.translate({ -30, 0, table_aabb[5] - patient_aabb[2] + 10 });
    blanket.setMaterial(lead, lead_dens);
    blanket.setSurfaceThickness(0.05);

    // Build the world
    world.build();
    return world;
}

void simulate(auto& world, auto& beam, std::uint32_t nThreads = 0)
{
    xraymc::Transport transport;
    transport.runConsole(world, beam, nThreads);
}

// Generate an image of the world using xraymc path tracing
void visualizeWorld(const auto& world, auto* beam = nullptr, const std::string& name = "world", double zoom = 1, double maxDose = 0, bool useLogColor = false)
{
    // Set the image resolution
    constexpr int resy = 1024 * 1;
    constexpr int resx = (resy * 3) / 2;

    // Create the visualization object
    xraymc::VisualizeWorld viz(world);
    viz.setLowestDoseColorToWhite(false);
    viz.setDistance(1000); // Setting camera view distance in cm

    // Create an image buffer
    auto buffer = viz.template createBuffer<double>(resx, resy);

    // setting some colors on objects
    std::array<std::uint8_t, 3> color;
    color = { 255, 195, 170 }; // skin
    viz.setColorOfItem(world.getItemPointerFromName("Patient"), color);
    viz.setColorOfItem(world.getItemPointerFromName("Nurse"), color);
    color = { 176, 108, 73 }; // dark skin
    viz.setColorOfItem(world.getItemPointerFromName("Doctor"), color);
    color = { 173, 232, 244 }; // baby blue
    viz.setColorOfItem(world.getItemPointerFromName("Floor"), color);
    color = { 242, 239, 223 }; // ivory white
    viz.setColorOfItem(world.getItemPointerFromName("Room"), color);
    viz.setColorOfItem(world.getItemPointerFromName("C-Arm"), color);
    color = { 55, 55, 55 }; // carbon
    viz.setColorOfItem(world.getItemPointerFromName("Table"), color);
    color = { 105, 155, 103 }; // white green
    viz.setColorOfItem(world.getItemPointerFromName("Table shield"), color);
    viz.setColorOfItem(world.getItemPointerFromName("Ceiling shield"), color);
    viz.setColorOfItem(world.getItemPointerFromName("Blanket"), color);

    // Set objects to be colored by the simulated dose values
    if (maxDose > 0) {
        viz.addColorByValueItem(world.getItemPointerFromName("Doctor"));
        viz.addColorByValueItem(world.getItemPointerFromName("Patient"));
        viz.addColorByValueItem(world.getItemPointerFromName("Nurse"));
        // Set min and max doses to include in color table
        if (useLogColor) {
            viz.setColorByValueMinMax(maxDose / 10000.0, maxDose);
            viz.setColorByValueLogScale(true); // Use log scale for dose coloring
        } else {
            viz.setColorByValueMinMax(0, maxDose);
        }
    }

    // Add lines illustrating the beam collimation
    if (beam) {
        viz.addLineProp(*beam, 95, 0.2);
    }

    // Pretty format int to string lambda function
    auto to_string = [](int n) -> std::string {
        std::stringstream ss;
        ss << std::setw(3) << std::setfill('0') << n;
        return ss.str();
    };

    // Set the camera azimuthal angle
    viz.setAzimuthalAngleDeg(60);

    // Iterate over polar angles to create images from different angles
    std::string message;
    constexpr int low_angle = 0;
    constexpr int high_angle = 360;
    constexpr int degStep = 30;
    for (int i = low_angle; i < high_angle; i = i + degStep) {
        std::string fn = name + to_string(i) + ".png";
        viz.setPolarAngleDeg(i);
        viz.suggestFOV(zoom); // Suggest camera FOV based on world size
        viz.generate(world, buffer); // Raytrace image
        viz.savePNG(fn, buffer); // Save image to png
    }
}

void loadAndShowSimulation(const std::string& filename)
{
    const auto buffer_exp = xraymc::Serializer::read(filename);
    if (!buffer_exp) {
        std::cout << "Could not read saved file: " << filename << " Exiting" << std::endl;
        return;
    }

    std::cout << "Loading simulation from " << filename << "...";

    // Getting the buffer
    auto buffer = std::span { buffer_exp.value() };

    // Loading the beam from file
    auto beam_buffer = xraymc::Serializer::getEmptyBuffer();
    auto beam_name = xraymc::Serializer::getNameIDTemplate();
    buffer = xraymc::Serializer::deserializeItem(beam_name, beam_buffer, buffer);
    // Create the beam object
    auto beam_opt = Beam::deserialize(beam_buffer);
    Beam* beam = nullptr;
    if (beam_opt)
        beam = &beam_opt.value();
    else
        beam = nullptr;

    // Loading the world from file
    auto world_buffer = xraymc::Serializer::getEmptyBuffer();
    auto world_name = xraymc::Serializer::getNameIDTemplate();
    buffer = xraymc::Serializer::deserializeItem(world_name, world_buffer, buffer);
    auto world_opt = World::deserialize(world_buffer);
    if (!world_opt) {
        std::cout << "Could not reconstruct world from: " << filename << " ...exiting" << std::endl;
        return;
    }
    auto& world = world_opt.value();
    std::cout << "Done\n";

    std::cout << "Generate some images of the setup...";
    visualizeWorld(world, beam, "show", 1.5);
    std::cout << "Done\n";

    // If the world contains a doctor, we find the max skin dose
    double max_dose = 0;
    auto doctor = std::get_if<TetMesh>(world.getItemPointerFromName("Doctor"));
    if (doctor) {
        const auto& outerTetIdx = doctor->outerContourTetrahedronIndices();
        for (std::size_t i = 0; i < outerTetIdx.size(); ++i)
            max_dose = std::max(max_dose, doctor->doseScored(outerTetIdx[i]).dose());
    }

    // Visualize the skin dose
    std::cout << "Generate some images of the scattered radiation...";
    visualizeWorld(world, beam, "doseLog", 1.5, max_dose / 10.0, true);
    visualizeWorld(world, beam, "dose", 1.5, max_dose / 10.0, false);
    std::cout << "Done\n";
}

void runSimulation(const std::string& filename)
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

    auto carm = std::get_if<CarmType>(world.getItemPointerFromName("C-Arm"));
    if (carm) {
        carm->setPrimaryAngleDeg(10);
        carm->setSecondaryAngleDeg(15);
        beam.setPosition(carm->beamSourcePos());
        beam.setDirectionCosines(carm->beamCosines());
    }
    // Rebuild after rotating the C-arm so the BVH acceleration structures
    // reflect the updated triangle positions.
    world.build();

    // Run the simulation
    std::cout << "Running simulation with " << beam.numberOfParticles() << " particles\n";
    xraymc::Transport transport;
    transport.runConsole(world, beam);

    // Save the simulation to a file
    std::cout << "Saving simulation in " << filename << "...";
    auto saveBuffer = xraymc::Serializer::getEmptyBuffer();
    xraymc::Serializer::serializeItem(beam, saveBuffer);
    xraymc::Serializer::serializeItem(world, saveBuffer);
    xraymc::Serializer::write(filename, saveBuffer);
    std::cout << "Done\n";
}

int main()
{
    // Lets test if we have already done a simulation
    // if not we run it, to force a new simulation set force_new_simulation to true

    const std::string savefile = "savefile.xr";

    bool force_new_simulation = false;
    if (!std::filesystem::exists(savefile) || force_new_simulation) {
        runSimulation(savefile);
    }

    loadAndShowSimulation(savefile);
    return EXIT_SUCCESS;
}