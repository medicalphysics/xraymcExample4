# xraymcExample4

## About

This example demonstrates how to implement a **user-defined world object** for the [XRayMClib](https://github.com/medicalphysics/XRayMClib) C++ library to simulate radiation transport in custom geometries.

## User-Defined World Objects

While XRayMC supports importing tetrahedral and triangulated meshes for arbitrary geometry, sometimes it's more convenient or efficient to implement custom geometry directly in code. User-defined world objects allow you to:

- Define geometric shapes that are difficult to represent as meshes (e.g., spheres, cylinders, toroidal shapes)
- Implement custom ray-geometry intersection logic
- Add specialized particle transport and scoring within the geometry
- Avoid the computational overhead of complex mesh representations

## The HollowSphere Example

This example implements a **HollowSphere** class that serves as a custom world object. The hollow sphere is a spherical shell made of concrete, featuring:

- **Geometric intersection**: Ray-sphere intersection calculations to track particle movement through the hollow shell
- **Particle transport**: Monte Carlo transport through the sphere material with energy deposition and dose scoring
- **Visualization support**: Surface normal calculations for rendering the geometry in 3D visualizations
- **Configurable parameters**: Adjustable outer radius, wall thickness, and material properties

### Key Components

- `HollowSphere` class (`include/hollowsphere.hpp`): Implements the custom world geometry. **Must satisfy the `xraymc::WorldItemType` concept**, which requires the following methods:
  | Method | Description |
  |---|---|
  | `translate(vec)` | Move the item by a displacement vector (cm). |
  | `center()` | Return the world-space centre as `array<double,3>`. |
  | `AABB()` | Return the axis-aligned bounding box as `array<double,6>` — `[xmin, ymin, zmin, xmax, ymax, zmax]`. |
  | `intersect(p)` | Compute the next ray–geometry intersection for particle `p`; returns a `WorldIntersectionResult` with the distance and a flag indicating whether the particle origin is inside the material. |
  | `intersectVisualization(p)` | Compute the next ray–geometry intersection for particle `p`; returns a templated `VisualizationIntersectionResult` with the distance and a flag indicating whether the particle origin is inside the material and surface normal vector. This method is not neccesary for Monte Carlo transport but allows the object to be visualized by the internal scene renderer in xraymc. |
  | `energyScored(index)` | Return the `EnergyScore` accumulator at the given tally index. |
  | `doseScored(index)` | Return the `DoseScore` accumulator at the given tally index. |
  | `clearEnergyScored()` | Reset all energy-score accumulators to zero. |
  | `clearDoseScored()` | Reset all dose-score accumulators to zero. |
  | `addEnergyScoredToDoseScore(f)` | Convert accumulated energy scores to dose using calibration factor `f` (e.g. keV/g → mGy scaling). |
  | `transport(p, state)` | Transport particle `p` through the item via Monte Carlo stepping until it exits or is absorbed, scoring imparted energy along the way. |
- `main.cpp`: Demonstrates how to create a world with the custom HollowSphere object and run a Monte Carlo X-ray simulation
- Integration with a lead-lined room to show multiple world objects working together

## Building and Running

Follow the standard CMake build procedure. The example will generate a PNG visualization of the geometry and run a Monte Carlo simulation with energy and dose scoring in the hollow sphere.



