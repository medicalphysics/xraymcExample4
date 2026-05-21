

#pragma once // Include guard

// Include the xraumc library
#include "xraymc/xraymc.hpp"

/**
 * @brief A hollow sphere geometry for Monte Carlo X-ray transport simulation.
 *
 * Models a spherical shell made of Portland concrete with a configurable wall
 * thickness. Supports particle transport, energy/dose scoring, and visualization
 * intersection queries as required by the xraymc world geometry interface.
 *
 * This class satisfies the xraymc::WorldItemType concept, which requires:
 *   - translate(vec)                  — move the item by a displacement vector
 *   - center()                        — return the world-space centre as array<double,3>
 *   - AABB()                          — return the axis-aligned bounding box as array<double,6>
 *   - intersect(p)                    — return a WorldIntersectionResult for particle p
 *   - energyScored(index)             — return the EnergyScore accumulator at index
 *   - doseScored(index)               — return the DoseScore accumulator at index
 *   - clearEnergyScored()             — reset all energy-score accumulators
 *   - clearDoseScored()               — reset all dose-score accumulators
 *   - addEnergyScoredToDoseScore(f)   — convert energy scores to dose with calibration factor f
 *   - transport(p, state)             — transport particle p through the item
 */
class HollowSphere {
public:
    /**
     * @brief Constructs a HollowSphere with Portland Concrete as the wall material.
     *
     * Default outer radius is 7 cm and wall thickness is 0.5 cm, centered at the origin.
     */
    HollowSphere()
        : m_material(xraymc::Material<>::byNistName("Concrete, Portland").value())
    {
        m_density = xraymc::NISTMaterials::density("Concrete, Portland");
        // Optionally the hollow halfsphere could be made on any material, example water:
        // m_material = xraymc::Material<>::byChemicalFormula("H2O").value();
    }

    /**
     * @brief Sets the wall thickness of the sphere. Optional
     * @param cm Wall thickness in centimetres. Must be positive and less than the outer radius.
     */
    void setWallThickness(double cm)
    {
        // Guard: thickness must be positive and less than the outer radius so the
        // inner radius (m_radius - cm) stays above zero.
        if (cm > 0 && cm < m_radius)
            m_thickness = cm;
    }

    /**
     * @brief Translates the sphere center by the given vector.
     * @param vector Displacement vector [dx, dy, dz] in centimetres.
     */
    void translate(const std::array<double, 3>& vector)
    {
        m_origin = xraymc::vectormath::add(m_origin, vector);
    }

    /**
     * @brief Returns the center position of the sphere.
     * @return Reference to the [x, y, z] origin array in centimetres.
     */
    const std::array<double, 3>& center() const
    {
        return m_origin;
    }

    /**
     * @brief Computes the axis-aligned bounding box of the sphere.
     * @return Array [xmin, ymin, zmin, xmax, ymax, zmax] in centimetres.
     */
    std::array<double, 6> AABB() const
    {
        // construct the Axis Aligned Bounding Box as an array [xmin, ymin, zmin, xmax, ymax, zmax]
        std::array<double, 6> aabb = {
            m_origin[0] - m_radius,
            m_origin[1] - m_radius,
            m_origin[2] - m_radius,
            m_origin[0] + m_radius,
            m_origin[1] + m_radius,
            m_origin[2] + m_radius
        };
        return aabb;
    }

    /**
     * @brief Computes the next ray-geometry intersection for a particle.
     *
     * Determines whether the particle is outside, inside the wall, or inside the
     * hollow cavity, and returns the corresponding intersection distance and
     * inside-item flag used by the transport loop.
     *
     * @param p Particle satisfying the xraymc::ParticleType concept.
     * @return WorldIntersectionResult with distance and containment state.
     */
    xraymc::WorldIntersectionResult intersect(const xraymc::ParticleType auto& p) const
    {
        // Cheap AABB pre-test: skip the more expensive sphere math if the ray
        // cannot possibly hit the bounding box.
        xraymc::WorldIntersectionResult intersection = xraymc::basicshape::AABB::intersect(p, AABB());
        if (intersection.valid()) {

            // Intersect both the outer shell surface and the inner (hollow) surface.
            // rayOriginIsInsideItem == true means the ray origin is inside that sphere.
            xraymc::WorldIntersectionResult sphere_outer = xraymc::basicshape::sphere::intersect(p, m_origin, m_radius);
            xraymc::WorldIntersectionResult sphere_inner = xraymc::basicshape::sphere::intersect(p, m_origin, m_radius - m_thickness);

            // Three spatial regions, mapped to the transport flag rayOriginIsInsideItem:
            //   Region 1 — outside the outer sphere:  next event = entry into outer shell, flag = false (not in material)
            //   Region 2 — inside outer, outside inner (inside the concrete wall): next event = inner surface, flag = true  (in material)
            //   Region 3 — inside both spheres (hollow cavity): next event = inner surface exit, flag = false (not in material)
            if (sphere_outer.rayOriginIsInsideItem == false) {
                // Region 1: particle is outside the object entirely.
                intersection = sphere_outer;
            } else {
                if (sphere_inner.rayOriginIsInsideItem == false) {
                    // Region 2: particle is inside the concrete wall — set flag so transport scores interactions.
                    intersection = sphere_inner;
                    intersection.rayOriginIsInsideItem = true;
                } else {
                    // Region 3: particle is in the hollow cavity — no material, so flag = false.
                    intersection = sphere_inner;
                    intersection.rayOriginIsInsideItem = false;
                }
            }
        }
        return intersection;
    }

    /**
     * @brief Computes an intersection result with surface normal for visualization.
     *
     * Calls intersect() and additionally calculates the outward surface normal at the
     * hit point, required by renderers to shade the geometry correctly.
     *
     * @tparam U Current item type used by the visualization backend.
     * @param p Particle satisfying the xraymc::ParticleType concept.
     * @return VisualizationIntersectionResult containing distance, inside flag, and normal.
     */
    template <typename U>
    xraymc::VisualizationIntersectionResult<U> intersectVisualization(const xraymc::ParticleType auto& p) const
    {
        xraymc::VisualizationIntersectionResult<U> vintersect;
        auto res = intersect(p);
        if (res.valid()) {
            vintersect.intersection = res.intersection;
            vintersect.rayOriginIsInsideItem = res.rayOriginIsInsideItem;

            // Outward surface normal: vector from sphere centre to hit point, normalised.
            // Works for both the outer and inner surfaces since both are centred at m_origin.
            auto hit_pos = xraymc::vectormath::add(p.pos, xraymc::vectormath::scale(p.dir, res.intersection));
            vintersect.normal = xraymc::vectormath::normalized(xraymc::vectormath::subtract(hit_pos, m_origin));
            vintersect.intersectionValid = true;
        }
        return vintersect;
    }

    /**
     * @brief Returns the accumulated energy score for a given tally index.
     * @param index Tally index (currently unused; the sphere has a single scorer).
     * @return Const reference to the EnergyScore accumulator.
     */
    const xraymc::EnergyScore& energyScored(std::uint64_t index = 0) const
    {
        return m_energyScored;
    }

    /**
     * @brief Returns the accumulated dose score for a given tally index.
     * @param index Tally index (currently unused; the sphere has a single scorer).
     * @return Const reference to the DoseScore accumulator.
     */
    const xraymc::DoseScore& doseScored(std::uint64_t index = 0) const
    {
        return m_doseScored;
    }

    /// Resets the energy scorer to zero.
    void clearEnergyScored()
    {
        m_energyScored.clear();
    }

    /// Resets the dose scorer to zero.
    void clearDoseScored()
    {
        m_doseScored.clear();
    }

    /**
     * @brief Converts accumulated energy to dose and adds it to the dose scorer.
     *
     * Computes the shell volume from the outer and inner radii, then calls
     * DoseScore::addScoredEnergy with the wall material density.
     *
     * @param factor Unit conversion factor (e.g. keV/g -> mGy scaling).
     */
    void addEnergyScoredToDoseScore(double factor)
    {
        // Shell volume = (4π/3)(R³ − r³), where R is the outer radius and r = R − thickness is the inner radius.
        // This is simply the volume of the full outer sphere minus the hollow interior.
        const double R3 = m_radius * m_radius * m_radius;
        const double r3 = (m_radius - m_thickness) * (m_radius - m_thickness) * (m_radius - m_thickness);
        const double volume = 4.0 * std::numbers::pi_v<double> / 3.0 * (R3 - r3);
        m_doseScored.addScoredEnergy(m_energyScored, volume, m_density, factor);
    }

    /**
     * @brief Transports a particle through the sphere wall until it exits or is absorbed.
     *
     * Assumes the particle has already entered the geometry. Repeatedly samples
     * a free path from an exponential distribution, performs photon interactions,
     * scores imparted energy, and advances the particle until it crosses the boundary
     * or is absorbed.
     *
     * @param p Particle satisfying the xraymc::ParticleType concept (modified in place).
     * @param state Random number generator state.
     */
    void transport(xraymc::ParticleType auto& p, xraymc::RandomState& state)
    {
        // The transport method always assume that the particle has just entered the geometry.
        // rayOriginIsInsideItem == true signals that the particle is inside the concrete wall
        // and must be stepped; false means it is in a void region or has left the object.

        xraymc::WorldIntersectionResult intersection = intersect(p);
        while (intersection.rayOriginIsInsideItem) {
            // Sample a free path from the exponential distribution using the inverse-CDF method:
            //   step = -ln(U) / (μ_total · ρ),  where U ~ Uniform(0,1).
            // This is the standard analog Monte Carlo path-length sampling under Beer-Lambert.
            xraymc::AttenuationValues attenuation = m_material.attenuationValues(p.energy);
            const double step = -std::log(state.randomUniform()) / (attenuation.sum() * m_density);

            if (step < intersection.intersection) {
                // Interaction occurs before the boundary: advance particle to the interaction site.
                p.translate(step);
                xraymc::interactions::InteractionResult res = xraymc::interactions::interact(attenuation, p, m_material, state);
                // Accumulate kinetic energy transferred to the medium (kerma proxy for dose scoring).
                m_energyScored.scoreEnergy(res.energyImparted);

                if (res.particleAlive) {
                    // Scatter changed direction and/or energy; recompute the next boundary distance.
                    intersection = intersect(p);
                } else {
                    // Photoelectric absorption: particle is gone, stop the loop.
                    intersection.rayOriginIsInsideItem = false;
                }
            } else {
                // Free path exceeds the distance to the boundary: move the particle exactly to the
                // surface (border_translate applies a small push to avoid floating-point re-entry).
                p.border_translate(intersection.intersection);
                intersection.rayOriginIsInsideItem = false;
            }
        }
    }

private:
    std::array<double, 3> m_origin = { 0, 0, 0 };
    double m_radius = 7.0; // outer radius of the sphere in cm
    double m_thickness = 0.5; // wall thickness in cm
    double m_density = 1.0;
    xraymc::Material<> m_material;
    xraymc::EnergyScore m_energyScored;
    xraymc::DoseScore m_doseScored;
};
