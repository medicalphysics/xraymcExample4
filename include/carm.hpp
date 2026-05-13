/*This file is part of xraymcCarm.

xraymcCarm is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

xraymcCarm is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with xraymcCarm. If not, see < https://www.gnu.org/licenses/>.

Copyright 2024 Erlend Andersen
*/

#pragma once

#include <xraymc/world/worlditems/triangulatedmesh.hpp>

#include <array>
#include <string>

// C-arm fluoroscopy unit modelled as a triangulated STL mesh.
// The STL file must have its X-ray source at the origin and the isocenter on
// the +Z axis when both gantry angles are zero.
template <int NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class Carm {
public:
    // C-arm frames are carbon-fibre composite; material is fixed to carbon regardless of the STL content.
    Carm(const std::string& path, std::size_t max_three_dept = 4)
        : m_mesh(path, max_three_dept)
    {
        const auto carbon = xraymc::Material<NMaterialShells>::byZ(6).value();
        const auto carbon_dens = xraymc::AtomHandler::Atom(6).standardDensity;
        m_mesh.setMaterial(carbon, carbon_dens);
    }

    // Deserialization constructor: accepts a pre-built mesh (e.g. loaded from a save file).
    Carm(const xraymc::TriangulatedMesh<NMaterialShells, LOWENERGYCORRECTION>& mesh)
        : m_mesh(mesh)
    {
        const auto carbon = xraymc::Material<NMaterialShells>::byZ(6).value();
        const auto carbon_dens = xraymc::AtomHandler::Atom(6).standardDensity;
        m_mesh.setMaterial(carbon, carbon_dens);
    }

    // The mesh does not store source/isocenter positions, so they must be
    // scaled explicitly to stay consistent with the geometry.
    void scale(double s)
    {
        m_mesh.scale(s);
        for (std::size_t i = 0; i < 3; ++i) {
            m_source_pos[i] *= s;
            m_isocenter[i] *= s;
        }
    }

    // World-space rotation (e.g. to tilt the whole unit into the room).
    // All five tracked vectors are rotated so that subsequent primary/secondary
    // angle calls remain correct relative to the new orientation.
    void rotate(const std::array<double, 3>& axis, double angle)
    {
        m_mesh.rotate(angle, axis);

        m_source_pos = xraymc::vectormath::rotate(m_source_pos, axis, angle);
        m_isocenter = xraymc::vectormath::rotate(m_isocenter, axis, angle);
        m_xAxis = xraymc::vectormath::rotate(m_xAxis, axis, angle);
        m_yAxis = xraymc::vectormath::rotate(m_yAxis, axis, angle);
        m_beam_cosines[0] = xraymc::vectormath::rotate(m_beam_cosines[0], axis, angle);
        m_beam_cosines[1] = xraymc::vectormath::rotate(m_beam_cosines[1], axis, angle);
    }

    // Beam cosines and isocenter travel with the mesh; source position follows.
    void translate(const std::array<double, 3>& translation)
    {
        m_mesh.translate(translation);
        m_source_pos = xraymc::vectormath::add(m_source_pos, translation);
        m_isocenter = xraymc::vectormath::add(m_isocenter, translation);
    }

    void setPrimaryAngleDeg(double angle) { setPrimaryAngle(angle * xraymc::DEG_TO_RAD<>()); }
    void setSecondaryAngleDeg(double angle) { setSecondaryAngle(angle * xraymc::DEG_TO_RAD<>()); }

    // Rotates the arm by the *delta* from the currently stored angle so that
    // repeated calls are additive and the mesh is never rotated from a stale pose.
    // The source orbits the isocenter: only the displacement vector is rotated,
    // keeping the source-to-isocenter distance constant.
    void setPrimaryAngle(double new_angle)
    {
        const double angle = new_angle - m_primary_angle;
        m_primary_angle = new_angle;

        m_mesh.rotate(angle, m_xAxis, m_isocenter);

        auto d = xraymc::vectormath::subtract(m_source_pos, m_isocenter);
        auto d_fin = xraymc::vectormath::rotate(d, m_xAxis, angle);
        m_source_pos = xraymc::vectormath::add(m_isocenter, d_fin);
        m_beam_cosines[0] = xraymc::vectormath::rotate(m_beam_cosines[0], m_xAxis, angle);
        m_beam_cosines[1] = xraymc::vectormath::rotate(m_beam_cosines[1], m_xAxis, angle);
    }
    void setSecondaryAngle(double new_angle)
    {
        const double angle = new_angle - m_secondary_angle;
        m_secondary_angle = new_angle;

        m_mesh.rotate(angle, m_yAxis, m_isocenter);

        auto d = xraymc::vectormath::subtract(m_source_pos, m_isocenter);
        auto d_fin = xraymc::vectormath::rotate(d, m_yAxis, angle);
        m_source_pos = xraymc::vectormath::add(m_isocenter, d_fin);
        m_beam_cosines[0] = xraymc::vectormath::rotate(m_beam_cosines[0], m_yAxis, angle);
        m_beam_cosines[1] = xraymc::vectormath::rotate(m_beam_cosines[1], m_yAxis, angle);
    }
    const std::array<double, 3>& center() const
    {
        return m_isocenter;
    }
    std::array<double, 6> AABB() const
    {
        return m_mesh.AABB();
    }
    xraymc::WorldIntersectionResult intersect(const xraymc::ParticleType auto& p) const
    {
        return m_mesh.intersect(p);
    }
    template <typename U>
    xraymc::VisualizationIntersectionResult<U> intersectVisualization(const xraymc::ParticleType auto& p) const
    {
        return m_mesh.template intersectVisualization<U>(p);
    }
    auto energyScored(std::size_t index = 0)
    {
        return m_mesh.energyScored(index);
    }
    auto doseScored(std::size_t index = 0) const
    {
        return m_mesh.doseScored(index);
    }
    void clearDoseScored()
    {
        m_mesh.clearDoseScored();
    }
    void clearEnergyScored()
    {
        m_mesh.clearEnergyScored();
    }
    void addEnergyScoredToDoseScore(double factor)
    {
        m_mesh.addEnergyScoredToDoseScore(factor);
    }
    void transport(xraymc::ParticleType auto& p, xraymc::RandomState& state)
    {
        m_mesh.transport(p, state);
    }

    std::array<double, 3> beamSourcePos() const
    {
        return m_source_pos;
    }
    // Two orthogonal unit vectors that define the detector plane orientation,
    // passed directly to DXBeam::setDirectionCosines().
    const std::array<std::array<double, 3>, 2>& beamCosines() const
    {
        return m_beam_cosines;
    }

    std::array<double, 3> isoCenter() const
    {
        return m_isocenter;
    }

    // 32-byte identifier written at the start of every serialized buffer.
    // Encodes the template parameters so deserialize() can reject a buffer that
    // was saved with a different physics configuration (different NMaterialShells
    // or LOWENERGYCORRECTION would produce incompatible cross-section tables).
    constexpr static std::array<char, 32> magicID()
    {
        std::string name = "CarmSTL1" + std::to_string(LOWENERGYCORRECTION) + std::to_string(NMaterialShells);
        name.resize(32, ' ');
        std::array<char, 32> k;
        std::copy(name.cbegin(), name.cend(), k.begin());
        return k;
    }

    static bool validMagicID(std::span<const char> data)
    {
        if (data.size() < 32)
            return false;
        const auto id = magicID();
        return std::search(data.cbegin(), data.cbegin() + 32, id.cbegin(), id.cend()) == data.cbegin();
    }

    std::vector<char> serialize() const
    {
        using Ser = xraymc::Serializer;

        auto buffer = Ser::getEmptyBuffer();
        Ser::serializeItem(m_mesh, buffer);

        Ser::serialize(m_source_pos, buffer);
        Ser::serialize(m_isocenter, buffer);
        Ser::serialize(m_xAxis, buffer);
        Ser::serialize(m_yAxis, buffer);
        Ser::serialize(m_beam_cosines[0], buffer);
        Ser::serialize(m_beam_cosines[1], buffer);
        Ser::serialize(m_primary_angle, buffer);
        Ser::serialize(m_secondary_angle, buffer);

        return buffer;
    }

    static std::optional<Carm<NMaterialShells, LOWENERGYCORRECTION>> deserialize(std::span<const char> buffer)
    {
        using Ser = xraymc::Serializer;

        auto mesh_buffer = Ser::getEmptyBuffer();
        auto name = Ser::getNameIDTemplate();
        buffer = Ser::deserializeItem(name, mesh_buffer, buffer);
        auto mesh_opt = xraymc::TriangulatedMesh<NMaterialShells, LOWENERGYCORRECTION>::deserialize(mesh_buffer);
        if (!mesh_opt)
            return std::nullopt;
        Carm<NMaterialShells, LOWENERGYCORRECTION> item(mesh_opt.value());

        buffer = Ser::deserialize(item.m_source_pos, buffer);
        buffer = Ser::deserialize(item.m_isocenter, buffer);
        buffer = Ser::deserialize(item.m_xAxis, buffer);
        buffer = Ser::deserialize(item.m_yAxis, buffer);
        buffer = Ser::deserialize(item.m_beam_cosines[0], buffer);
        buffer = Ser::deserialize(item.m_beam_cosines[1], buffer);
        buffer = Ser::deserialize(item.m_primary_angle, buffer);
        buffer = Ser::deserialize(item.m_secondary_angle, buffer);

        return std::make_optional(item);
    }

private:
    xraymc::TriangulatedMesh<NMaterialShells, LOWENERGYCORRECTION> m_mesh;
    std::array<double, 3> m_source_pos = { 0, 0, 0 };
    std::array<double, 3> m_isocenter = { 0, 0, 55 }; // 55 cm SID in the STL's local frame
    // Primary rotates around m_xAxis, secondary around m_yAxis; both axes are
    // updated by rotate() so they always match the mesh's current orientation.
    std::array<double, 3> m_xAxis = { 1, 0, 0 };
    std::array<double, 3> m_yAxis = { 0, 1, 0 };
    std::array<std::array<double, 3>, 2> m_beam_cosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };
    double m_primary_angle = 0;   // radians, current accumulated angle
    double m_secondary_angle = 0;
};
