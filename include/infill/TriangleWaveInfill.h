// Copyright (c) 2026 UltiMaker
// CuraEngine is released under the terms of the AGPLv3 or higher.

#ifndef INFILL_TRIANGLE_WAVE_INFILL_H
#define INFILL_TRIANGLE_WAVE_INFILL_H

#include "geometry/OpenLinesSet.h"
#include "utils/Coord_t.h"

namespace cura
{
class Shape;

class TriangleWaveInfill
{
public:
    /*!
     * Generate the triangle-wave infill pattern within a certain outline.
     *
     * The pattern consists of parallel continuous paths, each of which oscillates
     * perpendicular to the path direction as a triangle wave (a piecewise-linear
     * zigzag, see https://en.wikipedia.org/wiki/Triangle_wave). All waves share the
     * same phase, so adjacent paths stay at a constant distance from each other and
     * the material density matches that of the 'lines' pattern with the same
     * line_distance. Being path-based, each printed line is one long continuous
     * extrusion instead of a texture of disconnected segments.
     *
     * \param result_lines Output variable to store the resulting polyline segments in.
     * \param zig_zaggify Whether to connect the wave ends along the outline, forming
     * one single polyline or at least very few interruptions in the material flow.
     * \param line_distance Distance between adjacent waves, measured perpendicular to
     * the wave flanks. This determines the density of the pattern.
     * \param in_outline The outline in which to print the pattern.
     * \param fill_angle The angle (in degrees) of the general direction of the waves.
     */
    static void generateTotalTriangleWaveInfill(OpenLinesSet& result_lines, bool zig_zaggify, coord_t line_distance, const Shape& in_outline, const double fill_angle);
};

} // namespace cura

#endif // INFILL_TRIANGLE_WAVE_INFILL_H
