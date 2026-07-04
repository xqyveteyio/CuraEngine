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
     * The pattern is a single continuous triangle wave (a piecewise-linear zigzag,
     * see https://en.wikipedia.org/wiki/Triangle_wave) whose peaks and troughs span
     * the full extent of the region: the path oscillates diagonally between the two
     * opposite sides of the area, and the clipped tips are connected along the
     * boundary walls. This yields one single uninterrupted extrusion path per region,
     * which is the desired behaviour for e.g. concrete/cement printing.
     *
     * \param result_lines Output variable to store the resulting polyline segments in.
     * \param zig_zaggify Whether to connect the wave tips along the outline, forming
     * one single continuous path per region. Without it the wave falls apart into
     * separate V-shaped pieces wherever it is clipped by the boundary.
     * \param line_distance Horizontal distance between two successive flanks of the
     * wave, i.e. the spacing of the triangular teeth. This determines the density.
     * \param in_outline The outline in which to print the pattern.
     * \param fill_angle The angle (in degrees) of the general direction of the wave.
     */
    static void generateTotalTriangleWaveInfill(OpenLinesSet& result_lines, bool zig_zaggify, coord_t line_distance, const Shape& in_outline, const double fill_angle);
};

} // namespace cura

#endif // INFILL_TRIANGLE_WAVE_INFILL_H
