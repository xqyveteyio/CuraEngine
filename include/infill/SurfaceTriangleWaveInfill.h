// Copyright (c) 2026 UltiMaker
// CuraEngine is released under the terms of the AGPLv3 or higher.

#ifndef INFILL_SURFACE_TRIANGLE_WAVE_INFILL_H
#define INFILL_SURFACE_TRIANGLE_WAVE_INFILL_H

#include <vector>

#include "geometry/OpenLinesSet.h"
#include "utils/Coord_t.h"

namespace cura
{
class Shape;

/*!
 * Surface-conformal triangle wave infill ("曲面随形参数化三角折线填充").
 *
 * Designed for curved / freeform (architectural) members: arched walls, shells, members with a
 * varying cross section. A triangle wave is generated relative to the geometry of the part
 * instead of on a fixed XY grid:
 *
 *  1. Reference axis: per part, the medial axis (centerline) is computed and used as the
 *     longitudinal reference line. All wave periods and peak positions are parameterized by the
 *     arc length along this axis, so the wave follows the tangent direction of the curved
 *     surface instead of absolute coordinates.
 *  2. Modulus scaling: the axis is split into segments of equal modulus (fixed number of wave
 *     cells per segment). When the part narrows/widens or shortens/lengthens between layers,
 *     the peak COUNT stays fixed and the period/amplitude scale proportionally, so every peak
 *     keeps matching the corresponding groove of the layer below.
 *  3. Phase lock and mapping reuse: the wave of a layer is derived from the layer below. The
 *     axis orientation, segment boundaries and every single peak of the lower layer are
 *     projected onto the current axis, so the phase is inherited continuously; no global
 *     layer-to-layer lateral offset is ever applied. Peaks of consecutive layers therefore
 *     stack up along the curved surface (normal-direction correspondence) and interlock.
 *  4. Segmented modulus for strongly curved surfaces: the axis is cut where the accumulated
 *     turning angle becomes large; each segment carries its own fixed cell count with a smooth
 *     (continuous polyline) transition at the boundaries.
 */
class SurfaceWaveFillProvider
{
public:
    /*!
     * Sequentially build the waves of all layers, bottom up, each layer inheriting modulus and
     * phase from the layer below.
     * \param layer_outlines The infill areas of all layers (and all their parts).
     * \param line_distance The nominal distance along the axis between two wave cells (the wave
     *        period on the layer where a part first appears; afterwards it scales with the part).
     * \param line_width The width with which the pattern lines will be extruded.
     * \param wall_clearance Extra distance to keep from the boundary of the given outlines (for
     *        the extra infill walls and the minimum wall line width), constraining the wave.
     */
    SurfaceWaveFillProvider(const std::vector<Shape>& layer_outlines, coord_t line_distance, coord_t line_width, coord_t wall_clearance);

    /*!
     * Clip the pre-computed wave of one layer to the actual infill outline.
     * \param result_lines Output variable to store the resulting polylines in.
     * \param in_outline The infill area of this layer.
     * \param layer_idx The layer number, indexing the outlines given to the constructor.
     */
    void generate(OpenLinesSet& result_lines, const Shape& in_outline, size_t layer_idx) const;

    /*!
     * Stand-alone single-layer generation, without cross-layer inheritance. Used as fallback
     * when no provider is available (e.g. for support).
     */
    static void generateSingleLayer(OpenLinesSet& result_lines, const Shape& in_outline, coord_t line_distance, coord_t line_width);

private:
    coord_t line_width_;
    std::vector<OpenLinesSet> layer_waves_; //!< One pre-computed wave set per layer.
};

} // namespace cura

#endif // INFILL_SURFACE_TRIANGLE_WAVE_INFILL_H
