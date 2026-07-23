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
 *  1. Region division and reference axis: the fill regions are the same per-part infill areas
 *     the native grid-type patterns fill. Per region, the medial axis (centerline) is computed
 *     and used as the longitudinal reference line; all peak positions are parameterized by the
 *     arc length along this axis, so the wave follows the tangent direction of the curved
 *     surface instead of absolute coordinates.
 *  2. Wall contact with one apex angle: every wave peak/valley reaches from the axis
 *     perpendicularly (the perpendicular bisector at each apex is perpendicular to the medial
 *     axis) to the wall it faces and contacts it, with peaks and valleys alternating between
 *     the two sides. The apex spacing scales with the local cross-section width, so the apex
 *     ANGLE is the same everywhere ("等顶角"); the average spacing -- and thereby the angle --
 *     is set by the infill density. Per independent region, one single continuous open
 *     polyline is generated (no crossing grid, no self-intersection, no travel moves).
 *  3. Phase lock and mapping reuse: the wave of a layer is derived from the layer below. Every
 *     single peak of the lower layer is projected (along the surface normal) onto the current
 *     axis, so direction, fill angle and phase are inherited unchanged; no layer-to-layer
 *     lateral offset is ever applied. The wave only shifts along with the surface itself, so
 *     the lines of consecutive layers overlap and interlock even when the cross sections of a
 *     doubly-curved part drift sideways between layers.
 *  4. The first and last line of each path keep their own direction and are prolonged until
 *     they contact the wall, so both path ends are anchored on the innermost wall as well.
 *  5. The centerline is only kept inside the effective MAIN fill region: skeleton branches
 *     into small appendage areas attached to the boundary (rib pockets and the like) are
 *     removed, truncating the reference axis at the junction with the main region.
 *  6. A region contains nothing but the continuous triangle wave or its tip-truncated form:
 *     degenerate spikes (flanks folding onto each other) are cut off, and any stretch that
 *     would leave the region is replaced by the walk along the region boundary, connecting
 *     break points along the wall without crossing holes or non-printable areas.
 */
class SurfaceWaveFillProvider
{
public:
    /*!
     * Sequentially build the waves of all layers, bottom up, each layer inheriting modulus and
     * phase from the layer below.
     * \param layer_outlines The infill areas of all layers (and all their parts).
     * \param line_distance Half the wave period: the model-global fixed cell size is
     *        2 * line_distance, tying the wave's fold angle to the infill density.
     * \param line_width The width with which the pattern lines will be extruded.
     * \param wall_clearance Extra distance to keep from the boundary of the given outlines (for
     *        the extra infill walls that are subtracted from the areas at gcode time), so the
     *        wave is generated on the same region the pattern will actually fill.
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
