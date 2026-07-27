// Copyright (c) 2026 UltiMaker
// CuraEngine is released under the terms of the AGPLv3 or higher

#ifndef INFILL_MEDIAL_ZIGZAG_H
#define INFILL_MEDIAL_ZIGZAG_H

#include <vector>

#include "geometry/OpenLinesSet.h"
#include "utils/Coord_t.h"

namespace cura
{

class Shape;
class SliceMeshStorage;

/*!
 * Generate the medial zigzag lines for a single, stand-alone layer region (used e.g. for support,
 * where no cross-layer structure is available). The region itself acts as its own reference
 * cross-section: nodes are distributed along the walls on either side of the medial axis (uniformly
 * in each wall's own arc length, both walls sharing the same period count, half a period out of
 * phase) and connected into one alternating polyline per part. Parts without a usable medial axis
 * fall back to fixed cutting planes perpendicular to \p plane_angle. Lines are clipped against the
 * region, so holes are never crossed; a line leaving and re-entering the solid region continues as a
 * new polyline.
 */
void generateMedialZigzagLines(const Shape& region, coord_t line_distance, double plane_angle, coord_t plane_shift, OpenLinesSet& result);

/*!
 * Pre-computed medial zigzag infill for a whole mesh.
 *
 * The nodes are generated once, on the layer with the largest infill cross-section:
 *  - the medial axis of each part is computed there,
 *  - the common period count N follows from the axis length and the reference spacing (twice the
 *    line distance), rounded up so the spacing never exceeds the reference,
 *  - both walls receive nodes uniformly in their own arc length (a longer outer wall thus gets a
 *    larger actual spacing than a shorter inner wall), half a period out of phase in the normalized
 *    arc-length parameter, connected alternately into one polyline.
 *
 * Every node then defines a vertical reference plane, spanned by the vertical direction and the
 * perpendicular from the node to the medial axis. Going layer by layer outwards from the reference
 * cross-section, each node is mapped to the intersection of its plane with the next layer's wall,
 * considering only intersections on the node's side of the axis; the plane is then re-derived from
 * the node's new position and the local axis direction, so it follows the walls when the
 * cross-section drifts or twists along the height. On constant cross-sections the planes never
 * change, so the zigzag is vertically aligned by construction. When a node's plane doesn't
 * intersect the wall of a layer, the node is simply invalid there and skipped (never clamped to a
 * wall end or placed outside).
 */
class MedialZigzagGenerator
{
public:
    /*! Compute the zigzag lines of all layers of the mesh. */
    MedialZigzagGenerator(const SliceMeshStorage& mesh);

    /*! The pre-computed zigzag lines of the given layer (empty for out-of-range layers). */
    const OpenLinesSet& getLinesForLayer(size_t layer_nr) const;

private:
    std::vector<OpenLinesSet> lines_per_layer_;
};

} // namespace cura

#endif // INFILL_MEDIAL_ZIGZAG_H
