// Copyright (c) 2026 UltiMaker
// CuraEngine is released under the terms of the AGPLv3 or higher

#ifndef INFILL_MEDIAL_ZIGZAG_H
#define INFILL_MEDIAL_ZIGZAG_H

#include <vector>

#include "geometry/OpenLinesSet.h"
#include "geometry/Point2LL.h"
#include "utils/Coord_t.h"

namespace cura
{

class Shape;
class SliceMeshStorage;

/*!
 * Anchor information of one part's zigzag, kept from one layer to the next so that the nodes of a
 * layer can be generated based on the medial axis of that layer *and* the node positions of the
 * layer below.
 */
struct MedialZigzagAnchor
{
    bool is_ring{ false }; //!< Whether the part was generated as a ring (walls = outer wall + hole wall).
    size_t num_periods{ 0 }; //!< The number of node periods used, so the layer above can keep it (hysteresis).
    Point2LL match_point; //!< A point well inside the part, used to match parts between consecutive layers.
    Point2LL axis_front; //!< (ribbon) The ends of the medial axis, to keep the axis orientation stable across layers.
    Point2LL axis_back;
    Point2LL phase_point; //!< (ribbon) World position on the axis of a reference upper-side sample, to keep the node phase stable across layers.
    double axis_period{ 0.0 }; //!< (ribbon) The exact node spacing along the axis, kept across layers so jitter of the axis length can't re-space the grid.
};

/*!
 * Generate the medial zigzag lines for one layer region (which may contain multiple parts).
 *
 * Protruding details of the walls (e.g. rows of sawteeth) are removed by a morphological opening
 * before the nodes are placed, so that nodes stay on the main wall envelope instead of wandering
 * into the protrusions. The resulting lines are clipped against \p region, so lines never cross
 * holes or leave the region.
 *
 * \param region The infill region of this layer.
 * \param line_distance The reference line distance; the node period along the medial axis is twice this.
 * \param plane_angle The angle (degrees) of the cutting-plane fallback used for parts without a usable medial axis.
 * \param plane_shift The scanline shift of the cutting-plane fallback (from the infill origin).
 * \param previous_anchors The anchors generated for the layer below; empty when unknown.
 * \param[out] new_anchors The anchors of this layer, to pass to the layer above.
 * \param[out] result The generated zigzag polylines.
 */
void generateMedialZigzagLines(
    const Shape& region,
    coord_t line_distance,
    double plane_angle,
    coord_t plane_shift,
    const std::vector<MedialZigzagAnchor>& previous_anchors,
    std::vector<MedialZigzagAnchor>& new_anchors,
    OpenLinesSet& result);

/*!
 * Pre-computed medial zigzag infill for a whole mesh.
 *
 * The layers are generated sequentially from bottom to top, so that each layer's nodes are anchored
 * to the nodes of the layer below: parts are matched by position, and a matched part keeps the node
 * period count (with hysteresis), the axis orientation and the node phase of the part below it.
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
