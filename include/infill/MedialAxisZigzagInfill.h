// Copyright (c) 2026 UltiMaker
// CuraEngine is released under the terms of the AGPLv3 or higher.

#ifndef INFILL_MEDIAL_AXIS_ZIGZAG_INFILL_H
#define INFILL_MEDIAL_AXIS_ZIGZAG_INFILL_H

#include "geometry/OpenLinesSet.h"
#include "utils/Coord_t.h"

namespace cura
{
class Shape;

/*!
 * Medial-axis zigzag infill ("中轴波折填充").
 *
 * For each part of the outline the medial axis (centerline) is computed from the Voronoi
 * diagram of the boundary segments. Along this centerline a continuous triangular zigzag
 * (wave) path is generated, symmetric about the centerline. The wave amplitude locally
 * adapts to the part width, so the zigzag peaks just touch the boundary walls.
 *
 * Consecutive layers get a half-period phase shift so the triangles of one layer cross
 * those of the previous layer, interlocking the layers for better inter-layer bonding.
 * This makes the pattern well suited for large-scale (e.g. concrete) printing of wall-like
 * structures: each part is filled by (mostly) one single continuous path without material
 * interruptions.
 */
class MedialAxisZigzagInfill
{
public:
    /*!
     * Generate the medial-axis zigzag pattern within an outline.
     *
     * \param result_lines Output variable to store the resulting polylines in.
     * \param in_outline The area in which to print the pattern.
     * \param line_width The width with which the pattern lines will be extruded.
     * \param zag_spacing The distance along the centerline between two successive zigzag
     *        peaks (half the wave period). Taken from the infill line distance setting.
     * \param phase_flip Whether to shift the wave by half a period (alternated between
     *        consecutive layers to interlock them).
     */
    static void generateMedialAxisZigzagInfill(OpenLinesSet& result_lines, const Shape& in_outline, coord_t line_width, coord_t zag_spacing, bool phase_flip);
};

} // namespace cura

#endif // INFILL_MEDIAL_AXIS_ZIGZAG_INFILL_H
