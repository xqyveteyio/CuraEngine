// Copyright (c) 2026 UltiMaker
// CuraEngine is released under the terms of the AGPLv3 or higher.

#ifndef INFILL_TRIANGLE_WAVE_INFILL_H
#define INFILL_TRIANGLE_WAVE_INFILL_H

#include <vector>

#include "geometry/LinesSet.h"
#include "geometry/OpenLinesSet.h"
#include "geometry/PointMatrix.h"
#include "utils/Coord_t.h"

namespace cura
{
class Shape;

/*!
 * Pre-computed, model-global template for the triangle wave infill pattern.
 *
 * The pattern is a single continuous triangle wave (a piecewise-linear zigzag, see
 * https://en.wikipedia.org/wiki/Triangle_wave). To guarantee that the flanks of the wave stack
 * up into perfectly vertical walls across the layers, the wave is NOT generated per layer.
 * Instead, one fixed template wave is derived from all layers at once:
 *
 *  - The area is divided into vertical strips (columns), one tooth flank wide, on a fixed grid
 *    in absolute coordinates.
 *  - For every column, the highest (and lowest) extent of the infill area over ALL layers is
 *    determined; together these span the template region.
 *  - The complete triangle wave is drawn once inside this template region, with sharp apexes
 *    touching the template boundary.
 *
 * Each layer then prints this exact same wave, merely clipped to its own outline (only the
 * length of the lines changes). Where the wave is cut off by the walls, the path is interrupted
 * and the printer travels (G0) to the next piece instead of extruding along the boundary.
 */
class TriangleWaveFillProvider
{
public:
    /*!
     * Build the template wave from the infill areas of all layers.
     * \param layer_outlines The infill areas of all layers (and all their parts).
     * \param line_distance Horizontal distance between two successive flanks of the wave.
     * \param fill_angle The angle (in degrees) of the general direction of the wave. This must
     * be the same for all layers, otherwise the flanks cannot stack up between layers.
     */
    TriangleWaveFillProvider(const std::vector<Shape>& layer_outlines, coord_t line_distance, const double fill_angle);

    /*!
     * Clip the template wave to the outline of one layer.
     * \param result_lines Output variable to store the resulting polylines in.
     * \param in_outline The infill area of this layer.
     */
    void generate(OpenLinesSet& result_lines, const Shape& in_outline) const;

private:
    PointMatrix rotation_matrix_;
    OpenLinesSet template_wave_; // the fixed wave, in the rotated coordinate frame
};

/*!
 * Pre-computed, per-layer waves for the tracking ("advanced") triangle wave infill pattern.
 *
 * This variant is meant for models whose outline drifts sideways from layer to layer (leaning or
 * twisting shapes), where one fixed model-global template cannot follow the geometry. Only the
 * first layer (and parts newly appearing on higher layers) runs the full skeleton-driven
 * generation; every other layer is DERIVED tooth by tooth from the layer below, changing as
 * little as possible:
 *
 *  - Each tooth keeps its rib line (position and direction); only the tooth length is re-cut
 *    against the current outline, so teeth follow walls that moved and stay put everywhere else.
 *  - A wall end that suddenly disappears (e.g. the rib now looks through into another corridor)
 *    freezes the tooth instead of letting it jump there.
 *  - Teeth are added only where region newly appeared (at chain ends or in grown gaps), and
 *    removed only where their material vanished.
 *
 * Identical outlines therefore produce identical waves, and the layer-to-layer difference of the
 * infill is bounded by the layer-to-layer difference of the model itself, which keeps
 * consecutive layers in contact.
 */
class TriangleWaveTrackingProvider
{
public:
    /*!
     * Sequentially build the waves of all layers, bottom up.
     * \param layer_outlines The infill areas of all layers (and all their parts).
     * \param line_distance Horizontal distance between two successive flanks of the wave.
     * \param fill_angle The angle (in degrees) used for regions where no skeleton is available.
     */
    TriangleWaveTrackingProvider(const std::vector<Shape>& layer_outlines, coord_t line_distance, const double fill_angle);

    /*!
     * Clip the pre-computed wave of one layer to the actual infill outline.
     * \param result_lines Output variable to store the resulting polylines in.
     * \param in_outline The infill area of this layer.
     * \param layer_idx The layer number, indexing the outlines given to the constructor.
     */
    void generate(OpenLinesSet& result_lines, const Shape& in_outline, size_t layer_idx) const;

private:
    PointMatrix rotation_matrix_;
    std::vector<OpenLinesSet> layer_waves_; // one pre-computed wave per layer, in the rotated coordinate frame
};

class TriangleWaveInfill
{
public:
    /*!
     * Generate the triangle wave pattern from a single outline, without cross-layer template.
     *
     * This is the fallback used when no \ref TriangleWaveFillProvider is available (e.g. for
     * support). The wave apexes lie on the same fixed absolute grid, but the tooth length is
     * derived from this outline only, so flank slopes may vary between layers.
     *
     * \param result_lines Output variable to store the resulting polylines in.
     * \param line_distance Horizontal distance between two successive flanks of the wave.
     * \param in_outline The outline in which to print the pattern.
     * \param fill_angle The angle (in degrees) of the general direction of the wave.
     */
    static void generateTotalTriangleWaveInfill(OpenLinesSet& result_lines, coord_t line_distance, const Shape& in_outline, const double fill_angle);
};

} // namespace cura

#endif // INFILL_TRIANGLE_WAVE_INFILL_H
