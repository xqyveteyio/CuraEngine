// Copyright (c) 2026 UltiMaker
// CuraEngine is released under the terms of the AGPLv3 or higher.

#include "infill/TriangleWaveInfillSimple.h"

#include <algorithm>
#include <cmath>
#include <map>

#include "geometry/OpenPolyline.h"
#include "geometry/Shape.h"
#include "geometry/SingleShape.h"
#include "utils/AABB.h"

namespace cura
{

namespace
{

// Pull the tips inwards by a hair so that clipping against the outline never cuts the apex
// vertex itself and the tips stay sharp.
constexpr coord_t tip_inset = 20;

// The apex X position of column k. The columns lie on a fixed grid in absolute (model)
// coordinates, shifted by half a grid cell to avoid coinciding with typical model edges.
coord_t columnX(const int64_t k, const coord_t line_distance)
{
    return k * line_distance + line_distance / 2;
}

// For every column of the fixed absolute grid, find the lowest and highest boundary crossing of
// the given outline, and merge them into the extremes found so far.
void gatherColumnExtremes(const Shape& outline, const coord_t line_distance, std::map<int64_t, std::pair<coord_t, coord_t>>& extremes)
{
    const AABB aabb(outline);
    const Point2LL y_low(0, aabb.min_.Y - MM2INT(1));
    const Point2LL y_high(0, aabb.max_.Y + MM2INT(1));
    const double span = static_cast<double>(y_high.Y - y_low.Y);

    const int64_t k_start = static_cast<int64_t>(std::floor(static_cast<double>(aabb.min_.X - line_distance / 2) / line_distance));
    const int64_t k_end = static_cast<int64_t>(std::ceil(static_cast<double>(aabb.max_.X - line_distance / 2) / line_distance));

    for (int64_t k = k_start; k <= k_end; ++k)
    {
        const coord_t x = columnX(k, line_distance);
        const std::vector<float> crossings = outline.intersectionsWithSegment(Point2LL(x, y_low.Y), Point2LL(x, y_high.Y));
        if (crossings.size() < 2)
        {
            continue;
        }
        const auto [min_t, max_t] = std::minmax_element(crossings.begin(), crossings.end());
        const coord_t bottom = y_low.Y + static_cast<coord_t>(*min_t * span);
        const coord_t top = y_low.Y + static_cast<coord_t>(*max_t * span);

        auto [it, inserted] = extremes.try_emplace(k, bottom, top);
        if (! inserted)
        {
            it->second.first = std::min(it->second.first, bottom);
            it->second.second = std::max(it->second.second, top);
        }
    }
}

// Draw the complete triangle wave through the column extremes: troughs at even columns, peaks at
// odd columns, with sharp apexes touching the (template) boundary. Columns without material
// interrupt the wave.
OpenLinesSet buildWave(const std::map<int64_t, std::pair<coord_t, coord_t>>& extremes, const coord_t line_distance)
{
    OpenLinesSet wave;
    if (extremes.empty())
    {
        return wave;
    }

    std::vector<Point2LL> wave_points;
    auto flush_wave = [&wave, &wave_points]()
    {
        if (wave_points.size() >= 2)
        {
            wave.push_back(OpenPolyline{ wave_points });
        }
        wave_points.clear();
    };

    for (int64_t k = extremes.begin()->first; k <= extremes.rbegin()->first; ++k)
    {
        const auto it = extremes.find(k);
        if (it == extremes.end())
        {
            flush_wave();
            continue;
        }
        const bool is_peak = (((k % 2) + 2) % 2) == 1; // globally consistent up/down parity
        const coord_t apex_y = is_peak ? it->second.second - tip_inset : it->second.first + tip_inset;
        wave_points.emplace_back(columnX(k, line_distance), apex_y);
    }
    flush_wave();

    return wave;
}

// Build one independent wave per connected part of the given region, so that disjoint areas
// (multiple models or separate islands of one model) each get their own complete triangle wave
// instead of one wave spanning across the gaps between them.
OpenLinesSet buildWaves(const Shape& region, const coord_t line_distance)
{
    OpenLinesSet waves;
    for (const SingleShape& part : region.splitIntoParts())
    {
        std::map<int64_t, std::pair<coord_t, coord_t>> extremes;
        gatherColumnExtremes(part, line_distance, extremes);
        waves.push_back(buildWave(extremes, line_distance));
    }
    return waves;
}

// Clip a wave to the given outline (in the rotated frame) and rotate the result back to the
// original coordinate frame. Pieces cut off by the walls become separate polylines which get
// joined by travel moves later on.
OpenLinesSet clipWave(const OpenLinesSet& wave, const Shape& outline, const PointMatrix& rotation_matrix)
{
    constexpr bool restitch = true;
    OpenLinesSet clipped = outline.intersection(wave, restitch);

    const PointMatrix inverse_rotation = rotation_matrix.inverse();
    for (OpenPolyline& polyline : clipped)
    {
        polyline.applyMatrix(inverse_rotation);
    }
    return clipped;
}

} // namespace

TriangleWaveSimpleFillProvider::TriangleWaveSimpleFillProvider(const std::vector<Shape>& layer_outlines, coord_t line_distance, const double fill_angle)
    : rotation_matrix_(fill_angle)
{
    if (line_distance <= 0)
    {
        return;
    }

    // The template region is the union of the infill areas of all layers: for every column it
    // spans the highest and lowest extent that occurs anywhere in the model. Disjoint parts of
    // the union (separate models or islands) each get their own independent wave.
    Shape all_layers;
    for (const Shape& layer_outline : layer_outlines)
    {
        Shape rotated = layer_outline;
        rotated.applyMatrix(rotation_matrix_);
        all_layers.push_back(rotated);
    }
    all_layers = all_layers.unionPolygons();

    template_wave_ = buildWaves(all_layers, line_distance);
}

void TriangleWaveSimpleFillProvider::generate(OpenLinesSet& result_lines, const Shape& in_outline) const
{
    if (template_wave_.empty() || in_outline.empty())
    {
        return;
    }

    Shape outline = in_outline;
    outline.applyMatrix(rotation_matrix_);

    result_lines = clipWave(template_wave_, outline, rotation_matrix_);
}

void TriangleWaveInfillSimple::generateTotalTriangleWaveInfill(OpenLinesSet& result_lines, coord_t line_distance, const Shape& in_outline, const double fill_angle)
{
    if (line_distance <= 0 || in_outline.empty())
    {
        return;
    }

    const PointMatrix rotation_matrix(fill_angle);
    Shape outline = in_outline;
    outline.applyMatrix(rotation_matrix);

    result_lines = clipWave(buildWaves(outline, line_distance), outline, rotation_matrix);
}

} // namespace cura
