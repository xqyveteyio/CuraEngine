// Copyright (c) 2026 UltiMaker
// CuraEngine is released under the terms of the AGPLv3 or higher.

#include "infill/MedialAxisZigzagInfill.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <unordered_map>
#include <vector>

#include <boost/polygon/voronoi.hpp>

#include "BoostInterface.hpp"
#include "geometry/OpenPolyline.h"
#include "geometry/Shape.h"
#include "geometry/SingleShape.h"
#include "utils/VoronoiUtils.h"

namespace cura
{

namespace
{

using Segment = PolygonsSegmentIndex;
using vd_t = VoronoiUtils::vd_t;

/*!
 * A piece of the medial axis: the (discretized) polyline of one Voronoi edge, together with
 * the distance to the boundary (the medial axis "radius") at each of its points.
 */
struct MedialChain
{
    std::vector<Point2LL> points;
    std::vector<coord_t> radii;
    const vd_t::vertex_type* v0; //!< Graph node at points.front()
    const vd_t::vertex_type* v1; //!< Graph node at points.back()
    coord_t length{ 0 };
    bool removed{ false };
};

//! A centerline polyline assembled from medial chains, with the local boundary distance at each point.
struct CenterLine
{
    std::vector<Point2LL> points;
    std::vector<coord_t> radii;
};

//! The boundary vertices (poly_idx, vertex indices) that the source feature of a Voronoi cell touches.
struct CellFeature
{
    size_t poly;
    size_t vertex_a;
    size_t vertex_b; // Equal to vertex_a for point-cells.
};

CellFeature getCellFeature(const vd_t::cell_type& cell, const std::vector<Point2LL>& points, const std::vector<Segment>& segments)
{
    if (cell.contains_point())
    {
        const PolygonsPointIndex idx = VoronoiUtils::getSourcePointIndex(cell, points, segments);
        return { idx.poly_idx_, idx.point_idx_, idx.point_idx_ };
    }
    const Segment& segment = VoronoiUtils::getSourceSegment(cell, points, segments);
    const size_t poly_size = (*segment.polygons_)[segment.poly_idx_].size();
    return { segment.poly_idx_, segment.point_idx_, (segment.point_idx_ + 1) % poly_size };
}

//! Whether the source features of two cells share a boundary vertex (i.e. are adjacent on the outline).
bool featuresShareVertex(const CellFeature& a, const CellFeature& b)
{
    if (a.poly != b.poly)
    {
        return false;
    }
    return a.vertex_a == b.vertex_a || a.vertex_a == b.vertex_b || a.vertex_b == b.vertex_a || a.vertex_b == b.vertex_b;
}

/*!
 * Extract the pruned medial axis of one part as a set of chains.
 *
 * Voronoi edges between adjacent boundary features (segments/vertices which share an outline
 * vertex) are the "corner spokes" of the full medial axis; skipping them leaves the actual
 * centerline of the shape.
 */
std::vector<MedialChain> computeMedialChains(const SingleShape& part, const coord_t discretization_step)
{
    std::vector<MedialChain> chains;

    const std::vector<Point2LL> points; // Remains empty, the Voronoi diagram is built from segments only.
    std::vector<Segment> segments;
    for (size_t poly_idx = 0; poly_idx < part.size(); poly_idx++)
    {
        for (size_t point_idx = 0; point_idx < part[poly_idx].size(); point_idx++)
        {
            segments.emplace_back(&part, poly_idx, point_idx);
        }
    }
    if (segments.size() < 3)
    {
        return chains;
    }

    vd_t voronoi_diagram;
    construct_voronoi(segments.begin(), segments.end(), &voronoi_diagram);

    constexpr double transitioning_angle = std::numbers::pi / 4.0; // Only influences the discretization density of parabolic edges.

    for (const vd_t::edge_type& edge : voronoi_diagram.edges())
    {
        if (&edge > edge.twin()) // Each edge appears twice (as its own twin); process only one of the two.
        {
            continue;
        }
        if (! edge.is_finite() || ! edge.is_primary())
        {
            continue;
        }

        const vd_t::cell_type& cell_a = *edge.cell();
        const vd_t::cell_type& cell_b = *edge.twin()->cell();
        if (featuresShareVertex(getCellFeature(cell_a, points, segments), getCellFeature(cell_b, points, segments)))
        {
            continue; // Corner spoke emanating from a boundary vertex; not part of the pruned centerline.
        }

        const Point2LL p0 = VoronoiUtils::p(edge.vertex0());
        const Point2LL p1 = VoronoiUtils::p(edge.vertex1());
        if (p0 == p1)
        {
            continue;
        }
        if (! part.inside(p0, false) || ! part.inside(p1, false) || ! part.inside((p0 + p1) / 2, true))
        {
            continue; // Voronoi edges outside the part (e.g. across concavities or inside holes).
        }

        MedialChain chain;
        chain.v0 = edge.vertex0();
        chain.v1 = edge.vertex1();

        const bool a_is_point = cell_a.contains_point();
        const bool b_is_point = cell_b.contains_point();
        if (a_is_point != b_is_point)
        {
            // Parabolic edge between a reflex boundary vertex (focus) and a boundary segment (directrix).
            const vd_t::cell_type& point_cell = a_is_point ? cell_a : cell_b;
            const vd_t::cell_type& segment_cell = a_is_point ? cell_b : cell_a;
            const Point2LL focus = VoronoiUtils::getSourcePoint(point_cell, points, segments);
            const Segment& directrix = VoronoiUtils::getSourceSegment(segment_cell, points, segments);
            chain.points = VoronoiUtils::discretizeParabola(focus, directrix, p0, p1, discretization_step, transitioning_angle);
            chain.radii.reserve(chain.points.size());
            for (const Point2LL& p : chain.points)
            {
                chain.radii.push_back(vSize(p - focus)); // On the parabola, the distance to the focus equals the distance to the boundary.
            }
        }
        else if (a_is_point) // Both cells are points: straight edge, but the radius varies non-linearly along it.
        {
            const Point2LL focus = VoronoiUtils::getSourcePoint(cell_a, points, segments);
            const coord_t edge_length = vSize(p1 - p0);
            const size_t num_steps = std::max<size_t>(1, edge_length / std::max<coord_t>(discretization_step, 1));
            for (size_t step = 0; step <= num_steps; step++)
            {
                const Point2LL p = p0 + (p1 - p0) * static_cast<coord_t>(step) / static_cast<coord_t>(num_steps);
                chain.points.push_back(p);
                chain.radii.push_back(vSize(p - focus));
            }
        }
        else // Both cells are segments: straight edge along which the radius varies linearly.
        {
            chain.points = { p0, p1 };
            chain.radii = { VoronoiUtils::getDistance(p0, cell_a, points, segments), VoronoiUtils::getDistance(p1, cell_a, points, segments) };
        }

        for (size_t i = 1; i < chain.points.size(); i++)
        {
            chain.length += vSize(chain.points[i] - chain.points[i - 1]);
        }
        if (chain.length == 0)
        {
            continue;
        }
        chains.push_back(std::move(chain));
    }

    return chains;
}

/*!
 * Iteratively remove leaf chains which are just residual spokes towards (rounded/chamfered)
 * corners: a leaf branch shorter than the medial axis radius at its junction carries no
 * structural information about the centerline.
 */
void pruneMedialChains(std::vector<MedialChain>& chains, const coord_t line_width)
{
    std::unordered_map<const vd_t::vertex_type*, std::vector<size_t>> node_chains;
    for (size_t chain_idx = 0; chain_idx < chains.size(); chain_idx++)
    {
        node_chains[chains[chain_idx].v0].push_back(chain_idx);
        node_chains[chains[chain_idx].v1].push_back(chain_idx);
    }

    const auto degree = [&](const vd_t::vertex_type* node) -> size_t
    {
        size_t result = 0;
        for (const size_t chain_idx : node_chains[node])
        {
            if (! chains[chain_idx].removed)
            {
                result++;
            }
        }
        return result;
    };

    bool changed = true;
    while (changed)
    {
        changed = false;
        for (MedialChain& chain : chains)
        {
            if (chain.removed)
            {
                continue;
            }
            const size_t degree_v0 = degree(chain.v0);
            const size_t degree_v1 = degree(chain.v1);
            if (degree_v0 == 1 && degree_v1 > 1)
            {
                const coord_t junction_radius = chain.radii.back();
                if (chain.length < junction_radius * 3 / 2 + line_width)
                {
                    chain.removed = true;
                    changed = true;
                }
            }
            else if (degree_v1 == 1 && degree_v0 > 1)
            {
                const coord_t junction_radius = chain.radii.front();
                if (chain.length < junction_radius * 3 / 2 + line_width)
                {
                    chain.removed = true;
                    changed = true;
                }
            }
        }
    }
}

/*!
 * Assemble the remaining chains into as few and as long centerline polylines as possible.
 * At junctions the straightest continuation is preferred, so the main centerline keeps going
 * and side branches become separate polylines.
 */
std::vector<CenterLine> assembleCenterLines(std::vector<MedialChain>& chains)
{
    std::vector<CenterLine> center_lines;

    std::unordered_map<const vd_t::vertex_type*, std::vector<size_t>> node_chains;
    for (size_t chain_idx = 0; chain_idx < chains.size(); chain_idx++)
    {
        if (chains[chain_idx].removed)
        {
            continue;
        }
        node_chains[chains[chain_idx].v0].push_back(chain_idx);
        node_chains[chains[chain_idx].v1].push_back(chain_idx);
    }

    std::vector<bool> used(chains.size(), false);

    const auto appendChain = [&](CenterLine& line, const size_t chain_idx, const bool forward)
    {
        const MedialChain& chain = chains[chain_idx];
        const bool include_first = line.points.empty();
        if (forward)
        {
            for (size_t i = include_first ? 0 : 1; i < chain.points.size(); i++)
            {
                line.points.push_back(chain.points[i]);
                line.radii.push_back(chain.radii[i]);
            }
        }
        else
        {
            for (size_t i = include_first ? chain.points.size() : chain.points.size() - 1; i > 0; i--)
            {
                line.points.push_back(chain.points[i - 1]);
                line.radii.push_back(chain.radii[i - 1]);
            }
        }
    };

    const auto walkFrom = [&](const size_t start_chain, const vd_t::vertex_type* start_node)
    {
        CenterLine line;
        size_t chain_idx = start_chain;
        const vd_t::vertex_type* node = start_node;
        while (true)
        {
            used[chain_idx] = true;
            const MedialChain& chain = chains[chain_idx];
            const bool forward = (chain.v0 == node);
            appendChain(line, chain_idx, forward);
            node = forward ? chain.v1 : chain.v0;

            // Select the straightest unused continuation at this node.
            const Point2LL incoming_direction = line.points.back() - line.points[line.points.size() - 2];
            size_t best_chain = std::numeric_limits<size_t>::max();
            double best_alignment = -2.0;
            for (const size_t candidate_idx : node_chains[node])
            {
                const MedialChain& candidate = chains[candidate_idx];
                if (used[candidate_idx] || candidate.removed)
                {
                    continue;
                }
                const bool candidate_forward = (candidate.v0 == node);
                const Point2LL outgoing_direction = candidate_forward ? (candidate.points[1] - candidate.points[0])
                                                                      : (candidate.points[candidate.points.size() - 2] - candidate.points.back());
                const double lengths = vSize(incoming_direction) * vSize(outgoing_direction);
                if (lengths <= 0)
                {
                    continue;
                }
                const double alignment = static_cast<double>(dot(incoming_direction, outgoing_direction)) / lengths;
                if (alignment > best_alignment)
                {
                    best_alignment = alignment;
                    best_chain = candidate_idx;
                }
            }
            if (best_chain == std::numeric_limits<size_t>::max())
            {
                break;
            }
            chain_idx = best_chain;
        }
        if (line.points.size() >= 2)
        {
            center_lines.push_back(std::move(line));
        }
    };

    const auto degree = [&](const vd_t::vertex_type* node) -> size_t
    {
        size_t result = 0;
        for (const size_t chain_idx : node_chains[node])
        {
            if (! chains[chain_idx].removed && ! used[chain_idx])
            {
                result++;
            }
        }
        return result;
    };

    // First start walks from end-points of the medial axis, so paths span whole branches.
    for (size_t chain_idx = 0; chain_idx < chains.size(); chain_idx++)
    {
        const MedialChain& chain = chains[chain_idx];
        if (chain.removed || used[chain_idx])
        {
            continue;
        }
        if (degree(chain.v0) == 1)
        {
            walkFrom(chain_idx, chain.v0);
        }
        else if (degree(chain.v1) == 1)
        {
            walkFrom(chain_idx, chain.v1);
        }
    }
    // Then handle whatever remains (e.g. cyclic medial axes around holes).
    for (size_t chain_idx = 0; chain_idx < chains.size(); chain_idx++)
    {
        if (! chains[chain_idx].removed && ! used[chain_idx])
        {
            walkFrom(chain_idx, chains[chain_idx].v0);
        }
    }

    return center_lines;
}

/*!
 * Generate the triangular wave along one centerline. Peak positions are fixed on the arc
 * length of the centerline, so consecutive layers (with flipped phase) interlock.
 */
void generateWaveAlongCenterLine(const CenterLine& center_line, const coord_t half_period, const coord_t line_width, const bool phase_flip, OpenLinesSet& result_lines)
{
    std::vector<coord_t> cumulative_lengths;
    cumulative_lengths.reserve(center_line.points.size());
    cumulative_lengths.push_back(0);
    for (size_t i = 1; i < center_line.points.size(); i++)
    {
        cumulative_lengths.push_back(cumulative_lengths.back() + vSize(center_line.points[i] - center_line.points[i - 1]));
    }
    const coord_t total_length = cumulative_lengths.back();

    OpenPolyline wave;
    if (total_length < half_period || center_line.points.size() < 2)
    {
        // Too short for even a single zag; extrude along the centerline so the region still gets material.
        wave = OpenPolyline(std::initializer_list<Point2LL>{});
        for (const Point2LL& p : center_line.points)
        {
            wave.push_back(p);
        }
        if (wave.size() >= 2)
        {
            result_lines.push_back(std::move(wave));
        }
        return;
    }

    wave.push_back(center_line.points.front());
    // Flipping the starting side of a symmetric triangle wave equals shifting it by half a
    // period; alternating this per layer staggers the peaks so the layers interlock.
    int side = phase_flip ? -1 : 1;

    size_t segment_idx = 0;
    for (coord_t s = half_period / 2; s < total_length; s += half_period)
    {
        while (segment_idx + 2 < cumulative_lengths.size() && cumulative_lengths[segment_idx + 1] < s)
        {
            segment_idx++;
        }
        const coord_t segment_start = cumulative_lengths[segment_idx];
        const coord_t segment_length = cumulative_lengths[segment_idx + 1] - segment_start;
        const Point2LL segment_vector = center_line.points[segment_idx + 1] - center_line.points[segment_idx];
        const double along = segment_length > 0 ? static_cast<double>(s - segment_start) / static_cast<double>(segment_length) : 0.0;

        const Point2LL center = center_line.points[segment_idx] + Point2LL(std::llrint(segment_vector.X * along), std::llrint(segment_vector.Y * along));
        const coord_t radius = center_line.radii[segment_idx] + std::llrint((center_line.radii[segment_idx + 1] - center_line.radii[segment_idx]) * along);
        const coord_t amplitude = std::max<coord_t>(0, radius - line_width / 2);

        const Point2LL offset = normal(turn90CCW(segment_vector), amplitude);
        wave.push_back(side > 0 ? center + offset : center - offset);
        side = -side;
    }
    wave.push_back(center_line.points.back());
    result_lines.push_back(std::move(wave));
}

} // namespace

void MedialAxisZigzagInfill::generateMedialAxisZigzagInfill(
    OpenLinesSet& result_lines,
    const Shape& in_outline,
    const coord_t line_width,
    const coord_t zag_spacing,
    const bool phase_flip)
{
    const coord_t half_period = std::max<coord_t>(zag_spacing, line_width);
    const coord_t discretization_step = std::max<coord_t>(200, half_period / 4);

    for (const SingleShape& part : in_outline.splitIntoParts())
    {
        std::vector<MedialChain> chains = computeMedialChains(part, discretization_step);
        if (chains.empty())
        {
            continue;
        }
        pruneMedialChains(chains, line_width);

        OpenLinesSet part_lines;
        for (const CenterLine& center_line : assembleCenterLines(chains))
        {
            generateWaveAlongCenterLine(center_line, half_period, line_width, phase_flip, part_lines);
        }

        // Straight zags may cut across the boundary at sharp centerline turns; clip to be safe.
        result_lines.push_back(part.intersection(part_lines, true, line_width));
    }
}

} // namespace cura
