// Copyright (c) 2026 UltiMaker
// CuraEngine is released under the terms of the AGPLv3 or higher.

#include "infill/SurfaceTriangleWaveInfill.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <unordered_map>
#include <vector>

#ifdef SW_DEBUG
#include <cstdio>
#endif

#include <boost/polygon/voronoi.hpp>

#include "BoostInterface.hpp"
#include "geometry/OpenPolyline.h"
#include "geometry/Shape.h"
#include "geometry/SingleShape.h"
#include "utils/VoronoiUtils.h"
#include "utils/linearAlg2D.h"

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
 * centerline of the shape, which serves as the longitudinal reference axis of the (curved) part.
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
 * Iteratively remove leaf branches which are just residual spokes towards (rounded/chamfered)
 * corners or towards boundary noise: a leaf branch whose TOTAL length (from its tip up to the
 * junction, accumulated over already pruned chains) is shorter than the medial axis radius at
 * its junction carries no structural information about the centerline. Tracking the accumulated
 * length is essential: on noisy outlines the centerline is split into many short chains, and
 * removing them one by one without accumulation would cascade through the entire spine.
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

    // Length of the already pruned sub-branch hanging off a node.
    std::unordered_map<const vd_t::vertex_type*, coord_t> pruned_length;

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
            const vd_t::vertex_type* leaf = nullptr;
            const vd_t::vertex_type* junction = nullptr;
            coord_t junction_radius = 0;
            if (degree_v0 == 1 && degree_v1 > 1)
            {
                leaf = chain.v0;
                junction = chain.v1;
                junction_radius = chain.radii.back();
            }
            else if (degree_v1 == 1 && degree_v0 > 1)
            {
                leaf = chain.v1;
                junction = chain.v0;
                junction_radius = chain.radii.front();
            }
            else
            {
                continue;
            }
            const coord_t branch_length = chain.length + pruned_length[leaf];
            if (branch_length < junction_radius * 3 / 2 + line_width)
            {
                chain.removed = true;
                pruned_length[junction] = std::max(pruned_length[junction], branch_length);
                changed = true;
            }
        }
    }
}

/*!
 * Assemble the remaining chains into as few and as long centerline polylines as possible.
 * At junctions the straightest continuation is preferred, so the main reference axis keeps
 * going and side branches become separate polylines.
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

    // A continuation is only followed while it stays reasonably straight; side branches (e.g.
    // ribs standing on a main beam) then become their own axes instead of cutting the main
    // reference axis into arbitrary pieces (which would differ from layer to layer).
    constexpr double min_alignment = 0.5; // cos(60 degrees)
    // Both directions are measured over a window instead of a single (possibly tiny and noisy)
    // polyline segment, so boundary noise does not derail the walk at junctions.
    constexpr coord_t direction_window = 5000;

    // The direction of the last \p direction_window of the line, ending at its last point.
    const auto incomingDirection = [](const CenterLine& line) -> Point2LL
    {
        const Point2LL end = line.points.back();
        for (size_t i = line.points.size() - 1; i > 0; i--)
        {
            if (vSize(end - line.points[i - 1]) >= direction_window)
            {
                return end - line.points[i - 1];
            }
        }
        return end - line.points.front();
    };

    // The direction of the first \p direction_window of the chain, seen walking away from \p node.
    const auto outgoingDirection = [&chains](const size_t chain_idx, const bool forward) -> Point2LL
    {
        const MedialChain& chain = chains[chain_idx];
        const Point2LL start = forward ? chain.points.front() : chain.points.back();
        if (forward)
        {
            for (size_t i = 1; i < chain.points.size(); i++)
            {
                if (vSize(chain.points[i] - start) >= direction_window)
                {
                    return chain.points[i] - start;
                }
            }
            return chain.points.back() - start;
        }
        for (size_t i = chain.points.size() - 1; i > 0; i--)
        {
            if (vSize(chain.points[i - 1] - start) >= direction_window)
            {
                return chain.points[i - 1] - start;
            }
        }
        return chain.points.front() - start;
    };

    // Keep walking from \p node, always along the straightest unused continuation, and append
    // the traversed chains to \p line.
    const auto extend = [&](CenterLine& line, const vd_t::vertex_type* node)
    {
        while (true)
        {
            const Point2LL incoming_direction = incomingDirection(line);
            size_t best_chain = std::numeric_limits<size_t>::max();
            double best_alignment = min_alignment;
            for (const size_t candidate_idx : node_chains[node])
            {
                const MedialChain& candidate = chains[candidate_idx];
                if (used[candidate_idx] || candidate.removed)
                {
                    continue;
                }
                const bool candidate_forward = (candidate.v0 == node);
                const Point2LL outgoing_direction = outgoingDirection(candidate_idx, candidate_forward);
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
            used[best_chain] = true;
            const bool forward = (chains[best_chain].v0 == node);
            appendChain(line, best_chain, forward);
            node = forward ? chains[best_chain].v1 : chains[best_chain].v0;
        }
    };

    // Walk from a seed chain in BOTH directions, so the resulting axis always spans the whole
    // straight(ish) branch, independent of the enumeration order of the chains.
    const auto walkFrom = [&](const size_t start_chain)
    {
        CenterLine line;
        used[start_chain] = true;
        appendChain(line, start_chain, true);
        extend(line, chains[start_chain].v1);
        std::reverse(line.points.begin(), line.points.end());
        std::reverse(line.radii.begin(), line.radii.end());
        extend(line, chains[start_chain].v0);
        if (line.points.size() >= 2)
        {
            center_lines.push_back(std::move(line));
        }
    };

    // Seed with the longest chains first, so the main reference axis of the part is assembled
    // before its side branches consume any of its chains.
    std::vector<size_t> chain_order;
    for (size_t chain_idx = 0; chain_idx < chains.size(); chain_idx++)
    {
        if (! chains[chain_idx].removed)
        {
            chain_order.push_back(chain_idx);
        }
    }
    std::sort(
        chain_order.begin(),
        chain_order.end(),
        [&chains](const size_t a, const size_t b)
        {
            return chains[a].length > chains[b].length;
        });
    for (const size_t chain_idx : chain_order)
    {
        if (! used[chain_idx])
        {
            walkFrom(chain_idx);
        }
    }

    return center_lines;
}

/*!
 * The reference axis ("曲面基准轴线") of one branch of a part, on one layer, together with the
 * wave that is defined relative to it.
 *
 * The wave is fully parameterized by the arc length along this axis: every apex (wave peak) is
 * an arc-length position plus a side (sides strictly alternate). On the layer where a part
 * first appears, the apex positions are laid out with a fixed modulus per curvature segment; on
 * every following layer they are obtained by projecting the apexes of the layer below onto the
 * current axis, so the wave translates, bends and scales together with the surface while the
 * apex count and the relative phase stay locked.
 */
struct ReferenceAxis
{
    std::vector<Point2LL> points;
    std::vector<coord_t> radii; //!< Medial-axis radius (half the local cross-section width) at each point.
    std::vector<coord_t> cumulative; //!< Arc length from the start of the axis at each point.
    coord_t total_length{ 0 };

    // The wave, inherited layer over layer:
    std::vector<coord_t> apex_positions; //!< Arc-length positions of the wave apexes, strictly increasing.
    std::vector<bool> apex_sides; //!< For each apex, whether it lies on the left (positive) side of the axis; strictly alternating.
    coord_t amplitude{ 0 }; //!< The uniform amplitude of ALL apexes ("统一等幅"), limited by the narrowest cross section the wave passes.
};

//! The (windowed) direction of an axis at one of its ends, pointing away from the axis.
Point2LL axisOutwardDirection(const std::vector<Point2LL>& points, const bool at_back)
{
    constexpr coord_t window = 5000;
    const Point2LL end = at_back ? points.back() : points.front();
    if (at_back)
    {
        for (size_t i = points.size() - 1; i > 0; i--)
        {
            if (vSize(end - points[i - 1]) >= window)
            {
                return end - points[i - 1];
            }
        }
        return end - points.front();
    }
    for (size_t i = 1; i < points.size(); i++)
    {
        if (vSize(end - points[i]) >= window)
        {
            return end - points[i];
        }
    }
    return end - points.back();
}

/*!
 * Stitch centerline fragments into as few and as long axes as possible. The medial axis of a
 * noisy outline is often assembled into several collinear pieces (broken at rib junctions or
 * pruned spots); for a stable, layer-to-layer reproducible primary reference axis these pieces
 * are merged end-to-end whenever the gap is small and the directions continue each other.
 */
void stitchCenterLines(std::vector<CenterLine>& lines, const coord_t max_gap)
{
    bool merged_any = true;
    while (merged_any)
    {
        merged_any = false;
        // Always extend the longest line first, so the primary axis assembles greedily.
        std::sort(
            lines.begin(),
            lines.end(),
            [](const CenterLine& a, const CenterLine& b)
            {
                return a.points.size() > b.points.size();
            });
        for (size_t base_idx = 0; base_idx < lines.size() && ! merged_any; base_idx++)
        {
            CenterLine& base = lines[base_idx];
            if (base.points.size() < 2)
            {
                continue;
            }
            for (size_t other_idx = 0; other_idx < lines.size() && ! merged_any; other_idx++)
            {
                CenterLine& other = lines[other_idx];
                if (other_idx == base_idx || other.points.size() < 2)
                {
                    continue;
                }
                for (const bool base_back : { false, true })
                {
                    for (const bool other_front : { true, false })
                    {
                        const Point2LL base_end = base_back ? base.points.back() : base.points.front();
                        const Point2LL other_end = other_front ? other.points.front() : other.points.back();
                        if (vSize(base_end - other_end) > max_gap)
                        {
                            continue;
                        }
                        // The other line has to continue in roughly the direction the base ends with.
                        const Point2LL base_out = axisOutwardDirection(base.points, base_back);
                        const Point2LL other_in = axisOutwardDirection(other.points, ! other_front) * -1;
                        const double lengths = vSize(base_out) * vSize(other_in);
                        if (lengths <= 0 || static_cast<double>(dot(base_out, other_in)) / lengths < 0.3)
                        {
                            continue;
                        }

                        // Merge: base (oriented so the junction is at its back) + other (oriented
                        // so the junction is at its front).
                        if (! base_back)
                        {
                            std::reverse(base.points.begin(), base.points.end());
                            std::reverse(base.radii.begin(), base.radii.end());
                        }
                        if (! other_front)
                        {
                            std::reverse(other.points.begin(), other.points.end());
                            std::reverse(other.radii.begin(), other.radii.end());
                        }
                        const bool skip_duplicate = other.points.front() == base.points.back();
                        base.points.insert(base.points.end(), other.points.begin() + (skip_duplicate ? 1 : 0), other.points.end());
                        base.radii.insert(base.radii.end(), other.radii.begin() + (skip_duplicate ? 1 : 0), other.radii.end());
                        other.points.clear();
                        other.radii.clear();
                        merged_any = true;
                        break;
                    }
                    if (merged_any)
                    {
                        break;
                    }
                }
            }
        }
    }
    std::erase_if(
        lines,
        [](const CenterLine& line)
        {
            return line.points.size() < 2;
        });
}

void finalizeAxisGeometry(ReferenceAxis& axis)
{
    axis.cumulative.clear();
    axis.cumulative.reserve(axis.points.size());
    axis.cumulative.push_back(0);
    for (size_t i = 1; i < axis.points.size(); i++)
    {
        axis.cumulative.push_back(axis.cumulative.back() + vSize(axis.points[i] - axis.points[i - 1]));
    }
    axis.total_length = axis.cumulative.back();
}

//! Evaluate position, tangent and local radius of the axis at arc length \p s.
void sampleAxis(const ReferenceAxis& axis, const coord_t s, Point2LL& position, Point2LL& tangent, coord_t& radius)
{
    const auto it = std::upper_bound(axis.cumulative.begin(), axis.cumulative.end(), s);
    size_t segment_idx = (it == axis.cumulative.begin()) ? 0 : static_cast<size_t>(it - axis.cumulative.begin()) - 1;
    segment_idx = std::min(segment_idx, axis.points.size() - 2);

    const coord_t segment_start = axis.cumulative[segment_idx];
    const coord_t segment_length = axis.cumulative[segment_idx + 1] - segment_start;
    const Point2LL segment_vector = axis.points[segment_idx + 1] - axis.points[segment_idx];
    const double along = segment_length > 0 ? std::clamp(static_cast<double>(s - segment_start) / static_cast<double>(segment_length), 0.0, 1.0) : 0.0;

    position = axis.points[segment_idx] + Point2LL(std::llrint(segment_vector.X * along), std::llrint(segment_vector.Y * along));
    tangent = segment_vector;
    radius = axis.radii[segment_idx] + std::llrint((axis.radii[segment_idx + 1] - axis.radii[segment_idx]) * along);
}

/*!
 * Lay out a fresh wave on an axis that has no predecessor on the layer below, using segments of
 * equal modulus ("分段模数均分"): a segment boundary is inserted wherever the accumulated
 * turning angle of the axis exceeds 45 degrees, so that on strongly curved surfaces every
 * segment covers a stretch of roughly constant direction. Within one segment a FIXED number of
 * wave cells is placed (rounded from the nominal period); the actual period is then the segment
 * length divided by that count, uniform inside the segment and continuous over the boundary.
 */
void layoutFreshApexes(ReferenceAxis& axis, const coord_t cell_period)
{
    axis.apex_positions.clear();
    axis.apex_sides.clear();

    const coord_t min_segment_length = 2 * cell_period;
    std::vector<coord_t> cut_positions = { 0 };
    if (axis.points.size() > 2)
    {
        double accumulated_turn = 0.0;
        double previous_angle = std::atan2(static_cast<double>(axis.points[1].Y - axis.points[0].Y), static_cast<double>(axis.points[1].X - axis.points[0].X));
        coord_t last_cut = 0;
        for (size_t i = 1; i + 1 < axis.points.size(); i++)
        {
            const double angle
                = std::atan2(static_cast<double>(axis.points[i + 1].Y - axis.points[i].Y), static_cast<double>(axis.points[i + 1].X - axis.points[i].X));
            double turn = angle - previous_angle;
            while (turn > std::numbers::pi)
            {
                turn -= 2.0 * std::numbers::pi;
            }
            while (turn < -std::numbers::pi)
            {
                turn += 2.0 * std::numbers::pi;
            }
            previous_angle = angle;
            accumulated_turn += std::abs(turn);

            if (accumulated_turn > std::numbers::pi / 4.0 && axis.cumulative[i] - last_cut >= min_segment_length
                && axis.total_length - axis.cumulative[i] >= min_segment_length)
            {
                cut_positions.push_back(axis.cumulative[i]);
                last_cut = axis.cumulative[i];
                accumulated_turn = 0.0;
            }
        }
    }
    cut_positions.push_back(axis.total_length);

    bool side_positive = true;
    for (size_t segment_idx = 0; segment_idx + 1 < cut_positions.size(); segment_idx++)
    {
        const coord_t s0 = cut_positions[segment_idx];
        const coord_t segment_length = cut_positions[segment_idx + 1] - s0;
        const size_t cells = std::max<size_t>(1, static_cast<size_t>(std::llround(static_cast<double>(segment_length) / static_cast<double>(cell_period))));
        for (size_t cell = 0; cell < cells; cell++)
        {
            // Two apexes per cell, at 1/4 and 3/4 of the cell, on alternating sides of the axis.
            for (const double in_cell : { 0.25, 0.75 })
            {
                const coord_t s = s0 + std::llrint(static_cast<double>(segment_length) * (static_cast<double>(cell) + in_cell) / static_cast<double>(cells));
                axis.apex_positions.push_back(s);
                axis.apex_sides.push_back(side_positive);
                side_positive = ! side_positive;
            }
        }
    }
}

/*!
 * The uniform wave amplitude of an axis ("统一等幅"): limited by the narrowest cross section
 * along the span covered by the wave, minus half the extrusion width, so the equal-amplitude
 * wave fits inside the boundary everywhere and its lines never overlap the walls.
 */
void computeUniformAmplitude(ReferenceAxis& axis, const coord_t line_width)
{
    if (axis.apex_positions.empty())
    {
        axis.amplitude = 0;
        return;
    }
    const coord_t s_first = axis.apex_positions.front();
    const coord_t s_last = axis.apex_positions.back();
    coord_t min_radius = std::numeric_limits<coord_t>::max();
    for (size_t i = 0; i < axis.points.size(); i++)
    {
        if (axis.cumulative[i] >= s_first && axis.cumulative[i] <= s_last)
        {
            min_radius = std::min(min_radius, axis.radii[i]);
        }
    }
    // Also sample at the apex positions themselves (the range above may not contain any axis vertex).
    for (const coord_t s : { s_first, s_last })
    {
        Point2LL position;
        Point2LL tangent;
        coord_t radius;
        sampleAxis(axis, s, position, tangent, radius);
        min_radius = std::min(min_radius, radius);
    }
    axis.amplitude = std::max<coord_t>(0, min_radius - line_width / 2);
}

//! The 2D tip of one apex, computed on the axis it belongs to; the geometric feature that the next layer's wave binds to.
Point2LL apexTip(const ReferenceAxis& axis, const size_t apex_idx)
{
    Point2LL center;
    Point2LL tangent;
    coord_t radius;
    sampleAxis(axis, axis.apex_positions[apex_idx], center, tangent, radius);
    const Point2LL offset = normal(turn90CCW(tangent), axis.amplitude);
    return axis.apex_sides[apex_idx] ? center + offset : center - offset;
}

/*!
 * Remove self-intersections from a polyline by short-cutting: where two segments properly
 * cross (which can happen where apex tips fold over on the inside of a tightly bent axis),
 * the loop between them is cut out. This guarantees the infill lines never overlap each other.
 */
void removeSelfCrossings(std::vector<Point2LL>& points)
{
    const auto orientation = [](const Point2LL& a, const Point2LL& b, const Point2LL& c) -> int
    {
        const auto v = cross(b - a, c - a);
        return (v > 0) - (v < 0);
    };
    bool changed = true;
    while (changed)
    {
        changed = false;
        for (size_t i = 0; i + 1 < points.size() && ! changed; i++)
        {
            for (size_t j = i + 2; j + 1 < points.size() && ! changed; j++)
            {
                const Point2LL& p1 = points[i];
                const Point2LL& p2 = points[i + 1];
                const Point2LL& q1 = points[j];
                const Point2LL& q2 = points[j + 1];
                if (p1 == q1 || p1 == q2 || p2 == q1 || p2 == q2)
                {
                    continue;
                }
                if (orientation(p1, p2, q1) * orientation(p1, p2, q2) < 0 && orientation(q1, q2, p1) * orientation(q1, q2, p2) < 0)
                {
                    points.erase(points.begin() + static_cast<std::ptrdiff_t>(i) + 1, points.begin() + static_cast<std::ptrdiff_t>(j) + 1);
                    changed = true;
                }
            }
        }
    }
}

/*!
 * Generate the triangle wave that belongs to one reference axis: one single continuous,
 * unbroken polyline from axis start to axis end that visits every apex tip in order. All
 * apexes use the same amplitude, giving a uniform, symmetric equal-amplitude triangle wave.
 */
void generateWaveAlongAxis(const ReferenceAxis& axis, OpenLinesSet& result_lines)
{
    if (axis.points.size() < 2 || axis.apex_positions.size() < 2)
    {
        return;
    }

    std::vector<Point2LL> points;
    points.push_back(axis.points.front());
    for (size_t apex_idx = 0; apex_idx < axis.apex_positions.size(); apex_idx++)
    {
        points.push_back(apexTip(axis, apex_idx));
    }
    points.push_back(axis.points.back());
    removeSelfCrossings(points);

    OpenPolyline wave;
    for (const Point2LL& p : points)
    {
        wave.push_back(p);
    }
    result_lines.push_back(std::move(wave));
}

//! One apex of the layer below, as a geometric feature: where it sat on its axis and where its tip was.
struct ApexFeature
{
    Point2LL base;
    Point2LL tip;
};

//! All data of one part on one layer that the next layer derives its wave from.
struct PartRecord
{
    SingleShape part;
    std::vector<ApexFeature> apexes;
    Point2LL probe; //!< A point well inside the part (midpoint of its longest axis), for part-to-part matching.
};

Point2LL axisMidpoint(const ReferenceAxis& axis)
{
    Point2LL position;
    Point2LL tangent;
    coord_t radius;
    sampleAxis(axis, axis.total_length / 2, position, tangent, radius);
    return position;
}

/*!
 * Find the part of the layer below that the given part stands on, by testing whether the
 * midpoint of the current reference axis (guaranteed to lie inside the current part) lies
 * inside a part of the layer below, or vice versa.
 */
const PartRecord* findPartBelow(const std::vector<PartRecord>& records_below, const Point2LL& current_probe, const SingleShape& current_part)
{
    for (const PartRecord& record : records_below)
    {
        if (record.part.inside(current_probe, true))
        {
            return &record;
        }
    }
    // The parts may have drifted sideways more than their overlap; also try the reverse test.
    for (const PartRecord& record : records_below)
    {
        if (current_part.inside(record.probe, true))
        {
            return &record;
        }
    }
    return nullptr;
}

//! Project point \p p onto the axis polyline; returns the arc length of the closest point and the distance to it.
void projectOntoAxis(const ReferenceAxis& axis, const Point2LL& p, coord_t& arc_length, coord_t& distance)
{
    distance = std::numeric_limits<coord_t>::max();
    arc_length = 0;
    for (size_t i = 0; i + 1 < axis.points.size(); i++)
    {
        const Point2LL closest = LinearAlg2D::getClosestOnLineSegment(p, axis.points[i], axis.points[i + 1]);
        const coord_t d = vSize(closest - p);
        if (d < distance)
        {
            distance = d;
            arc_length = axis.cumulative[i] + vSize(closest - axis.points[i]);
        }
    }
}

/*!
 * Derive the wave of the given axes from the apex features of the part below ("分层路径映射复用"):
 * every apex of the layer below is projected (along the surface normal, i.e. onto the nearest
 * point of the new axis) into the current layer, keeping its side. The apex count and order --
 * and thereby the phase -- are inherited unchanged; the positions only move as far as the
 * surface itself moved, so consecutive layers' peaks stay locked onto each other. Where the
 * axis extends beyond the inherited apexes (newly appeared region), the wave is continued with
 * the nominal period.
 */
void inheritApexesFromBelow(std::vector<ReferenceAxis>& axes, const std::vector<ApexFeature>& below_apexes, const coord_t cell_period, const coord_t line_width)
{
    struct ProjectedApex
    {
        coord_t arc_length;
        bool side_positive;
    };
    std::vector<std::vector<ProjectedApex>> projected(axes.size());

    for (const ApexFeature& apex : below_apexes)
    {
        // Bind the apex to the nearest current axis.
        size_t best_axis = std::numeric_limits<size_t>::max();
        coord_t best_distance = std::numeric_limits<coord_t>::max();
        coord_t best_arc_length = 0;
        for (size_t axis_idx = 0; axis_idx < axes.size(); axis_idx++)
        {
            coord_t arc_length;
            coord_t distance;
            projectOntoAxis(axes[axis_idx], apex.base, arc_length, distance);
            if (distance < best_distance)
            {
                best_distance = distance;
                best_axis = axis_idx;
                best_arc_length = arc_length;
            }
        }
        if (best_axis == std::numeric_limits<size_t>::max() || best_distance > cell_period)
        {
            continue; // This piece of the part below has no counterpart on this layer.
        }

        // Re-derive the side geometrically from the tip, so it is independent of the (arbitrary)
        // orientation of either axis polyline.
        const ReferenceAxis& axis = axes[best_axis];
        Point2LL center;
        Point2LL tangent;
        coord_t radius;
        sampleAxis(axis, best_arc_length, center, tangent, radius);
        const bool side_positive = cross(tangent, apex.tip - center) > 0;
        projected[best_axis].push_back({ best_arc_length, side_positive });
    }

    for (size_t axis_idx = 0; axis_idx < axes.size(); axis_idx++)
    {
        ReferenceAxis& axis = axes[axis_idx];
        std::vector<ProjectedApex>& apexes = projected[axis_idx];
        if (apexes.empty())
        {
            layoutFreshApexes(axis, cell_period); // Newly appeared branch: no wave below to inherit from.
            continue;
        }
        std::sort(
            apexes.begin(),
            apexes.end(),
            [](const ProjectedApex& a, const ProjectedApex& b)
            {
                return a.arc_length < b.arc_length;
            });

        axis.apex_positions.clear();
        axis.apex_sides.clear();
        for (const ProjectedApex& apex : apexes)
        {
            // Apexes may collapse where the axis got shorter (their positions clamp against the
            // axis ends); keep only the first of such a bunch, the shape has no room for more.
            if (! axis.apex_positions.empty() && apex.arc_length - axis.apex_positions.back() < line_width)
            {
                continue;
            }
            // The wave must stay a strictly alternating (symmetric, non-self-overlapping)
            // zigzag; if projection produced two consecutive same-side apexes (e.g. because an
            // in-between apex was dropped), skip the later one.
            if (! axis.apex_sides.empty() && apex.side_positive == axis.apex_sides.back())
            {
                continue;
            }
            axis.apex_positions.push_back(apex.arc_length);
            axis.apex_sides.push_back(apex.side_positive);
        }

        if (axis.apex_positions.size() < 2)
        {
            // Not enough inherited apexes for even one period; lay the wave out afresh so this
            // axis still gets at least one full triangle wave period.
            layoutFreshApexes(axis, cell_period);
            continue;
        }

        // Continue the wave with the nominal period into regions the axis newly grew into.
        const coord_t half_period = cell_period / 2;
        while (axis.apex_positions.front() > half_period * 3 / 2)
        {
            axis.apex_positions.insert(axis.apex_positions.begin(), axis.apex_positions.front() - half_period);
            axis.apex_sides.insert(axis.apex_sides.begin(), ! axis.apex_sides.front());
        }
        while (axis.total_length - axis.apex_positions.back() > half_period * 3 / 2)
        {
            axis.apex_positions.push_back(axis.apex_positions.back() + half_period);
            axis.apex_sides.push_back(! axis.apex_sides.back());
        }
    }
}

/*!
 * Build the waves of one layer. For a part that also exists on the layer below, the wave is
 * derived from the apexes of the layer below by normal projection ("分层路径映射复用"), so no
 * additional lateral phase offset is ever introduced ("锁定填充相位"). Newly appearing parts
 * and branches get a freshly laid out wave with segmented, fixed modulus.
 */
void buildLayerWaves(
    const Shape& outline,
    const coord_t cell_period,
    const coord_t line_width,
    const std::vector<PartRecord>& records_below,
    std::vector<PartRecord>& records_out,
    OpenLinesSet& lines_out)
{
    const coord_t discretization_step = std::max<coord_t>(200, cell_period / 8);

    // The longest of the parts whose wave was cancelled for being too short; if the whole layer
    // ends up empty, this one still gets a single compressed wave period ("至少每层生成一个周期").
    ReferenceAxis best_skipped_axis;
    SingleShape best_skipped_part;

    for (SingleShape& part : outline.splitIntoParts())
    {
        std::vector<MedialChain> chains = computeMedialChains(part, discretization_step);
        if (chains.empty())
        {
            continue;
        }
#ifdef SW_DEBUG
        {
            coord_t total = 0;
            for (const MedialChain& chain : chains)
            {
                total += chain.length;
            }
            std::fprintf(stderr, "part: %zu chains, total len=%lld\n", chains.size(), static_cast<long long>(total));
        }
#endif
        pruneMedialChains(chains, line_width);
#ifdef SW_DEBUG
        {
            coord_t total = 0;
            size_t kept = 0;
            for (const MedialChain& chain : chains)
            {
                if (! chain.removed)
                {
                    total += chain.length;
                    kept++;
                }
            }
            std::fprintf(stderr, "part: %zu chains kept after prune, total len=%lld\n", kept, static_cast<long long>(total));
        }
#endif

        std::vector<CenterLine> center_lines = assembleCenterLines(chains);
        stitchCenterLines(center_lines, cell_period);

        // Keep ONLY the primary (longest) reference axis of the part: one axis means one single
        // continuous, unbroken wave polyline per part and no side-branch waves that could cross
        // or overlap the main wave ("单条不间断折线, 无交叉网格, 线条不重叠").
        ReferenceAxis primary_axis;
        for (CenterLine& center_line : center_lines)
        {
            ReferenceAxis axis;
            axis.points = std::move(center_line.points);
            axis.radii = std::move(center_line.radii);
            finalizeAxisGeometry(axis);
            if (axis.total_length > primary_axis.total_length)
            {
                primary_axis = std::move(axis);
            }
        }
        // Cancel waves on parts too short for even one full period ("较短的填充线条可取消");
        // one period needs at least half a nominal cell (it is then compressed to fit).
        if (primary_axis.total_length < cell_period / 2 || primary_axis.points.size() < 2)
        {
            if (primary_axis.points.size() >= 2 && primary_axis.total_length > 2 * line_width && primary_axis.total_length > best_skipped_axis.total_length)
            {
                best_skipped_axis = std::move(primary_axis);
                best_skipped_part = part;
            }
            continue;
        }
        std::vector<ReferenceAxis> axes;
        axes.push_back(std::move(primary_axis));
        const Point2LL probe = axisMidpoint(axes.front());

        const PartRecord* below = findPartBelow(records_below, probe, part);
        if (below != nullptr && ! below->apexes.empty())
        {
            inheritApexesFromBelow(axes, below->apexes, cell_period, line_width);
        }
        else
        {
            layoutFreshApexes(axes.front(), cell_period);
        }
        computeUniformAmplitude(axes.front(), line_width);

#ifdef SW_DEBUG
        std::fprintf(
            stderr,
            "  axis len=%lld apexes=%zu amplitude=%lld inherited=%d\n",
            static_cast<long long>(axes.front().total_length),
            axes.front().apex_positions.size(),
            static_cast<long long>(axes.front().amplitude),
            below != nullptr && ! below->apexes.empty());
#endif

        OpenLinesSet part_lines;
        PartRecord record;
        {
            const ReferenceAxis& axis = axes.front();
            generateWaveAlongAxis(axis, part_lines);
            for (size_t apex_idx = 0; apex_idx < axis.apex_positions.size(); apex_idx++)
            {
                Point2LL base;
                Point2LL tangent;
                coord_t radius;
                sampleAxis(axis, axis.apex_positions[apex_idx], base, tangent, radius);
                record.apexes.push_back({ base, apexTip(axis, apex_idx) });
            }
        }

        // Straight wave flanks may cut across the boundary at sharp axis turns; clip to be safe.
        lines_out.push_back(part.intersection(part_lines, true, line_width));

        record.part = std::move(part);
        record.probe = probe;
        records_out.push_back(std::move(record));
    }

    // Nothing on this layer was long enough for a regular wave: place one single compressed
    // period on the longest cancelled part, so every layer carries at least one triangle wave.
    if (lines_out.empty() && best_skipped_axis.points.size() >= 2)
    {
        ReferenceAxis& axis = best_skipped_axis;
        axis.apex_positions = { axis.total_length / 4, axis.total_length * 3 / 4 };
        axis.apex_sides = { true, false };
        computeUniformAmplitude(axis, line_width);

        OpenLinesSet part_lines;
        generateWaveAlongAxis(axis, part_lines);
        lines_out.push_back(best_skipped_part.intersection(part_lines, true, line_width));

        PartRecord record;
        for (size_t apex_idx = 0; apex_idx < axis.apex_positions.size(); apex_idx++)
        {
            Point2LL base;
            Point2LL tangent;
            coord_t radius;
            sampleAxis(axis, axis.apex_positions[apex_idx], base, tangent, radius);
            record.apexes.push_back({ base, apexTip(axis, apex_idx) });
        }
        record.probe = axisMidpoint(axis);
        record.part = std::move(best_skipped_part);
        records_out.push_back(std::move(record));
    }
}

} // namespace

SurfaceWaveFillProvider::SurfaceWaveFillProvider(const std::vector<Shape>& layer_outlines, const coord_t line_distance, const coord_t line_width, const coord_t wall_clearance)
    : line_width_(line_width)
{
    // One wave cell (full period) spans two flank spacings along the axis.
    const coord_t cell_period = std::max<coord_t>(2 * line_distance, 2 * line_width);

    layer_waves_.resize(layer_outlines.size());
    std::vector<PartRecord> records_below;
    for (size_t layer_idx = 0; layer_idx < layer_outlines.size(); layer_idx++)
    {
        // Constrain the wave by the walls: keep the wall clearance (extra infill walls plus
        // minimum wall line width) away from the boundary of the available area.
        const Shape constrained_outline = wall_clearance > 0 ? layer_outlines[layer_idx].offset(-wall_clearance) : layer_outlines[layer_idx];
        std::vector<PartRecord> records;
        buildLayerWaves(constrained_outline, cell_period, line_width_, records_below, records, layer_waves_[layer_idx]);
        records_below = std::move(records);
    }
}

void SurfaceWaveFillProvider::generate(OpenLinesSet& result_lines, const Shape& in_outline, const size_t layer_idx) const
{
    if (layer_idx >= layer_waves_.size())
    {
        return;
    }
    result_lines.push_back(in_outline.intersection(layer_waves_[layer_idx], true, line_width_));
}

void SurfaceWaveFillProvider::generateSingleLayer(OpenLinesSet& result_lines, const Shape& in_outline, const coord_t line_distance, const coord_t line_width)
{
    const coord_t cell_period = std::max<coord_t>(2 * line_distance, 2 * line_width);
    std::vector<PartRecord> no_records_below;
    std::vector<PartRecord> records;
    buildLayerWaves(in_outline, cell_period, line_width, no_records_below, records, result_lines);
}

} // namespace cura
