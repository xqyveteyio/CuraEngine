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
 *
 * Additionally, branches into small appendage regions attached to the boundary (rib pockets,
 * bosses -- any lobe whose effective span stays below half a wave cell) are removed as well:
 * the centerline is only kept inside the effective MAIN fill region ("中轴线仅在有效主填充区域内
 * 生成"), and where such a branch was produced by the skeleton it is deleted, truncating the
 * main axis at the junction with the main region ("删除该短支, 主中轴线截断于连接处").
 */
void pruneMedialChains(std::vector<MedialChain>& chains, const coord_t line_width, const coord_t cell_period)
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
            if (branch_length < std::max(junction_radius * 3 / 2 + line_width, cell_period / 2))
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
 * first appears, the apex positions are laid out with the global fixed period; on every
 * following layer they are obtained by projecting the apexes of the layer below onto the
 * current axis, so the wave translates and bends together with the surface while the apex
 * count and the relative phase stay locked. The amplitude of each apex is not stored: every
 * peak/valley reaches out from the axis to the wall it faces ("波峰波谷与壁接触"), so the wave
 * amplitude follows the local cross-section width and outline shape automatically.
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
 * The width-normalized parameter profile of an axis ("归一化弧长参数"): cumulative
 * tau(s) = integral of ds / width(s) at every axis vertex, with the width taken as both medial
 * radii summed. This is the stable per-region parameter domain in which the wave phase lives:
 * uniform spacing in tau equals one constant apex angle, and a position expressed as a tau
 * FRACTION can be transferred between the cross sections of consecutive layers without using
 * any spatial coordinates.
 */
std::vector<double> tauProfile(const ReferenceAxis& axis)
{
    std::vector<double> profile(axis.points.size(), 0.0);
    for (size_t i = 0; i + 1 < axis.points.size(); i++)
    {
        const double width = std::max<double>(static_cast<double>(axis.radii[i] + axis.radii[i + 1]), 1.0);
        profile[i + 1] = profile[i] + static_cast<double>(axis.cumulative[i + 1] - axis.cumulative[i]) / width;
    }
    return profile;
}

//! The width-normalized total length of an axis (the end value of the tau profile).
double tauTotal(const ReferenceAxis& axis)
{
    return tauProfile(axis).back();
}

//! Convert an (absolute) tau value into the arc length along the axis, by linear interpolation of the profile.
coord_t tauToArc(const ReferenceAxis& axis, const std::vector<double>& profile, const double tau)
{
    if (tau <= 0.0)
    {
        return 0;
    }
    size_t i = 0;
    while (i + 2 < profile.size() && profile[i + 1] < tau)
    {
        i++;
    }
    if (tau >= profile.back())
    {
        return axis.total_length;
    }
    const double along = (tau - profile[i]) / std::max(profile[i + 1] - profile[i], 1e-12);
    return axis.cumulative[i] + std::llrint(std::clamp(along, 0.0, 1.0) * static_cast<double>(axis.cumulative[i + 1] - axis.cumulative[i]));
}

//! Convert an arc length along the axis into the (absolute) tau value, by linear interpolation of the profile.
double arcToTau(const ReferenceAxis& axis, const std::vector<double>& profile, const coord_t s)
{
    if (s <= 0)
    {
        return 0.0;
    }
    size_t i = 0;
    while (i + 2 < axis.cumulative.size() && axis.cumulative[i + 1] < s)
    {
        i++;
    }
    if (s >= axis.total_length)
    {
        return profile.back();
    }
    const double along = static_cast<double>(s - axis.cumulative[i]) / std::max<double>(static_cast<double>(axis.cumulative[i + 1] - axis.cumulative[i]), 1.0);
    return profile[i] + std::clamp(along, 0.0, 1.0) * (profile[i + 1] - profile[i]);
}

//! Distribute \p apex_count apexes uniformly in the width-normalized parameter tau (given at every axis vertex), alternating sides.
void distributeApexesByTau(ReferenceAxis& axis, const std::vector<double>& tau, const size_t apex_count)
{
    axis.apex_positions.clear();
    axis.apex_sides.clear();
    const double delta_tau = tau.back() / static_cast<double>(apex_count);
    bool side_positive = true;
    size_t segment_idx = 0;
    for (size_t apex_idx = 0; apex_idx < apex_count; apex_idx++)
    {
        const double target = (static_cast<double>(apex_idx) + 0.5) * delta_tau;
        while (segment_idx + 2 < tau.size() && tau[segment_idx + 1] < target)
        {
            segment_idx++;
        }
        const double along = (target - tau[segment_idx]) / std::max(tau[segment_idx + 1] - tau[segment_idx], 1e-12);
        const coord_t s
            = axis.cumulative[segment_idx] + std::llrint(std::clamp(along, 0.0, 1.0) * static_cast<double>(axis.cumulative[segment_idx + 1] - axis.cumulative[segment_idx]));
        axis.apex_positions.push_back(s);
        axis.apex_sides.push_back(side_positive);
        side_positive = ! side_positive;
    }
}

Point2LL apexTip(const ReferenceAxis& axis, size_t apex_idx, const SingleShape& part); // Defined below.

/*!
 * Lay out a fresh wave on an axis that has no predecessor on the layer below, with ONE apex
 * angle for the whole region ("三角波顶角角度大小一致"). Since every apex tip contacts the wall
 * across the local width, a constant apex angle requires the spacing between consecutive
 * apexes to be proportional to the local cross-section width ("根据区域宽度自动改变节距，但顶角
 * 角度保持一致"). The apexes are therefore distributed uniformly in the width-normalized arc
 * length tau(s) = integral of ds / width(s), which makes delta-tau -- and thereby the tangent
 * of the half apex angle -- one constant everywhere; a second pass re-measures the width from
 * the actual wall-contact tips and redistributes. The apex count is chosen such that the
 * AVERAGE spacing matches the nominal half period from the infill density, so the density
 * setting is what sets the apex angle.
 */
/*!
 * Redistribute a FIXED number of apexes over the axis such that ALL apex angles are equal
 * ("所有尖角角度不变"): uniform in the width-normalized parameter, then refined twice with the
 * actually measured wall-contact distances. NOTE: this resets the sides to alternate starting
 * positive; the caller re-aligns the parity when a phase has to be preserved.
 */
void redistributeApexesEqualAngle(ReferenceAxis& axis, const SingleShape& part, const size_t apex_count)
{
    // First pass: width-normalized cumulative parameter at every axis vertex, with the width
    // (the local tip-to-tip span of the wave) estimated as both medial radii summed.
    std::vector<double> tau(axis.points.size(), 0.0);
    for (size_t i = 0; i + 1 < axis.points.size(); i++)
    {
        const double width = std::max<double>(static_cast<double>(axis.radii[i] + axis.radii[i + 1]), 1.0);
        tau[i + 1] = tau[i] + static_cast<double>(axis.cumulative[i + 1] - axis.cumulative[i]) / width;
    }
    distributeApexesByTau(axis, tau, apex_count);

    // Refinement: the actual wall-contact distance of a tip can deviate from the medial radius
    // (skewed cross sections, bumpy outlines). Measure the real tip distances and redistribute
    // with the measured widths, homogenizing the apex angle further.
    for (size_t iteration = 0; iteration < 2; iteration++)
    {
        std::vector<coord_t> tip_distance(axis.apex_positions.size());
        for (size_t apex_idx = 0; apex_idx < axis.apex_positions.size(); apex_idx++)
        {
            Point2LL center;
            Point2LL tangent;
            coord_t radius;
            sampleAxis(axis, axis.apex_positions[apex_idx], center, tangent, radius);
            tip_distance[apex_idx] = vSize(apexTip(axis, apex_idx, part) - center);
        }
        // Effective width per axis vertex: interpolated from the measured spans of the
        // enclosing flank (the transverse extent between two consecutive opposite-side tips).
        std::vector<double> refined_tau(axis.points.size(), 0.0);
        size_t apex_idx = 0;
        for (size_t i = 0; i + 1 < axis.points.size(); i++)
        {
            const coord_t s_mid = (axis.cumulative[i] + axis.cumulative[i + 1]) / 2;
            while (apex_idx + 2 < axis.apex_positions.size() && axis.apex_positions[apex_idx + 1] < s_mid)
            {
                apex_idx++;
            }
            const double width = std::max<double>(static_cast<double>(tip_distance[apex_idx] + tip_distance[apex_idx + 1]), 1.0);
            refined_tau[i + 1] = refined_tau[i] + static_cast<double>(axis.cumulative[i + 1] - axis.cumulative[i]) / width;
        }
        distributeApexesByTau(axis, refined_tau, apex_count);
    }
}

void layoutFreshApexes(ReferenceAxis& axis, const coord_t cell_period, const SingleShape& part)
{
    axis.apex_positions.clear();
    axis.apex_sides.clear();

    const coord_t half_period = cell_period / 2;
    if (axis.total_length < cell_period || half_period <= 0)
    {
        axis.apex_positions = { axis.total_length / 4, axis.total_length * 3 / 4 };
        axis.apex_sides = { true, false };
        return;
    }

    const size_t apex_count = std::max<size_t>(2, static_cast<size_t>(axis.total_length / half_period));
    redistributeApexesEqualAngle(axis, part, apex_count);
}

/*!
 * Cast a ray from \p from in direction \p direction and return the nearest intersection with
 * the boundary of \p part (outer wall or hole), pulled back inside by a tiny margin so the
 * point survives later clipping. Returns \p fallback when nothing is hit within \p max_length.
 */
Point2LL castToBoundary(const SingleShape& part, const Point2LL& from, const Point2LL& direction, const coord_t max_length, const Point2LL& fallback)
{
    const double dx = static_cast<double>(direction.X);
    const double dy = static_cast<double>(direction.Y);
    const double dir_length = std::hypot(dx, dy);
    if (dir_length <= 0.0)
    {
        return fallback;
    }
    double best_t = -1.0;
    for (const Polygon& poly : part)
    {
        for (size_t i = 0; i < poly.size(); i++)
        {
            const Point2LL& a = poly[i];
            const Point2LL& b = poly[(i + 1) % poly.size()];
            const double ex = static_cast<double>(b.X - a.X);
            const double ey = static_cast<double>(b.Y - a.Y);
            const double det = ex * dy - dx * ey;
            if (det == 0.0)
            {
                continue; // Ray parallel to this boundary segment.
            }
            const double fx = static_cast<double>(a.X - from.X);
            const double fy = static_cast<double>(a.Y - from.Y);
            const double t = (ex * fy - fx * ey) / det; // Distance along the ray, in units of |direction|.
            const double u = (dx * fy - dy * fx) / det; // Position within the segment.
            if (t > 0.0 && u >= 0.0 && u <= 1.0 && (best_t < 0.0 || t < best_t))
            {
                best_t = t;
            }
        }
    }
    if (best_t < 0.0 || best_t * dir_length > static_cast<double>(max_length))
    {
        return fallback;
    }
    // Stay a hair inside the boundary, so the contact point is not eaten by the clip.
    const double pullback = std::min(10.0 / dir_length, best_t);
    best_t -= pullback;
    return from + Point2LL(std::llrint(dx * best_t), std::llrint(dy * best_t));
}

//! One intersection of a path segment with the region boundary.
struct BoundaryCrossing
{
    double t; //!< Position along the path segment (0..1).
    size_t poly_idx; //!< Which polygon of the part (outer boundary or hole) is crossed.
    size_t edge_idx; //!< Which edge of that polygon.
    Point2LL point; //!< The intersection point.
};

//! All crossings of segment \p p - \p q with the boundary of \p part, ordered along the segment.
std::vector<BoundaryCrossing> collectBoundaryCrossings(const SingleShape& part, const Point2LL& p, const Point2LL& q)
{
    std::vector<BoundaryCrossing> crossings;
    const double dx = static_cast<double>(q.X - p.X);
    const double dy = static_cast<double>(q.Y - p.Y);
    for (size_t poly_idx = 0; poly_idx < part.size(); poly_idx++)
    {
        const Polygon& poly = part[poly_idx];
        for (size_t i = 0; i < poly.size(); i++)
        {
            const Point2LL& a = poly[i];
            const Point2LL& b = poly[(i + 1) % poly.size()];
            const double ex = static_cast<double>(b.X - a.X);
            const double ey = static_cast<double>(b.Y - a.Y);
            const double det = ex * dy - dx * ey;
            if (det == 0.0)
            {
                continue;
            }
            const double fx = static_cast<double>(a.X - p.X);
            const double fy = static_cast<double>(a.Y - p.Y);
            const double t = (ex * fy - fx * ey) / det;
            const double u = (dx * fy - dy * fx) / det;
            if (t > 0.0 && t < 1.0 && u >= 0.0 && u < 1.0)
            {
                crossings.push_back({ t, poly_idx, i, p + Point2LL(std::llrint(dx * t), std::llrint(dy * t)) });
            }
        }
    }
    std::sort(
        crossings.begin(),
        crossings.end(),
        [](const BoundaryCrossing& a, const BoundaryCrossing& b)
        {
            return a.t < b.t;
        });
    return crossings;
}

/*!
 * Ensure a path point lies STRICTLY inside the region: a point on or just outside the boundary
 * (an apex tip bound exactly to the wall, or a plateau end sticking out past a curved wall) is
 * moved onto the nearest boundary edge and pulled a hair to its inside, so the final clip can
 * never cut the path apart at such a point.
 */
Point2LL pullInside(const SingleShape& part, const Point2LL& p)
{
    Point2LL result = p;
    for (size_t attempt = 0; attempt < 3; attempt++)
    {
        if (part.inside(result, false))
        {
            return result;
        }
        coord_t best_distance = std::numeric_limits<coord_t>::max();
        Point2LL moved = result;
        for (const Polygon& poly : part)
        {
            for (size_t i = 0; i < poly.size(); i++)
            {
                const Point2LL closest = LinearAlg2D::getClosestOnLineSegment(result, poly[i], poly[(i + 1) % poly.size()]);
                const coord_t d = vSize(closest - result);
                if (d < best_distance)
                {
                    best_distance = d;
                    // The material is to the left of the edge direction, for outer boundaries and holes alike.
                    moved = closest + normal(turn90CCW(poly[(i + 1) % poly.size()] - poly[i]), 25 * (static_cast<coord_t>(attempt) + 1));
                }
            }
        }
        result = moved;
    }
    return result;
}

//! A point on a boundary edge, pulled slightly to the inside of the region (the material is to the left of the edge direction, for outer boundaries and holes alike).
Point2LL insetOntoEdge(const Polygon& poly, const size_t edge_idx, const Point2LL& point)
{
    const Point2LL direction = poly[(edge_idx + 1) % poly.size()] - poly[edge_idx];
    return point + normal(turn90CCW(direction), 15);
}

//! A polygon vertex, pulled slightly to the inside of the region along its angle bisector.
Point2LL insetVertex(const Polygon& poly, const size_t vertex_idx)
{
    const size_t n = poly.size();
    const Point2LL previous = poly[(vertex_idx + n - 1) % n];
    const Point2LL vertex = poly[vertex_idx];
    const Point2LL next = poly[(vertex_idx + 1) % n];
    const Point2LL inward = normal(turn90CCW(vertex - previous), 1000) + normal(turn90CCW(next - vertex), 1000);
    return vertex + normal(inward, 15);
}

/*!
 * Append the walk along polygon \p poly from crossing \p from to crossing \p to (both on that
 * polygon), taking the shorter of the two possible directions, to \p out. Every emitted point
 * is pulled slightly inside, so the resulting path hugs the wall without leaving the region.
 */
void appendBoundaryWalk(const Polygon& poly, const BoundaryCrossing& from, const BoundaryCrossing& to, std::vector<Point2LL>& out)
{
    const size_t n = poly.size();

    // Cumulative arc length of every vertex around the polygon.
    std::vector<coord_t> cumulative(n + 1, 0);
    for (size_t i = 0; i < n; i++)
    {
        cumulative[i + 1] = cumulative[i] + vSize(poly[(i + 1) % n] - poly[i]);
    }
    const coord_t perimeter = cumulative[n];
    if (perimeter <= 0)
    {
        return;
    }
    const coord_t pos_from = cumulative[from.edge_idx] + vSize(from.point - poly[from.edge_idx]);
    const coord_t pos_to = cumulative[to.edge_idx] + vSize(to.point - poly[to.edge_idx]);
    const coord_t forward_distance = (pos_to - pos_from + perimeter) % perimeter;
    const bool forward = forward_distance <= perimeter - forward_distance;

    out.push_back(insetOntoEdge(poly, from.edge_idx, from.point));
    if (forward)
    {
        if (from.edge_idx != to.edge_idx)
        {
            // Forward, the vertices passed are from.edge_idx+1 up to and including to.edge_idx.
            size_t edge = (from.edge_idx + 1) % n;
            while (true)
            {
                out.push_back(insetVertex(poly, edge));
                if (edge == to.edge_idx)
                {
                    break;
                }
                edge = (edge + 1) % n;
            }
        }
    }
    else
    {
        // Backward, the vertices passed are from.edge_idx down to and including to.edge_idx+1.
        for (size_t edge = from.edge_idx; edge != to.edge_idx;)
        {
            out.push_back(insetVertex(poly, edge));
            edge = (edge + n - 1) % n;
        }
    }
    out.push_back(insetOntoEdge(poly, to.edge_idx, to.point));
}

/*!
 * The 2D tip of one apex: the point where the wave peak/valley CONTACTS the wall it faces
 * ("波峰波谷与壁接触"). Found by casting a ray from the axis perpendicularly (the perpendicular
 * bisector of the apex triangle is perpendicular to the medial axis, "中垂线⊥中轴线") to the
 * nearest boundary of the part, so the amplitude follows the local cross-section width and
 * outline shape automatically ("随局部截面宽度和轮廓形状自动变幅").
 */
Point2LL apexTip(const ReferenceAxis& axis, const size_t apex_idx, const SingleShape& part)
{
    Point2LL center;
    Point2LL tangent;
    coord_t radius;
    sampleAxis(axis, axis.apex_positions[apex_idx], center, tangent, radius);
    const Point2LL perpendicular = turn90CCW(tangent);
    const Point2LL offset_direction = axis.apex_sides[apex_idx] ? perpendicular : perpendicular * -1;
    const Point2LL fallback = center + normal(offset_direction, radius);
    // The reach is limited to twice the local half-width: enough to contact the facing wall
    // even on skewed cross sections, but not enough to shoot through the narrow mouth of a
    // boundary pocket (e.g. between ribs) and bind to a wall far inside it.
    return castToBoundary(part, center, offset_direction, 2 * radius + 1000, fallback);
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
 * unbroken polyline that visits every apex tip in order. Every tip touches the wall it faces,
 * with peaks and valleys alternating between the two sides of the axis. The first and the last
 * line of the path keep their own direction and are prolonged until they contact the wall
 * ("首、末线条按首段/末段线条方向延伸至壁"), so both path ends are connected to the innermost
 * wall without changing the wave angle or adding a period.
 *
 * The result contains ONLY the continuous triangle wave itself and, where geometry forces it,
 * its tip-truncated form ("只有正交变换的连续三角波, 或者被削顶的连续三角波"): degenerate spikes
 * are cut off, and any stretch that would leave the region is replaced by the walk along the
 * region boundary, connecting the break points along the wall without crossing holes.
 */
void generateWaveAlongAxis(const ReferenceAxis& axis, const SingleShape& part, const coord_t line_width, OpenLinesSet& result_lines)
{
    if (axis.points.size() < 2 || axis.apex_positions.size() < 2)
    {
        return;
    }

    std::vector<Point2LL> tips;
    for (size_t apex_idx = 0; apex_idx < axis.apex_positions.size(); apex_idx++)
    {
        tips.push_back(apexTip(axis, apex_idx, part));
    }

    std::vector<Point2LL> points;
    {
        // First line of the path: from the first apex tip through the axis start point,
        // extended in that same direction until it contacts the wall at the region's end.
        const Point2LL anchor = axis.points.front();
        const coord_t reach = 4 * axis.radii.front() + 4 * vSize(anchor - tips.front()) + 1000;
        points.push_back(castToBoundary(part, anchor, anchor - tips.front(), reach, anchor));
    }
    for (size_t apex_idx = 0; apex_idx < tips.size(); apex_idx++)
    {
        Point2LL center;
        Point2LL tangent;
        coord_t radius;
        sampleAxis(axis, axis.apex_positions[apex_idx], center, tangent, radius);

        // Effective overlap instead of exact coincidence ("有效搭接"): each peak gets a SHORT
        // PLATEAU along the wall instead of an ideal sharp point, so a small lateral phase
        // shift between consecutive layers still leaves the plateau areas overlapping.
        const coord_t prev_gap = (apex_idx == 0) ? axis.apex_positions[apex_idx] : axis.apex_positions[apex_idx] - axis.apex_positions[apex_idx - 1];
        const coord_t next_gap
            = (apex_idx + 1 == tips.size()) ? axis.total_length - axis.apex_positions[apex_idx] : axis.apex_positions[apex_idx + 1] - axis.apex_positions[apex_idx];
        const coord_t plateau_half = line_width / 2;
        if (std::min(prev_gap, next_gap) > 4 * line_width && plateau_half > 0)
        {
            const Point2LL along_wall = normal(tangent, plateau_half);
            points.push_back(tips[apex_idx] - along_wall);
            points.push_back(tips[apex_idx] + along_wall);
        }
        else
        {
            points.push_back(tips[apex_idx]);
        }
    }
    {
        // Last line of the path, extended likewise beyond the axis end point.
        const Point2LL anchor = axis.points.back();
        const coord_t reach = 4 * axis.radii.back() + 4 * vSize(anchor - tips.back()) + 1000;
        points.push_back(castToBoundary(part, anchor, anchor - tips.back(), reach, anchor));
    }

    // Truncate degenerate spikes ("削掉三角波尖点"): where the two flank lines of an apex fold
    // onto each other (within a 5% overlap error, i.e. the apex angle is below ~2*atan(0.05)),
    // the tip would be a doubled-up line squeezed into a pinch. The tip is cut off where the
    // flanks are one line width apart; the resulting cut runs across the pinch, and if it lies
    // outside the region it is re-routed along the boundary below.
    constexpr double spike_apex_angle = 0.1; // Radians; 2*atan(0.05) for the 5% overlap criterion.
    for (size_t i = 1; i + 1 < points.size(); i++)
    {
        const Point2LL to_prev = points[i - 1] - points[i];
        const Point2LL to_next = points[i + 1] - points[i];
        const double len_prev = std::hypot(static_cast<double>(to_prev.X), static_cast<double>(to_prev.Y));
        const double len_next = std::hypot(static_cast<double>(to_next.X), static_cast<double>(to_next.Y));
        if (len_prev <= 0.0 || len_next <= 0.0)
        {
            continue;
        }
        const double cosine = std::clamp((static_cast<double>(to_prev.X) * to_next.X + static_cast<double>(to_prev.Y) * to_next.Y) / (len_prev * len_next), -1.0, 1.0);
        const double angle = std::acos(cosine);
        if (angle >= spike_apex_angle)
        {
            continue;
        }
        const coord_t cut = std::min<coord_t>(
            static_cast<coord_t>(static_cast<double>(line_width) / std::max(angle, 0.01)),
            static_cast<coord_t>(std::min(len_prev, len_next) * 0.6)); // At most 60% of the shorter flank.
        if (cut < line_width)
        {
            continue;
        }
        const Point2LL cut_a = points[i] + normal(to_prev, cut);
        const Point2LL cut_b = points[i] + normal(to_next, cut);
        points[i] = cut_a;
        points.insert(points.begin() + static_cast<std::ptrdiff_t>(i) + 1, cut_b);
        i++; // Skip the just-inserted point.
    }

    // Wall contact means the tips sit ON the boundary (and plateau ends can stick out past a
    // curved wall); pull every path point strictly inside, so no point of the single continuous
    // path can be nibbled off by the final clip.
    for (Point2LL& p : points)
    {
        p = pullInside(part, p);
    }

    // Edge break points connect along the wall ("边缘断点沿壁连接, 不穿越孔洞"): every stretch of
    // the path that would leave the region (a flank cutting across a boundary concavity such as
    // a rib mouth, a truncated tip, or a stretch crossing a hole) is replaced by the walk along
    // the region boundary between the exit and the re-entry point. The path then consists solely
    // of triangle wave flanks and wall-following connections, stays in one piece, and never
    // enters holes or non-printable areas.
    {
        std::vector<Point2LL> conformed;
        conformed.reserve(points.size());
        conformed.push_back(points.front());
        for (size_t i = 0; i + 1 < points.size(); i++)
        {
            const std::vector<BoundaryCrossing> crossings = collectBoundaryCrossings(part, points[i], points[i + 1]);
            for (size_t k = 0; k + 1 < crossings.size(); k++)
            {
                const Point2LL middle = (crossings[k].point + crossings[k + 1].point) / 2;
                if (part.inside(middle, true))
                {
                    continue; // This stretch stays inside; only the stretches outside are re-routed.
                }
                if (crossings[k].poly_idx == crossings[k + 1].poly_idx)
                {
                    appendBoundaryWalk(part[crossings[k].poly_idx], crossings[k], crossings[k + 1], conformed);
                }
            }
            conformed.push_back(points[i + 1]);
        }
        points = std::move(conformed);
    }
#ifdef SW_DEBUG
    for (size_t i = 0; i < points.size(); i++)
    {
        if (! part.inside(points[i], true))
        {
            std::fprintf(stderr, "POST-CONFORM POINT OUTSIDE: %zu/%zu (%lld,%lld)\n", i, points.size(), static_cast<long long>(points[i].X), static_cast<long long>(points[i].Y));
        }
    }
    for (size_t i = 0; i + 1 < points.size(); i++)
    {
        const std::vector<BoundaryCrossing> crossings = collectBoundaryCrossings(part, points[i], points[i + 1]);
        for (size_t k = 0; k + 1 < crossings.size(); k++)
        {
            const Point2LL middle = (crossings[k].point + crossings[k + 1].point) / 2;
            if (! part.inside(middle, true))
            {
                std::fprintf(
                    stderr,
                    "POST-CONFORM OUTSIDE: seg %zu (%lld,%lld)-(%lld,%lld) out (%lld,%lld)-(%lld,%lld) polys %zu/%zu\n",
                    i,
                    static_cast<long long>(points[i].X),
                    static_cast<long long>(points[i].Y),
                    static_cast<long long>(points[i + 1].X),
                    static_cast<long long>(points[i + 1].Y),
                    static_cast<long long>(crossings[k].point.X),
                    static_cast<long long>(crossings[k].point.Y),
                    static_cast<long long>(crossings[k + 1].point.X),
                    static_cast<long long>(crossings[k + 1].point.Y),
                    crossings[k].poly_idx,
                    crossings[k + 1].poly_idx);
            }
        }
    }
#endif
    removeSelfCrossings(points);
#ifdef SW_DEBUG
    for (size_t i = 0; i + 1 < points.size(); i++)
    {
        const std::vector<BoundaryCrossing> crossings = collectBoundaryCrossings(part, points[i], points[i + 1]);
        for (size_t k = 0; k + 1 < crossings.size(); k++)
        {
            const Point2LL middle = (crossings[k].point + crossings[k + 1].point) / 2;
            if (! part.inside(middle, true))
            {
                std::fprintf(
                    stderr,
                    "POST-DECROSS OUTSIDE: seg %zu out (%lld,%lld)-(%lld,%lld)\n",
                    i,
                    static_cast<long long>(crossings[k].point.X),
                    static_cast<long long>(crossings[k].point.Y),
                    static_cast<long long>(crossings[k + 1].point.X),
                    static_cast<long long>(crossings[k + 1].point.Y));
            }
        }
    }
#endif

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
    double delta_tau{ 0.0 }; //!< Width-normalized apex spacing; uniquely determines the target apex angle. Fixed at part birth and inherited unchanged.

    // The continuous phase field of this part ("三维连续相位场"), one slice of it: where each
    // apex sits as a FRACTION of the width-normalized parameter domain, plus the geometric
    // orientation of that domain. The next layer transfers the phase through this parameter
    // mapping instead of projecting spatial coordinates.
    Point2LL axis_front;
    Point2LL axis_back;
    std::vector<double> apex_tau_fractions; //!< Normalized (0..1) tau position of every apex.
    std::vector<bool> apex_sides; //!< Side of every apex, relative to the front->back direction of the recorded axis.
};

Point2LL axisMidpoint(const ReferenceAxis& axis)
{
    Point2LL position;
    Point2LL tangent;
    coord_t radius;
    sampleAxis(axis, axis.total_length / 2, position, tangent, radius);
    return position;
}

//! The minimum distance from a point to the boundary of a shape.
coord_t distanceToShapeBoundary(const SingleShape& shape, const Point2LL& p)
{
    coord_t best = std::numeric_limits<coord_t>::max();
    for (const Polygon& poly : shape)
    {
        for (size_t i = 0; i < poly.size(); i++)
        {
            const Point2LL closest = LinearAlg2D::getClosestOnLineSegment(p, poly[i], poly[(i + 1) % poly.size()]);
            best = std::min(best, vSize(closest - p));
        }
    }
    return best;
}

/*!
 * Find the part of the layer below that the given part derives its phase from ("父区域"):
 *  - Normally the part it stands on (probe containment, both directions).
 *  - When several parts of the layer below merged into this part, the DOMINANT parent (most
 *    apexes) provides the phase; the other parents' phases simply end where they meet it, which
 *    places the unavoidable phase seam of a merge at the junction, a single low-impact defect
 *    instead of a global re-layout ("区域合并: 选择低影响位置设置相位缝").
 *  - A part that stands on nothing inherits from the geometrically NEAREST part of the layer
 *    below within \p search_radius, extrapolating that parent's phase field sideways instead of
 *    initializing an unrelated phase next to it ("新区域从父区域的相位场外推").
 */
const PartRecord* findPartBelow(const std::vector<PartRecord>& records_below, const Point2LL& current_probe, const SingleShape& current_part, const coord_t search_radius)
{
    const PartRecord* dominant = nullptr;
    for (const PartRecord& record : records_below)
    {
        if (record.part.inside(current_probe, true) || current_part.inside(record.probe, true))
        {
            if (dominant == nullptr || record.apexes.size() > dominant->apexes.size())
            {
                dominant = &record;
            }
        }
    }
    if (dominant != nullptr)
    {
        return dominant;
    }
    coord_t best_distance = search_radius;
    for (const PartRecord& record : records_below)
    {
        const coord_t d = distanceToShapeBoundary(record.part, current_probe);
        if (d < best_distance)
        {
            best_distance = d;
            dominant = &record;
        }
    }
    return dominant;
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
void inheritApexesFromBelow(
    std::vector<ReferenceAxis>& axes,
    const std::vector<ApexFeature>& below_apexes,
    const coord_t cell_period,
    const coord_t line_width,
    const SingleShape& part)
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
            layoutFreshApexes(axis, cell_period, part); // Newly appeared branch: no wave below to inherit from.
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
            layoutFreshApexes(axis, cell_period, part);
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
 * Transfer the phase of the layer below through the PARAMETER MAPPING between the two cross
 * sections ("截面间参数映射"): every apex of the layer below is carried over as its normalized
 * width-normalized-arc-length fraction, NOT as a spatial coordinate. This keeps the apex order
 * and the equal-angle property intact by construction, and is insensitive to the per-layer
 * jitter of the extracted skeleton. The only spatial information used is the pairing of the two
 * axes' endpoints, which fixes the orientation of the parameter domain.
 *
 * Only applicable when the two cross sections correspond one to one; a strong length mismatch
 * (region split or merge) makes the caller fall back to spatial projection for this layer.
 */
bool transferApexesByParameter(ReferenceAxis& axis, const PartRecord& below)
{
    if (below.apex_tau_fractions.size() < 2 || below.apex_tau_fractions.size() != below.apex_sides.size() || axis.points.size() < 2)
    {
        return false;
    }

    // Orient the current parameter domain like the recorded one, by pairing the axis endpoints.
    const coord_t same = vSize(axis.points.front() - below.axis_front) + vSize(axis.points.back() - below.axis_back);
    const coord_t flipped = vSize(axis.points.front() - below.axis_back) + vSize(axis.points.back() - below.axis_front);
    const bool reversed = flipped < same;

    // The endpoints must roughly correspond; if they moved by a large part of the axis length,
    // the region topology changed and the parameter domains are not comparable.
    if (std::min(same, flipped) > axis.total_length)
    {
        return false;
    }

    const std::vector<double> profile = tauProfile(axis);
    axis.apex_positions.clear();
    axis.apex_sides.clear();
    for (size_t i = 0; i < below.apex_tau_fractions.size(); i++)
    {
        const size_t src = reversed ? below.apex_tau_fractions.size() - 1 - i : i;
        const double fraction = reversed ? 1.0 - below.apex_tau_fractions[src] : below.apex_tau_fractions[src];
        // Reversing the axis direction also mirrors left and right.
        const bool side = reversed ? ! below.apex_sides[src] : below.apex_sides[src];
        axis.apex_positions.push_back(tauToArc(axis, profile, fraction * profile.back()));
        axis.apex_sides.push_back(side);
    }
    return true;
}

/*!
 * Distributed phase correction ("分布式相位校正"): the practical, linearized form of the global
 * optimization -- inter-layer correspondence is kept (the inherited positions are the starting
 * point and may only move a bounded distance per layer, preserving the effective overlap), and
 * the equal-apex-angle layout acts as a soft attractor of which only a small share of the
 * residual error is corrected per layer, so deviations bleed away over many layers instead of
 * causing one sudden lateral jump ("每层只修正剩余误差的一部分, 避免漂移积累后突然重置").
 *
 * The apex COUNT is controlled here as well: when the region's parameter length demands a
 * different count (perimeter growing/shrinking with height), at most ONE full period (a pair of
 * apexes, keeping the strict side alternation) is inserted into the widest interior gap -- a
 * peak gradually splitting in two over the following layers ("受控分叉") -- or the closest
 * interior pair is removed ("受控合并"), with a hysteresis of one apex so the count does not
 * oscillate. End gaps are padded/trimmed separately where the axis grew or shrank.
 */
void applyDistributedPhaseCorrection(ReferenceAxis& axis, const double delta_tau, const coord_t line_width)
{
    if (axis.apex_positions.size() < 2 || delta_tau <= 0.0)
    {
        return;
    }
    const std::vector<double> profile = tauProfile(axis);
    const double tau_total = profile.back();

    // Work in the parameter domain.
    std::vector<double> taus(axis.apex_positions.size());
    for (size_t i = 0; i < taus.size(); i++)
    {
        taus[i] = arcToTau(axis, profile, axis.apex_positions[i]);
    }
    std::vector<bool> sides = axis.apex_sides;

    // Pad the ends where the axis grew into new territory, with the nominal spacing.
    while (taus.front() > 1.5 * delta_tau)
    {
        taus.insert(taus.begin(), taus.front() - delta_tau);
        sides.insert(sides.begin(), ! sides.front());
    }
    while (tau_total - taus.back() > 1.5 * delta_tau)
    {
        taus.push_back(taus.back() + delta_tau);
        sides.push_back(! sides.back());
    }

    // Controlled split / merge: at most one apex PAIR per layer, in the region interior.
    const size_t target_count = std::max<size_t>(2, static_cast<size_t>(std::llround(tau_total / delta_tau)));
    if (target_count >= taus.size() + 2 && taus.size() >= 2)
    {
        size_t widest = 0;
        for (size_t i = 1; i + 1 < taus.size(); i++) // Interior gaps only.
        {
            if (i + 1 < taus.size() && taus[i + 1] - taus[i] > taus[widest + 1] - taus[widest])
            {
                widest = i;
            }
        }
        const double gap = taus[widest + 1] - taus[widest];
        taus.insert(taus.begin() + static_cast<std::ptrdiff_t>(widest) + 1, { taus[widest] + gap / 3.0, taus[widest] + gap * 2.0 / 3.0 });
        sides.insert(sides.begin() + static_cast<std::ptrdiff_t>(widest) + 1, { ! sides[widest], sides[widest] });
    }
    else if (target_count + 2 <= taus.size() && taus.size() >= 4)
    {
        size_t tightest = 1;
        for (size_t i = 1; i + 2 < taus.size(); i++) // A pair (i, i+1), both interior.
        {
            if (taus[i + 1] - taus[i] < taus[tightest + 1] - taus[tightest])
            {
                tightest = i;
            }
        }
        taus.erase(taus.begin() + static_cast<std::ptrdiff_t>(tightest), taus.begin() + static_cast<std::ptrdiff_t>(tightest) + 2);
        sides.erase(sides.begin() + static_cast<std::ptrdiff_t>(tightest), sides.begin() + static_cast<std::ptrdiff_t>(tightest) + 2);
    }

    // Equal-angle relaxation: pull every apex a SMALL share towards the uniform layout, with the
    // per-layer displacement capped so consecutive layers always keep overlapping.
    constexpr double correction_gain = 0.15; // 15% of the residual error per layer.
    const double max_shift_tau
        = arcToTau(axis, profile, std::min<coord_t>(axis.total_length, 2 * line_width)); // Displacement cap, expressed in tau near the axis start (approximation).
    const double uniform_spacing = tau_total / static_cast<double>(taus.size());
    std::vector<double> corrected(taus.size());
    for (size_t i = 0; i < taus.size(); i++)
    {
        const double target = (static_cast<double>(i) + 0.5) * uniform_spacing;
        const double shift = std::clamp(correction_gain * (target - taus[i]), -max_shift_tau, max_shift_tau);
        corrected[i] = taus[i] + shift;
    }

    // Order constraint: apexes must stay strictly ordered with a minimum separation, so the
    // zigzag can never fold back over itself.
    const double min_separation = std::min(0.25 * delta_tau, tau_total / static_cast<double>(corrected.size()));
    for (size_t i = 0; i < corrected.size(); i++)
    {
        const double lower = (i == 0) ? 0.0 : corrected[i - 1] + min_separation;
        corrected[i] = std::clamp(corrected[i], lower, tau_total);
    }

    axis.apex_positions.clear();
    axis.apex_sides = std::move(sides);
    for (const double tau : corrected)
    {
        axis.apex_positions.push_back(tauToArc(axis, profile, tau));
    }
}

/*!
 * Build the waves of one layer. For a part that also exists on the layer below, the phase is
 * transferred through the cross-section parameter mapping (spatial projection only as fallback
 * at topology changes), then corrected by a small share towards the equal-angle layout. Newly
 * appearing parts extrapolate the phase of their nearest parent, or get a fresh equal-angle
 * layout when fully isolated.
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
        pruneMedialChains(chains, line_width, cell_period);
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

        const PartRecord* below = findPartBelow(records_below, probe, part, 2 * cell_period);
        double delta_tau = 0.0;
        bool inherited = false;
        if (below != nullptr && ! below->apexes.empty() && below->delta_tau > 0.0)
        {
            delta_tau = below->delta_tau;

            // Normal case: the phase is transferred through the cross-section parameter mapping
            // ("截面间参数映射"), which preserves the apex order and the equal-angle property by
            // construction. At topology changes (region split/merge, strong drift) the mapping
            // is not valid; then this one layer falls back to spatial projection, which handles
            // partial correspondence: children of a split both start from the parent's phase,
            // a merge continues the dominant parent's phase.
            inherited = transferApexesByParameter(axes.front(), *below);
            if (! inherited)
            {
                inheritApexesFromBelow(axes, below->apexes, cell_period, line_width, part);
                inherited = axes.front().apex_positions.size() >= 2;
            }
            if (inherited)
            {
                applyDistributedPhaseCorrection(axes.front(), delta_tau, line_width);
            }
        }
        if (! inherited || axes.front().apex_positions.size() < 2)
        {
            // Fully isolated new region: only here the phase starts fresh, zeroed on the
            // region's own skeleton with the equal-angle layout.
            layoutFreshApexes(axes.front(), cell_period, part);
            delta_tau = tauTotal(axes.front()) / static_cast<double>(std::max<size_t>(1, axes.front().apex_positions.size()));
        }
#ifdef SW_DEBUG
        std::fprintf(
            stderr,
            "  axis len=%lld apexes=%zu inherited=%d\n",
            static_cast<long long>(axes.front().total_length),
            axes.front().apex_positions.size(),
            below != nullptr && ! below->apexes.empty());
#endif

        OpenLinesSet part_lines;
        PartRecord record;
        {
            const ReferenceAxis& axis = axes.front();
            generateWaveAlongAxis(axis, part, line_width, part_lines);
            const std::vector<double> profile = tauProfile(axis);
            for (size_t apex_idx = 0; apex_idx < axis.apex_positions.size(); apex_idx++)
            {
                Point2LL base;
                Point2LL tangent;
                coord_t radius;
                sampleAxis(axis, axis.apex_positions[apex_idx], base, tangent, radius);
                record.apexes.push_back({ base, apexTip(axis, apex_idx, part) });
                record.apex_tau_fractions.push_back(arcToTau(axis, profile, axis.apex_positions[apex_idx]) / std::max(profile.back(), 1e-12));
                record.apex_sides.push_back(axis.apex_sides[apex_idx]);
            }
            record.axis_front = axis.points.front();
            record.axis_back = axis.points.back();
            record.delta_tau = delta_tau;
        }

        // Straight wave flanks may still get cut where the boundary is degenerate (e.g. an apex
        // reaching through the pinched mouth of a rib pocket); clip to the part and cancel the
        // short cut-off stubs ("较短的填充线条可取消"), so that one single continuous polyline
        // remains per region and no separate fragment needs a travel move.
        OpenLinesSet clipped = part.intersection(part_lines, true, line_width);
        if (clipped.size() > 1)
        {
            coord_t longest = 0;
            for (const OpenPolyline& polyline : clipped)
            {
                longest = std::max(longest, polyline.length());
            }
            OpenLinesSet kept;
            for (OpenPolyline& polyline : clipped)
            {
                if (polyline.length() >= cell_period || polyline.length() == longest)
                {
                    kept.push_back(std::move(polyline));
                }
            }
            clipped = std::move(kept);
        }
        lines_out.push_back(std::move(clipped));

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

        OpenLinesSet part_lines;
        generateWaveAlongAxis(axis, best_skipped_part, line_width, part_lines);
        lines_out.push_back(best_skipped_part.intersection(part_lines, true, line_width));

        PartRecord record;
        for (size_t apex_idx = 0; apex_idx < axis.apex_positions.size(); apex_idx++)
        {
            Point2LL base;
            Point2LL tangent;
            coord_t radius;
            sampleAxis(axis, axis.apex_positions[apex_idx], base, tangent, radius);
            record.apexes.push_back({ base, apexTip(axis, apex_idx, best_skipped_part) });
            record.apex_sides.push_back(axis.apex_sides[apex_idx]);
        }
        record.apex_tau_fractions = { 0.25, 0.75 };
        record.axis_front = axis.points.front();
        record.axis_back = axis.points.back();
        record.delta_tau = tauTotal(axis) / 2.0;
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
        // Same region division as the native grid-type patterns: the given infill areas, minus
        // the extra infill walls (like Infill::generateWallToolPaths does at gcode time). The
        // wave peaks then contact the boundary of the region the pattern actually fills.
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
