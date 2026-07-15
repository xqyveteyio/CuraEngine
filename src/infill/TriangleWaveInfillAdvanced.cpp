// Copyright (c) 2026 UltiMaker
// CuraEngine is released under the terms of the AGPLv3 or higher.

#include "infill/TriangleWaveInfillAdvanced.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

#include "BoostInterface.hpp" // needed for the boost voronoi traits of PolygonsSegmentIndex
#include "geometry/OpenPolyline.h"
#include "geometry/Shape.h"
#include "geometry/SingleShape.h"
#include "utils/AABB.h"
#include "utils/PolygonsSegmentIndex.h"
#include "utils/VoronoiUtils.h"

namespace cura
{

namespace
{

// Pull the tips inwards by a hair so that clipping against the outline never cuts the apex
// vertex itself and the tips stay sharp.
constexpr coord_t tip_inset = 20;

// ================================ straight wave (fallback & junction patches) ================================

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

// Straight axis-aligned triangle wave over the given region (used as fallback for regions without
// a usable skeleton, and to fill the junction patches).
OpenLinesSet buildStraightWave(const Shape& region, const coord_t line_distance)
{
    std::map<int64_t, std::pair<coord_t, coord_t>> extremes;
    gatherColumnExtremes(region, line_distance, extremes);
    return buildWave(extremes, line_distance);
}

// ================================ medial axis (skeleton) extraction ================================

// Undirected graph of the (approximate) medial axis of a region. Each node knows its clearance:
// the distance from the node to the region boundary (the local half-width).
struct MedialAxisGraph
{
    struct Node
    {
        Point2LL p;
        coord_t clearance{ 0 };
        std::vector<size_t> adj; // indices of the neighboring nodes
    };
    std::vector<Node> nodes;

    void removeEdge(const size_t a, const size_t b)
    {
        auto& adj_a = nodes[a].adj;
        auto& adj_b = nodes[b].adj;
        adj_a.erase(std::find(adj_a.begin(), adj_a.end(), b));
        adj_b.erase(std::find(adj_b.begin(), adj_b.end(), a));
    }
};

// A maximal skeleton path between two terminal nodes (endpoints or junctions), or a closed loop.
struct SkeletonBranch
{
    std::vector<Point2LL> points;
    size_t start_node{ 0 };
    size_t end_node{ 0 };
    bool closed{ false };
    bool start_at_junction{ false };
    bool end_at_junction{ false };
};

// Approximate the medial axis of the part with the internal edges of the voronoi diagram of its
// boundary segments. Parabolic edges are approximated by their chord; that is accurate enough
// here because the wave ribs re-measure the exact wall positions by ray casting later on.
MedialAxisGraph extractMedialAxis(const SingleShape& part)
{
    MedialAxisGraph graph;

    const std::vector<Point2LL> points; // remains empty; all sites are segments
    std::vector<PolygonsSegmentIndex> segments;
    for (size_t poly_idx = 0; poly_idx < part.size(); poly_idx++)
    {
        for (size_t point_idx = 0; point_idx < part[poly_idx].size(); point_idx++)
        {
            segments.emplace_back(&part, poly_idx, point_idx);
        }
    }
    if (segments.empty())
    {
        return graph;
    }

    VoronoiUtils::vd_t voronoi_diagram;
    boost::polygon::construct_voronoi(segments.begin(), segments.end(), &voronoi_diagram);

    std::map<std::pair<coord_t, coord_t>, size_t> node_ids;
    auto node_for = [&](const Point2LL& p, const VoronoiUtils::vd_t::cell_type& cell) -> size_t
    {
        const auto [it, inserted] = node_ids.try_emplace(std::make_pair(p.X, p.Y), graph.nodes.size());
        if (inserted)
        {
            graph.nodes.push_back({ p, VoronoiUtils::getDistance(p, cell, points, segments), {} });
        }
        return it->second;
    };

    for (const auto& edge : voronoi_diagram.edges())
    {
        if (! edge.is_primary() || ! edge.is_finite() || &edge > edge.twin())
        {
            continue; // handle every twin pair only once
        }
        const Point2LL v0 = VoronoiUtils::p(edge.vertex0());
        const Point2LL v1 = VoronoiUtils::p(edge.vertex1());
        if (v0 == v1 || ! part.inside(v0, false) || ! part.inside(v1, false))
        {
            continue; // keep only the skeleton edges strictly inside the region
        }
        const size_t n0 = node_for(v0, *edge.cell());
        const size_t n1 = node_for(v1, *edge.cell());
        if (n0 == n1 || std::find(graph.nodes[n0].adj.begin(), graph.nodes[n0].adj.end(), n1) != graph.nodes[n0].adj.end())
        {
            continue;
        }
        graph.nodes[n0].adj.push_back(n1);
        graph.nodes[n1].adj.push_back(n0);
    }

    return graph;
}

// Remove the spurs which the medial axis grows into every convex boundary corner: leaf edges that
// dive steeply towards the boundary, and leaf ends of regions too narrow to hold a wave anyway.
void pruneSpurs(MedialAxisGraph& graph, const coord_t line_distance)
{
    const coord_t min_useful_clearance = line_distance / 2;
    bool changed = true;
    while (changed)
    {
        changed = false;
        for (size_t leaf = 0; leaf < graph.nodes.size(); ++leaf)
        {
            if (graph.nodes[leaf].adj.size() != 1)
            {
                continue;
            }
            const size_t neighbor = graph.nodes[leaf].adj.front();
            const coord_t rise = graph.nodes[neighbor].clearance - graph.nodes[leaf].clearance;
            const coord_t length = vSize(graph.nodes[neighbor].p - graph.nodes[leaf].p);
            // Corner spurs ascend with a slope of at least cos(corner_angle / 2) (~0.7 for right
            // angles); genuinely tapering limbs ascend much slower and are kept.
            const bool is_corner_spur = rise > (length * 3) / 5;
            const bool too_narrow = graph.nodes[leaf].clearance < min_useful_clearance && rise > 0;
            if (is_corner_spur || too_narrow)
            {
                graph.removeEdge(leaf, neighbor);
                changed = true;
            }
        }
    }
}

// Split the pruned skeleton into maximal branches: open paths running between terminals (nodes of
// degree 1 or >= 3, the latter being the junctions), plus closed loops of pure degree-2 nodes
// (e.g. the skeleton of an annulus).
std::vector<SkeletonBranch> extractBranches(const MedialAxisGraph& graph, std::vector<size_t>& junctions)
{
    std::vector<SkeletonBranch> branches;
    std::set<std::pair<size_t, size_t>> visited;
    auto edge_key = [](const size_t a, const size_t b)
    {
        return std::make_pair(std::min(a, b), std::max(a, b));
    };

    auto walk = [&](const size_t start, const size_t first) -> SkeletonBranch
    {
        SkeletonBranch branch;
        branch.start_node = start;
        branch.points.push_back(graph.nodes[start].p);
        size_t prev = start;
        size_t current = first;
        visited.insert(edge_key(prev, current));
        while (true)
        {
            branch.points.push_back(graph.nodes[current].p);
            if (graph.nodes[current].adj.size() != 2 || current == start)
            {
                break;
            }
            const auto& adj = graph.nodes[current].adj;
            const size_t next = (adj[0] == prev) ? adj[1] : adj[0];
            visited.insert(edge_key(current, next));
            prev = current;
            current = next;
        }
        branch.end_node = current;
        branch.closed = (current == start);
        branch.start_at_junction = graph.nodes[start].adj.size() >= 3;
        branch.end_at_junction = graph.nodes[current].adj.size() >= 3;
        return branch;
    };

    for (size_t i = 0; i < graph.nodes.size(); ++i)
    {
        const size_t degree = graph.nodes[i].adj.size();
        if (degree >= 3)
        {
            junctions.push_back(i);
        }
        if (degree == 2) // interior of a branch, or part of a loop; handled below
        {
            continue;
        }
        for (const size_t neighbor : graph.nodes[i].adj)
        {
            if (! visited.contains(edge_key(i, neighbor)))
            {
                branches.push_back(walk(i, neighbor));
            }
        }
    }

    // Whatever is left unvisited consists of closed loops of degree-2 nodes.
    for (size_t i = 0; i < graph.nodes.size(); ++i)
    {
        if (graph.nodes[i].adj.size() != 2)
        {
            continue;
        }
        for (const size_t neighbor : graph.nodes[i].adj)
        {
            if (! visited.contains(edge_key(i, neighbor)))
            {
                SkeletonBranch loop = walk(i, neighbor);
                loop.points.pop_back(); // drop the repeated start point
                branches.push_back(loop);
            }
        }
    }

    // Canonical direction, so the left/right tooth parity of the wave is deterministic.
    for (SkeletonBranch& branch : branches)
    {
        if (branch.closed)
        {
            continue;
        }
        const Point2LL& a = branch.points.front();
        const Point2LL& b = branch.points.back();
        if (std::make_pair(b.X, b.Y) < std::make_pair(a.X, a.Y))
        {
            std::reverse(branch.points.begin(), branch.points.end());
            std::swap(branch.start_node, branch.end_node);
            std::swap(branch.start_at_junction, branch.end_at_junction);
        }
    }

    return branches;
}

// ================================ skeleton driven wave ================================

struct Vec2d
{
    double x;
    double y;
};

Point2LL roundedPoint(const double x, const double y)
{
    return { std::llrint(x), std::llrint(y) };
}

// Arc-length parametrization of a branch polyline, with linear interpolation and central
// difference tangents.
class PathSampler
{
public:
    PathSampler(const std::vector<Point2LL>& points, const bool closed)
        : points_(points)
        , closed_(closed)
    {
        cumulative_.push_back(0.0);
        for (size_t i = 1; i < points_.size(); ++i)
        {
            cumulative_.push_back(cumulative_.back() + vSizeMM(points_[i] - points_[i - 1]) * 1000.0);
        }
        if (closed_ && points_.size() >= 2)
        {
            cumulative_.push_back(cumulative_.back() + vSizeMM(points_.front() - points_.back()) * 1000.0);
        }
    }

    double length() const
    {
        return cumulative_.back();
    }

    Point2LL at(double s) const
    {
        if (closed_)
        {
            s = std::fmod(std::fmod(s, length()) + length(), length());
        }
        else
        {
            s = std::clamp(s, 0.0, length());
        }
        const auto it = std::upper_bound(cumulative_.begin(), cumulative_.end(), s);
        const size_t segment = std::min<size_t>(std::distance(cumulative_.begin(), it), cumulative_.size() - 1) - 1;
        const double segment_length = cumulative_[segment + 1] - cumulative_[segment];
        const double t = (segment_length > 0.0) ? (s - cumulative_[segment]) / segment_length : 0.0;
        const Point2LL& a = points_[segment];
        const Point2LL& b = points_[(segment + 1) % points_.size()];
        return roundedPoint(a.X + (b.X - a.X) * t, a.Y + (b.Y - a.Y) * t);
    }

    Vec2d tangentAt(const double s) const
    {
        constexpr double h = 100.0; // central difference step [μm]
        const Point2LL before = at(closed_ ? s - h : std::max(s - h, 0.0));
        const Point2LL after = at(closed_ ? s + h : std::min(s + h, length()));
        const double dx = static_cast<double>(after.X - before.X);
        const double dy = static_cast<double>(after.Y - before.Y);
        const double len = std::hypot(dx, dy);
        if (len <= 0.0)
        {
            return { 1.0, 0.0 };
        }
        return { dx / len, dy / len };
    }

private:
    std::vector<Point2LL> points_;
    bool closed_;
    std::vector<double> cumulative_;
};

// Distance from p to the region wall along the given direction, found by ray casting. Returns
// ray_length when no wall is hit.
double wallDistance(const Shape& region, const Point2LL& p, const Vec2d& direction, const double ray_length)
{
    const Point2LL far = roundedPoint(p.X + direction.x * ray_length, p.Y + direction.y * ray_length);
    const std::vector<float> crossings = region.intersectionsWithSegment(p, far);
    double nearest = 1.0;
    bool found = false;
    for (const float t : crossings)
    {
        if (t > 0.0f && t < nearest)
        {
            nearest = t;
            found = true;
        }
    }
    return found ? nearest * ray_length : ray_length;
}

// Build the triangle wave for one limb: sample rib positions at an even spacing (close to the
// requested line distance) along the skeleton path, and put the apexes alternately on the left
// and right wall, perpendicular to the local skeleton direction.
OpenLinesSet
    buildRibWave(const PathSampler& sampler, const double s_begin, const double s_end, const bool closed, const Shape& region, const coord_t line_distance, const double ray_length)
{
    OpenLinesSet wave;
    const double span = closed ? sampler.length() : s_end - s_begin;
    if (span < static_cast<double>(line_distance) / 2.0)
    {
        return wave;
    }

    size_t rib_count;
    double spacing;
    if (closed)
    {
        // A closed loop needs an even tooth count for the alternating wave to close onto itself.
        rib_count = 2 * std::max<size_t>(1, static_cast<size_t>(std::llround(span / (2.0 * line_distance))));
        spacing = span / static_cast<double>(rib_count);
    }
    else
    {
        rib_count = std::max<size_t>(2, static_cast<size_t>(std::llround(span / line_distance)));
        spacing = span / static_cast<double>(rib_count);
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

    for (size_t i = 0; i < rib_count; ++i)
    {
        const double s = closed ? s_begin + spacing * static_cast<double>(i) : s_begin + spacing * (static_cast<double>(i) + 0.5);
        const Point2LL p = sampler.at(s);
        const Vec2d tangent = sampler.tangentAt(s);
        const double side = (i % 2 == 0) ? 1.0 : -1.0;
        const Vec2d normal{ -tangent.y * side, tangent.x * side };

        if (! region.inside(p, true))
        {
            flush_wave(); // rib foot outside the limb (sharp curvature artifact): interrupt the wave
            continue;
        }
        const double wall = wallDistance(region, p, normal, ray_length);
        const double apex_distance = std::max(0.0, wall - static_cast<double>(tip_inset));
        wave_points.push_back(roundedPoint(p.X + normal.x * apex_distance, p.Y + normal.y * apex_distance));
    }

    if (closed && wave_points.size() >= 3)
    {
        wave_points.push_back(wave_points.front()); // close the loop
    }
    flush_wave();

    return wave;
}

// Thin quad which severs the region perpendicular to the skeleton at the given point, so that
// splitIntoParts() separates the limbs from the junction patch. The quad is guaranteed to poke
// through both walls.
Polygon cutBand(const Shape& region, const Point2LL& p, const Vec2d& tangent, const coord_t line_distance, const double ray_length)
{
    constexpr double half_thickness = 25.0;
    const Vec2d normal{ -tangent.y, tangent.x };
    const double d_plus = wallDistance(region, p, normal, ray_length) + static_cast<double>(line_distance);
    const double d_minus = wallDistance(region, p, { -normal.x, -normal.y }, ray_length) + static_cast<double>(line_distance);

    const double tx = tangent.x * half_thickness;
    const double ty = tangent.y * half_thickness;
    Polygon band;
    band.push_back(roundedPoint(p.X + normal.x * d_plus - tx, p.Y + normal.y * d_plus - ty));
    band.push_back(roundedPoint(p.X + normal.x * d_plus + tx, p.Y + normal.y * d_plus + ty));
    band.push_back(roundedPoint(p.X - normal.x * d_minus + tx, p.Y - normal.y * d_minus + ty));
    band.push_back(roundedPoint(p.X - normal.x * d_minus - tx, p.Y - normal.y * d_minus - ty));
    return band;
}

// Skeleton driven triangle wave for one connected part of the template region:
// - extract and prune the medial axis,
// - sever the region at every skeleton junction ("block off" the branch point), which splits it
//   into limb regions and junction patches,
// - run one wave along the skeleton of each limb,
// - fill the junction patches with a separate small straight wave.
// Returns an empty set when the skeleton degenerates (caller falls back to the straight wave).
OpenLinesSet buildSkeletonWaves(const SingleShape& part, const coord_t line_distance)
{
    OpenLinesSet waves;

    MedialAxisGraph graph = extractMedialAxis(part);
    pruneSpurs(graph, line_distance);
    std::vector<size_t> junctions;
    std::vector<SkeletonBranch> branches = extractBranches(graph, junctions);
    if (branches.empty())
    {
        return waves;
    }

    const AABB part_box(part);
    const double ray_length = vSizeMM(part_box.max_ - part_box.min_) * 1000.0 + static_cast<double>(line_distance);

    // Decide the cut positions: on every branch, at the junction clearance distance away from the
    // junction node. Branches too short to be severed are absorbed into the junction patch.
    struct BranchPlan
    {
        PathSampler sampler;
        double s_begin;
        double s_end;
        bool absorbed{ false };
    };
    std::vector<BranchPlan> plans;
    Shape cut_bands;
    const double min_limb_length = static_cast<double>(line_distance) / 2.0;

    for (const SkeletonBranch& branch : branches)
    {
        PathSampler sampler(branch.points, branch.closed);
        BranchPlan plan{ std::move(sampler), 0.0, 0.0, false };
        plan.s_end = plan.sampler.length();

        if (! branch.closed)
        {
            const double cut_start = branch.start_at_junction ? static_cast<double>(graph.nodes[branch.start_node].clearance) : 0.0;
            const double cut_end = branch.end_at_junction ? static_cast<double>(graph.nodes[branch.end_node].clearance) : 0.0;
            if (plan.sampler.length() < cut_start + cut_end + min_limb_length)
            {
                // Too short to sever: the whole branch becomes part of the junction patch.
                plan.absorbed = branch.start_at_junction || branch.end_at_junction;
            }
            else
            {
                if (branch.start_at_junction)
                {
                    plan.s_begin = cut_start;
                    cut_bands.push_back(cutBand(part, plan.sampler.at(cut_start), plan.sampler.tangentAt(cut_start), line_distance, ray_length));
                }
                if (branch.end_at_junction)
                {
                    plan.s_end = plan.sampler.length() - cut_end;
                    cut_bands.push_back(cutBand(part, plan.sampler.at(plan.s_end), plan.sampler.tangentAt(plan.s_end), line_distance, ray_length));
                }
            }
        }
        plans.push_back(std::move(plan));
    }

    // Sever the region and sort the pieces into limbs and junction patches.
    const Shape remaining = cut_bands.empty() ? Shape(part) : part.difference(cut_bands.unionPolygons());
    std::vector<SingleShape> sub_parts = remaining.splitIntoParts();
    std::vector<bool> is_patch(sub_parts.size(), false);
    for (size_t i = 0; i < sub_parts.size(); ++i)
    {
        for (const size_t junction : junctions)
        {
            if (sub_parts[i].inside(graph.nodes[junction].p, true))
            {
                is_patch[i] = true;
                break;
            }
        }
    }

    auto find_part = [&sub_parts](const Point2LL& p) -> int
    {
        for (size_t i = 0; i < sub_parts.size(); ++i)
        {
            if (sub_parts[i].inside(p, true))
            {
                return static_cast<int>(i);
            }
        }
        return -1;
    };

    // One wave per limb, running along its skeleton branch.
    for (size_t branch_idx = 0; branch_idx < branches.size(); ++branch_idx)
    {
        const SkeletonBranch& branch = branches[branch_idx];
        BranchPlan& plan = plans[branch_idx];
        if (plan.absorbed)
        {
            continue; // covered by the junction patch fill
        }

        // Extend the free (non-junction) ends of the branch to the walls, so the wave also covers
        // the end caps of the limb; the medial axis itself always stops half a width short.
        if (! branch.closed)
        {
            std::vector<Point2LL> extended = branch.points;
            bool changed_points = false;
            if (! branch.start_at_junction)
            {
                const Vec2d tangent = plan.sampler.tangentAt(0.0);
                const double d = wallDistance(part, extended.front(), { -tangent.x, -tangent.y }, ray_length);
                const double extension = std::max(0.0, d - tip_inset);
                if (extension > 0.0)
                {
                    extended.insert(extended.begin(), roundedPoint(extended.front().X - tangent.x * extension, extended.front().Y - tangent.y * extension));
                    plan.s_end += extension; // all arc positions shift by the length prepended at the front
                    changed_points = true;
                }
            }
            if (! branch.end_at_junction)
            {
                const Vec2d tangent = plan.sampler.tangentAt(plan.sampler.length());
                const double d = wallDistance(part, extended.back(), tangent, ray_length);
                const double extension = std::max(0.0, d - tip_inset);
                if (extension > 0.0)
                {
                    extended.push_back(roundedPoint(extended.back().X + tangent.x * extension, extended.back().Y + tangent.y * extension));
                    plan.s_end += extension;
                    changed_points = true;
                }
            }
            if (changed_points)
            {
                plan.sampler = PathSampler(extended, false);
                plan.s_end = std::min(plan.s_end, plan.sampler.length());
            }
        }

        const double s_mid = branch.closed ? 0.0 : (plan.s_begin + plan.s_end) / 2.0;
        const int part_idx = find_part(plan.sampler.at(s_mid));
        if (part_idx < 0 || is_patch[part_idx])
        {
            continue; // limb collapsed into a patch after cutting; the patch fill covers it
        }

        waves.push_back(buildRibWave(plan.sampler, plan.s_begin, plan.s_end, branch.closed, sub_parts[part_idx], line_distance, ray_length));
    }

    // Finally fill the blocked-off junction patches with their own separate wave.
    for (size_t i = 0; i < sub_parts.size(); ++i)
    {
        if (is_patch[i])
        {
            waves.push_back(buildStraightWave(sub_parts[i], line_distance));
        }
    }

    return waves;
}

// ================================ assembly ================================

// Build one independent wave per connected part of the given region, so that disjoint areas
// (multiple models or separate islands of one model) each get their own complete triangle wave
// instead of one wave spanning across the gaps between them. Each part uses the skeleton driven
// wave when possible, and the straight axis-aligned wave otherwise.
OpenLinesSet buildWaves(const Shape& region, const coord_t line_distance)
{
    OpenLinesSet waves;
    for (const SingleShape& part : region.splitIntoParts())
    {
        OpenLinesSet part_waves = buildSkeletonWaves(part, line_distance);
        if (part_waves.empty())
        {
            part_waves = buildStraightWave(part, line_distance);
        }
        waves.push_back(part_waves);
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

TriangleWaveAdvancedFillProvider::TriangleWaveAdvancedFillProvider(const std::vector<Shape>& layer_outlines, coord_t line_distance, const double fill_angle)
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

void TriangleWaveAdvancedFillProvider::generate(OpenLinesSet& result_lines, const Shape& in_outline) const
{
    if (template_wave_.empty() || in_outline.empty())
    {
        return;
    }

    Shape outline = in_outline;
    outline.applyMatrix(rotation_matrix_);

    result_lines = clipWave(template_wave_, outline, rotation_matrix_);
}

void TriangleWaveInfillAdvanced::generateTotalTriangleWaveInfill(OpenLinesSet& result_lines, coord_t line_distance, const Shape& in_outline, const double fill_angle)
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
