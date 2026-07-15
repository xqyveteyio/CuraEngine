// Copyright (c) 2026 UltiMaker
// CuraEngine is released under the terms of the AGPLv3 or higher.

#include "infill/TriangleWaveInfillEpic.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

#ifdef TW_EPIC_BENCH
#include <chrono>
#include <cstdio>
#endif

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

// (Definitions live in the tracking section further down; declared here because the straight
// wave builder can also record its teeth as chains for the tracking mode.)
struct ToothChain;
OpenLinesSet buildWave(const std::map<int64_t, std::pair<coord_t, coord_t>>& extremes, const coord_t line_distance, std::vector<ToothChain>* chains_out);

// Straight axis-aligned triangle wave over the given region (used as fallback for regions without
// a usable skeleton, and to fill the junction patches).
OpenLinesSet buildStraightWave(const Shape& region, const coord_t line_distance, std::vector<ToothChain>* chains_out = nullptr)
{
    std::map<int64_t, std::pair<coord_t, coord_t>> extremes;
    gatherColumnExtremes(region, line_distance, extremes);
    return buildWave(extremes, line_distance, chains_out);
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

// One tooth of a tracked wave. The tooth is defined by its rib: the line through 'foot' along
// 'dir'. The apex (the wave vertex) lies on that line towards dir * side, the base is the chord
// end on the opposite wall. In tracking mode each layer derives its teeth 1:1 from the layer
// below: the rib line stays where it is, only the chord against the current outline is re-cut.
struct TrackedTooth
{
    Point2LL foot; // center of the rib chord (in the rotated frame)
    Vec2d dir; // unit rib direction (unsigned)
    double side; // the apex lies along dir * side
    Point2LL apex; // the wave vertex (wall minus tip inset)
    Point2LL base; // the opposite chord end
};

// An ordered run of alternating teeth; the wave polyline connects the apexes in order.
struct ToothChain
{
    std::vector<TrackedTooth> teeth;
    bool closed{ false };
};

// Draw the complete triangle wave through the column extremes: troughs at even columns, peaks at
// odd columns, with sharp apexes touching the (template) boundary. Columns without material
// interrupt the wave. When \p chains_out is given, the teeth are recorded for the tracking mode.
OpenLinesSet buildWave(const std::map<int64_t, std::pair<coord_t, coord_t>>& extremes, const coord_t line_distance, std::vector<ToothChain>* chains_out)
{
    OpenLinesSet wave;
    if (extremes.empty())
    {
        return wave;
    }

    std::vector<Point2LL> wave_points;
    ToothChain chain;
    auto flush_wave = [&]()
    {
        if (wave_points.size() >= 2)
        {
            wave.push_back(OpenPolyline{ wave_points });
        }
        wave_points.clear();
        if (chains_out != nullptr && chain.teeth.size() >= 2)
        {
            chains_out->push_back(chain);
        }
        chain = ToothChain{};
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
        const coord_t x = columnX(k, line_distance);
        const coord_t apex_y = is_peak ? it->second.second - tip_inset : it->second.first + tip_inset;
        wave_points.emplace_back(x, apex_y);

        if (chains_out != nullptr)
        {
            const coord_t base_y = is_peak ? it->second.first + tip_inset : it->second.second - tip_inset;
            const Point2LL apex(x, apex_y);
            const Point2LL base(x, base_y);
            const Point2LL foot(x, (apex_y + base_y) / 2);
            chain.teeth.push_back({ foot, Vec2d{ 0.0, 1.0 }, is_peak ? 1.0 : -1.0, apex, base });
        }
    }
    flush_wave();

    return wave;
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

// Build the triangle wave for one limb, putting the apexes alternately on the left and right
// wall, perpendicular to the local skeleton direction. The rib positions are distributed evenly
// (with a spacing close to the requested line distance) over the skeleton path. When \p
// chains_out is given, the teeth are also recorded as chains for the tracking mode: the next
// layer then derives its teeth from these instead of regenerating.
OpenLinesSet buildRibWave(
    const PathSampler& sampler,
    const double s_begin,
    const double s_end,
    const bool closed,
    const Shape& region,
    const coord_t line_distance,
    const double ray_length,
    std::vector<ToothChain>* chains_out)
{
    OpenLinesSet wave;
    const double span = closed ? sampler.length() : s_end - s_begin;
    if (span < static_cast<double>(line_distance) / 2.0)
    {
        return wave;
    }

    // Collect the rib positions and their sides (+1: "left" of the path direction).
    std::vector<std::pair<double, double>> ribs; // (arc position, side)
    if (closed)
    {
        // A closed loop needs an even tooth count for the alternating wave to close onto itself.
        const size_t rib_count = 2 * std::max<size_t>(1, static_cast<size_t>(std::llround(span / (2.0 * line_distance))));
        const double spacing = span / static_cast<double>(rib_count);
        for (size_t i = 0; i < rib_count; ++i)
        {
            ribs.emplace_back(spacing * static_cast<double>(i), (i % 2 == 0) ? 1.0 : -1.0);
        }
    }
    else
    {
        const size_t rib_count = std::max<size_t>(2, static_cast<size_t>(std::llround(span / line_distance)));
        const double spacing = span / static_cast<double>(rib_count);
        for (size_t i = 0; i < rib_count; ++i)
        {
            ribs.emplace_back(s_begin + spacing * (static_cast<double>(i) + 0.5), (i % 2 == 0) ? 1.0 : -1.0);
        }
    }

    std::vector<Point2LL> wave_points;
    ToothChain chain;
    bool interrupted = false;
    auto flush_wave = [&]()
    {
        if (wave_points.size() >= 2)
        {
            wave.push_back(OpenPolyline{ wave_points });
        }
        wave_points.clear();
        if (chains_out != nullptr && chain.teeth.size() >= 2)
        {
            chains_out->push_back(chain);
        }
        chain = ToothChain{};
    };

    for (const auto& [s, side] : ribs)
    {
        const Point2LL p = sampler.at(s);
        const Vec2d tangent = sampler.tangentAt(s);
        const Vec2d normal{ -tangent.y, tangent.x };
        const Vec2d apex_dir{ normal.x * side, normal.y * side };

        if (! region.inside(p, true))
        {
            flush_wave(); // rib foot outside the limb (sharp curvature artifact): interrupt the wave
            interrupted = true;
            continue;
        }
        const double apex_wall = wallDistance(region, p, apex_dir, ray_length);
        const double apex_distance = std::max(0.0, apex_wall - static_cast<double>(tip_inset));
        const Point2LL apex = roundedPoint(p.X + apex_dir.x * apex_distance, p.Y + apex_dir.y * apex_distance);
        wave_points.push_back(apex);

        if (chains_out != nullptr)
        {
            const double base_wall = wallDistance(region, p, { -apex_dir.x, -apex_dir.y }, ray_length);
            const double base_distance = std::max(0.0, base_wall - static_cast<double>(tip_inset));
            const Point2LL base = roundedPoint(p.X - apex_dir.x * base_distance, p.Y - apex_dir.y * base_distance);
            const Point2LL foot = roundedPoint((apex.X + base.X) / 2.0, (apex.Y + base.Y) / 2.0);
            chain.teeth.push_back({ foot, normal, side, apex, base });
        }
    }

    if (closed && wave_points.size() >= 3)
    {
        wave_points.push_back(wave_points.front()); // close the loop
        if (! interrupted)
        {
            chain.closed = true;
        }
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
OpenLinesSet buildSkeletonWaves(const SingleShape& part, const coord_t line_distance, std::vector<ToothChain>* chains_out = nullptr)
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

        waves.push_back(buildRibWave(plan.sampler, plan.s_begin, plan.s_end, branch.closed, sub_parts[part_idx], line_distance, ray_length, chains_out));
    }

    // Finally fill the blocked-off junction patches with their own separate wave.
    for (size_t i = 0; i < sub_parts.size(); ++i)
    {
        if (is_patch[i])
        {
            waves.push_back(buildStraightWave(sub_parts[i], line_distance, chains_out));
        }
    }

    return waves;
}

// ================================ tracking: per-tooth derivation ================================

// The maximal interval of the rib line inside the region, as signed distances from the foot
// (positive towards dir). When the line crosses the region several times (holes, other
// corridors), the interval closest to the foot is chosen; 'valid' is false when there is none
// within max_foot_distance.
struct RibChord
{
    double s_low{ 0.0 };
    double s_high{ 0.0 };
    bool valid{ false };
};

RibChord ribChordNear(const Shape& region, const Point2LL& foot, const Vec2d& dir, const double ray_length, const double max_foot_distance)
{
    const Point2LL p1 = roundedPoint(foot.X - dir.x * ray_length, foot.Y - dir.y * ray_length);
    const Point2LL p2 = roundedPoint(foot.X + dir.x * ray_length, foot.Y + dir.y * ray_length);
    std::vector<float> ts = region.intersectionsWithSegment(p1, p2);
    if (ts.size() < 2 || ts.size() % 2 != 0)
    {
        return {}; // degenerate (tangential) crossing; report failure so the tooth freezes
    }
    std::sort(ts.begin(), ts.end());

    RibChord best;
    double best_distance = max_foot_distance;
    for (size_t i = 0; i + 1 < ts.size(); i += 2) // p1 lies far outside, so the even intervals are the inside ones
    {
        const double s_low = (static_cast<double>(ts[i]) - 0.5) * 2.0 * ray_length;
        const double s_high = (static_cast<double>(ts[i + 1]) - 0.5) * 2.0 * ray_length;
        const double distance = (s_low > 0.0) ? s_low : ((s_high < 0.0) ? -s_high : 0.0);
        if (distance < best_distance)
        {
            best_distance = distance;
            best = { s_low, s_high, true };
        }
    }
    return best;
}

// Derive one tooth from its previous-layer state: the rib line stays where it is, only the chord
// against the current outline is re-cut and both ends follow their wall. An end chases a receding
// wall without limit only while the movement per layer stays small; when the wall suddenly
// vanishes (e.g. the rib now looks through a junction opening into another corridor), the end
// freezes in place instead of jumping there. Walls moving inwards are always followed exactly.
// Returns false when the tooth has no material left to live in.
bool updateTooth(TrackedTooth& tooth, const Shape& region, const coord_t line_distance, const double ray_length)
{
    const RibChord chord = ribChordNear(region, tooth.foot, tooth.dir, ray_length, static_cast<double>(line_distance));
    if (! chord.valid)
    {
        return region.inside(tooth.foot, true); // freeze entirely, or die when the material is gone
    }
    if (chord.s_high - chord.s_low < 4.0 * tip_inset)
    {
        return false; // sliver, too narrow to hold a tooth
    }

    const double jump_threshold = static_cast<double>(line_distance) / 2.0;
    const Vec2d apex_dir{ tooth.dir.x * tooth.side, tooth.dir.y * tooth.side };

    // Work in signed distances from the foot: 'a' towards the apex wall, 'b' towards the base wall.
    const double a_wall = (tooth.side > 0.0) ? chord.s_high : -chord.s_low;
    const double b_wall = (tooth.side > 0.0) ? -chord.s_low : chord.s_high;

    auto derive_end = [jump_threshold](const double candidate, const double previous, const double other_wall_limit)
    {
        double result = candidate;
        if (candidate - previous > jump_threshold)
        {
            result = previous; // wall vanished: freeze instead of chasing it
        }
        return std::clamp(result, -other_wall_limit, candidate);
    };

    const double a_prev = (tooth.apex.X - tooth.foot.X) * apex_dir.x + (tooth.apex.Y - tooth.foot.Y) * apex_dir.y;
    const double b_prev = -((tooth.base.X - tooth.foot.X) * apex_dir.x + (tooth.base.Y - tooth.foot.Y) * apex_dir.y);
    const double a_new = derive_end(a_wall - tip_inset, a_prev, b_wall - tip_inset);
    const double b_new = derive_end(b_wall - tip_inset, b_prev, a_wall - tip_inset);

    tooth.apex = roundedPoint(tooth.foot.X + apex_dir.x * a_new, tooth.foot.Y + apex_dir.y * a_new);
    tooth.base = roundedPoint(tooth.foot.X - apex_dir.x * b_new, tooth.foot.Y - apex_dir.y * b_new);
    tooth.foot = roundedPoint((tooth.apex.X + tooth.base.X) / 2.0, (tooth.apex.Y + tooth.base.Y) / 2.0);
    return true;
}

// Create a brand new tooth (gap fill or chain extension) at the given foot/dir/side.
bool makeTooth(TrackedTooth& tooth, const Shape& region, const coord_t line_distance, const double ray_length)
{
    const Point2LL requested_foot = tooth.foot;
    const RibChord chord = ribChordNear(region, tooth.foot, tooth.dir, ray_length, 0.75 * static_cast<double>(line_distance));
    if (! chord.valid || chord.s_high - chord.s_low < 4.0 * tip_inset)
    {
        return false;
    }
    const Vec2d apex_dir{ tooth.dir.x * tooth.side, tooth.dir.y * tooth.side };
    const double a = ((tooth.side > 0.0) ? chord.s_high : -chord.s_low) - tip_inset;
    const double b = ((tooth.side > 0.0) ? -chord.s_low : chord.s_high) - tip_inset;
    tooth.apex = roundedPoint(tooth.foot.X + apex_dir.x * a, tooth.foot.Y + apex_dir.y * a);
    tooth.base = roundedPoint(tooth.foot.X - apex_dir.x * b, tooth.foot.Y - apex_dir.y * b);
    tooth.foot = roundedPoint((tooth.apex.X + tooth.base.X) / 2.0, (tooth.apex.Y + tooth.base.Y) / 2.0);

    // Reject the tooth when re-centering onto the chord has thrown its foot far away from the
    // requested position (the rib looked through an opening into a wide area, e.g. across the
    // center of a star-shaped region). Such runaway teeth do not close the gap they were created
    // for, so the same gap would spawn ever more teeth on every following layer, which used to
    // let the tooth count (and the slicing time per layer) grow exponentially.
    const double displacement = std::hypot(static_cast<double>(tooth.foot.X - requested_foot.X), static_cast<double>(tooth.foot.Y - requested_foot.Y));
    if (displacement > 0.75 * static_cast<double>(line_distance))
    {
        return false;
    }
    return true;
}

// Drop teeth whose feet have drifted closer than half a pitch to their predecessor. Without
// this, teeth accumulate without bound: the per-layer re-centering of the feet compresses the
// spacing on one side of a gap, the gap re-opens on the other side and is re-filled with new
// teeth on every layer, so the tooth count (and with it the time per layer) keeps growing.
void pruneCrowdedTeeth(ToothChain& chain, const coord_t line_distance)
{
    const double min_spacing = 0.5 * static_cast<double>(line_distance);
    std::vector<TrackedTooth> kept;
    for (const TrackedTooth& tooth : chain.teeth)
    {
        if (! kept.empty())
        {
            const double distance = std::hypot(static_cast<double>(tooth.foot.X - kept.back().foot.X), static_cast<double>(tooth.foot.Y - kept.back().foot.Y));
            if (distance < min_spacing)
            {
                continue;
            }
        }
        kept.push_back(tooth);
    }
    if (chain.closed && kept.size() >= 2)
    {
        const double wrap_distance = std::hypot(static_cast<double>(kept.back().foot.X - kept.front().foot.X), static_cast<double>(kept.back().foot.Y - kept.front().foot.Y));
        if (wrap_distance < min_spacing)
        {
            kept.pop_back();
        }
    }
    chain.teeth = std::move(kept);
}

// Drop teeth which break the left/right alternation (this can happen when a tooth in between
// died); closed chains additionally need an even count for the wave to close onto itself.
void fixParity(ToothChain& chain)
{
    std::vector<TrackedTooth> kept;
    for (const TrackedTooth& tooth : chain.teeth)
    {
        if (kept.empty() || kept.back().side != tooth.side)
        {
            kept.push_back(tooth);
        }
    }
    if (chain.closed)
    {
        while (kept.size() >= 2 && (kept.size() % 2 != 0 || kept.front().side == kept.back().side))
        {
            kept.pop_back();
        }
    }
    chain.teeth = std::move(kept);
}

// Slowly re-align the rib directions with the run of the chain (perpendicular to the line through
// the neighboring feet), limited to a small rotation per layer so that a slowly twisting model is
// followed without ever changing the pattern abruptly.
void smoothDirections(ToothChain& chain)
{
    constexpr double max_rotation = 0.06; // [rad] per layer
    const size_t n = chain.teeth.size();
    if (n < 3)
    {
        return;
    }
    const std::vector<TrackedTooth> before = chain.teeth; // read the tangents from the unmodified state
    for (size_t i = 0; i < n; ++i)
    {
        const TrackedTooth& prev = before[(i == 0) ? (chain.closed ? n - 1 : 0) : i - 1];
        const TrackedTooth& next = before[(i + 1 >= n) ? (chain.closed ? 0 : n - 1) : i + 1];
        const double tx = static_cast<double>(next.foot.X - prev.foot.X);
        const double ty = static_cast<double>(next.foot.Y - prev.foot.Y);
        const double len = std::hypot(tx, ty);
        if (len <= 0.0)
        {
            continue;
        }
        Vec2d target{ -ty / len, tx / len }; // perpendicular to the chain
        const Vec2d dir = chain.teeth[i].dir;
        if (target.x * dir.x + target.y * dir.y < 0.0)
        {
            target = { -target.x, -target.y }; // keep the sign continuous
        }
        const double angle = std::atan2(dir.x * target.y - dir.y * target.x, dir.x * target.x + dir.y * target.y);
        const double clamped = std::clamp(angle, -max_rotation, max_rotation);
        const double c = std::cos(clamped);
        const double s = std::sin(clamped);
        chain.teeth[i].dir = { dir.x * c - dir.y * s, dir.x * s + dir.y * c };
    }
}

// Insert new teeth where the feet of neighboring teeth have moved more than ~1.5 teeth apart
// (the region grew in the middle of a chain). The number of inserted teeth respects the
// alternation parity of the two neighbors.
void fillGaps(ToothChain& chain, const Shape& region, const coord_t line_distance, const double ray_length)
{
    const double pitch = static_cast<double>(line_distance);
    const size_t n = chain.teeth.size();
    if (n < 2)
    {
        return;
    }
    std::vector<TrackedTooth> result;
    const size_t pair_count = chain.closed ? n : n - 1;
    for (size_t i = 0; i < pair_count; ++i)
    {
        const TrackedTooth& current = chain.teeth[i];
        const TrackedTooth& following = chain.teeth[(i + 1) % n];
        result.push_back(current);

        const double dx = static_cast<double>(following.foot.X - current.foot.X);
        const double dy = static_cast<double>(following.foot.Y - current.foot.Y);
        const double distance = std::hypot(dx, dy);
        if (distance < 1.6 * pitch)
        {
            continue;
        }
        int64_t count = std::llround(distance / pitch) - 1;
        const bool need_even = (current.side != following.side); // normal alternation
        if ((count % 2 == 0) != need_even)
        {
            count += 1;
        }
        for (int64_t k = 1; k <= count; ++k)
        {
            const double t = static_cast<double>(k) / static_cast<double>(count + 1);
            Vec2d dir{ current.dir.x + (following.dir.x - current.dir.x) * t, current.dir.y + (following.dir.y - current.dir.y) * t };
            const double dir_len = std::hypot(dir.x, dir.y);
            if (dir_len <= 0.0)
            {
                continue;
            }
            dir = { dir.x / dir_len, dir.y / dir_len };
            TrackedTooth tooth{ roundedPoint(current.foot.X + dx * t, current.foot.Y + dy * t), dir, (k % 2 == 1) ? -current.side : current.side, {}, {} };
            if (makeTooth(tooth, region, line_distance, ray_length))
            {
                result.push_back(tooth);
            }
        }
    }
    if (! chain.closed)
    {
        result.push_back(chain.teeth.back());
    }
    chain.teeth = std::move(result);
}

// Grow an open chain at its two ends, one pitch at a time, into region which is not covered yet
// (the model grew at the end of a corridor).
void extendEnds(ToothChain& chain, const Shape& region, const coord_t line_distance, const double ray_length, std::vector<Point2LL>& all_feet)
{
    if (chain.closed || chain.teeth.size() < 2)
    {
        return;
    }
    const double pitch = static_cast<double>(line_distance);
    const double near_distance = 0.8 * pitch;
    auto near_existing = [&all_feet, near_distance](const Point2LL& p)
    {
        for (const Point2LL& foot : all_feet)
        {
            if (std::hypot(static_cast<double>(p.X - foot.X), static_cast<double>(p.Y - foot.Y)) < near_distance)
            {
                return true;
            }
        }
        return false;
    };

    for (int end = 0; end < 2; ++end)
    {
        for (int guard = 0; guard < 256; ++guard)
        {
            const size_t n = chain.teeth.size();
            const TrackedTooth& last = (end == 0) ? chain.teeth.front() : chain.teeth.back();
            const TrackedTooth& before_last = (end == 0) ? chain.teeth[1] : chain.teeth[n - 2];
            const double sx = static_cast<double>(last.foot.X - before_last.foot.X);
            const double sy = static_cast<double>(last.foot.Y - before_last.foot.Y);
            const double len = std::hypot(sx, sy);
            if (len <= 0.0)
            {
                break;
            }
            const Point2LL candidate = roundedPoint(last.foot.X + sx / len * pitch, last.foot.Y + sy / len * pitch);
            if (! region.inside(candidate, false) || near_existing(candidate))
            {
                break;
            }
            TrackedTooth tooth{ candidate, last.dir, -last.side, {}, {} };
            if (! makeTooth(tooth, region, line_distance, ray_length))
            {
                break;
            }
            all_feet.push_back(tooth.foot);
            if (end == 0)
            {
                chain.teeth.insert(chain.teeth.begin(), tooth);
            }
            else
            {
                chain.teeth.push_back(tooth);
            }
        }
    }
}

// Derive the chains of one layer from the previous layer's chains: update every tooth in place
// (dead teeth split their chain into fragments), then re-align directions, fill gaps and grow the
// open ends into newly appeared region.
std::vector<ToothChain> deriveChains(const std::vector<ToothChain>& previous, const Shape& region, const coord_t line_distance, const double ray_length)
{
    std::vector<ToothChain> derived;

    for (const ToothChain& chain : previous)
    {
        std::vector<std::pair<TrackedTooth, bool>> updated;
        updated.reserve(chain.teeth.size());
        bool all_alive = true;
        for (const TrackedTooth& tooth : chain.teeth)
        {
            TrackedTooth copy = tooth;
            const bool alive = updateTooth(copy, region, line_distance, ray_length);
            all_alive &= alive;
            updated.emplace_back(copy, alive);
        }

        if (chain.closed && all_alive)
        {
            ToothChain intact;
            intact.closed = true;
            for (const auto& [tooth, alive] : updated)
            {
                intact.teeth.push_back(tooth);
            }
            derived.push_back(std::move(intact));
            continue;
        }

        // Split into open fragments at the dead teeth. For a (broken) closed chain, start the
        // walk behind a dead tooth so that the wrap-around run stays in one piece.
        const size_t n = updated.size();
        size_t start = 0;
        if (chain.closed)
        {
            while (start < n && updated[start].second)
            {
                ++start;
            }
            ++start; // first index behind the first dead tooth
        }
        ToothChain fragment;
        for (size_t k = 0; k < n; ++k)
        {
            const auto& [tooth, alive] = updated[(start + k) % n];
            if (alive)
            {
                fragment.teeth.push_back(tooth);
            }
            else if (! fragment.teeth.empty())
            {
                derived.push_back(std::move(fragment));
                fragment = ToothChain{};
            }
        }
        if (! fragment.teeth.empty())
        {
            derived.push_back(std::move(fragment));
        }
    }

    for (ToothChain& chain : derived)
    {
        pruneCrowdedTeeth(chain, line_distance);
        fixParity(chain);
        smoothDirections(chain);
        fillGaps(chain, region, line_distance, ray_length);
    }

    std::vector<Point2LL> all_feet;
    for (const ToothChain& chain : derived)
    {
        for (const TrackedTooth& tooth : chain.teeth)
        {
            all_feet.push_back(tooth.foot);
        }
    }
    for (ToothChain& chain : derived)
    {
        extendEnds(chain, region, line_distance, ray_length, all_feet);
    }

    std::erase_if(
        derived,
        [](const ToothChain& chain)
        {
            return chain.teeth.size() < 2;
        });
    return derived;
}

// The printed wave of a set of chains: the polyline through the apexes of each chain.
OpenLinesSet chainsToWaves(const std::vector<ToothChain>& chains)
{
    OpenLinesSet waves;
    for (const ToothChain& chain : chains)
    {
        if (chain.teeth.size() < 2)
        {
            continue;
        }
        std::vector<Point2LL> points;
        points.reserve(chain.teeth.size() + 1);
        for (const TrackedTooth& tooth : chain.teeth)
        {
            points.push_back(tooth.apex);
        }
        if (chain.closed)
        {
            points.push_back(points.front());
        }
        waves.push_back(OpenPolyline{ points });
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

TriangleWaveEpicFillProvider::TriangleWaveEpicFillProvider(const std::vector<Shape>& layer_outlines, coord_t line_distance, const double fill_angle)
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

void TriangleWaveEpicFillProvider::generate(OpenLinesSet& result_lines, const Shape& in_outline) const
{
    if (template_wave_.empty() || in_outline.empty())
    {
        return;
    }

    Shape outline = in_outline;
    outline.applyMatrix(rotation_matrix_);

    result_lines = clipWave(template_wave_, outline, rotation_matrix_);
}

TriangleWaveEpicTrackingProvider::TriangleWaveEpicTrackingProvider(const std::vector<Shape>& layer_outlines, coord_t line_distance, const double fill_angle)
    : rotation_matrix_(fill_angle)
{
    if (line_distance <= 0)
    {
        return;
    }

    layer_waves_.resize(layer_outlines.size());

#ifdef TW_EPIC_BENCH
    const auto bench_start = std::chrono::steady_clock::now();
    auto bench_last = bench_start;
#endif

    // Build the layers bottom up. Only the first layer (and any part which newly appears on a
    // higher layer) runs the full skeleton pipeline; every other layer is DERIVED tooth by tooth
    // from the layer below, changing as little as possible: rib lines stay in place, only the
    // tooth lengths follow the walls, and teeth are only added/removed where the region actually
    // appeared or disappeared. Identical outlines therefore yield identical waves.
    std::vector<ToothChain> chains;
    for (size_t layer_idx = 0; layer_idx < layer_outlines.size(); ++layer_idx)
    {
        Shape rotated = layer_outlines[layer_idx];
        rotated.applyMatrix(rotation_matrix_);
        rotated = rotated.unionPolygons();

        const AABB region_box(rotated);
        const double ray_length = vSizeMM(region_box.max_ - region_box.min_) * 1000.0 + static_cast<double>(line_distance);

        std::vector<ToothChain> layer_chains = deriveChains(chains, rotated, line_distance, ray_length);
        OpenLinesSet waves = chainsToWaves(layer_chains);

        // Parts not reached by any derived tooth (new islands, or the very first layer) get a
        // fresh wave from the full skeleton pipeline.
        for (const SingleShape& part : rotated.splitIntoParts())
        {
            const bool covered = std::any_of(
                layer_chains.begin(),
                layer_chains.end(),
                [&part](const ToothChain& chain)
                {
                    return std::any_of(
                        chain.teeth.begin(),
                        chain.teeth.end(),
                        [&part](const TrackedTooth& tooth)
                        {
                            return part.inside(tooth.foot, true);
                        });
                });
            if (covered)
            {
                continue;
            }
            std::vector<ToothChain> fresh_chains;
            OpenLinesSet part_waves = buildSkeletonWaves(part, line_distance, &fresh_chains);
            if (part_waves.empty())
            {
                part_waves = buildStraightWave(part, line_distance, &fresh_chains);
            }
            waves.push_back(part_waves);
            layer_chains.insert(layer_chains.end(), fresh_chains.begin(), fresh_chains.end());
        }

        layer_waves_[layer_idx] = std::move(waves);
        chains = std::move(layer_chains);

#ifdef TW_EPIC_BENCH
        const auto bench_now = std::chrono::steady_clock::now();
        size_t teeth_count = 0;
        for (const ToothChain& c : chains)
        {
            teeth_count += c.teeth.size();
        }
        std::fprintf(
            stderr,
            "layer %zu: %.1f ms (total %.1f s), chains %zu, teeth %zu\n",
            layer_idx,
            std::chrono::duration<double, std::milli>(bench_now - bench_last).count(),
            std::chrono::duration<double>(bench_now - bench_start).count(),
            chains.size(),
            teeth_count);
        bench_last = bench_now;
#endif
    }
}

void TriangleWaveEpicTrackingProvider::generate(OpenLinesSet& result_lines, const Shape& in_outline, size_t layer_idx) const
{
    if (layer_idx >= layer_waves_.size() || layer_waves_[layer_idx].empty() || in_outline.empty())
    {
        return;
    }

    Shape outline = in_outline;
    outline.applyMatrix(rotation_matrix_);

    result_lines = clipWave(layer_waves_[layer_idx], outline, rotation_matrix_);
}

void TriangleWaveInfillEpic::generateTotalTriangleWaveInfill(OpenLinesSet& result_lines, coord_t line_distance, const Shape& in_outline, const double fill_angle)
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
