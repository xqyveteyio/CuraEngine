// Copyright (c) 2026 UltiMaker
// CuraEngine is released under the terms of the AGPLv3 or higher

#include "infill/MedialZigzag.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <queue>
#include <vector>

#include "geometry/OpenPolyline.h"
#include "geometry/PointMatrix.h"
#include "geometry/Polygon.h"
#include "geometry/Shape.h"
#include "geometry/SingleShape.h"
#include "settings/types/Angle.h"
#include "sliceDataStorage.h"
#include "utils/AABB.h"
#include "utils/PolygonsSegmentIndex.h"
#include "utils/Simplify.h"
#include "utils/VoronoiUtils.h"

#include "BoostInterface.hpp"

namespace cura
{

namespace
{

double distanceF(const Point2LL& a, const Point2LL& b)
{
    return std::hypot(double(a.X - b.X), double(a.Y - b.Y));
}

/*!
 * A wall (or the medial axis) as a chain of points with precomputed cumulative arc lengths,
 * so that points can be sampled at arbitrary arc-length positions along the chain.
 */
struct WallChain
{
    std::vector<Point2LL> points_;
    std::vector<double> cumulative_; //!< Arc length from the start of the chain up to each point.
    double total_ = 0.0;
    bool closed_ = false;

    void finish()
    {
        cumulative_.assign(points_.size(), 0.0);
        for (size_t point_idx = 1; point_idx < points_.size(); point_idx++)
        {
            cumulative_[point_idx] = cumulative_[point_idx - 1] + distanceF(points_[point_idx], points_[point_idx - 1]);
        }
        total_ = cumulative_.empty() ? 0.0 : cumulative_.back();
        if (closed_ && points_.size() >= 2)
        {
            total_ += distanceF(points_.front(), points_.back());
        }
    }

    double segmentLength(const size_t start_idx) const
    {
        return (start_idx + 1 < points_.size()) ? cumulative_[start_idx + 1] - cumulative_[start_idx] : total_ - cumulative_[start_idx];
    }

    Point2LL pointAt(double s) const
    {
        if (points_.empty())
        {
            return Point2LL();
        }
        if (points_.size() == 1 || total_ <= 0.0)
        {
            return points_.front();
        }
        if (closed_)
        {
            s = std::fmod(s, total_);
            if (s < 0.0)
            {
                s += total_;
            }
        }
        else
        {
            s = std::clamp(s, 0.0, total_);
        }
        const auto it = std::upper_bound(cumulative_.begin(), cumulative_.end(), s);
        const size_t idx = (it == cumulative_.begin()) ? 0 : std::distance(cumulative_.begin(), it) - 1;
        const Point2LL& p0 = points_[idx];
        const Point2LL& p1 = points_[(idx + 1) % points_.size()];
        const double segment_length = segmentLength(idx);
        const double t = (segment_length <= 0.0) ? 0.0 : std::clamp((s - cumulative_[idx]) / segment_length, 0.0, 1.0);
        return Point2LL(std::llround(p0.X + t * (p1.X - p0.X)), std::llround(p0.Y + t * (p1.Y - p0.Y)));
    }
};

/*!
 * Make a closed wall chain from a polygon, oriented counter-clockwise so that all walls of a part
 * advance in the same rotational direction.
 */
WallChain closedWall(const Polygon& polygon)
{
    WallChain wall;
    wall.closed_ = true;
    wall.points_.assign(polygon.begin(), polygon.end());
    if (polygon.area() < 0)
    {
        std::reverse(wall.points_.begin(), wall.points_.end());
    }
    wall.finish();
    return wall;
}

/*!
 * The arc position along an open chain of the point closest to \p point. Used to project the phase
 * anchor of the layer below onto the medial axis of the current layer.
 */
double nearestArcPosition(const WallChain& chain, const Point2LL& point)
{
    double best_distance2 = std::numeric_limits<double>::max();
    double best_s = 0.0;
    for (size_t point_idx = 0; point_idx + 1 < chain.points_.size(); point_idx++)
    {
        const Point2LL& p0 = chain.points_[point_idx];
        const Point2LL& p1 = chain.points_[point_idx + 1];
        const double dx = double(p1.X - p0.X);
        const double dy = double(p1.Y - p0.Y);
        const double length2 = dx * dx + dy * dy;
        const double t = (length2 <= 0.0) ? 0.0 : std::clamp((double(point.X - p0.X) * dx + double(point.Y - p0.Y) * dy) / length2, 0.0, 1.0);
        const double proj_x = p0.X + t * dx;
        const double proj_y = p0.Y + t * dy;
        const double distance2 = (point.X - proj_x) * (point.X - proj_x) + (point.Y - proj_y) * (point.Y - proj_y);
        if (distance2 < best_distance2)
        {
            best_distance2 = distance2;
            best_s = chain.cumulative_[point_idx] + t * std::sqrt(length2);
        }
    }
    return best_s;
}

/*!
 * The arc position where the horizontal ray from \p origin towards +X crosses the wall, taking the
 * outermost (highest X) crossing. Used to give the walls of a ring a common, layer-stable parameter origin.
 */
double rayCrossingParameter(const WallChain& wall, const Point2LL& origin)
{
    double best_x = std::numeric_limits<double>::lowest();
    double best_s = 0.0;
    const size_t point_count = wall.points_.size();
    for (size_t point_idx = 0; point_idx < point_count; point_idx++)
    {
        const Point2LL& p0 = wall.points_[point_idx];
        const Point2LL& p1 = wall.points_[(point_idx + 1) % point_count];
        if ((p0.Y <= origin.Y) == (p1.Y <= origin.Y))
        {
            continue;
        }
        const double t = double(origin.Y - p0.Y) / double(p1.Y - p0.Y);
        const double x = p0.X + t * double(p1.X - p0.X);
        if (x >= origin.X && x > best_x)
        {
            best_x = x;
            best_s = wall.cumulative_[point_idx] + t * wall.segmentLength(point_idx);
        }
    }
    return best_s;
}

/*!
 * The arc position of the first crossing of the ray (\p origin, direction \p direction) with the wall,
 * or nullopt when the ray doesn't hit the wall.
 */
std::optional<double> rayBoundaryArcPosition(const WallChain& wall, const Point2LL& origin, const Point2LL& direction)
{
    const auto cross = [](const double ax, const double ay, const double bx, const double by)
    {
        return ax * by - ay * bx;
    };
    double best_t = std::numeric_limits<double>::max();
    std::optional<double> best_s;
    const size_t point_count = wall.points_.size();
    for (size_t point_idx = 0; point_idx < point_count; point_idx++)
    {
        const Point2LL& p0 = wall.points_[point_idx];
        const Point2LL& p1 = wall.points_[(point_idx + 1) % point_count];
        const double denominator = cross(direction.X, direction.Y, p1.X - p0.X, p1.Y - p0.Y);
        if (std::abs(denominator) < 1e-9)
        {
            continue;
        }
        const double t = cross(p0.X - origin.X, p0.Y - origin.Y, p1.X - p0.X, p1.Y - p0.Y) / denominator;
        const double u = cross(p0.X - origin.X, p0.Y - origin.Y, direction.X, direction.Y) / denominator;
        if (t > 1e-9 && u >= 0.0 && u < 1.0 && t < best_t)
        {
            best_t = t;
            best_s = wall.cumulative_[point_idx] + u * wall.segmentLength(point_idx);
        }
    }
    return best_s;
}

/*!
 * Approximate the medial axis of a part without holes: build the voronoi diagram of the boundary
 * segments, keep the internal skeleton edges (pruning the spurs that connect to boundary corners),
 * and take the longest path through the remaining skeleton graph.
 *
 * \return The medial axis polyline, or nullopt when there is no usable axis (e.g. round-ish parts).
 */
std::optional<OpenPolyline> computeMedialAxisPath(const SingleShape& part)
{
    using vd_t = boost::polygon::voronoi_diagram<double>;

    Shape voronoi_input;
    voronoi_input.push_back(part.outerPolygon());
    const Polygon& outline = voronoi_input[0];
    const size_t num_segments = outline.size();
    if (num_segments < 3)
    {
        return std::nullopt;
    }

    std::vector<Point2LL> points; // Remains empty; the input consists purely of boundary segments.
    std::vector<PolygonsSegmentIndex> segments;
    segments.reserve(num_segments);
    for (size_t segment_idx = 0; segment_idx < num_segments; segment_idx++)
    {
        segments.emplace_back(&voronoi_input, 0, segment_idx);
    }

    vd_t voronoi_diagram;
    boost::polygon::construct_voronoi(segments.begin(), segments.end(), &voronoi_diagram);

    // The boundary vertices touched by the source feature (segment or segment endpoint) of a voronoi cell.
    const auto touched_vertices = [&](const vd_t::cell_type& cell) -> std::pair<size_t, size_t>
    {
        const size_t segment_idx = cell.source_index();
        if (cell.contains_segment())
        {
            return { segment_idx, (segment_idx + 1) % num_segments };
        }
        const Point2LL source_point = VoronoiUtils::getSourcePoint(cell, points, segments);
        const size_t vertex_idx = (source_point == outline[segment_idx]) ? segment_idx : (segment_idx + 1) % num_segments;
        return { vertex_idx, vertex_idx };
    };

    std::map<std::pair<coord_t, coord_t>, size_t> node_ids;
    std::vector<Point2LL> nodes;
    std::vector<std::vector<std::pair<size_t, double>>> adjacency;
    const auto node_id = [&](const Point2LL& p)
    {
        const auto [it, inserted] = node_ids.emplace(std::make_pair(p.X, p.Y), nodes.size());
        if (inserted)
        {
            nodes.push_back(p);
            adjacency.emplace_back();
        }
        return it->second;
    };

    size_t start_node = std::numeric_limits<size_t>::max();
    double longest_edge = -1.0;
    for (const vd_t::edge_type& edge : voronoi_diagram.edges())
    {
        if (! edge.is_finite() || ! edge.is_primary() || &edge > edge.twin())
        {
            continue; // Handle each edge pair once.
        }
        const Point2LL p0 = VoronoiUtils::p(edge.vertex0());
        const Point2LL p1 = VoronoiUtils::p(edge.vertex1());
        if (p0 == p1)
        {
            continue;
        }
        const auto [a0, a1] = touched_vertices(*edge.cell());
        const auto [b0, b1] = touched_vertices(*edge.twin()->cell());
        if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1)
        {
            continue; // Spur between two adjacent boundary features (runs into a corner); not part of the medial axis.
        }
        if (! part.inside(Point2LL((p0.X + p1.X) / 2, (p0.Y + p1.Y) / 2), false))
        {
            continue; // Voronoi edge outside of the part.
        }
        const double length = distanceF(p0, p1);
        const size_t node_0 = node_id(p0);
        const size_t node_1 = node_id(p1);
        adjacency[node_0].emplace_back(node_1, length);
        adjacency[node_1].emplace_back(node_0, length);
        if (length > longest_edge)
        {
            longest_edge = length;
            start_node = node_0;
        }
    }
    if (nodes.size() < 2 || start_node == std::numeric_limits<size_t>::max())
    {
        return std::nullopt;
    }

    // The medial axis is approximated by the longest path through the skeleton: run two farthest-point
    // (Dijkstra) searches, starting in the component that contains the longest skeleton edge.
    const auto farthest_from = [&](const size_t start) -> std::pair<std::vector<double>, std::vector<size_t>>
    {
        std::vector<double> distance(nodes.size(), std::numeric_limits<double>::infinity());
        std::vector<size_t> parent(nodes.size(), std::numeric_limits<size_t>::max());
        using QueueEntry = std::pair<double, size_t>;
        std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
        distance[start] = 0.0;
        queue.emplace(0.0, start);
        while (! queue.empty())
        {
            const auto [d, u] = queue.top();
            queue.pop();
            if (d > distance[u])
            {
                continue;
            }
            for (const auto& [v, w] : adjacency[u])
            {
                if (distance[u] + w < distance[v])
                {
                    distance[v] = distance[u] + w;
                    parent[v] = u;
                    queue.emplace(distance[v], v);
                }
            }
        }
        return { std::move(distance), std::move(parent) };
    };
    const auto farthest_reachable = [](const std::vector<double>& distance)
    {
        size_t best = 0;
        double best_distance = -1.0;
        for (size_t node = 0; node < distance.size(); node++)
        {
            if (std::isfinite(distance[node]) && distance[node] > best_distance)
            {
                best_distance = distance[node];
                best = node;
            }
        }
        return best;
    };

    const size_t end_a = farthest_reachable(farthest_from(start_node).first);
    const auto [distance_b, parent_b] = farthest_from(end_a);
    const size_t end_b = farthest_reachable(distance_b);

    OpenPolyline axis;
    for (size_t node = end_b; node != std::numeric_limits<size_t>::max(); node = parent_b[node])
    {
        axis.push_back(nodes[node]);
    }
    if (axis.size() < 2)
    {
        return std::nullopt;
    }
    return axis;
}

/*!
 * The anchor of the layer below whose match point lies inside \p part, or nullptr when the part has
 * no counterpart in the layer below.
 */
const MedialZigzagAnchor* findPreviousAnchor(const std::vector<MedialZigzagAnchor>& previous_anchors, const SingleShape& part)
{
    for (const MedialZigzagAnchor& anchor : previous_anchors)
    {
        if (part.inside(anchor.match_point, true))
        {
            return &anchor;
        }
    }
    return nullptr;
}

/*!
 * The number of node periods to use: normally the rounded ratio of axis length and reference period,
 * but when the layer below used a period count that is still within tolerance, that count is kept
 * (hysteresis), so the node grid doesn't re-space on the layer where the ratio crosses a rounding boundary.
 */
size_t choosePeriodCount(const double axis_length, const coord_t period, const MedialZigzagAnchor* previous, const bool ring)
{
    const double target = axis_length / period;
    const size_t minimum = ring ? 3 : 1;
    if (previous != nullptr && previous->is_ring == ring && previous->num_periods >= minimum && std::abs(target - double(previous->num_periods)) <= 0.75)
    {
        return previous->num_periods;
    }
    return std::max<size_t>(minimum, std::llround(target));
}

/*!
 * Generate the zigzag for a ring-like part (a part with a dominant hole, e.g. an annulus): the two
 * walls are the outer wall and the hole wall. Both walls get the same number of node periods, each
 * distributed uniformly along its own arc length, half a period out of phase, connected into a
 * closed alternating zigzag loop.
 */
bool generateRingZigzag(const SingleShape& part, const coord_t period, const MedialZigzagAnchor* previous, MedialZigzagAnchor& anchor, OpenLinesSet& zigzag_lines)
{
    // Find the dominant hole.
    size_t hole_idx = 0;
    double hole_length = 0.0;
    for (size_t poly_idx = 1; poly_idx < part.size(); poly_idx++)
    {
        const double length = double(part[poly_idx].length());
        if (length > hole_length)
        {
            hole_length = length;
            hole_idx = poly_idx;
        }
    }
    const double outer_length = double(part.outerPolygon().length());
    if (hole_idx == 0 || hole_length < 0.25 * outer_length)
    {
        return false; // The hole is a small cutout rather than the inner wall of a ring; handle the part differently.
    }

    const WallChain outer_wall = closedWall(part.outerPolygon());
    const WallChain inner_wall = closedWall(part[hole_idx]);
    if (outer_wall.total_ <= 0.0 || inner_wall.total_ <= 0.0)
    {
        return false;
    }

    // The medial axis of a ring is approximated by the average of its two walls. The common number of
    // node periods follows from its length and the reference node period, kept from the layer below
    // while it is still within tolerance.
    const double axis_length = (outer_wall.total_ + inner_wall.total_) / 2.0;
    const size_t num_periods = choosePeriodCount(axis_length, period, previous, true);

    // Give both walls a common parameter origin (on the +X ray from the hole center) so that the phase
    // relation between the walls is well-defined and stable across layers.
    const Point2LL origin = AABB(part[hole_idx]).getMiddle();
    const double anchor_outer = rayCrossingParameter(outer_wall, origin);
    const double anchor_inner = rayCrossingParameter(inner_wall, origin);

    OpenPolyline zigzag;
    for (size_t node_idx = 0; node_idx < num_periods; node_idx++)
    {
        zigzag.push_back(outer_wall.pointAt(anchor_outer + double(node_idx) / num_periods * outer_wall.total_));
        zigzag.push_back(inner_wall.pointAt(anchor_inner + (node_idx + 0.5) / num_periods * inner_wall.total_));
    }
    zigzag.push_back(outer_wall.pointAt(anchor_outer)); // Close the loop.
    zigzag_lines.push_back(std::move(zigzag));

    anchor.is_ring = true;
    anchor.num_periods = num_periods;
    const Point2LL on_outer = outer_wall.pointAt(anchor_outer);
    const Point2LL on_inner = inner_wall.pointAt(anchor_inner);
    anchor.match_point = Point2LL((on_outer.X + on_inner.X) / 2, (on_outer.Y + on_inner.Y) / 2); // Between the walls: inside the ring body.
    return true;
}

/*!
 * Generate the zigzag for an elongated part without holes (slab, L-shape, C-shape, ...): compute the
 * medial axis and sample node positions on a regular arc-length grid along it, alternating between the
 * two sides. At each sample position the node is placed where the local normal of the axis (the local
 * 'node cutting plane') first hits the wall on that side.
 *
 * When the layer below produced a matching part, the grid of this layer is anchored to it: the axis
 * keeps the same orientation, the period count is kept while within tolerance, and the grid phase is
 * chosen such that the samples land where the samples of the layer below were (the reference sample of
 * the layer below is projected onto the current axis). Without a layer below, the grid is centered on
 * the middle of the axis and the orientation is chosen canonically, both of which are stable choices too.
 */
bool generateRibbonZigzag(const SingleShape& part, const coord_t period, const MedialZigzagAnchor* previous, MedialZigzagAnchor& anchor, OpenLinesSet& zigzag_lines)
{
    std::optional<OpenPolyline> axis = computeMedialAxisPath(part);
    if (! axis.has_value())
    {
        return false;
    }

    // Keep the axis orientation stable across layers: match the end points to the layer below, or
    // without a layer below orient by coordinates. The orientation determines which side is 'upper'.
    const auto oriented_backwards = [&]() -> bool
    {
        const Point2LL& front = axis->front();
        const Point2LL& back = axis->back();
        if (previous != nullptr && ! previous->is_ring)
        {
            const double keep = distanceF(front, previous->axis_front) + distanceF(back, previous->axis_back);
            const double flip = distanceF(front, previous->axis_back) + distanceF(back, previous->axis_front);
            return flip < keep;
        }
        return std::make_pair(back.X, back.Y) < std::make_pair(front.X, front.Y);
    };
    if (oriented_backwards())
    {
        std::reverse(axis->begin(), axis->end());
    }

    WallChain axis_chain;
    axis_chain.points_.assign(axis->begin(), axis->end());
    axis_chain.finish();
    const double axis_length = axis_chain.total_;
    if (axis_length < period / 2.0)
    {
        return false; // Too short to give a meaningful direction; handle the part differently.
    }

    const WallChain wall = closedWall(part.outerPolygon());

    // The number of periods and the node spacing along the axis. The spacing is kept *exactly* from the
    // layer below while the part still fits it (hysteresis): the measured axis length jitters a little
    // from layer to layer (especially at its ends), and re-deriving the spacing from it every layer
    // would wobble all node positions far from the grid anchor.
    const double target = axis_length / period;
    size_t num_periods;
    double axis_period;
    if (previous != nullptr && ! previous->is_ring && previous->axis_period > 0.0 && previous->num_periods >= 1
        && std::abs(target - double(previous->num_periods)) <= 0.75)
    {
        num_periods = previous->num_periods;
        axis_period = previous->axis_period;
    }
    else
    {
        num_periods = std::max<size_t>(1, std::llround(target));
        axis_period = axis_length / double(num_periods);
    }

    // The phase of the sample grid: the arc position of one upper-side sample. Anchored to the layer
    // below when available, otherwise to the middle of the axis.
    const double grid_anchor = (previous != nullptr && ! previous->is_ring) ? nearestArcPosition(axis_chain, previous->phase_point) : axis_length / 2.0;

    // Enumerate the grid: upper samples at grid_anchor + k * axis_period, lower samples half a period
    // further, keeping a margin to the axis ends (where the axis direction is poorly defined).
    const double margin = 0.2 * axis_period;
    const long step_min = static_cast<long>(std::ceil((margin - grid_anchor) / (0.5 * axis_period)));
    const long step_max = static_cast<long>(std::floor((axis_length - margin - grid_anchor) / (0.5 * axis_period)));

    bool produced = false;
    OpenPolyline zigzag;
    const auto flush_zigzag = [&zigzag, &zigzag_lines, &produced]()
    {
        if (zigzag.size() >= 2)
        {
            zigzag_lines.push_back(std::move(zigzag));
            produced = true;
        }
        zigzag.clear();
    };

    for (long step = step_min; step <= step_max; step++)
    {
        const double s = grid_anchor + 0.5 * axis_period * step;
        const Point2LL origin = axis_chain.pointAt(s);

        // The local axis direction, averaged over half a period to be insensitive to skeleton noise.
        const double window = std::min(axis_length / 2.0, period / 2.0);
        const Point2LL ahead = axis_chain.pointAt(std::min(s + window, axis_length));
        const Point2LL behind = axis_chain.pointAt(std::max(s - window, 0.0));
        const double direction_x = double(ahead.X - behind.X);
        const double direction_y = double(ahead.Y - behind.Y);
        const double direction_length = std::hypot(direction_x, direction_y);
        if (direction_length <= 0.0)
        {
            flush_zigzag();
            continue;
        }

        // The local normal, pointing towards the side this node belongs to. Even grid steps are the
        // upper-side samples, odd steps the lower-side samples.
        const bool upper = ((step % 2) + 2) % 2 == 0;
        const double side = upper ? 1.0 : -1.0;
        constexpr double direction_scale = 65536.0; // The ray direction is only used for its direction, but has integer coordinates.
        const Point2LL normal(std::llround(-direction_y / direction_length * direction_scale * side), std::llround(direction_x / direction_length * direction_scale * side));

        // The node lies where the normal ray first hits the wall of this side.
        const std::optional<double> s_wall = rayBoundaryArcPosition(wall, origin, normal);
        if (! s_wall.has_value())
        {
            flush_zigzag();
            continue;
        }
        zigzag.push_back(wall.pointAt(*s_wall));
    }
    flush_zigzag();

    if (produced)
    {
        anchor.is_ring = false;
        anchor.num_periods = num_periods;
        anchor.axis_period = axis_period;
        anchor.axis_front = axis_chain.points_.front();
        anchor.axis_back = axis_chain.points_.back();
        anchor.match_point = axis_chain.pointAt(axis_length / 2.0);
        // The reference sample for the layer above: the upper-side grid sample closest to the axis middle.
        const double s_reference = std::clamp(grid_anchor + std::round((axis_length / 2.0 - grid_anchor) / axis_period) * axis_period, 0.0, axis_length);
        anchor.phase_point = axis_chain.pointAt(s_reference);
    }
    return produced;
}

int computeScanSegmentIdx(const int x, const int line_width)
{
    if (x < 0)
    {
        return (x + 1) / line_width - 1;
    }
    return x / line_width;
}

/*!
 * Fallback for parts without a usable medial axis (e.g. round-ish blobs): place the nodes where fixed,
 * equally spaced cutting planes (perpendicular to the infill angle, anchored to the infill origin)
 * cross the upper resp. lower wall of the part. This scheme is anchored in world space, so it is
 * inherently stable across layers.
 */
void generatePlaneZigzag(const Polygon& outer_polygon, const PointMatrix& rotation_matrix, const coord_t shift, const coord_t plane_distance, OpenLinesSet& zigzag_lines)
{
    Polygon outline = outer_polygon;
    outline.applyMatrix(rotation_matrix);
    if (outline.size() < 3)
    {
        return;
    }

    coord_t x_min = std::numeric_limits<coord_t>::max();
    coord_t x_max = std::numeric_limits<coord_t>::lowest();
    for (const Point2LL& point : outline)
    {
        x_min = std::min(x_min, point.X);
        x_max = std::max(x_max, point.X);
    }
    const int first_plane_idx = computeScanSegmentIdx(x_min - shift, plane_distance) + 1;
    const int last_plane_idx = computeScanSegmentIdx(x_max - shift, plane_distance);

    OpenPolyline zigzag;
    const auto flush_zigzag = [&zigzag, &zigzag_lines]()
    {
        if (zigzag.size() >= 2)
        {
            zigzag_lines.push_back(std::move(zigzag));
        }
        zigzag.clear();
    };

    for (int plane_idx = first_plane_idx; plane_idx <= last_plane_idx; ++plane_idx)
    {
        const coord_t x = plane_idx * plane_distance + shift;

        // Find where this cutting plane crosses the wall.
        coord_t y_min = std::numeric_limits<coord_t>::max();
        coord_t y_max = std::numeric_limits<coord_t>::lowest();
        bool crossed = false;
        for (size_t point_idx = 0; point_idx < outline.size(); point_idx++)
        {
            const Point2LL& p0 = outline[point_idx];
            const Point2LL& p1 = outline[(point_idx + 1) % outline.size()];
            if ((p0.X <= x && p1.X > x) || (p1.X <= x && p0.X > x))
            {
                const coord_t y = p0.Y + (p1.Y - p0.Y) * (x - p0.X) / (p1.X - p0.X);
                y_min = std::min(y_min, y);
                y_max = std::max(y_max, y);
                crossed = true;
            }
        }
        if (! crossed)
        {
            flush_zigzag();
            continue;
        }

        // Even planes place their node on the upper wall, odd planes on the lower wall. The parity is
        // based on the absolute plane index, so it is the same on every layer and in every part.
        const bool upper_wall = (plane_idx % 2) == 0;
        zigzag.push_back(rotation_matrix.unapply(Point2LL(x, upper_wall ? y_max : y_min)));
    }
    flush_zigzag();
}

} // namespace

void generateMedialZigzagLines(
    const Shape& region,
    const coord_t line_distance,
    const double plane_angle,
    coord_t plane_shift,
    const std::vector<MedialZigzagAnchor>& previous_anchors,
    std::vector<MedialZigzagAnchor>& new_anchors,
    OpenLinesSet& result)
{
    if (line_distance <= 0 || region.empty())
    {
        return;
    }

    // Distance between two nodes on the same wall, measured along the wall; consecutive zigzag nodes
    // (which alternate between the two walls) are half of this period apart.
    const coord_t period = 2 * line_distance;

    // Parameters for the cutting-plane fallback used for parts without a usable medial axis.
    const PointMatrix rotation_matrix(plane_angle);
    plane_shift = ((plane_shift % line_distance) + line_distance) % line_distance;

    // Simplify the walls before analyzing them; the voronoi diagram and the node placement don't need
    // (and shouldn't be disturbed by) micron-sized boundary details.
    const Simplify simplifier(100, 25, 0);

    const coord_t opening_radius = line_distance / 2;

    OpenLinesSet zigzag_lines;
    for (const SingleShape& raw_part : region.splitIntoParts())
    {
        SingleShape part{ simplifier.polygon(raw_part) };
        if (part.empty() || part.outerPolygon().size() < 3)
        {
            continue;
        }

        // Remove protruding details (e.g. rows of sawteeth) from the walls with a morphological opening,
        // so the nodes are placed on the main wall envelope and never wander into a protrusion. The
        // opening is only used when it doesn't shatter or shrink the part (which happens when the part
        // itself is barely wider than the opening); the lines are clipped against the real region at the
        // end either way, so they always stay inside the solid region.
        const std::vector<SingleShape> opened_parts = part.offset(-opening_radius).offset(opening_radius).splitIntoParts();
        if (opened_parts.size() == 1 && opened_parts.front().area() >= 0.7 * part.area())
        {
            part = SingleShape{ simplifier.polygon(opened_parts.front()) };
            if (part.empty() || part.outerPolygon().size() < 3)
            {
                continue;
            }
        }
        const MedialZigzagAnchor* previous = findPreviousAnchor(previous_anchors, part);

        MedialZigzagAnchor anchor;
        bool generated = false;
        if (part.size() > 1)
        {
            // Ring-like part: the two walls are the outer wall and the wall of the dominant hole.
            generated = generateRingZigzag(part, period, previous, anchor, zigzag_lines);
        }
        if (! generated)
        {
            // Elongated part: the two walls are the boundary chains on either side of the medial axis.
            generated = generateRibbonZigzag(part, period, previous, anchor, zigzag_lines);
        }
        if (generated)
        {
            new_anchors.push_back(anchor);
        }
        else
        {
            generatePlaneZigzag(part.outerPolygon(), rotation_matrix, plane_shift, line_distance, zigzag_lines);
        }
    }

    // Clip against the actual region (including its holes and details): segments crossing a hole or
    // leaving the part through a concavity are cut at the boundary and continue where they re-enter.
    result.push_back(region.intersection(zigzag_lines));
}

MedialZigzagGenerator::MedialZigzagGenerator(const SliceMeshStorage& mesh)
{
    const coord_t line_distance = mesh.settings.get<coord_t>("infill_line_distance");
    const auto infill_wall_line_count = static_cast<coord_t>(mesh.settings.get<size_t>("infill_wall_line_count"));
    const coord_t infill_line_width = mesh.settings.get<coord_t>("infill_line_width");
    const coord_t infill_overlap = mesh.settings.get<coord_t>("infill_overlap_mm");
    // Match the contour that Infill will actually use for this pattern (infill walls, overlap, and the
    // half line width inset since the pattern prints along the walls), minus a small epsilon, so the
    // pre-computed lines survive Infill's clip without getting trimmed at every node.
    const coord_t region_offset = -infill_wall_line_count * infill_line_width + infill_overlap - infill_line_width / 2 - 10;

    // Parameters of the cutting-plane fallback, replicating what FffGcodeWriter passes for this pattern:
    // angle 0 unless overridden by infill_angles, origin in the middle of the mesh.
    const auto angles = mesh.settings.get<std::vector<AngleDegrees>>("infill_angles");
    const double plane_angle = angles.empty() ? 0.0 : double(angles.front());
    const Point3LL mesh_middle = mesh.bounding_box.getMiddle();
    const Point2LL origin(mesh_middle.x_ + mesh.settings.get<coord_t>("infill_offset_x"), mesh_middle.y_ + mesh.settings.get<coord_t>("infill_offset_y"));
    coord_t plane_shift = 0;
    if (origin.X != 0 || origin.Y != 0)
    {
        const double rotation_rads = plane_angle * std::numbers::pi / 180;
        plane_shift = origin.X * std::cos(rotation_rads) - origin.Y * std::sin(rotation_rads);
    }

    // Generate the layers bottom-up, anchoring each layer to the anchors of the layer below.
    lines_per_layer_.resize(mesh.layers.size());
    std::vector<MedialZigzagAnchor> previous_anchors;
    std::vector<MedialZigzagAnchor> current_anchors;
    for (size_t layer_nr = 0; layer_nr < mesh.layers.size(); layer_nr++)
    {
        Shape region;
        for (const SliceLayerPart& part : mesh.layers[layer_nr].parts)
        {
            region.push_back(part.getOwnInfillArea().offset(region_offset));
        }
        region = region.unionPolygons();
        if (region.empty())
        {
            continue; // Keep the previous anchors, so the pattern continues unchanged above a gap.
        }
        current_anchors.clear();
        generateMedialZigzagLines(region, line_distance, plane_angle, plane_shift, previous_anchors, current_anchors, lines_per_layer_[layer_nr]);
        std::swap(previous_anchors, current_anchors);
    }
}

const OpenLinesSet& MedialZigzagGenerator::getLinesForLayer(const size_t layer_nr) const
{
    static const OpenLinesSet empty;
    return layer_nr < lines_per_layer_.size() ? lines_per_layer_[layer_nr] : empty;
}

} // namespace cura
