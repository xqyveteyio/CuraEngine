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
 * The arc position along a chain of the point closest to \p point.
 */
double nearestArcPosition(const WallChain& chain, const Point2LL& point)
{
    const size_t point_count = chain.points_.size();
    const size_t segment_count = chain.closed_ ? point_count : (point_count > 0 ? point_count - 1 : 0);
    double best_distance2 = std::numeric_limits<double>::max();
    double best_s = 0.0;
    for (size_t segment_idx = 0; segment_idx < segment_count; segment_idx++)
    {
        const Point2LL& p0 = chain.points_[segment_idx];
        const Point2LL& p1 = chain.points_[(segment_idx + 1) % point_count];
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
            best_s = chain.cumulative_[segment_idx] + t * chain.segmentLength(segment_idx);
        }
    }
    return best_s;
}

/*!
 * The portion of a closed wall from arc position \p s0 forward (in wall direction) to \p s1, as an
 * open chain.
 */
WallChain subChain(const WallChain& wall, const double s0, const double s1)
{
    WallChain chain;
    chain.closed_ = false;
    double span = s1 - s0;
    if (span < 0.0)
    {
        span += wall.total_;
    }
    chain.points_.push_back(wall.pointAt(s0));
    // Collect the wall vertices that lie strictly between s0 and s1 (walking forward from s0), in order.
    std::vector<std::pair<double, size_t>> in_range;
    for (size_t point_idx = 0; point_idx < wall.points_.size(); point_idx++)
    {
        double delta = wall.cumulative_[point_idx] - s0;
        if (delta < 0.0)
        {
            delta += wall.total_;
        }
        if (delta > 1.0 && delta < span - 1.0) // 1 micron margin against duplicating the end points.
        {
            in_range.emplace_back(delta, point_idx);
        }
    }
    std::sort(in_range.begin(), in_range.end());
    for (const auto& [delta, point_idx] : in_range)
    {
        chain.points_.push_back(wall.points_[point_idx]);
    }
    chain.points_.push_back(wall.pointAt(s0 + span));
    chain.finish();
    return chain;
}

/*!
 * The arc position where the ray from \p origin in the given direction first crosses the wall, or
 * nullopt when it never does.
 */
std::optional<double> rayArcPosition(const WallChain& wall, const Point2LL& origin, const double direction_x, const double direction_y)
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
        const double denominator = cross(direction_x, direction_y, p1.X - p0.X, p1.Y - p0.Y);
        if (std::abs(denominator) < 1e-9)
        {
            continue;
        }
        const double t = cross(p0.X - origin.X, p0.Y - origin.Y, p1.X - p0.X, p1.Y - p0.Y) / denominator;
        const double u = cross(p0.X - origin.X, p0.Y - origin.Y, direction_x, direction_y) / denominator;
        if (t > 1e-6 && u >= 0.0 && u < 1.0 && t < best_t)
        {
            best_t = t;
            best_s = wall.cumulative_[point_idx] + u * wall.segmentLength(point_idx);
        }
    }
    return best_s;
}

/*!
 * The arc position where the horizontal ray from \p origin towards +X crosses the wall, taking the
 * outermost (highest X) crossing. Used to give the walls of a ring a common, deterministic parameter origin.
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

/*!
 * One node of the zigzag. The node position together with the foot of its perpendicular on the
 * medial axis spans the node's vertical reference plane (in 2D: the ray from the axis foot through
 * the node), which is used to map the node onto the next layer. When the cross-section drifts or
 * rotates from layer to layer, the plane is re-derived each layer from the node's new position and
 * the local medial axis, so it keeps following the walls; on constant cross-sections the re-derived
 * plane is identical every layer, so the nodes stay on fixed vertical planes.
 */
struct TrackedNode
{
    Point2LL axis_foot_; //!< The foot of the perpendicular from the node to the medial axis.
    Point2LL position_; //!< The node position on the wall.
};

/*!
 * The zigzag structure of one connected part: its nodes in zigzag order (alternating between the
 * two walls), plus the part outline used to find the corresponding part on the next layer.
 */
struct TrackedPart
{
    bool closed_ = false; //!< Ring parts close their zigzag into a loop.
    Shape shape_; //!< The outline of the part on the last processed layer, used to find the corresponding part on the next layer by overlap.
    std::vector<TrackedNode> nodes_;
};

/*!
 * Remove protruding details (e.g. rows of sawteeth) from a part's walls with a morphological opening,
 * so that nodes are placed on (or mapped to) the main wall envelope and never wander into a
 * protrusion. The opening is only used when it doesn't shatter or shrink the part, which happens when
 * the part itself is barely wider than the opening; the generated lines are clipped against the real
 * region at the end either way, so they always stay inside the solid region.
 */
SingleShape openedPart(const SingleShape& part, const coord_t opening_radius, const Simplify& simplifier)
{
    const std::vector<SingleShape> opened_parts = part.offset(-opening_radius).offset(opening_radius).splitIntoParts();
    if (opened_parts.size() == 1 && opened_parts.front().area() >= 0.7 * part.area())
    {
        SingleShape opened{ simplifier.polygon(opened_parts.front()) };
        if (! opened.empty() && opened.outerPolygon().size() >= 3)
        {
            return opened;
        }
    }
    return part;
}

/*!
 * Build the reference nodes for a ring-like part (a part with a dominant hole, e.g. an annulus): the
 * two walls are the outer wall and the hole wall. Both walls get the same number of node periods,
 * distributed uniformly along each wall's own arc length, half a period out of phase in the
 * normalized parameter, connected alternately into a closed loop. The medial axis of a ring is the
 * midway curve between the two walls.
 */
bool buildRingReference(const SingleShape& part, const coord_t period, TrackedPart& reference)
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

    // The number of node periods follows from the medial axis length (approximated by the average of
    // the two walls) and the reference spacing, rounded up so the spacing never exceeds the reference.
    const double axis_length = (outer_wall.total_ + inner_wall.total_) / 2.0;
    const size_t num_periods = std::max<size_t>(3, static_cast<size_t>(std::ceil(axis_length / period)));

    // Give both walls a common parameter origin (on the +X ray from the hole center) so that the phase
    // relation between the walls is well-defined and deterministic.
    const Point2LL origin = AABB(part[hole_idx]).getMiddle();
    const double anchor_outer = rayCrossingParameter(outer_wall, origin);
    const double anchor_inner = rayCrossingParameter(inner_wall, origin);

    // The ring's medial axis: the curve midway between the two walls, sampled densely.
    WallChain axis;
    axis.closed_ = true;
    const size_t axis_samples = std::max<size_t>(64, 8 * num_periods);
    for (size_t sample_idx = 0; sample_idx < axis_samples; sample_idx++)
    {
        const double fraction = double(sample_idx) / double(axis_samples);
        const Point2LL on_outer = outer_wall.pointAt(anchor_outer + fraction * outer_wall.total_);
        const Point2LL on_inner = inner_wall.pointAt(anchor_inner + fraction * inner_wall.total_);
        axis.points_.emplace_back((on_outer.X + on_inner.X) / 2, (on_outer.Y + on_inner.Y) / 2);
    }
    axis.finish();

    const auto add_node = [&reference, &axis](const Point2LL& position) -> bool
    {
        const Point2LL foot = axis.pointAt(nearestArcPosition(axis, position));
        if (foot == position)
        {
            return false; // Degenerate: no perpendicular direction, so no reference plane. Skip this node.
        }
        reference.nodes_.push_back({ foot, position });
        return true;
    };

    for (size_t node_idx = 0; node_idx < num_periods; node_idx++)
    {
        add_node(outer_wall.pointAt(anchor_outer + double(node_idx) / num_periods * outer_wall.total_));
        add_node(inner_wall.pointAt(anchor_inner + (node_idx + 0.5) / num_periods * inner_wall.total_));
    }
    if (reference.nodes_.size() < 3)
    {
        return false;
    }
    reference.closed_ = true;
    return true;
}

/*!
 * Build the reference nodes for an elongated part without a dominant hole (slab, L-shape, C-shape...):
 * compute the medial axis and split the boundary at the two axis ends into the two side walls. The
 * common period count N = ceil(axis length / reference spacing) is shared by both walls; each wall
 * distributes its nodes uniformly along its own arc length (the longer wall thus gets a larger actual
 * spacing), with the two walls half a period out of phase in the normalized arc-length parameter:
 * the longer wall's nodes lie at i/N (starting right at the axis start), the shorter wall's nodes at
 * (i+1/2)/N. The nodes are connected alternately into one open polyline.
 */
bool buildRibbonReference(const SingleShape& part, const coord_t period, TrackedPart& reference)
{
    std::optional<OpenPolyline> axis = computeMedialAxisPath(part);
    if (! axis.has_value())
    {
        return false;
    }
    // Deterministic axis orientation (which end is the start doesn't matter, but it must be stable).
    if (std::make_pair(axis->back().X, axis->back().Y) < std::make_pair(axis->front().X, axis->front().Y))
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
    const size_t num_periods = std::max<size_t>(1, static_cast<size_t>(std::ceil(axis_length / period)));

    // Split the boundary into the two side walls where the extensions of the medial axis (following
    // its end tangents) hit the wall. This is well-defined even when the axis end is equidistant to
    // several walls (e.g. in a rectangle, where the nearest wall point would be ambiguous).
    const WallChain wall = closedWall(part.outerPolygon());
    if (wall.total_ <= 0.0)
    {
        return false;
    }
    const auto split_position = [&wall, &axis_chain](const Point2LL& end, const Point2LL& inward) -> double
    {
        const std::optional<double> extended = rayArcPosition(wall, end, double(end.X - inward.X), double(end.Y - inward.Y));
        return extended.has_value() ? *extended : nearestArcPosition(wall, end);
    };
    const double split_front = split_position(axis_chain.points_.front(), axis_chain.points_[1]);
    const double split_back = split_position(axis_chain.points_.back(), axis_chain.points_[axis_chain.points_.size() - 2]);
    WallChain side_a = subChain(wall, split_front, split_back);
    WallChain side_b = subChain(wall, split_back, split_front);
    // side_b runs from the axis end back to the axis start; reverse it so both sides are parameterized
    // from the axis start towards the axis end.
    std::reverse(side_b.points_.begin(), side_b.points_.end());
    side_b.finish();
    if (side_a.total_ <= 0.0 || side_b.total_ <= 0.0)
    {
        return false;
    }
    const WallChain& long_wall = (side_a.total_ >= side_b.total_) ? side_a : side_b;
    const WallChain& short_wall = (side_a.total_ >= side_b.total_) ? side_b : side_a;

    const auto add_node = [&reference, &axis_chain](const WallChain& on_wall, const double u) -> bool
    {
        const Point2LL position = on_wall.pointAt(u * on_wall.total_);
        const Point2LL foot = axis_chain.pointAt(nearestArcPosition(axis_chain, position));
        if (foot == position)
        {
            return false; // Degenerate: no perpendicular direction, so no reference plane. Skip this node.
        }
        reference.nodes_.push_back({ foot, position });
        return true;
    };

    // Alternate between the walls: the longer wall at i/N (starting at the axis start), the shorter
    // wall at (i+1/2)/N, so the two walls are half a period out of phase in the normalized parameter.
    for (size_t node_idx = 0; node_idx <= num_periods; node_idx++)
    {
        add_node(long_wall, double(node_idx) / num_periods);
        if (node_idx < num_periods)
        {
            add_node(short_wall, (node_idx + 0.5) / num_periods);
        }
    }
    if (reference.nodes_.size() < 2)
    {
        return false;
    }
    reference.closed_ = false;
    return true;
}

/*!
 * Build the reference zigzag structures for all parts of the reference cross-section.
 */
std::vector<TrackedPart> buildTrackedParts(const Shape& region, const coord_t line_distance, const Simplify& simplifier)
{
    // Distance between two nodes on the same wall, measured along the wall; consecutive zigzag nodes
    // (which alternate between the two walls) are half of this period apart.
    const coord_t period = 2 * line_distance;
    const coord_t opening_radius = line_distance / 2;

    std::vector<TrackedPart> references;
    for (const SingleShape& raw_part : region.splitIntoParts())
    {
        SingleShape part{ simplifier.polygon(raw_part) };
        if (part.empty() || part.outerPolygon().size() < 3)
        {
            continue;
        }
        part = openedPart(part, opening_radius, simplifier);

        TrackedPart reference;
        bool built = false;
        if (part.size() > 1)
        {
            built = buildRingReference(part, period, reference);
        }
        if (! built)
        {
            built = buildRibbonReference(part, period, reference);
        }
        if (built)
        {
            reference.shape_ = part;
            references.push_back(std::move(reference));
        }
        // Parts without a usable reference structure get the cutting-plane fallback on every layer.
    }
    return references;
}

/*!
 * All crossings of a node's reference plane (the line through the axis foot and the node position)
 * with the walls of a part. The crossing parameter t is 0 at the axis foot and 1 at the node
 * position, so t > 0 is "the node's side of the axis".
 */
struct PlaneCrossing
{
    double t_;
    Point2LL point_;
};

void collectPlaneCrossings(const SingleShape& part, const Point2LL& origin, const double direction_x, const double direction_y, std::vector<PlaneCrossing>& crossings)
{
    const auto cross = [](const double ax, const double ay, const double bx, const double by)
    {
        return ax * by - ay * bx;
    };
    for (const Polygon& polygon : part)
    {
        const size_t point_count = polygon.size();
        for (size_t point_idx = 0; point_idx < point_count; point_idx++)
        {
            const Point2LL& p0 = polygon[point_idx];
            const Point2LL& p1 = polygon[(point_idx + 1) % point_count];
            const double denominator = cross(direction_x, direction_y, p1.X - p0.X, p1.Y - p0.Y);
            if (std::abs(denominator) < 1e-9)
            {
                continue;
            }
            const double t = cross(p0.X - origin.X, p0.Y - origin.Y, p1.X - p0.X, p1.Y - p0.Y) / denominator;
            const double u = cross(p0.X - origin.X, p0.Y - origin.Y, direction_x, direction_y) / denominator;
            if (u < 0.0 || u >= 1.0)
            {
                continue;
            }
            crossings.push_back({ t, Point2LL(std::llround(origin.X + t * direction_x), std::llround(origin.Y + t * direction_y)) });
        }
    }
}

/*!
 * Map the nodes of the tracked parts onto the walls of this layer and update the tracked state, so
 * that the next layer can be mapped from this one. For every node, its reference plane (spanned by
 * the vertical direction and the perpendicular from the node to the medial axis, as recorded in the
 * tracked state) is intersected with the walls of the corresponding part; only crossings on the
 * node's side of the axis count, and the crossing closest to the node's previous position is taken.
 * Nodes whose plane has no such crossing are invalid on this layer: they are skipped (never clamped
 * to a wall end or placed outside the walls) and keep their previous plane, so they can become valid
 * again on a later layer.
 *
 * After mapping, each node's reference plane is re-derived from its new position: the chord of the
 * plane through the part gives a fresh axis point (the chord midpoint), and the perpendicular to the
 * local axis direction (estimated from the neighboring nodes' chord midpoints) gives the plane
 * orientation. On constant cross-sections this reproduces exactly the same plane, keeping the nodes
 * on fixed vertical planes; on drifting or twisting cross-sections the planes follow the geometry.
 *
 * Layer parts without a corresponding tracked part get the cutting-plane fallback. The lines are
 * clipped against the real region, so no lines cross holes; a line leaving the solid region and
 * re-entering it continues as a new polyline.
 */
OpenLinesSet propagateLayer(
    std::vector<TrackedPart>& tracks,
    const Shape& region,
    const coord_t line_distance,
    const PointMatrix& rotation_matrix,
    const coord_t plane_shift,
    const Simplify& simplifier)
{
    const coord_t opening_radius = line_distance / 2;

    // Match every part of this layer to the tracked part it overlaps most with. Overlap (rather than
    // a point-containment test) keeps the correspondence intact when the cross-section drifts
    // sideways from layer to layer.
    std::vector<SingleShape> parts;
    for (const SingleShape& raw_part : region.splitIntoParts())
    {
        SingleShape part{ simplifier.polygon(raw_part) };
        if (! part.empty() && part.outerPolygon().size() >= 3)
        {
            parts.push_back(std::move(part));
        }
    }
    std::vector<std::vector<size_t>> track_parts(tracks.size());
    std::vector<size_t> fallback_parts;
    for (size_t part_idx = 0; part_idx < parts.size(); part_idx++)
    {
        size_t best_track = std::numeric_limits<size_t>::max();
        double best_overlap = 0.0;
        for (size_t track_idx = 0; track_idx < tracks.size(); track_idx++)
        {
            const double overlap = parts[part_idx].intersection(tracks[track_idx].shape_).area();
            if (overlap > best_overlap)
            {
                best_overlap = overlap;
                best_track = track_idx;
            }
        }
        if (best_track == std::numeric_limits<size_t>::max())
        {
            fallback_parts.push_back(part_idx);
        }
        else
        {
            track_parts[best_track].push_back(part_idx);
        }
    }

    OpenLinesSet zigzag_lines;
    for (size_t track_idx = 0; track_idx < tracks.size(); track_idx++)
    {
        TrackedPart& track = tracks[track_idx];
        if (track_parts[track_idx].empty())
        {
            continue; // The part vanished on this layer; the nodes keep their planes in case it reappears.
        }

        // Node placement uses the morphologically opened walls, so nodes stay on the main wall
        // envelope instead of being captured by protruding details (e.g. rows of sawteeth).
        std::vector<SingleShape> walls;
        Shape combined_shape;
        for (const size_t part_idx : track_parts[track_idx])
        {
            walls.push_back(openedPart(parts[part_idx], opening_radius, simplifier));
            combined_shape.push_back(parts[part_idx]);
        }

        // First pass: intersect each node's plane with the walls. The node moves to the crossing on
        // its side of the axis that is closest to its previous position; the chord midpoint (between
        // this crossing and the one across the part) samples the local medial axis.
        struct MappedNode
        {
            bool valid_ = false;
            Point2LL position_;
            std::optional<Point2LL> chord_midpoint_;
        };
        std::vector<MappedNode> mapped(track.nodes_.size());
        std::vector<PlaneCrossing> crossings;
        for (size_t node_idx = 0; node_idx < track.nodes_.size(); node_idx++)
        {
            const TrackedNode& node = track.nodes_[node_idx];
            const double direction_x = double(node.position_.X - node.axis_foot_.X);
            const double direction_y = double(node.position_.Y - node.axis_foot_.Y);
            if (direction_x == 0.0 && direction_y == 0.0)
            {
                continue;
            }
            crossings.clear();
            for (const SingleShape& wall_part : walls)
            {
                collectPlaneCrossings(wall_part, node.axis_foot_, direction_x, direction_y, crossings);
            }

            // The node's new position: the crossing on the node's side (t > 0), closest to the node's
            // previous position (t == 1). No crossing there means the node is invalid on this layer.
            double t_node = 0.0;
            double best_error = std::numeric_limits<double>::max();
            for (const PlaneCrossing& crossing : crossings)
            {
                const double error = std::abs(crossing.t_ - 1.0);
                if (crossing.t_ > 1e-6 && error < best_error)
                {
                    best_error = error;
                    t_node = crossing.t_;
                    mapped[node_idx].position_ = crossing.point_;
                    mapped[node_idx].valid_ = true;
                }
            }
            if (! mapped[node_idx].valid_)
            {
                continue;
            }

            // The crossing across the part (the other end of the plane's chord through the solid
            // material): the nearest crossing before the node whose chord midpoint lies inside.
            double t_partner = std::numeric_limits<double>::lowest();
            for (const PlaneCrossing& crossing : crossings)
            {
                if (crossing.t_ < t_node - 1e-6 && crossing.t_ > t_partner)
                {
                    t_partner = crossing.t_;
                }
            }
            if (t_partner > std::numeric_limits<double>::lowest())
            {
                const double t_mid = (t_node + t_partner) / 2.0;
                const Point2LL midpoint(std::llround(node.axis_foot_.X + t_mid * direction_x), std::llround(node.axis_foot_.Y + t_mid * direction_y));
                for (const SingleShape& wall_part : walls)
                {
                    if (wall_part.inside(midpoint, true))
                    {
                        mapped[node_idx].chord_midpoint_ = midpoint;
                        break;
                    }
                }
            }
        }

        // Second pass: re-derive each node's plane from its new position. The local axis direction is
        // estimated from the neighboring chord midpoints; the new plane is perpendicular to it.
        std::vector<size_t> midpoint_nodes;
        for (size_t node_idx = 0; node_idx < track.nodes_.size(); node_idx++)
        {
            if (! track.closed_ && (node_idx == 0 || node_idx + 1 == track.nodes_.size()))
            {
                continue; // End node chords run along the axis instead of across it; not usable as axis samples.
            }
            if (mapped[node_idx].valid_ && mapped[node_idx].chord_midpoint_.has_value())
            {
                midpoint_nodes.push_back(node_idx);
            }
        }
        for (size_t node_idx = 0; node_idx < track.nodes_.size(); node_idx++)
        {
            if (! mapped[node_idx].valid_)
            {
                continue; // Invalid on this layer: the node keeps its previous plane.
            }
            TrackedNode& node = track.nodes_[node_idx];
            // The end nodes of an open chain sit at the axis ends, where their planes contain the
            // axis direction rather than being perpendicular to it, so the perpendicular re-derivation
            // below doesn't apply to them: their planes only translate along with the node.
            const bool is_end_node = ! track.closed_ && (node_idx == 0 || node_idx + 1 == track.nodes_.size());
            if (is_end_node || ! mapped[node_idx].chord_midpoint_.has_value())
            {
                // Translate the plane with the node, keeping its orientation.
                node.axis_foot_ += mapped[node_idx].position_ - node.position_;
                node.position_ = mapped[node_idx].position_;
                continue;
            }
            const Point2LL midpoint = *mapped[node_idx].chord_midpoint_;

            // The neighboring chord midpoints around this node (wrapping around for rings).
            const auto neighbor_it = std::lower_bound(midpoint_nodes.begin(), midpoint_nodes.end(), node_idx);
            const size_t rank = std::distance(midpoint_nodes.begin(), neighbor_it);
            const size_t count = midpoint_nodes.size();
            std::optional<Point2LL> before;
            std::optional<Point2LL> after;
            if (count >= 2)
            {
                if (track.closed_)
                {
                    before = *mapped[midpoint_nodes[(rank + count - 1) % count]].chord_midpoint_;
                    after = *mapped[midpoint_nodes[(rank + 1) % count]].chord_midpoint_;
                }
                else
                {
                    before = *mapped[midpoint_nodes[rank == 0 ? 0 : rank - 1]].chord_midpoint_;
                    after = *mapped[midpoint_nodes[std::min(rank + 1, count - 1)]].chord_midpoint_;
                }
            }
            const double chord_length = distanceF(mapped[node_idx].position_, midpoint);
            double axis_dx = before.has_value() ? double(after->X - before->X) : 0.0;
            double axis_dy = before.has_value() ? double(after->Y - before->Y) : 0.0;
            const double axis_length = std::hypot(axis_dx, axis_dy);
            bool rotated = false;
            if (axis_length > 1.0 && chord_length > 10.0)
            {
                // The new plane: through the node, perpendicular to the local axis direction, with the
                // axis foot at the chord's distance so the crossing search stays calibrated (t == 1 at
                // the node).
                double normal_x = -axis_dy / axis_length;
                double normal_y = axis_dx / axis_length;
                if (normal_x * (mapped[node_idx].position_.X - midpoint.X) + normal_y * (mapped[node_idx].position_.Y - midpoint.Y) < 0.0)
                {
                    normal_x = -normal_x;
                    normal_y = -normal_y;
                }
                // Sanity-cap the rotation: consecutive layers only twist the cross-section a little, so
                // a wildly different plane orientation indicates a degenerate tangent estimate (e.g. from
                // a chord caught in a boundary detail) and is ignored, keeping the previous orientation.
                const double old_length = distanceF(node.position_, node.axis_foot_);
                const double along_old = old_length <= 0.0 ? 1.0 : (normal_x * (node.position_.X - node.axis_foot_.X) + normal_y * (node.position_.Y - node.axis_foot_.Y)) / old_length;
                if (along_old > std::numbers::sqrt2 / 2.0) // Less than 45 degrees away from the previous plane.
                {
                    node.axis_foot_ = mapped[node_idx].position_ - Point2LL(std::llround(normal_x * chord_length), std::llround(normal_y * chord_length));
                    rotated = true;
                }
            }
            if (! rotated)
            {
                node.axis_foot_ += mapped[node_idx].position_ - node.position_;
            }
            node.position_ = mapped[node_idx].position_;
        }
        track.shape_ = combined_shape;

        // The zigzag polyline of this layer: the valid nodes in zigzag order.
        OpenPolyline zigzag;
        for (const MappedNode& node : mapped)
        {
            if (node.valid_)
            {
                zigzag.push_back(node.position_);
            }
        }
        if (track.closed_ && zigzag.size() >= 3)
        {
            zigzag.push_back(zigzag.front()); // Close the ring loop.
        }
        if (zigzag.size() >= 2)
        {
            zigzag_lines.push_back(std::move(zigzag));
        }
        else
        {
            for (const size_t part_idx : track_parts[track_idx])
            {
                generatePlaneZigzag(parts[part_idx].outerPolygon(), rotation_matrix, plane_shift, line_distance, zigzag_lines);
            }
        }
    }
    for (const size_t part_idx : fallback_parts)
    {
        generatePlaneZigzag(parts[part_idx].outerPolygon(), rotation_matrix, plane_shift, line_distance, zigzag_lines);
    }

    // Clip against the actual region (including its holes and details): no lines inside holes; segments
    // crossing a hole or leaving the part through a concavity are cut at the boundary and continue as a
    // new polyline where they re-enter the solid region.
    return region.intersection(zigzag_lines);
}

} // namespace

void generateMedialZigzagLines(const Shape& region, const coord_t line_distance, const double plane_angle, coord_t plane_shift, OpenLinesSet& result)
{
    if (line_distance <= 0 || region.empty())
    {
        return;
    }
    // Simplify the walls before analyzing them; the voronoi diagram and the node placement don't need
    // (and shouldn't be disturbed by) micron-sized boundary details.
    const Simplify simplifier(100, 25, 0);
    const PointMatrix rotation_matrix(plane_angle);
    plane_shift = ((plane_shift % line_distance) + line_distance) % line_distance;

    // The stand-alone region acts as its own reference cross-section; mapping it onto itself
    // reproduces the reference nodes.
    std::vector<TrackedPart> tracks = buildTrackedParts(region, line_distance, simplifier);
    result.push_back(propagateLayer(tracks, region, line_distance, rotation_matrix, plane_shift, simplifier));
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
    if (line_distance > 0)
    {
        plane_shift = ((plane_shift % line_distance) + line_distance) % line_distance;
    }

    lines_per_layer_.resize(mesh.layers.size());
    if (line_distance <= 0 || mesh.layers.empty())
    {
        return;
    }
    const Simplify simplifier(100, 25, 0);
    const PointMatrix rotation_matrix(plane_angle);

    // The infill regions of all layers, and the reference cross-section: the layer with the largest
    // infill area. The reference nodes generated there define the reference planes through which the
    // nodes are mapped from layer to layer, upwards and downwards from the reference cross-section.
    std::vector<Shape> regions(mesh.layers.size());
    size_t reference_layer = 0;
    double largest_area = -1.0;
    for (size_t layer_nr = 0; layer_nr < mesh.layers.size(); layer_nr++)
    {
        for (const SliceLayerPart& part : mesh.layers[layer_nr].parts)
        {
            regions[layer_nr].push_back(part.getOwnInfillArea().offset(region_offset));
        }
        regions[layer_nr] = regions[layer_nr].unionPolygons();
        const double area = regions[layer_nr].area();
        if (area > largest_area)
        {
            largest_area = area;
            reference_layer = layer_nr;
        }
    }
    if (largest_area <= 0.0)
    {
        return;
    }
    std::vector<TrackedPart> tracks = buildTrackedParts(regions[reference_layer], line_distance, simplifier);
    lines_per_layer_[reference_layer] = propagateLayer(tracks, regions[reference_layer], line_distance, rotation_matrix, plane_shift, simplifier);

    std::vector<TrackedPart> tracks_down = tracks; // Both directions start from the reference cross-section's state.
    for (size_t layer_nr = reference_layer + 1; layer_nr < mesh.layers.size(); layer_nr++)
    {
        if (! regions[layer_nr].empty())
        {
            lines_per_layer_[layer_nr] = propagateLayer(tracks, regions[layer_nr], line_distance, rotation_matrix, plane_shift, simplifier);
        }
    }
    for (size_t layer_nr = reference_layer; layer_nr-- > 0;)
    {
        if (! regions[layer_nr].empty())
        {
            lines_per_layer_[layer_nr] = propagateLayer(tracks_down, regions[layer_nr], line_distance, rotation_matrix, plane_shift, simplifier);
        }
    }
}

const OpenLinesSet& MedialZigzagGenerator::getLinesForLayer(const size_t layer_nr) const
{
    static const OpenLinesSet empty;
    return layer_nr < lines_per_layer_.size() ? lines_per_layer_[layer_nr] : empty;
}

} // namespace cura
