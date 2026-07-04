// Copyright (c) 2026 UltiMaker
// CuraEngine is released under the terms of the AGPLv3 or higher.

#include "infill/TriangleWaveInfill.h"

#include <cmath>
#include <limits>
#include <numbers>

#include "geometry/OpenPolyline.h"
#include "geometry/PointMatrix.h"
#include "geometry/Polygon.h"
#include "geometry/Shape.h"
#include "utils/AABB.h"
#include "utils/linearAlg2D.h"

namespace cura
{

void TriangleWaveInfill::generateTotalTriangleWaveInfill(OpenLinesSet& result_lines, bool zig_zaggify, coord_t line_distance, const Shape& in_outline, const double fill_angle)
{
    if (line_distance <= 0 || in_outline.empty())
    {
        return;
    }

    // Work in a rotated coordinate frame so that the waves always run along the X axis.
    // The results are rotated back at the very end.
    const PointMatrix rotation_matrix(fill_angle);
    Shape outline = in_outline;
    outline.applyMatrix(rotation_matrix);

    const AABB aabb(outline);

    // The wave flanks rise/fall at 45 degrees, so the peak-to-peak amplitude equals half the period.
    // Rows are spaced line_distance * sqrt(2) apart vertically which makes the perpendicular
    // distance between the flanks of adjacent waves exactly line_distance, producing the same
    // material density as the 'lines' pattern.
    const coord_t pitch = line_distance * std::numbers::sqrt2;

    // Subdivide the flanks so that boundary crossings are detected on short segments,
    // similar to what the gyroid pattern does.
    int num_steps = 1;
    coord_t step = line_distance * 2;
    while (step > 500 && num_steps < 8)
    {
        num_steps *= 2;
        step = line_distance * 2 / num_steps;
    }
    const coord_t half_period = step * num_steps; // recalculate to avoid rounding errors
    const coord_t amplitude = half_period; // peak-to-peak, gives 45 degree flanks

    OpenLinesSet result;
    std::vector<Point2LL> chains[2]; // [start_points[], end_points[]]
    std::vector<unsigned> connected_to[2]; // [chain_indices[], chain_indices[]]
    std::vector<int> line_numbers; // which row a chain is part of

    unsigned num_rows = 0;
    for (coord_t y = (std::floor(static_cast<double>(aabb.min_.Y) / pitch) - 1) * pitch; y <= aabb.max_.Y + amplitude; y += pitch)
    {
        bool is_first_point = true;
        Point2LL last;
        bool last_inside = false;
        unsigned chain_end_index = 0;
        Point2LL chain_end[2];
        const coord_t x_min = (std::floor(static_cast<double>(aabb.min_.X) / (2 * half_period)) - 1) * 2 * half_period;
        unsigned sample_index = 0;
        for (coord_t x = x_min; x <= aabb.max_.X + 2 * half_period; x += step, ++sample_index)
        {
            // Triangle wave: rises from -amplitude/2 to +amplitude/2 over the first half period, then falls back.
            const unsigned phase = sample_index % (2 * num_steps);
            const coord_t y_offset = (phase <= static_cast<unsigned>(num_steps)) ? -amplitude / 2 + amplitude * static_cast<coord_t>(phase) / num_steps
                                                                                 : amplitude / 2 - amplitude * static_cast<coord_t>(phase - num_steps) / num_steps;
            const Point2LL current(x, y + y_offset);
            const bool current_inside = outline.inside(current, true);
            if (! is_first_point)
            {
                if (last_inside && current_inside)
                {
                    // segment doesn't hit the boundary, add it wholly
                    result.addSegment(last, current);
                }
                else if (last_inside != current_inside)
                {
                    // segment hits the boundary, add the part that's inside the boundary
                    OpenLinesSet line;
                    line.addSegment(last, current);
                    constexpr bool restitch = false; // only a single line doesn't need stitching
                    line = outline.intersection(line, restitch);
                    if (line.size() > 0)
                    {
                        // some of the segment is inside the boundary
                        result.addSegment(line[0][0], line[0][1]);
                        if (zig_zaggify)
                        {
                            chain_end[chain_end_index] = line[0][(line[0][0] != last && line[0][0] != current) ? 0 : 1];
                            if (++chain_end_index == 2)
                            {
                                chains[0].push_back(chain_end[0]);
                                chains[1].push_back(chain_end[1]);
                                chain_end_index = 0;
                                connected_to[0].push_back(std::numeric_limits<unsigned>::max());
                                connected_to[1].push_back(std::numeric_limits<unsigned>::max());
                                line_numbers.push_back(num_rows);
                            }
                        }
                    }
                    else
                    {
                        // none of the segment is inside the boundary so the point that's actually on the boundary
                        // is the chain end
                        if (zig_zaggify)
                        {
                            chain_end[chain_end_index] = (last_inside) ? last : current;
                            if (++chain_end_index == 2)
                            {
                                chains[0].push_back(chain_end[0]);
                                chains[1].push_back(chain_end[1]);
                                chain_end_index = 0;
                                connected_to[0].push_back(std::numeric_limits<unsigned>::max());
                                connected_to[1].push_back(std::numeric_limits<unsigned>::max());
                                line_numbers.push_back(num_rows);
                            }
                        }
                    }
                }
            }
            last = current;
            last_inside = current_inside;
            is_first_point = false;
        }
        ++num_rows;
    }

    if (zig_zaggify && chains[0].size() > 0)
    {
        // zig-zaggification consists of joining alternate chain ends to make a chain of chains
        // the basic algorithm is that we follow the infill area boundary and as we progress we are either drawing a connector or not
        // whenever we come across the end of a chain we toggle the connector drawing state
        // things are made more complicated by the fact that we want to avoid generating loops and so we need to keep track
        // of the identity of the first chain in a connected sequence

        int chain_ends_remaining = chains[0].size() * 2;

        for (const Polygon& outline_poly : outline)
        {
            std::vector<Point2LL> connector_points; // the points that make up a connector line

            // we need to remember the first chain processed and the path to it from the first outline point
            // so that later we can possibly connect to it from the last chain processed
            unsigned first_chain_chain_index = std::numeric_limits<unsigned>::max();
            std::vector<Point2LL> path_to_first_chain;

            bool drawing = false; // true when a connector line is being (potentially) created

            // keep track of the chain+point that a connector line started at
            unsigned connector_start_chain_index = std::numeric_limits<unsigned>::max();
            unsigned connector_start_point_index = std::numeric_limits<unsigned>::max();

            Point2LL cur_point; // current point of interest - either an outline point or a chain end

            // go round all of the region's outline and find the chain ends that meet it
            // quit the loop early if we have seen all the chain ends and are not currently drawing a connector
            for (unsigned outline_point_index = 0; (chain_ends_remaining > 0 || drawing) && outline_point_index < outline_poly.size(); ++outline_point_index)
            {
                Point2LL op0 = outline_poly[outline_point_index];
                Point2LL op1 = outline_poly[(outline_point_index + 1) % outline_poly.size()];
                std::vector<unsigned> points_on_outline_chain_index;
                std::vector<unsigned> points_on_outline_point_index;

                // collect the chain ends that meet this segment of the outline
                for (unsigned chain_index = 0; chain_index < chains[0].size(); ++chain_index)
                {
                    for (unsigned point_index = 0; point_index < 2; ++point_index)
                    {
                        // don't include chain ends that are close to the segment but are beyond the segment ends
                        short beyond = 0;
                        if (LinearAlg2D::getDist2FromLineSegment(op0, chains[point_index][chain_index], op1, &beyond) < 10 && ! beyond)
                        {
                            points_on_outline_point_index.push_back(point_index);
                            points_on_outline_chain_index.push_back(chain_index);
                        }
                    }
                }

                if (outline_point_index == 0 || vSize2(op0 - cur_point) > MM2INT(0.1))
                {
                    // this is either the first outline point or it is another outline point that is not too close to cur_point

                    if (first_chain_chain_index == std::numeric_limits<unsigned>::max())
                    {
                        // include the outline point in the path to the first chain
                        path_to_first_chain.push_back(op0);
                    }

                    cur_point = op0;
                    if (drawing)
                    {
                        // include the start point of this outline segment in the connector
                        connector_points.push_back(op0);
                    }
                }

                // iterate through each of the chain ends that meet the current outline segment
                while (points_on_outline_chain_index.size() > 0)
                {
                    // find the nearest chain end to the current point
                    unsigned nearest_point_index = 0;
                    double nearest_point_dist2 = std::numeric_limits<float>::infinity();
                    for (unsigned pi = 0; pi < points_on_outline_chain_index.size(); ++pi)
                    {
                        double dist2 = vSize2f(chains[points_on_outline_point_index[pi]][points_on_outline_chain_index[pi]] - cur_point);
                        if (dist2 < nearest_point_dist2)
                        {
                            nearest_point_dist2 = dist2;
                            nearest_point_index = pi;
                        }
                    }
                    const unsigned point_index = points_on_outline_point_index[nearest_point_index];
                    const unsigned chain_index = points_on_outline_chain_index[nearest_point_index];

                    // make the chain end the current point and add it to the connector line
                    cur_point = chains[point_index][chain_index];

                    if (drawing && connector_points.size() > 0 && vSize2(cur_point - connector_points.back()) < MM2INT(0.1))
                    {
                        // this chain end will be too close to the last connector point so throw away the last connector point
                        connector_points.pop_back();
                    }
                    connector_points.push_back(cur_point);

                    if (first_chain_chain_index == std::numeric_limits<unsigned>::max())
                    {
                        // this is the first chain to be processed, remember it
                        first_chain_chain_index = chain_index;
                        path_to_first_chain.push_back(cur_point);
                    }

                    if (drawing)
                    {
                        // add the connector line segments but only if
                        //  1 - the start/end points are not the opposite ends of the same chain
                        //  2 - the other end of the current chain is not connected to the chain the connector line is coming from

                        if (chain_index != connector_start_chain_index && connected_to[(point_index + 1) % 2][chain_index] != connector_start_chain_index)
                        {
                            result.push_back(OpenPolyline{ connector_points });
                            drawing = false;
                            connector_points.clear();
                            // remember the connection
                            connected_to[point_index][chain_index] = connector_start_chain_index;
                            connected_to[connector_start_point_index][connector_start_chain_index] = chain_index;
                        }
                        else
                        {
                            // start a new connector from the current location
                            connector_points.clear();
                            connector_points.push_back(cur_point);

                            // remember the chain+point that the connector started from
                            connector_start_chain_index = chain_index;
                            connector_start_point_index = point_index;
                        }
                    }
                    else
                    {
                        // we have just jumped a gap so now we want to start drawing again
                        drawing = true;

                        // if this connector is the first to be created or we are not connecting chains from the same row,
                        // remember the chain+point that this connector is starting from
                        if (connector_start_chain_index == std::numeric_limits<unsigned>::max() || line_numbers[chain_index] != line_numbers[connector_start_chain_index])
                        {
                            connector_start_chain_index = chain_index;
                            connector_start_point_index = point_index;
                        }
                    }

                    // done with this chain end
                    points_on_outline_chain_index.erase(points_on_outline_chain_index.begin() + nearest_point_index);
                    points_on_outline_point_index.erase(points_on_outline_point_index.begin() + nearest_point_index);

                    // decrement total amount of work to do
                    --chain_ends_remaining;
                }
            }

            // we have now visited all the points in the outline, if a connector was (potentially) being drawn
            // check whether the first chain is already connected to the last chain and, if not, draw the
            // connector between
            if (drawing && first_chain_chain_index != std::numeric_limits<unsigned>::max() && first_chain_chain_index != connector_start_chain_index
                && connected_to[0][first_chain_chain_index] != connector_start_chain_index && connected_to[1][first_chain_chain_index] != connector_start_chain_index)
            {
                // output the connector line segments from the last chain to the first point in the outline
                connector_points.push_back(outline_poly[0]);
                result.push_back(OpenPolyline{ connector_points });
                // output the connector line segments from the first point in the outline to the first chain
                result.push_back(OpenPolyline{ path_to_first_chain });
            }

            if (chain_ends_remaining < 1)
            {
                break;
            }
        }
    }

    // rotate the resulting paths back to the original coordinate frame
    const PointMatrix inverse_rotation = rotation_matrix.inverse();
    for (OpenPolyline& polyline : result)
    {
        polyline.applyMatrix(inverse_rotation);
    }

    result_lines = result;
}

} // namespace cura
