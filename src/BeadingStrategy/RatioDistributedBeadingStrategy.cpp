//Copyright (c) 2021 Ultimaker B.V.
//CuraEngine is released under the terms of the AGPLv3 or higher.

#include "RatioDistributedBeadingStrategy.h"

#include <algorithm>
#include <numeric>

namespace cura
{
RatioDistributedBeadingStrategy::RatioDistributedBeadingStrategy(const std::vector<coord_t>& optimal_width_values,
                                                                 const coord_t minimum_line_width,
                                                                 const coord_t default_transition_length,
                                                                 const AngleRadians transitioning_angle,
                                                                 const Ratio wall_split_middle_threshold,
                                                                 const Ratio wall_add_middle_threshold,
                                                                 const int distribution_radius)
    : BeadingStrategy(optimal_width_values[0], default_transition_length, transitioning_angle)
    , optimal_width_values(optimal_width_values)
    , minimum_line_width(minimum_line_width)
    , wall_split_middle_threshold(wall_split_middle_threshold)
    , wall_add_middle_threshold(wall_add_middle_threshold)
{
    if(distribution_radius >= 2)
    {
        one_over_distribution_radius_squared = 1.0f / (distribution_radius - 1) * 1.0f / (distribution_radius - 1);
    }
    else
    {
        one_over_distribution_radius_squared = 1.0f / 1 * 1.0f / 1;
    }

    name = "RatioDistributedBeadingStrategy";
}

coord_t RatioDistributedBeadingStrategy::getOptimalThickness(coord_t bead_count) const
{
    // return optimal_width * bead_count;
    std::vector<coord_t> full_beads_width = getFixedOptimalWidthValues(bead_count);

    // if (optimal_width_values.size() == 2)
    // {
    //     logDebug("Optimal thickness for %d = %d.\n", bead_count, std::accumulate(full_beads_width.cbegin(), full_beads_width.cend(), static_cast<coord_t>(0)));
    // }
    return std::accumulate(full_beads_width.cbegin(), full_beads_width.cend(), static_cast<coord_t>(0));
}

coord_t RatioDistributedBeadingStrategy::getTransitionThickness(coord_t lower_bead_count) const
{
    // return lower_bead_count * optimal_width + optimal_width * (lower_bead_count % 2 == 1 ? wall_split_middle_threshold : wall_add_middle_threshold);
    std::vector<coord_t> full_beads_width = getFixedOptimalWidthValues(lower_bead_count + 1);

    const coord_t transition_thickness = getOptimalThickness(lower_bead_count) + full_beads_width[lower_bead_count / 2] * (lower_bead_count % 2 == 1 ? wall_split_middle_threshold : wall_add_middle_threshold);

    // if (optimal_width_values.size() == 2)
    // {
    //     logDebug("Optimal transition thickness for %d = %d.\n", lower_bead_count, transition_thickness);
    // }
    return transition_thickness;
}

coord_t RatioDistributedBeadingStrategy::getOptimalBeadCount(coord_t thickness) const
{
    // return (thickness + optimal_width / 2) / optimal_width;
    const coord_t max_width = std::accumulate(optimal_width_values.cbegin(), optimal_width_values.cend(), static_cast<coord_t>(0));
    if (thickness >= max_width)
    {
        const coord_t width_diff = thickness - max_width;
        const coord_t optimal_remaining_width = optimal_width_values[optimal_width_values.size() / 2];
        // if (optimal_width_values.size() == 2)
        // {
        //     logDebug("Optimal bead count for %d = %d.\n", thickness, optimal_width_values.size() + ((width_diff + optimal_remaining_width / 2) / optimal_remaining_width));
        // }
        return optimal_width_values.size() + ((width_diff + optimal_remaining_width / 2) / optimal_remaining_width);
    }
    else
    {
        coord_t current_thickness = 0;
        for (size_t i = 0; i < optimal_width_values.size() / 2; i++)
        {
            current_thickness += 2 * optimal_width_values[i];
            if (current_thickness >= thickness)
            {
                return (i + 1) * 2;
            }
        }
        // if (optimal_width_values.size() == 2)
        // {
        //     logDebug("Optimal bead count for %d = %d.\n", thickness, optimal_width_values.size());
        // }
        return optimal_width_values.size();
    }
}

std::vector<coord_t> RatioDistributedBeadingStrategy::getFixedOptimalWidthValues(coord_t bead_count) const
{
    coord_t beads_diff = bead_count - optimal_width_values.size();
    std::vector<coord_t> full_beads_widths = optimal_width_values;
    if (beads_diff > 0)
    {
        full_beads_widths.insert(full_beads_widths.begin() + full_beads_widths.size() / 2,
                                 beads_diff, optimal_width_values[full_beads_widths.size() / 2]);
    }
    else
    {
        beads_diff = -beads_diff;
        auto left_bound = full_beads_widths.begin() + full_beads_widths.size() / 2 - beads_diff / 2;
        auto right_bound = full_beads_widths.begin() + full_beads_widths.size() / 2 + beads_diff / 2 + ((beads_diff % 2 == 1) ? 1 : 0);
        full_beads_widths.erase(left_bound, right_bound);
    }
    // if (optimal_width_values.size() == 2)
    // {
    //     logDebug("Bead count after fixing: wanted = %d, actual %d.\n", bead_count, full_beads_widths.size());
    // }
    return full_beads_widths;
}

BeadingStrategy::Beading RatioDistributedBeadingStrategy::compute(coord_t thickness, coord_t bead_count, coord_t distance_to_source) const
{
    Beading ret;

    ret.total_thickness = thickness;
    std::vector<coord_t> full_beads_width = getFixedOptimalWidthValues(bead_count);
    if (bead_count > 2)
    {
        const coord_t to_be_divided = thickness - getOptimalThickness(bead_count);// getOptimalThickness(bead_count);
        // const coord_t to_be_divided = thickness - bead_count * optimal_width;// getOptimalThickness(bead_count);
        const float middle = static_cast<float>(bead_count - 1) / 2;

        const auto getWeight = [middle, this](coord_t bead_idx)
        {
            const float dev_from_middle = bead_idx - middle;
            return std::max(0.0f, 1.0f - one_over_distribution_radius_squared * dev_from_middle * dev_from_middle);
        };

        std::vector<float> weights;
        weights.resize(bead_count);
        for (coord_t bead_idx = 0; bead_idx < bead_count; bead_idx++)
        {
            weights[bead_idx] = getWeight(bead_idx);
        }

        const float total_weight = std::accumulate(weights.cbegin(), weights.cend(), 0.f);
        for (coord_t bead_idx = 0; bead_idx < bead_count; bead_idx++)
        {
            const float weight_fraction = weights[bead_idx] / total_weight;
            const coord_t splitup_left_over_weight = to_be_divided * weight_fraction;
            // const coord_t width = optimal_width + splitup_left_over_weight;
            const coord_t width = full_beads_width[bead_idx] + splitup_left_over_weight;
            // if (optimal_width_values.size() == 2)
            // {
            //     logDebug("Bead width %d = %d.\n", bead_idx, width);
            // }
            if (bead_idx == 0)
            {
                ret.toolpath_locations.emplace_back(width / 2);
            }
            else
            {
                ret.toolpath_locations.emplace_back(ret.toolpath_locations.back() + (ret.bead_widths.back() + width) / 2);
            }
            ret.bead_widths.emplace_back(width);
        }
        ret.left_over = 0;
    }
    else if (bead_count == 2)
    {
        const coord_t outer_width = thickness / 2;
        ret.bead_widths.emplace_back(outer_width);
        ret.bead_widths.emplace_back(outer_width);
        ret.toolpath_locations.emplace_back(outer_width / 2);
        ret.toolpath_locations.emplace_back(thickness - outer_width / 2);
        ret.left_over = 0;
    }
    else if (bead_count == 1)
    {
        const coord_t outer_width = thickness;
        ret.bead_widths.emplace_back(outer_width);
        ret.toolpath_locations.emplace_back(outer_width / 2);
        ret.left_over = 0;
    }
    else
    {
        ret.left_over = thickness;
    }

    return ret;
}
} // namespace cura
