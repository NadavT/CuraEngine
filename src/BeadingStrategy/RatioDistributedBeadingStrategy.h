//Copyright (c) 2020 Ultimaker B.V.
//CuraEngine is released under the terms of the AGPLv3 or higher.

#ifndef RATIO_DISTRIBUTED_BEADING_STRATEGY_H
#define RATIO_DISTRIBUTED_BEADING_STRATEGY_H

#include "BeadingStrategy.h"

#include "../settings/types/Ratio.h" // For the wall transition threshold.

#include <vector>

namespace cura
{
/*!
     * A meta-beading-strategy that takes outer and inner wall widths into account.
     *
     * The outer wall will try to keep a constant width by only applying the beading strategy on the inner walls. This
     * ensures that this outer wall doesn't react to changes happening to inner walls. It will limit print artifacts on
     * the surface of the print. Although this strategy technically deviates from the original philosophy of the paper.
     * It will generally results in better prints because of a smoother motion and less variation in extrusion width in
     * the outer walls.
     *
     * If the thickness of the model is less then two times the optimal outer wall width and once the minimum inner wall
     * width it will keep the minimum inner wall at a minimum constant and vary the outer wall widths symmetrical. Until
     * The thickness of the model is that of at least twice the optimal outer wall width it will then use two
     * symmetrical outer walls only. Until it transitions into a single outer wall. These last scenario's are always
     * symmetrical in nature, disregarding the user specified strategy.
     */
class RatioDistributedBeadingStrategy : public BeadingStrategy
{
  public:
    RatioDistributedBeadingStrategy(const std::vector<coord_t>& optimal_width_values,
                                    const std::vector<Ratio>& optimal_width_ratios,
                                    const coord_t minimum_line_width,
                                    const coord_t maximum_line_width,
                                    const coord_t default_transition_length,
                                    const AngleRadians transitioning_angle,
                                    const Ratio wall_split_middle_threshold,
                                    const Ratio wall_add_middle_threshold,
                                    const int distribution_radius);

    virtual ~RatioDistributedBeadingStrategy() override = default;

    Beading compute(coord_t thickness, coord_t bead_count, coord_t distance_to_source) const override;

    coord_t getOptimalThickness(coord_t bead_count) const override;
    coord_t getTransitionThickness(coord_t lower_bead_count) const override;
    coord_t getOptimalBeadCount(coord_t thickness) const override;

  protected:
    std::vector<coord_t> getFixedOptimalWidthValues(coord_t bead_count) const;

  protected:
    std::vector<coord_t> optimal_width_values;
    std::vector<Ratio> optimal_width_ratios;
    const coord_t minimum_line_width;
    const coord_t maximum_line_width;

    // For uneven numbers of lines: Minimum factor of the optimal width for which the middle line will be split into two lines.
    Ratio wall_split_middle_threshold;

    // For even numbers of lines: Minimum factor of the optimal width for which a new middle line will be added between the two innermost lines.
    Ratio wall_add_middle_threshold;

    float one_over_distribution_radius_squared; // (1 / distribution_radius)^2
};

} // namespace cura
#endif // RATIO_DISTRIBUTED_BEADING_STRATEGY_H
