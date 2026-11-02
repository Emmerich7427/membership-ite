################################################################################
# Module: 01_interval_utils.R
# Purpose: Utilities for interval unions and prediction-set summaries
################################################################################

empty_interval_union <- function() {
    matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("lower", "upper")))
}


normalize_interval_union <- function(union, merge_touching = TRUE) {
    if (is.null(union) || length(union) == 0) {
        return(empty_interval_union())
    }

    union <- as.matrix(union)
    if (ncol(union) != 2) {
        if (length(union) == 2) {
            union <- matrix(as.numeric(union), nrow = 1)
        } else {
            stop("interval union must have exactly two columns.")
        }
    }
    if (nrow(union) == 0) {
        return(empty_interval_union())
    }

    storage.mode(union) <- "numeric"
    union <- union[!is.na(union[, 1]) & !is.na(union[, 2]), , drop = FALSE]
    if (nrow(union) == 0) {
        return(empty_interval_union())
    }

    lower <- union[, 1]
    upper <- union[, 2]
    union <- cbind(lower, upper)
    union <- union[upper >= lower, , drop = FALSE]
    if (nrow(union) == 0) {
        return(empty_interval_union())
    }

    union <- union[order(union[, 1], union[, 2]), , drop = FALSE]
    merged <- matrix(union[1, ], nrow = 1)
    for (k in seq_len(nrow(union))[-1]) {
        last <- nrow(merged)
        overlaps <- if (merge_touching) {
            union[k, 1] <= merged[last, 2]
        } else {
            union[k, 1] < merged[last, 2]
        }

        if (overlaps) {
            merged[last, 2] <- max(merged[last, 2], union[k, 2])
        } else {
            merged <- rbind(merged, union[k, ])
        }
    }

    colnames(merged) <- c("lower", "upper")
    merged
}


validate_linear_grid <- function(grid, values) {
    grid <- as.numeric(grid)
    values <- as.numeric(values)

    if (length(grid) != length(values)) {
        stop("grid and values must have the same length.")
    }
    if (length(grid) < 2) {
        stop("grid and values must contain at least two points.")
    }
    if (any(!is.finite(grid))) {
        stop("grid values must be finite.")
    }
    if (any(diff(grid) <= 0)) {
        stop("grid values must be strictly increasing.")
    }
    if (any(is.na(values))) {
        stop("function values must not contain NA or NaN.")
    }

    list(grid = grid, values = values)
}


piecewise_linear_eval <- function(grid,
                                  values,
                                  x,
                                  left_value = values[1],
                                  right_value = values[length(values)]) {
    checked <- validate_linear_grid(grid, values)
    grid <- checked$grid
    values <- checked$values
    x <- as.numeric(x)
    out <- rep(NA_real_, length(x))

    if (length(x) == 0) {
        return(out)
    }

    left_value <- as.numeric(left_value)[1]
    right_value <- as.numeric(right_value)[1]
    if (is.na(left_value) || is.na(right_value)) {
        stop("left_value and right_value must not be NA or NaN.")
    }

    valid_x <- !is.na(x)
    left <- valid_x & x < grid[1]
    right <- valid_x & x > grid[length(grid)]
    inside <- valid_x & !left & !right

    out[left] <- left_value
    out[right] <- right_value
    if (any(inside)) {
        out[inside] <- stats::approx(
            x = grid,
            y = values,
            xout = x[inside],
            rule = 1,
            ties = "ordered"
        )$y
    }

    out
}


linear_superlevel_intervals <- function(grid,
                                        values,
                                        threshold,
                                        left_value = -Inf,
                                        right_value = -Inf,
                                        tol = sqrt(.Machine$double.eps)) {
    checked <- validate_linear_grid(grid, values)
    grid <- checked$grid
    values <- checked$values
    threshold <- as.numeric(threshold)[1]
    left_value <- as.numeric(left_value)[1]
    right_value <- as.numeric(right_value)[1]

    if (is.na(threshold)) {
        stop("threshold must not be NA or NaN.")
    }
    if (is.na(left_value) || is.na(right_value)) {
        stop("left_value and right_value must not be NA or NaN.")
    }
    if (!is.finite(tol) || tol < 0) {
        stop("tol must be a nonnegative finite number.")
    }

    if (threshold == -Inf) {
        return(normalize_interval_union(matrix(c(-Inf, Inf), nrow = 1)))
    }

    if (threshold == Inf) {
        pieces <- matrix(numeric(0), ncol = 2)
        if (left_value == Inf) {
            pieces <- rbind(pieces, c(-Inf, grid[1]))
        }
        inf_nodes <- which(values == Inf)
        if (length(inf_nodes) > 0) {
            pieces <- rbind(pieces, cbind(grid[inf_nodes], grid[inf_nodes]))
        }
        if (length(values) > 1) {
            for (j in seq_len(length(values) - 1)) {
                if (values[j] == Inf && values[j + 1] == Inf) {
                    pieces <- rbind(pieces, c(grid[j], grid[j + 1]))
                }
            }
        }
        if (right_value == Inf) {
            pieces <- rbind(pieces, c(grid[length(grid)], Inf))
        }
        return(normalize_interval_union(pieces))
    }

    if (any(!is.finite(values))) {
        stop("function values must be finite when threshold is finite.")
    }

    pieces <- matrix(numeric(0), ncol = 2)
    pass <- values >= threshold - tol

    if (left_value >= threshold - tol) {
        if (pass[1]) {
            pieces <- rbind(pieces, c(-Inf, grid[1]))
        } else {
            pieces <- rbind(pieces, c(-Inf, grid[1]))
        }
    }

    for (j in seq_len(length(grid) - 1)) {
        x0 <- grid[j]
        x1 <- grid[j + 1]
        y0 <- values[j]
        y1 <- values[j + 1]
        p0 <- pass[j]
        p1 <- pass[j + 1]

        if (p0 && p1) {
            pieces <- rbind(pieces, c(x0, x1))
        } else if (p0 != p1) {
            if (abs(y1 - y0) <= tol) {
                crossing <- if (p0) x0 else x1
            } else {
                crossing <- x0 + (threshold - y0) / (y1 - y0) * (x1 - x0)
                crossing <- min(max(crossing, x0), x1)
            }

            if (p0) {
                pieces <- rbind(pieces, c(x0, crossing))
            } else {
                pieces <- rbind(pieces, c(crossing, x1))
            }
        }
    }

    if (right_value >= threshold - tol) {
        if (pass[length(pass)]) {
            pieces <- rbind(pieces, c(grid[length(grid)], Inf))
        } else {
            pieces <- rbind(pieces, c(grid[length(grid)], Inf))
        }
    }

    normalize_interval_union(pieces)
}


linear_min_over_union <- function(grid,
                                  values,
                                  interval_union,
                                  region = range(grid),
                                  outside_value = 0) {
    checked <- validate_linear_grid(grid, values)
    grid <- checked$grid
    values <- checked$values
    union <- normalize_interval_union(interval_union)

    if (nrow(union) == 0) {
        return(Inf)
    }

    region <- as.numeric(region)
    if (length(region) != 2 || any(is.na(region)) || any(!is.finite(region))) {
        stop("region must be a finite numeric vector of length two.")
    }
    region <- sort(region)
    if (region[1] >= region[2]) {
        stop("region must have positive width.")
    }

    outside_value <- as.numeric(outside_value)[1]
    if (is.na(outside_value)) {
        stop("outside_value must not be NA or NaN.")
    }

    candidates <- numeric(0)
    has_outside <- any(union[, 1] < region[1] | union[, 2] > region[2])
    if (has_outside) {
        candidates <- c(candidates, outside_value)
    }

    for (k in seq_len(nrow(union))) {
        lower <- max(union[k, 1], region[1])
        upper <- min(union[k, 2], region[2])
        if (lower <= upper) {
            locs <- c(lower, upper, grid[grid >= lower & grid <= upper])
            locs <- sort(unique(locs))
            candidates <- c(
                candidates,
                piecewise_linear_eval(
                    grid = grid,
                    values = values,
                    x = locs,
                    left_value = outside_value,
                    right_value = outside_value
                )
            )
        }
    }

    if (length(candidates) == 0) {
        return(outside_value)
    }

    min(candidates, na.rm = TRUE)
}


linear_quantile_over_union <- function(grid,
                                       values,
                                       interval_union,
                                       gamma,
                                       region = range(grid),
                                       outside_value = 0,
                                       tol = sqrt(.Machine$double.eps)) {
    if (length(gamma) != 1 || is.na(gamma) || gamma < 0 || gamma >= 1) {
        stop("gamma must be in [0, 1).")
    }
    if (!is.finite(tol) || tol < 0) {
        stop("tol must be a nonnegative finite number.")
    }
    if (gamma == 0) {
        return(linear_min_over_union(
            grid = grid,
            values = values,
            interval_union = interval_union,
            region = region,
            outside_value = outside_value
        ))
    }

    checked <- validate_linear_grid(grid, values)
    grid <- checked$grid
    values <- checked$values
    if (any(!is.finite(values))) {
        stop("function values must be finite for length-quantile calculation.")
    }

    region <- as.numeric(region)
    if (length(region) != 2 || any(is.na(region)) || any(!is.finite(region))) {
        stop("region must be a finite numeric vector of length two.")
    }
    region <- sort(region)
    if (region[1] >= region[2]) {
        stop("region must have positive width.")
    }

    outside_value <- as.numeric(outside_value)[1]
    if (is.na(outside_value) || !is.finite(outside_value)) {
        stop("outside_value must be finite and not NA or NaN.")
    }

    union <- normalize_interval_union(interval_union)
    if (nrow(union) == 0) {
        return(Inf)
    }
    if (any(!is.finite(union[, 1]) | !is.finite(union[, 2]))) {
        stop("Cannot compute a length quantile over an unbounded interval union.")
    }

    atoms <- data.frame(value = numeric(0), mass = numeric(0))
    ramps <- data.frame(low = numeric(0), high = numeric(0), density = numeric(0))

    add_atom <- function(value, mass) {
        if (mass > tol) {
            atoms <<- rbind(atoms, data.frame(value = value, mass = mass))
        }
    }

    add_linear_piece <- function(left, right) {
        if (right - left <= tol) {
            return(NULL)
        }

        endpoints <- piecewise_linear_eval(
            grid = grid,
            values = values,
            x = c(left, right),
            left_value = outside_value,
            right_value = outside_value
        )
        value_left <- endpoints[1]
        value_right <- endpoints[2]
        piece_length <- right - left

        if (abs(value_right - value_left) <= tol) {
            add_atom((value_left + value_right) / 2, piece_length)
        } else {
            low <- min(value_left, value_right)
            high <- max(value_left, value_right)
            ramps <<- rbind(
                ramps,
                data.frame(
                    low = low,
                    high = high,
                    density = piece_length / (high - low)
                )
            )
        }

        NULL
    }

    for (k in seq_len(nrow(union))) {
        lower <- union[k, 1]
        upper <- union[k, 2]

        if (lower < region[1]) {
            outside_upper <- min(upper, region[1])
            if (outside_upper > lower) {
                add_atom(outside_value, outside_upper - lower)
            }
        }
        if (upper > region[2]) {
            outside_lower <- max(lower, region[2])
            if (upper > outside_lower) {
                add_atom(outside_value, upper - outside_lower)
            }
        }

        inside_lower <- max(lower, region[1])
        inside_upper <- min(upper, region[2])
        if (inside_upper > inside_lower) {
            locs <- sort(unique(c(
                inside_lower,
                inside_upper,
                grid[grid > inside_lower & grid < inside_upper]
            )))
            for (j in seq_len(length(locs) - 1)) {
                add_linear_piece(locs[j], locs[j + 1])
            }
        }
    }

    total_mass <- sum(atoms$mass) + sum((ramps$high - ramps$low) * ramps$density)
    if (!is.finite(total_mass)) {
        stop("Cannot compute a length quantile with non-finite total length.")
    }
    if (total_mass <= tol) {
        return(linear_min_over_union(
            grid = grid,
            values = values,
            interval_union = union,
            region = region,
            outside_value = outside_value
        ))
    }

    breakpoints <- c(atoms$value, ramps$low, ramps$high)
    breakpoints <- sort(unique(breakpoints[is.finite(breakpoints)]))
    if (length(breakpoints) == 0) {
        stop("No finite function values available for length-quantile calculation.")
    }

    target <- gamma * total_mass
    cumulative <- 0

    for (idx in seq_along(breakpoints)) {
        point <- breakpoints[idx]
        if (nrow(atoms) > 0) {
            atom_mass <- sum(atoms$mass[atoms$value == point])
            cumulative <- cumulative + atom_mass
            if (cumulative + tol >= target) {
                return(point)
            }
        }

        if (idx < length(breakpoints)) {
            next_point <- breakpoints[idx + 1]
            if (next_point > point) {
                active <- if (nrow(ramps) == 0) {
                    rep(FALSE, 0)
                } else {
                    ramps$low < next_point - tol & ramps$high > point + tol
                }
                slope <- sum(ramps$density[active])
                if (slope > 0) {
                    interval_mass <- slope * (next_point - point)
                    if (cumulative + interval_mass + tol >= target) {
                        return(point + (target - cumulative) / slope)
                    }
                    cumulative <- cumulative + interval_mass
                }
            }
        }
    }

    breakpoints[length(breakpoints)]
}


interval_union_is_subset <- function(union_a, union_b, tol = 1e-10) {
    if (!is.finite(tol) || tol < 0) {
        stop("tol must be a nonnegative finite number.")
    }

    union_a <- normalize_interval_union(union_a)
    union_b <- normalize_interval_union(union_b)

    if (nrow(union_a) == 0) {
        return(TRUE)
    }
    if (nrow(union_b) == 0) {
        return(FALSE)
    }

    for (i in seq_len(nrow(union_a))) {
        covered <- any(
            union_b[, 1] <= union_a[i, 1] + tol &
                union_b[, 2] >= union_a[i, 2] - tol
        )
        if (!covered) {
            return(FALSE)
        }
    }

    TRUE
}


make_regular_grid <- function(region,
                              n_grid = 500,
                              grid_resolution = 0.05,
                              max_grid_size = 2000) {
    region <- as.numeric(region)
    if (length(region) != 2 || any(is.na(region)) || any(!is.finite(region))) {
        stop("region must be a finite numeric vector of length two.")
    }
    region <- sort(region)
    if (region[1] >= region[2]) {
        stop("region must have positive width.")
    }
    if (length(n_grid) != 1 || is.na(n_grid) || n_grid < 2) {
        stop("n_grid must be at least 2.")
    }
    if (length(grid_resolution) != 1 || is.na(grid_resolution) || grid_resolution <= 0) {
        stop("grid_resolution must be positive.")
    }
    if (length(max_grid_size) != 1 || is.na(max_grid_size) || max_grid_size < 2) {
        stop("max_grid_size must be at least 2.")
    }

    n_grid <- as.integer(n_grid)
    max_grid_size <- as.integer(max_grid_size)
    region_width <- region[2] - region[1]
    n_grid_actual <- min(
        max(n_grid, ceiling(region_width / grid_resolution) + 1L),
        max_grid_size
    )

    seq(region[1], region[2], length.out = n_grid_actual)
}


consecutive_points_to_intervals <- function(grid_points, tolerance = NULL) {
    if (length(grid_points) == 0) {
        return(empty_interval_union())
    }

    grid_points <- sort(grid_points)
    if (is.null(tolerance)) {
        tolerance <- if (length(grid_points) > 1) median(diff(grid_points)) * 1.5 else 0
    }

    if (length(grid_points) == 1) {
        intervals <- matrix(c(grid_points[1], grid_points[1]), nrow = 1)
    } else {
        gaps <- diff(grid_points)
        is_gap <- gaps > tolerance
        starts <- c(1, which(is_gap) + 1)
        ends <- c(which(is_gap), length(grid_points))
        intervals <- cbind(grid_points[starts], grid_points[ends])
    }

    normalize_interval_union(intervals)
}


union_total_length <- function(union) {
    union <- normalize_interval_union(union)
    if (nrow(union) == 0) {
        return(0)
    }
    sum(pmax(union[, 2] - union[, 1], 0))
}


point_in_union <- function(point, union) {
    union <- normalize_interval_union(union)
    if (nrow(union) == 0) {
        return(FALSE)
    }
    any(point >= union[, 1] & point <= union[, 2])
}


points_in_union <- function(points, union) {
    union <- normalize_interval_union(union)
    if (nrow(union) == 0) {
        return(rep(FALSE, length(points)))
    }

    inside <- rep(FALSE, length(points))
    for (k in seq_len(nrow(union))) {
        inside <- inside | (points >= union[k, 1] & points <= union[k, 2])
    }
    inside
}


sample_from_union <- function(union, n_samples) {
    union <- normalize_interval_union(union)
    if (n_samples <= 0) {
        return(numeric(0))
    }
    if (nrow(union) == 0 || union_total_length(union) <= 0) {
        stop("Cannot sample from an empty interval union.")
    }

    lengths <- union[, 2] - union[, 1]
    if (any(!is.finite(lengths))) {
        stop("Cannot sample from an unbounded interval union.")
    }
    counts <- rmultinom(1, n_samples, lengths / sum(lengths))[, 1]
    samples <- numeric(0)
    for (k in seq_len(nrow(union))) {
        if (counts[k] > 0) {
            samples <- c(samples, runif(counts[k], union[k, 1], union[k, 2]))
        }
    }
    samples
}


difference_interval_unions <- function(union_a, union_b) {
    union_a <- normalize_interval_union(union_a)
    union_b <- normalize_interval_union(union_b)
    if (nrow(union_a) == 0 || nrow(union_b) == 0) {
        return(empty_interval_union())
    }

    pieces <- matrix(numeric(0), ncol = 2)
    for (i in seq_len(nrow(union_a))) {
        for (j in seq_len(nrow(union_b))) {
            pieces <- rbind(
                pieces,
                c(union_a[i, 1] - union_b[j, 2], union_a[i, 2] - union_b[j, 1])
            )
        }
    }
    normalize_interval_union(pieces)
}


summarize_prediction_sets <- function(true_values, sets) {
    if (is.data.frame(sets) && all(c("lower", "upper") %in% names(sets))) {
        sets <- lapply(seq_len(nrow(sets)), function(i) {
            matrix(c(sets$lower[i], sets$upper[i]), nrow = 1)
        })
    }

    n <- length(true_values)
    coverage <- logical(n)
    lengths <- numeric(n)
    n_intervals <- integer(n)
    for (i in seq_len(n)) {
        set_i <- normalize_interval_union(sets[[i]])
        coverage[i] <- point_in_union(true_values[i], set_i)
        lengths[i] <- union_total_length(set_i)
        n_intervals[i] <- nrow(set_i)
    }

    list(
        coverage = mean(coverage),
        mean_length = mean(lengths),
        se_length = stats::sd(lengths) / sqrt(length(lengths)),
        min_length = min(lengths),
        max_length = max(lengths),
        mean_n_intervals = mean(n_intervals),
        coverage_vec = coverage,
        lengths = lengths,
        n_intervals = n_intervals
    )
}
