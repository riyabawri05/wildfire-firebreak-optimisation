# ============================================================
# Fire Simulator — Context A: Wildfire Firebreak Planning
# MA324 Mathematical Modelling and Simulation
# ============================================================
#
# Simulates fire spread on a grid using a probabilistic
# cellular automaton based on Alexandridis et al. (2008)
# DOI: 10.1016/j.amc.2008.06.046
#
# Cell states:
#   0 = NO_FUEL   (bare ground / firebreak — cannot burn)
#   1 = UNBURNED  (vegetated — can catch fire)
#   2 = BURNING   (on fire — spreads to neighbours this step)
#   3 = BURNED    (fire exhausted — cannot burn again)
#
# Spread probability from burning cell (i,j) to unburned
# neighbour (i',j'):
#
#   p_burn = p0 * (1 + p_veg) * (1 + p_den) * p_wind
#
# clamped to [0, 1], where:
#
#   p0    = 0.25            base ignition probability
#   p_veg = vegetation factor of receiving cell
#           grassland: -0.3, shrubland: 0.0, forest: +0.4
#   p_den = density factor of receiving cell
#           sparse: -0.4, normal: 0.0, dense: +0.3
#   p_wind = exp(c1 * V) * exp(V * c2 * (cos(theta) - 1))
#           c1 = 0.045, c2 = 0.05
#           V = wind speed (km/h)
#           theta = angle between fire push direction and
#                   direction from burning cell to neighbour
#
# Wind convention:
#   wind_dir = direction wind comes FROM (compass: 0 = N,
#   pi/2 = E, pi = S, 3pi/2 = W). Fire is pushed toward
#   wind_dir + pi.
#
# Neighbourhood: Moore (8 adjacent cells), no wrapping.
# ============================================================


# === Constants ===============================================
P0 = 0.25
C1 = 0.045
C2 = 0.05

P_VEG = c(-0.3, 0.0, 0.4)   # grass, shrub, forest
P_DEN = c(-0.4, 0.0, 0.3)   # sparse, normal, dense

# Neighbour offsets and compass angles (0 = N, pi/2 = E, ...)
NEIGHBOUR_DR  = c(-1, -1,  0,  1,  1,  1,  0, -1)
NEIGHBOUR_DC  = c( 0,  1,  1,  1,  0, -1, -1, -1)
NEIGHBOUR_PHI = c(0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4)


# === Internal helpers ========================================

parse_landscape = function(landscape) {
    list(veg_type = landscape %/% 10L, density = landscape %% 10L)
}

compute_p_burn = function(veg_j, den_j, wind_speed, wind_dir, phi_n) {
    push_dir = wind_dir + pi
    theta = push_dir - phi_n
    p_wind = exp(C1 * wind_speed) * exp(wind_speed * C2 * (cos(theta) - 1))
    p = P0 * (1 + P_VEG[veg_j]) * (1 + P_DEN[den_j]) * p_wind
    min(max(p, 0), 1)
}


# === Main simulator ==========================================

#' Simulate fire spread on a landscape grid
#'
#' Probabilistic cellular automaton: fire spreads from burning
#' cells to unburned Moore neighbours with probability p_burn.
#' Bare ground (code 0) blocks fire entirely.
#'
#' @param landscape integer matrix (grid_size x grid_size),
#'   cells encoded as veg_type * 10 + density (see header)
#' @param ignition length-2 vector c(row, col) or matrix of
#'   ignition points (one per row)
#' @param wind_speed wind speed in km/h (default 0)
#' @param wind_dir direction wind comes FROM in radians,
#'   compass convention: 0=N, pi/2=E, pi=S, 3pi/2=W (default 0)
#' @return list with burned (logical matrix), burn_time
#'   (integer matrix, NA if not burned), total_burned (integer)
simulate_fire = function(landscape, ignition, wind_speed = 0, wind_dir = 0) {
    nr = nrow(landscape); nc = ncol(landscape)
    parsed = parse_landscape(landscape)

    # State: 0=NO_FUEL, 1=UNBURNED, 2=BURNING, 3=BURNED
    state = matrix(1L, nr, nc)
    state[landscape == 0L] = 0L

    if (is.vector(ignition) && length(ignition) == 2)
        ignition = matrix(ignition, nrow = 1)
    for (k in 1:nrow(ignition)) {
        r = ignition[k, 1]; cc = ignition[k, 2]
        if (state[r, cc] == 0L) stop("Cannot ignite bare ground at (", r, ",", cc, ")")
        state[r, cc] = 2L
    }

    burn_time = matrix(NA_integer_, nr, nc)
    burn_time[state == 2L] = 0L
    time = 0L

    repeat {
        burning = which(state == 2L, arr.ind = TRUE)
        if (nrow(burning) == 0) break
        time = time + 1L
        will_ignite = matrix(FALSE, nr, nc)

        for (b in 1:nrow(burning)) {
            ri = burning[b, 1]; ci = burning[b, 2]
            for (d in 1:8) {
                rj = ri + NEIGHBOUR_DR[d]
                cj = ci + NEIGHBOUR_DC[d]
                if (rj < 1 || rj > nr || cj < 1 || cj > nc) next
                if (state[rj, cj] != 1L) next
                p = compute_p_burn(
                    parsed$veg_type[rj, cj], parsed$density[rj, cj],
                    wind_speed, wind_dir, NEIGHBOUR_PHI[d]
                )
                if (runif(1) < p) will_ignite[rj, cj] = TRUE
            }
        }

        state[state == 2L] = 3L
        state[will_ignite] = 2L
        burn_time[will_ignite] = time
    }

    burned = (state == 3L)
    list(burned = burned, burn_time = burn_time, total_burned = sum(burned))
}


# === Helpers =================================================

#' Set cells to bare ground (firebreaks)
#'
#' @param landscape integer matrix
#' @param cells matrix with columns (row, col), or length-2 vector
#' @return modified landscape
set_firebreaks = function(landscape, cells) {
    out = landscape
    if (is.vector(cells) && length(cells) == 2)
        cells = matrix(cells, nrow = 1)
    if (nrow(cells) == 0) return(out)
    for (k in 1:nrow(cells))
        out[cells[k, 1], cells[k, 2]] = 0L
    out
}

#' Compute weighted damage from a burn matrix
#'
#' @param burned logical matrix
#' @param targets data frame with columns row, col, weight
#' @return numeric total damage
compute_damage = function(burned, targets) {
    sum(targets$weight * burned[cbind(targets$row, targets$col)])
}

#' Compute spread probabilities for all adjacent pairs
#'
#' @param landscape integer matrix
#' @param wind_speed wind speed in km/h (default 0)
#' @param wind_dir wind direction in radians (default 0)
#' @return data frame with from_row, from_col, to_row, to_col, p_burn
spread_probabilities = function(landscape, wind_speed = 0, wind_dir = 0) {
    nr = nrow(landscape); nc = ncol(landscape)
    parsed = parse_landscape(landscape)
    rows = list()
    idx = 0

    for (i in 1:nr) for (j in 1:nc) {
        if (landscape[i, j] == 0L) next
        for (d in 1:8) {
            ri = i + NEIGHBOUR_DR[d]
            ci = j + NEIGHBOUR_DC[d]
            if (ri < 1 || ri > nr || ci < 1 || ci > nc) next
            if (landscape[ri, ci] == 0L) next
            p = compute_p_burn(
                parsed$veg_type[ri, ci], parsed$density[ri, ci],
                wind_speed, wind_dir, NEIGHBOUR_PHI[d]
            )
            idx = idx + 1
            rows[[idx]] = data.frame(
                from_row = i, from_col = j,
                to_row = ri, to_col = ci,
                p_burn = p
            )
        }
    }
    do.call(rbind, rows)
}

#' Plot landscape with optional burn and firebreak overlays
#'
#' @param landscape integer matrix
#' @param burned logical matrix (optional)
#' @param firebreaks matrix with columns (row, col) (optional)
#' @param targets data frame with row, col, weight (optional)
#' @param main plot title
plot_fire = function(landscape, burned = NULL, firebreaks = NULL,
                     targets = NULL, main = "") {
    nr = nrow(landscape); nc = ncol(landscape)

    # --- Vegetation colour map ---
    veg_codes = c(0, 11, 12, 13, 21, 22, 23, 31, 32, 33)
    veg_cols  = c(
        "grey85",                           # bare ground
        "#D5E8D4", "#97D077", "#56A439",    # grassland s/n/d
        "#FFF2CC", "#FFD966", "#D6A21E",    # shrubland s/n/d
        "#B5CEA8", "#548235", "#2E5E0F"     # forest s/n/d
    )
    col_map = setNames(veg_cols, as.character(veg_codes))

    # --- Canvas ---
    par(mar = c(1, 1, 2.5, 1))
    plot(NULL, xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5),
         asp = 1, axes = FALSE, xlab = "", ylab = "", main = main)

    # --- Layer 1: vegetation base ---
    for (i in 1:nr) for (j in 1:nc)
        rect(j - 0.5, nr - i + 0.5, j + 0.5, nr - i + 1.5,
             col = col_map[as.character(landscape[i, j])],
             border = "white", lwd = 0.3)

    # --- Layer 2: burned cells (semi-transparent terracotta) ---
    if (NROW(burned) > 0)
        for (i in 1:nr) for (j in 1:nc)
            if (burned[i, j])
                rect(j - 0.5, nr - i + 0.5, j + 0.5, nr - i + 1.5,
                     col = "#D6604D88", border = NA)

    # --- Layer 3: firebreak borders (black) ---
    if (NROW(firebreaks) > 0)
        for (k in 1:nrow(firebreaks))
            rect(firebreaks[k, 2] - 0.5, nr - firebreaks[k, 1] + 0.5,
                 firebreaks[k, 2] + 0.5, nr - firebreaks[k, 1] + 1.5,
                 col = NA, border = "black", lwd = 3.5)

    # --- Layer 4: target borders (blue=wetland, terracotta=settlement) ---
    if (NROW(targets) > 0)
        for (k in 1:nrow(targets)) {
            col = "#2166AC"
            rect(targets$col[k] - 0.5, nr - targets$row[k] + 0.5,
                 targets$col[k] + 0.5, nr - targets$row[k] + 1.5,
                 col = NA, border = col, lwd = 2.5)
        }
}


# === Tests (run with: Rscript simulate_fire.r) ===============
if (sys.nframe() == 0) {
    library(unittest, quietly = TRUE)

    L = as.matrix(read.table("landscape.csv", sep = ","))
    targets = read.csv("targets.csv")
    n = nrow(L)
    mid = ceiling(n / 2)
    K = 40  # replications per stochastic test

    # --- Aggregate behaviour ---

    ok(local({
        fracs = replicate(K, {
            r = simulate_fire(L, c(mid, mid))
            r$total_burned / (n * n)
        })
        cat("  avg burn fraction:", round(100 * mean(fracs)), "%\n")
        mean(fracs) > 0.3 && mean(fracs) < 0.7
    }), "centre, no wind: avg 30-70% burned")

    ok(local({
        L_cor = L; L_cor[mid, ] = 0L
        norths = replicate(K, {
            r = simulate_fire(L_cor, c(n, mid))
            sum(r$burned[1:(mid - 1), ])
        })
        cat("  max north burned:", max(norths), "\n")
        all(norths == 0)
    }), "corridor blocks all S to N fire")

    ok(local({
        east = 0; west = 0
        for (s in 1:K) {
            set.seed(s)
            r = simulate_fire(L, c(mid, mid), wind_speed = 10, wind_dir = 3 * pi / 2)
            east = east + sum(r$burned[, (mid + 1):n])
            west = west + sum(r$burned[, 1:(mid - 1)])
        }
        cat("  avg E burned:", east / K, " avg W burned:", west / K, "\n")
        east > west * 1.6
    }), "wind from W: east burns >1.6x west")

    ok(local({
        front = cbind(rep(n, n), 1:n)
        fracs = replicate(K, {
            r = simulate_fire(L, front)
            r$total_burned / (n * n)
        })
        cat("  avg burn fraction:", round(100 * mean(fracs)), "%\n")
        mean(fracs) > 0.55
    }), "south front, no wind: avg >55% burned")

    ok(local({
        L_bare = matrix(0L, n, n)
        L_bare[mid, mid] = 22L
        r = simulate_fire(L_bare, c(mid, mid))
        r$total_burned == 1
    }), "bare ground: only ignition cell burns")

    ok(local({
        set.seed(324)
        r1 = simulate_fire(L, c(mid, mid))
        set.seed(324)
        r2 = simulate_fire(L, c(mid, mid))
        identical(r1$burned, r2$burned)
    }), "same seed gives identical result")

    # --- Formula correctness ---

    ok(local({
        p = compute_p_burn(2, 2, 0, 0, 0)
        cat("  shrubland normal, no wind:", p, "\n")
        abs(p - P0) < 1e-10
    }), "p_burn: shrubland normal, no wind = p0")

    ok(local({
        p = compute_p_burn(3, 3, 0, 0, 0)
        expected = P0 * (1 + 0.4) * (1 + 0.3)
        cat("  dense forest, no wind:", p, " expected:", expected, "\n")
        abs(p - expected) < 1e-10
    }), "p_burn: dense forest, no wind = p0 * 1.4 * 1.3")

    # --- Wind direction convention ---

    ok(local({
        south = 0; north = 0
        for (s in 1:K) {
            set.seed(s)
            r = simulate_fire(L, c(mid, mid), wind_speed = 10, wind_dir = 0)
            south = south + sum(r$burned[(mid + 1):n, ])
            north = north + sum(r$burned[1:(mid - 1), ])
        }
        cat("  avg S burned:", south / K, " avg N burned:", north / K, "\n")
        south > north
    }), "wind from N: south burns more than north")

    ok(local({
        south = 0; north = 0
        for (s in 1:K) {
            set.seed(s)
            r = simulate_fire(L, c(mid, mid), wind_speed = 15, wind_dir = pi)
            south = south + sum(r$burned[(mid + 1):n, ])
            north = north + sum(r$burned[1:(mid - 1), ])
        }
        cat("  avg S burned:", south / K, " avg N burned:", north / K, "\n")
        north > south
    }), "wind from S: north burns more than south")

    # --- Vegetation differentiation ---

    ok(local({
        L_forest = matrix(33L, n, n)
        L_grass  = matrix(11L, n, n)
        f_fracs = replicate(K, {
            r = simulate_fire(L_forest, c(mid, mid))
            r$total_burned / (n * n)
        })
        g_fracs = replicate(K, {
            r = simulate_fire(L_grass, c(mid, mid))
            r$total_burned / (n * n)
        })
        cat("  avg forest:", round(100 * mean(f_fracs)), "%",
            " avg grass:", round(100 * mean(g_fracs)), "%\n")
        mean(f_fracs) > mean(g_fracs)
    }), "dense forest burns more than sparse grassland")

    # --- burn_time ---

    ok(local({
        set.seed(324)
        L_uni = matrix(22L, 5, 5)
        r = simulate_fire(L_uni, c(3, 3))
        t0 = r$burn_time[3, 3]
        neighbours_t = r$burn_time[2:4, 2:4]
        neighbours_t = neighbours_t[!is.na(neighbours_t) & neighbours_t > 0]
        cat("  ignition t:", t0, " neighbour min t:", min(neighbours_t), "\n")
        t0 == 0 && all(neighbours_t >= 1)
    }), "burn_time: ignition=0, neighbours>=1")

    # --- Synchronous update ---

    ok(local({
        L_strip = matrix(0L, 1, 5)
        L_strip[1, ] = 22L
        caught_at_1 = FALSE
        for (s in 1:200) {
            set.seed(s)
            r = simulate_fire(L_strip, c(1, 1))
            if (!is.na(r$burn_time[1, 3]) && r$burn_time[1, 3] < 2)
                caught_at_1 = TRUE
        }
        cat("  cell 3 burned at t<2:", caught_at_1, "\n")
        !caught_at_1
    }), "synchronous: cell 3 never burns before t=2")

    # --- Helpers ---

    ok(local({
        L3 = matrix(22L, 3, 3)
        edges = spread_probabilities(L3, 0, 0)
        cat("  edges:", nrow(edges), "\n")
        nrow(edges) == 40 && all(abs(edges$p_burn - P0) < 1e-10)
    }), "spread_probabilities: 3x3 shrubland = 40 edges, all p0")

    ok(local({
        L3 = matrix(22L, 3, 3)
        L3[2, 2] = 0L
        edges = spread_probabilities(L3, 0, 0)
        has_centre = any(
            (edges$from_row == 2 & edges$from_col == 2) |
            (edges$to_row == 2 & edges$to_col == 2)
        )
        cat("  edges involving centre:", has_centre, "\n")
        !has_centre
    }), "spread_probabilities: bare ground excluded from edges")

    ok(local({
        burned = matrix(FALSE, n, n)
        burned[3, 3] = TRUE
        burned[3, 16] = TRUE
        d = compute_damage(burned, targets)
        cat("  damage:", d, " expected:", 5 + 10, "\n")
        d == 15
    }), "compute_damage: correct for known burned cells")

    ok(local({
        L_test = matrix(22L, 3, 3)
        L_out = set_firebreaks(L_test, matrix(c(1, 1, 2, 2), nrow = 2, byrow = TRUE))
        L_out[1, 1] == 0L && L_out[2, 2] == 0L && L_out[3, 3] == 22L
    }), "set_firebreaks: clears specified cells only")

    # --- Visual check ---

    r1 = simulate_fire(L, c(mid, mid))
    L_cor = L; L_cor[mid, ] = 0L
    r2 = simulate_fire(L_cor, c(n, mid))
    r3 = simulate_fire(L, c(mid, mid), wind_speed = 10, wind_dir = 3 * pi / 2)
    r4 = simulate_fire(L, cbind(rep(n, n), 1:n))

    par(mfrow = c(2, 2), bg = "white")
    plot_fire(L, r1$burned, targets = targets, main = "Centre, no wind")
    plot_fire(L_cor, r2$burned, targets = targets, main = "S + corridor")
    plot_fire(L, r3$burned, targets = targets, main = "Wind from W")
    plot_fire(L, r4$burned, targets = targets, main = "South front")
}