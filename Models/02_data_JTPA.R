################################################################################
# Module: 02_data_JTPA.R
# Purpose: Prepare National JTPA Study real data in the same data_simul
#          structure used by S1-S4
################################################################################


################################################################################
# Helpers
################################################################################

to_numeric_clean <- function(x) {
    x <- trimws(as.character(x))
    x[x == ""] <- NA_character_
    suppressWarnings(as.numeric(x))
}


################################################################################
# Load and preprocess the fixed JTPA empirical pool
################################################################################

load_jtpa_pool <- function(expbif_path,
                           earns_path,
                           sample_group = c("full", "adult"),
                           outcome = c("30m", "12m", "6m", "month30"),
                           y_scale = 1000) {

    sample_group <- match.arg(sample_group)
    outcome <- match.arg(outcome)

    if (!requireNamespace("haven", quietly = TRUE)) {
        stop("Package `haven` is required to read JTPA Stata files.")
    }

    if (!file.exists(expbif_path)) {
        stop("Cannot find expbif.dta: ", expbif_path)
    }

    if (!file.exists(earns_path)) {
        stop("Cannot find earns.dta: ", earns_path)
    }

    bif <- haven::read_dta(expbif_path)
    earns <- haven::read_dta(earns_path)

    # --------------------------------------------------------------------------
    # Match the earnings sample to the experimental baseline file.
    # Keep the ordering of earns.dta.
    # --------------------------------------------------------------------------

    idx <- match(as.character(earns$recid),
                 as.character(bif$recid))

    if (anyNA(idx)) {
        stop(
            sum(is.na(idx)),
            " earnings observations cannot be matched to expbif.dta."
        )
    }

    b <- bif[idx, , drop = FALSE]

    # --------------------------------------------------------------------------
    # Optional adult-sample restriction
    # --------------------------------------------------------------------------

    keep_group <- rep(TRUE, nrow(earns))

    if (sample_group == "adult") {
        adultm <- to_numeric_clean(b$adultm)
        adultf <- to_numeric_clean(b$adultf)

        keep_group <- adultm == 1 | adultf == 1
    }

    b <- b[keep_group, , drop = FALSE]
    earns <- earns[keep_group, , drop = FALSE]

    # --------------------------------------------------------------------------
    # Treatment
    #
    # treatmt is the randomized JTPA assignment indicator:
    #     1 = treatment
    #     0 = control
    # --------------------------------------------------------------------------

    T <- to_numeric_clean(b$treatmt)

    if (!all(T[!is.na(T)] %in% c(0, 1))) {
        stop("Unexpected coding in JTPA `treatmt`.")
    }

    # --------------------------------------------------------------------------
    # Baseline covariates
    #
    # Seven-covariate specification:
    #   age
    #   female
    #   black
    #   hispanic
    #   years of education
    #   married
    #   previous-year earnings
    # --------------------------------------------------------------------------

    age <- to_numeric_clean(b$age)
    female <- to_numeric_clean(b$female)
    black <- to_numeric_clean(b$black)
    hispanic <- to_numeric_clean(b$hispanic)

    education <- to_numeric_clean(b$bfeduca)
    education[education == 99] <- NA_real_

    married <- to_numeric_clean(b$mstatmar)
    married[married == 9] <- NA_real_

    pre_earnings <- to_numeric_clean(b$yearearn)
    pre_earnings[pre_earnings == 99999] <- NA_real_

    # Put previous earnings in $1,000 units as well.
    pre_earnings <- pre_earnings / y_scale

    X <- cbind(
        age = age,
        female = female,
        black = black,
        hispanic = hispanic,
        education = education,
        married = married,
        pre_earnings = pre_earnings
    )

    storage.mode(X) <- "double"

    # --------------------------------------------------------------------------
    # Outcome
    #
    # NEWERN01,...,NEWERN30 are monthly post-assignment earnings.
    # --------------------------------------------------------------------------

    earn_names <- sprintf("newern%02d", 1:30)

    missing_earn_vars <- setdiff(earn_names, names(earns))

    if (length(missing_earn_vars) > 0L) {
        stop(
            "Missing earnings variables: ",
            paste(missing_earn_vars, collapse = ", ")
        )
    }

    earn_mat <- sapply(
        earns[, earn_names, drop = FALSE],
        as.numeric
    )

    if (!is.matrix(earn_mat)) {
        earn_mat <- as.matrix(earn_mat)
    }

    Y <- switch(
        outcome,
        "30m" =
            rowSums(earn_mat[, 1:30, drop = FALSE]),

        "12m" =
            rowSums(earn_mat[, 19:30, drop = FALSE]),

        "6m" =
            rowSums(earn_mat[, 25:30, drop = FALSE]),

        "month30" =
            earn_mat[, 30]
    )

    # Work in thousands of dollars.
    Y <- Y / y_scale

    # --------------------------------------------------------------------------
    # Complete-case empirical pool
    # --------------------------------------------------------------------------

    keep <- complete.cases(X) &
            !is.na(T) &
            !is.na(Y) &
            is.finite(Y)

    X <- X[keep, , drop = FALSE]
    T <- T[keep]
    Y <- Y[keep]

    if (any(!is.finite(X))) {
        stop("Non-finite covariates remain after JTPA preprocessing.")
    }

    if (!all(T %in% c(0, 1))) {
        stop("Treatment must be binary after preprocessing.")
    }

    list(
        X = X,
        T = as.numeric(T),
        Y = as.numeric(Y),
        ids = as.character(earns$recid)[keep]
    )
}


################################################################################
# Build the same list structure used by build_data_simul()
################################################################################

build_data_real <- function(X,
                            T,
                            Y,
                            ntrain,
                            ntest,
                            seed,
                            ps_value = NULL) {

    X <- as.matrix(X)
    T <- as.numeric(T)
    Y <- as.numeric(Y)

    n <- nrow(X)

    if (length(T) != n || length(Y) != n) {
        stop("X, T, and Y have incompatible dimensions.")
    }

    n_required <- ntrain + ntest

    if (n_required > n) {
        stop(
            "Requested ntrain + ntest = ", n_required,
            ", but the JTPA empirical pool contains only ", n,
            " usable observations."
        )
    }

    if (ntrain < 1L || ntest < 1L) {
        stop("ntrain and ntest must be positive.")
    }

    # --------------------------------------------------------------------------
    # Fixed real dataset; seed changes only the train/test split.
    # --------------------------------------------------------------------------

    set.seed(seed)

    ids <- sample.int(n, n_required, replace = FALSE)

    trainid <- ids[seq_len(ntrain)]
    testid <- ids[ntrain + seq_len(ntest)]

    # --------------------------------------------------------------------------
    # Randomized trial propensity.
    #
    # If no design probability is supplied, use one constant empirical
    # assignment probability. Since treatment is randomized, e(X) is constant.
    # --------------------------------------------------------------------------

    if (is.null(ps_value)) {
        ps_value <- mean(T)
    }

    if (!is.finite(ps_value) ||
        ps_value <= 0 ||
        ps_value >= 1) {
        stop("ps_value must lie strictly between 0 and 1.")
    }

    pstrain <- rep(ps_value, ntrain)
    pstest <- rep(ps_value, ntest)

    # --------------------------------------------------------------------------
    # Oracle quantities do not exist in real data.
    # Preserve the simulation interface with numeric NA placeholders.
    # --------------------------------------------------------------------------

    oracle_na <- rep(NA_real_, ntest)

    list(
        Xtrain = X[trainid, , drop = FALSE],
        Ttrain = T[trainid],
        Ytrain = Y[trainid],

        Xtest = X[testid, , drop = FALSE],
        Ttest = T[testid],
        Ytest = Y[testid],

        Y0test = oracle_na,
        Y1test = oracle_na,
        tautest = oracle_na,
        CATEtest = oracle_na,
        stdtest = oracle_na,

        pstrain = pstrain,
        pstest = pstest
    )
}


################################################################################
# Public interface
################################################################################

generate_data_JTPA <- function(expbif_path,
                               earns_path,
                               ntrain,
                               ntest,
                               seed = 0,
                               sample_group = "full",
                               outcome = "30m",
                               y_scale = 1000,
                               ps_value = NULL) {

    pool <- load_jtpa_pool(
        expbif_path = expbif_path,
        earns_path = earns_path,
        sample_group = sample_group,
        outcome = outcome,
        y_scale = y_scale
    )

    build_data_real(
        X = pool$X,
        T = pool$T,
        Y = pool$Y,
        ntrain = ntrain,
        ntest = ntest,
        seed = seed,
        ps_value = ps_value
    )
}
