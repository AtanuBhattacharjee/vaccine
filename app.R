## =============================================================================
##  Which optimal biological dose?
##  Interactive explorer for estimand ambiguity under informative cluster size.
##
##  Companion to: Bhattacharjee A. "Which optimal biological dose? Estimand
##  ambiguity under informative cluster size, with application to multi-epitope
##  cancer vaccines."
##
##  Deploy to Posit Connect Cloud -- see DEPLOY.md.
##  Depends only on shiny + shinydashboard; all numerical work is base R.
## =============================================================================

library(shiny)
library(shinydashboard)

## ---------------------------------------------------------------- palette ---
COL_B   <- "#1f6fb4"   # breadth objective  Q_B
COL_A   <- "#d1495b"   # cluster-averaged   Q_A
COL_POP <- "#c9ced6"
COL_GRID<- "#e8ebef"
COL_INK <- "#1b2733"

## ============================ MODEL ==========================================
## Generative model of Section 4.1, evaluated exactly -- no simulation.
##
##   stimulus  l_ij = log A + b_i + log w_ij + u_i
##   logit pi  = th0 + th1 l + th2 l^2 - phi c log(m_i) (u_i - log d_ref)
##
## Conditional on HLA supertype S the shift b_i + kappa*tau*S is normal and the
## panel size m_i is independent of it, so both objectives are finite sums over
## two strata. At tau = 0 with equal quality the strata coincide and the
## separation is exactly zero, which is Proposition 1.
## -----------------------------------------------------------------------------
TH0 <- 0.20; TH1 <- 0.00; TH2 <- -0.50
SIGMA_W <- 0.5; D_REF <- 400
DOSES <- c(100, 250, 600, 1500)
MLO <- 4; MHI <- 30

## Gauss-Hermite (probabilists'), Golub-Welsch
gh_nodes <- function(n) {
  J <- matrix(0, n, n)
  if (n > 1) {
    k <- seq_len(n - 1)
    J[cbind(k, k + 1)] <- sqrt(k); J[cbind(k + 1, k)] <- sqrt(k)
  }
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = (e$vectors[1, o])^2)
}
GH  <- gh_nodes(24)
GHq <- gh_nodes(10)     # coarser quality quadrature on the competition path

## F(z): patient-level mean response, integrated over epitope quality
.ZG <- seq(-14, 14, length.out = 8001)
.FT <- {
  ell <- outer(.ZG, SIGMA_W * GH$x, "+")
  as.vector(plogis(TH0 + TH1 * ell + TH2 * ell^2) %*% GH$w)
}
Ffun <- function(z) approx(.ZG, .FT, xout = z, rule = 2)$y

## exact first and second derivatives of F (do NOT difference the table above --
## its knots are 0.0035 apart and a small step collapses the second difference)
Fderiv <- function(z) {
  ell <- outer(z, SIGMA_W * GH$x, "+")
  pr  <- plogis(TH0 + TH1 * ell + TH2 * ell^2)
  q1  <- TH1 + 2 * TH2 * ell
  v   <- pr * (1 - pr)
  list(d1 = as.vector((v * q1) %*% GH$w),
       d2 = as.vector((v * (2 * TH2 + (1 - 2 * pr) * q1^2)) %*% GH$w))
}

## clipped Poisson, as in pmin(pmax(rpois(lam), 4), 30)
pois_pmf <- function(lam) {
  raw <- dpois(0:(MHI + 60), lam)
  p <- numeric(MHI - MLO + 1)
  for (m in MLO:MHI) {
    i <- m - MLO + 1
    p[i] <- if (m == MLO) sum(raw[1:(MLO + 1)])
            else if (m == MHI) sum(raw[(MHI + 1):length(raw)])
            else raw[m + 1]
  }
  p / sum(p)
}

strata <- function(p) {
  lapply(list(c(1, p$p_common, p$lam_common), c(0, 1 - p$p_common, p$lam_rare)),
         function(v) {
           pm <- pois_pmf(v[3])
           list(S = v[1], pr = v[2], pm = pm, mgrid = MLO:MHI,
                mbar = sum((MLO:MHI) * pm),
                shift = p$tau * v[1] + p$kappa * p$tau * v[1])
         })
}

objectives <- function(u, st, logA, phi, cval, sigma_b, breadth_only = FALSE) {
  QB <- numeric(length(u)); QA <- numeric(length(u))
  for (s in st) {
    if (phi == 0) {
      z <- outer(u, logA + s$shift + sigma_b * GH$x, "+")
      acc <- as.vector(matrix(Ffun(z), nrow = length(u)) %*% GH$w)
      QA <- QA + s$pr * acc; QB <- QB + s$pr * s$mbar * acc
    } else {
      keep <- which(s$pm > 1e-4)
      a0 <- numeric(length(u)); a1 <- numeric(length(u))
      for (k in keep) {
        m <- s$mgrid[k]
        comp <- phi * cval * log(m) * (u - log(D_REF))
        acc <- numeric(length(u))
        for (g in seq_along(GH$x)) {
          z <- u + logA + s$shift + sigma_b * GH$x[g]
          for (q in seq_along(GHq$x)) {
            ell <- z + SIGMA_W * GHq$x[q]
            acc <- acc + GH$w[g] * GHq$w[q] *
              plogis(TH0 + TH1 * ell + TH2 * ell^2 - comp)
          }
        }
        a0 <- a0 + s$pm[k] * acc; a1 <- a1 + s$pm[k] * m * acc
      }
      QA <- QA + s$pr * a0; QB <- QB + s$pr * a1
    }
  }
  if (breadth_only) QB else list(QB = QB, QA = QA)
}

argmax_fine <- function(u, Q) {
  k <- which.max(Q)
  if (k > 1 && k < length(Q)) {
    y0 <- Q[k - 1]; y1 <- Q[k]; y2 <- Q[k + 1]
    den <- y0 - 2 * y1 + y2
    if (den != 0) return(u[k] + 0.5 * (u[2] - u[1]) * (y0 - y2) / den)
  }
  u[k]
}

## peak of Q_B, located on a coarse grid then refined locally
peak_B <- function(st, logA, phi, cval, sigma_b, centre = log(400), half = 1.2,
                   n = 90) {
  u  <- seq(centre - half, centre + half, length.out = n)
  u1 <- argmax_fine(u, objectives(u, st, logA, phi, cval, sigma_b, TRUE))
  st1 <- u[2] - u[1]
  uu <- seq(u1 - 2 * st1, u1 + 2 * st1, length.out = 41)
  argmax_fine(uu, objectives(uu, st, logA, phi, cval, sigma_b, TRUE))
}

## log A such that d*_B sits at the target
calibrate <- function(st, phi, cval, sigma_b, target = 396) {
  u0 <- seq(-4, 4, length.out = 900)
  logA <- argmax_fine(u0, objectives(u0, st, 0, 0, cval, sigma_b, TRUE)) -
    log(target)
  if (phi == 0) return(logA)
  f <- function(a) peak_B(st, a, phi, cval, sigma_b) - log(target)
  out <- try(uniroot(f, c(logA - 3.5, logA + 3.5), tol = 1e-7)$root,
             silent = TRUE)
  if (inherits(out, "try-error")) logA else out
}

solve_model <- function(p) {
  st   <- strata(p)
  logA <- calibrate(st, p$phi, p$cval, p$sigma_b)
  u    <- seq(log(60), log(3000), length.out = if (p$phi > 0) 240 else 420)
  Q    <- objectives(u, st, logA, p$phi, p$cval, p$sigma_b)
  step <- u[2] - u[1]
  refine <- function(u0, which) {
    uu <- seq(u0 - 2 * step, u0 + 2 * step, length.out = 41)
    argmax_fine(uu, objectives(uu, st, logA, p$phi, p$cval, p$sigma_b)[[which]])
  }
  uB <- refine(argmax_fine(u, Q$QB), "QB")
  uA <- refine(argmax_fine(u, Q$QA), "QA")

  ## first-order prediction, eq. (3), over the same strata
  Em <- Eg <- Emg <- Eh <- 0
  for (s in st) {
    z <- uB + logA + s$shift + p$sigma_b * GH$x
    d <- Fderiv(z)
    w <- s$pr * GH$w
    Em  <- Em  + sum(w) * s$mbar
    Eg  <- Eg  + sum(w * d$d1)
    Emg <- Emg + s$mbar * sum(w * d$d1)
    Eh  <- Eh  + sum(w * d$d2)
  }
  list(u = u, QB = Q$QB, QA = Q$QA, uB = uB, uA = uA, logA = logA, st = st,
       sigma_b = p$sigma_b, pred = (Emg - Em * Eg) / (Em * Eh),
       dB = exp(uB), dA = exp(uA), gap = 100 * (exp(uA - uB) - 1))
}

nearest_level <- function(d) DOSES[which.min(abs(log(DOSES) - log(d)))]

## ---- hand-drawn explainer ---------------------------------------------------
## Two patients, two doses, the same responses summarised two ways. Drawn with
## deliberately wobbly strokes so it reads as a sketch on a napkin rather than
## as a result.
doodle_svg <- function() {
  ## a row of epitope circles, the first `filled` of them shaded
  dots <- function(x, y, n, filled, per_row = 10, r = 7.2, gap = 17) {
    out <- character(0)
    for (i in seq_len(n)) {
      col <- (i - 1) %% per_row; row <- (i - 1) %/% per_row
      cx <- x + col * gap; cy <- y + row * (gap + 1)
      wob <- ((i * 37) %% 5 - 2) * 0.35            # deterministic jitter
      out <- c(out, sprintf(
        '<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" stroke="#1b2733" stroke-width="1.7"/>',
        cx + wob, cy - wob, r, if (i <= filled) "#f5a524" else "#ffffff"))
    }
    paste(out, collapse = "")
  }
  ## a stick figure
  person <- function(x, y, s = 1) sprintf(
    '<g stroke="#1b2733" stroke-width="2.4" fill="none" stroke-linecap="round"
        transform="translate(%.1f,%.1f) scale(%.2f)">
       <circle cx="0" cy="0" r="9" fill="#ffffff"/>
       <path d="M0 9 L0 30 M0 15 L-11 24 M0 15 L11 24 M0 30 L-9 46 M0 30 L9 46"/>
     </g>', x, y, s)
  ## a wobbly rounded box
  card <- function(x, y, w, h, fill = "#ffffff", stroke = "#1b2733") sprintf(
    '<path d="M%.0f %.0f q%.0f -6 %.0f 0 q6 %.0f 0 %.0f q-%.0f 6 -%.0f 0 q-6 -%.0f 0 -%.0f z"
        fill="%s" stroke="%s" stroke-width="2.4" stroke-linejoin="round"/>',
    x, y, w / 2, w, h / 2, h, w / 2, w, h / 2, h, fill, stroke)
  tick <- function(x, y, col = "#2e8b57") sprintf(
    '<path d="M%.0f %.0f l7 8 l14 -18" stroke="%s" stroke-width="4.5"
        fill="none" stroke-linecap="round" stroke-linejoin="round"/>', x, y, col)

  paste0(
'<svg viewBox="0 0 980 400" width="100%" style="max-width:980px"
      xmlns="http://www.w3.org/2000/svg"
      font-family="Segoe Print, Bradley Hand, Chalkboard SE, Comic Sans MS, cursive">',

## ---------------- left panel: the lower dose --------------------------------
card(30, 22, 430, 250, "#f4f8fd", "#1f6fb4"),
'<text x="245" y="52" text-anchor="middle" font-size="23" font-weight="bold"
       fill="#1f6fb4">at 250 &#181;g</text>',

person(80, 100),
'<text x="80" y="172" text-anchor="middle" font-size="15" fill="#1b2733">20 peptides</text>',
dots(120, 92, 20, 10),
'<text x="330" y="172" text-anchor="middle" font-size="16" font-weight="bold"
       fill="#1b2733">10 responded</text>',

person(80, 208),
'<text x="80" y="256" text-anchor="middle" font-size="15" fill="#1b2733">8 peptides</text>',
dots(120, 205, 8, 3),
'<text x="330" y="256" text-anchor="middle" font-size="16" font-weight="bold"
       fill="#1b2733">3 responded</text>',

## ---------------- right panel: the higher dose ------------------------------
card(520, 22, 430, 250, "#fdf5f6", "#d1495b"),
'<text x="735" y="52" text-anchor="middle" font-size="23" font-weight="bold"
       fill="#d1495b">at 600 &#181;g</text>',

person(570, 100),
'<text x="570" y="172" text-anchor="middle" font-size="15" fill="#1b2733">20 peptides</text>',
dots(610, 92, 20, 8),
'<text x="820" y="172" text-anchor="middle" font-size="16" font-weight="bold"
       fill="#1b2733">8 responded</text>',

person(570, 208),
'<text x="570" y="256" text-anchor="middle" font-size="15" fill="#1b2733">8 peptides</text>',
dots(610, 205, 8, 4),
'<text x="820" y="256" text-anchor="middle" font-size="16" font-weight="bold"
       fill="#1b2733">4 responded</text>',

## ---------------- the two ways of adding up ---------------------------------
card(30, 292, 430, 92, "#ffffff", "#1f6fb4"),
'<text x="52" y="322" font-size="19" font-weight="bold" fill="#1f6fb4">Count every epitope</text>',
'<text x="52" y="352" font-size="19" fill="#1b2733">10 + 3 = <tspan font-weight="bold">13</tspan>',
'   vs   8 + 4 = 12</text>',
tick(392, 336),
'<text x="52" y="374" font-size="16" fill="#555">&#8594; the lower dose wins</text>',

card(520, 292, 430, 92, "#ffffff", "#d1495b"),
'<text x="542" y="322" font-size="19" font-weight="bold" fill="#d1495b">Average each patient</text>',
'<text x="542" y="352" font-size="19" fill="#1b2733">0.44   vs   <tspan font-weight="bold">0.45</tspan></text>',
tick(882, 336, "#d1495b"),
'<text x="542" y="374" font-size="16" fill="#555">&#8594; the higher dose wins</text>',

'</svg>')
}

## ===================== CLINICIAN-FACING LAYER ================================
## The model parameters tau, kappa and sigma_b are latent. No investigator can
## supply them. What an investigator CAN supply, from a completed dose-finding
## cohort, is:
##
##   * the dose levels administered
##   * which dose gave the largest total number of responding epitopes
##   * the average panel size in each HLA group, and the split between groups
##   * the per-epitope response rate observed in each HLA group at that dose
##
## Only the latent SHIFT between the two HLA groups matters for the separation,
## so kappa is folded into tau and a single quantity, delta, is fitted to the
## observed difference in per-epitope response rate. sigma_b is left at its
## reference value because the separation is nearly insensitive to it.
## -----------------------------------------------------------------------------

## per-epitope response rate implied for each HLA group at dose u
group_rates <- function(p, logA, u) {
  st <- strata(p)
  vapply(st, function(s) {
    mean(sum(GH$w * Ffun(u + logA + s$shift + p$sigma_b * GH$x)))
  }, 0)
}

## Fit the model to two observable quantities. Only the latent SHIFT between
## HLA groups (delta) and the spread of immune capacity (sigma_b) are free:
##   delta   is set by the observed GAP in per-epitope response rate
##   sigma_b is set by the observed LEVEL of response, pooled
## sigma_b is a legitimate free parameter here precisely because the separation
## between the estimands is almost insensitive to it (Section 4.6), so using it
## to match the response level does not distort the answer.
fit_model <- function(anchor_dose, rate_common, rate_rare, lam_common,
                      lam_rare, p_common) {
  gap_obs  <- rate_common - rate_rare
  mean_obs <- p_common * rate_common + (1 - p_common) * rate_rare

  rates_at <- function(delta, sigma_b) {
    p  <- list(tau = delta, kappa = 0, sigma_b = sigma_b,
               lam_common = lam_common, lam_rare = lam_rare,
               p_common = p_common, phi = 0, cval = 0.35)
    st <- strata(p)
    a  <- calibrate(st, 0, p$cval, sigma_b, target = anchor_dose)
    list(p = p, logA = a, r = group_rates(p, a, log(anchor_dose)))
  }
  solve_delta <- function(sigma_b) {
    if (gap_obs <= 0) return(0)
    f <- function(d) { r <- rates_at(d, sigma_b)$r; r[1] - r[2] - gap_obs }
    if (f(2.5) < 0) return(2.5)
    out <- try(uniroot(f, c(0, 2.5), tol = 1e-5)$root, silent = TRUE)
    if (inherits(out, "try-error")) 0 else out
  }
  g <- function(sigma_b) {
    d <- solve_delta(sigma_b)
    r <- rates_at(d, sigma_b)$r
    p_common * r[1] + (1 - p_common) * r[2] - mean_obs
  }
  lo <- 0.15; hi <- 1.6
  clamp <- NULL
  if (g(lo) < 0) { sig <- lo; clamp <- "high" }
  else if (g(hi) > 0) { sig <- hi; clamp <- "low" }
  else {
    out <- try(uniroot(g, c(lo, hi), tol = 1e-4)$root, silent = TRUE)
    sig <- if (inherits(out, "try-error")) 0.8 else out
  }
  delta <- solve_delta(sig)
  fin <- rates_at(delta, sig)
  list(delta = delta, sigma_b = sig, clamp = clamp,
       rate_common = fin$r[1], rate_rare = fin$r[2],
       obs_common = rate_common, obs_rare = rate_rare,
       gap_obs = gap_obs, mean_obs = mean_obs)
}

## everything the Clinic tab needs, from investigator-supplied quantities
clinic_solve <- function(inp) {
  levels <- sort(unique(inp$levels))
  fit <- fit_model(inp$anchor, inp$rate_common, inp$rate_rare,
                   inp$lam_common, inp$lam_rare, inp$p_common)
  p <- list(tau = fit$delta, kappa = 0, sigma_b = fit$sigma_b,
            lam_common = inp$lam_common, lam_rare = inp$lam_rare,
            p_common = inp$p_common, phi = 0, cval = 0.35)
  st   <- strata(p)
  logA <- calibrate(st, 0, p$cval, p$sigma_b, target = inp$anchor)
  u    <- seq(log(min(levels) / 3), log(max(levels) * 3), length.out = 420)
  Q    <- objectives(u, st, logA, 0, p$cval, p$sigma_b)
  step <- u[2] - u[1]
  refine <- function(u0, which) {
    uu <- seq(u0 - 2 * step, u0 + 2 * step, length.out = 41)
    argmax_fine(uu, objectives(uu, st, logA, 0, p$cval, p$sigma_b)[[which]])
  }
  uB <- refine(argmax_fine(u, Q$QB), "QB")
  uA <- refine(argmax_fine(u, Q$QA), "QA")
  dB <- exp(uB); dA <- exp(uA)
  pick <- function(d) levels[which.min(abs(log(levels) - log(d)))]
  list(p = p, st = st, logA = logA, u = u, QB = Q$QB, QA = Q$QA,
       dB = dB, dA = dA, gap = 100 * (dA / dB - 1),
       levB = pick(dB), levA = pick(dA), levels = levels, fit = fit,
       rate_common = fit$rate_common, rate_rare = fit$rate_rare)
}

## the paragraph an investigator can paste into a protocol
protocol_text <- function(cl, choice, npat) {
  target  <- if (choice == "breadth") "Q_B" else "Q_A"
  wording <- if (choice == "breadth")
    paste("the dose maximising the expected total number of epitopes to which",
          "a response is induced across the trial population (the breadth",
          "estimand)")
  else
    paste("the dose maximising the expected per-epitope response probability",
          "for a typical epitope in a typical patient (the cluster-averaged",
          "estimand)")
  analysis <- if (choice == "breadth")
    paste("a binomial regression of the patient-level response count k_i on",
          "the number of epitopes administered m_i, with linear and quadratic",
          "terms in log dose. This weights each patient by m_i and therefore",
          "targets the stated estimand. An exchangeable generalized estimating",
          "equation will NOT be used, because the estimated intracluster",
          "correlation would silently move the target between the two",
          "estimands; if a marginal model is preferred, an independence working",
          "correlation will be specified, which recovers the count weighting",
          "exactly.")
  else
    paste("a regression of the per-patient mean response rate on linear and",
          "quadratic terms in log dose, giving every patient weight one, or",
          "equivalently a cluster-weighted analysis weighting epitope j of",
          "patient i by 1/m_i. This targets the stated estimand. An",
          "exchangeable generalized estimating equation will NOT be used, for",
          "the reason given above.")
  chosen <- if (choice == "breadth") cl$levB else cl$levA
  other  <- if (choice == "breadth") cl$levA else cl$levB

  paste0(
"ESTIMAND (dose optimisation)\n\n",
"The optimal biological dose is defined as ", wording, ", denoted ", target,
". This definition is stated because the efficacy outcome is observed once per\n",
"epitope and the number of epitopes administered is fixed by the patient's HLA\n",
"type, so it is informative; under these conditions the two natural population\n",
"summaries of the same patient-level dose-response are maximised at different\n",
"doses, and the reported optimal dose is not interpretable without stating\n",
"which summary is intended.\n\n",
"ANALYSIS\n\n",
"The primary dose-optimisation analysis will be ", analysis, "\n\n",
"The recommended dose will be the administered level maximising the fitted\n",
"curve. Under the projection used at design stage this is ", chosen,
" micrograms per peptide",
if (cl$levB != cl$levA)
  paste0("; the alternative estimand would have selected ", other,
         " micrograms, which is why the choice above is pre-specified.\n\n")
else
  paste0(". On the dose grid used here the alternative estimand would have\n",
         "selected the same level, but this is a property of this grid and not\n",
         "a general one.\n\n"),
"PRECISION\n\n",
"With ", npat, " patients the confidence interval for the optimal dose is\n",
"expected to span a wide range, and to contain the alternative estimand in a\n",
"large majority of replicates. The trial is therefore not expected to\n",
"adjudicate between the two estimands empirically, which is the reason the\n",
"choice is made here rather than at analysis.\n")
}

## ============================== UI ===========================================
BIG_CSS <- HTML("
  body, .content-wrapper, .main-sidebar { font-size: 16px; }
  .main-header .logo { font-size: 20px; font-weight: 700; letter-spacing: .2px; }
  .box-title { font-size: 19px !important; font-weight: 700 !important; }
  .small-box h3 { font-size: 40px !important; font-weight: 800 !important; }
  .small-box p  { font-size: 16px !important; font-weight: 600 !important; }
  .irs-grid-text { font-size: 12px; }
  .control-label { font-size: 15px !important; font-weight: 700 !important; }
  .lead-note { font-size: 17px; line-height: 1.6; }
  .guide h3 { font-weight: 800; margin-top: 26px; font-size: 22px; }
  .guide p, .guide li { font-size: 17px; line-height: 1.65; }
  .guide strong { font-weight: 800; }
  .callout { font-size: 17px; }
  table.gt { font-size: 16px; }
")

slider_box <- box(
  title = "Population and product", status = "primary", solidHeader = TRUE,
  width = 12,
  h4(strong("Informativeness of cluster size")),
  sliderInput("tau", HTML("&tau; &mdash; capacity gain in the common supertype"),
              0, 1, 0.5, 0.02, width = "100%"),
  sliderInput("kappa", HTML("&kappa; &mdash; supertype &rarr; epitope quality coupling"),
              0, 1.6, 0.67, 0.02, width = "100%"),
  hr(),
  h4(strong("Population")),
  sliderInput("sigma_b", HTML("&sigma;<sub>b</sub> &mdash; spread of immune capacity"),
              0.2, 1.4, 0.8, 0.05, width = "100%"),
  fluidRow(
    column(6, sliderInput("lam_common", "Mean panel, common supertype",
                          12, 28, 20, 1, width = "100%")),
    column(6, sliderInput("lam_rare", "Mean panel, rare supertype",
                          4, 20, 11, 1, width = "100%"))),
  sliderInput("p_common", "Prevalence of the common supertype",
              0.2, 0.9, 0.6, 0.02, width = "100%"),
  hr(),
  h4(strong("Antigenic competition")),
  sliderInput("phi", HTML("&phi; &mdash; competition on / off"),
              0, 1, 0, 0.25, width = "100%"),
  sliderInput("cval", "c — competition scale", 0, 1, 0.35, 0.05, width = "100%"),
  actionButton("reset", "Reset to the paper's reference configuration",
               icon = icon("rotate-left"), class = "btn-primary btn-lg")
)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Which optimal biological dose?", titleWidth = 380),
  dashboardSidebar(
    width = 240,
    sidebarMenu(
      menuItem("Clinic", tabName = "clinic", icon = icon("stethoscope")),
      menuItem("Explore", tabName = "explore", icon = icon("sliders")),
      menuItem("Sensitivity", tabName = "sens", icon = icon("table")),
      menuItem("How to use this", tabName = "guide", icon = icon("book-open")),
      menuItem("About", tabName = "about", icon = icon("circle-info"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(BIG_CSS)),
    tabItems(

      ## ------------------------------------------------------------- CLINIC --
      tabItem(
        "clinic",
        fluidRow(box(
          width = 12, status = "primary", solidHeader = TRUE,
          title = "The problem, in one picture",
          div(style = "text-align:center;", HTML(doodle_svg())),
          p(class = "lead-note", style = "margin-top:12px;",
            "Same two patients. Same responses. ",
            strong("Counting epitopes picks the lower dose; averaging patients
                    picks the higher one."),
            " Neither is wrong. Your protocol has to say which one you mean —
             otherwise the dose you carry forward depends on who runs the
             analysis.")
        )),

        fluidRow(column(
          width = 12,
          div(style = "max-width:900px;margin:0 auto;",

            box(width = 12, status = "success", solidHeader = TRUE,
                title = "Five short steps",
                p(class = "lead-note",
                  "One question per screen. Each asks for something you can read
                   off a cohort you have already run — no statistical
                   parameters. At the end you get a written summary to take to
                   your trial statistician, so the two of you are deciding the
                   same thing on purpose rather than discovering it afterwards."),
                p(class = "lead-note",
                  strong("In a hurry?"), " Load an example and change only what
                   differs from your trial."),
                div(
                  actionButton("preset_neo", "Personalised neoantigen (≈20 peptides)",
                               icon = icon("dna"), class = "btn-primary"),
                  actionButton("preset_panel", "Fixed panel, HLA-selected",
                               icon = icon("list-check"), class = "btn-primary"),
                  actionButton("preset_small", "Small pilot, narrow dose range",
                               icon = icon("flask"), class = "btn-primary"))),

            box(width = 12, status = "warning", solidHeader = TRUE,
                title = uiOutput("w_title", inline = TRUE),
                uiOutput("w_progress"),

                ## ---- step 1 ------------------------------------------------
                conditionalPanel(
                  "output.step == 1",
                  h3(style = "font-weight:800;margin-top:4px;",
                     "Which dose do you actually want?"),
                  p(class = "lead-note",
                    "There are two defensible answers and they are not the same
                     dose. Pick the one that matches what the product is for."),
                  radioButtons(
                    "clin_choice", NULL, width = "100%",
                    choiceNames = list(
                      HTML("<b style='font-size:19px'>Cover the most epitopes.</b>
                            <br><span style='font-size:16px'>A patient given 20
                            peptides counts 20 times as much as one given 1.
                            Choose this if breadth of coverage across the
                            population is the goal.</span>"),
                      HTML("<b style='font-size:19px'>Best response per epitope.</b>
                            <br><span style='font-size:16px'>Every patient counts
                            equally, whatever their panel size. Choose this if
                            the question is how well a typical epitope performs
                            in a typical patient.</span>")),
                    choiceValues = c("breadth", "average"))),

                ## ---- step 2 ------------------------------------------------
                conditionalPanel(
                  "output.step == 2",
                  h3(style = "font-weight:800;margin-top:4px;",
                     "Which doses did you give?"),
                  p(class = "lead-note",
                    "The dose levels administered in your dose-finding cohort,
                     in micrograms per peptide, separated by commas. Two or
                     more."),
                  textInput("clin_levels", NULL, "100, 250, 600, 1500",
                            width = "100%"),
                  p(style = "font-size:15px;color:#666",
                    "Example: 100, 250, 600, 1500"),
                  br(),
                  sliderInput("clin_n", "How many patients were in that cohort?",
                              10, 300, 60, 5, width = "100%")),

                ## ---- step 3 ------------------------------------------------
                conditionalPanel(
                  "output.step == 3",
                  h3(style = "font-weight:800;margin-top:4px;",
                     "Which of those doses produced the most responding
                      epitopes in total?"),
                  p(class = "lead-note",
                    "Add up the responding epitopes across all patients at each
                     dose and pick the winner. Total count, not a rate — a
                     patient with many peptides contributes many. This anchors
                     the model to your data."),
                  selectInput("clin_anchor_sel", NULL,
                              choices = c(100, 250, 600, 1500), selected = 600,
                              width = "100%"),
                  p(style = "font-size:15px;color:#666",
                    "If two doses tie, pick the lower one.")),

                ## ---- step 4 ------------------------------------------------
                conditionalPanel(
                  "output.step == 4",
                  h3(style = "font-weight:800;margin-top:4px;",
                     "How many peptides did each patient get?"),
                  p(class = "lead-note",
                    "Split your patients into the commoner HLA group and the
                     rest. If you did not record HLA, use whatever divides your
                     patients into those who received more peptides and those
                     who received fewer."),
                  fluidRow(
                    column(6, numericInput(
                      "clin_lam_c", "Average peptides — commoner group",
                      20, 4, 30, width = "100%")),
                    column(6, numericInput(
                      "clin_lam_r", "Average peptides — other group",
                      11, 4, 30, width = "100%"))),
                  sliderInput("clin_pc",
                              "What share of your patients were in the commoner group?",
                              0.2, 0.9, 0.6, 0.05, width = "100%")),

                ## ---- step 5 ------------------------------------------------
                conditionalPanel(
                  "output.step == 5",
                  h3(style = "font-weight:800;margin-top:4px;",
                     "How well did each group respond?"),
                  p(class = "lead-note",
                    "At the dose you chose in step 3. If you have the numbers,
                     enter them; if not, describe it and the dashboard will use
                     a reasonable value."),
                  radioButtons("clin_mode", NULL, width = "100%", inline = TRUE,
                               choiceNames = list("I have the numbers",
                                                  "I'll describe it"),
                               choiceValues = c("rates", "rough"),
                               selected = "rough"),
                  conditionalPanel(
                    "input.clin_mode == 'rates'",
                    p(style = "font-size:16px;color:#555",
                      "Responding epitopes divided by epitopes given, averaged
                       over the patients in each group. A number between 0 and 1."),
                    fluidRow(
                      column(6, numericInput("clin_rate_c", "Commoner group",
                                             0.52, 0.02, 0.95, 0.01,
                                             width = "100%")),
                      column(6, numericInput("clin_rate_r", "Other group",
                                             0.44, 0.02, 0.95, 0.01,
                                             width = "100%")))),
                  conditionalPanel(
                    "input.clin_mode == 'rough'",
                    sliderInput("clin_rate_all",
                                "Roughly what fraction of epitopes responded, overall?",
                                0.05, 0.90, 0.48, 0.01, width = "100%"),
                    radioButtons(
                      "clin_gap_desc",
                      "Did the commoner HLA group respond better than the other?",
                      width = "100%",
                      choiceNames = list("No difference worth mentioning",
                                         "A little better",
                                         "Clearly better",
                                         "Much better"),
                      choiceValues = c("none", "small", "moderate", "large"),
                      selected = "moderate"))),

                hr(),
                fluidRow(
                  column(4, conditionalPanel(
                    "output.step > 1",
                    actionButton("w_back", "Back", icon = icon("arrow-left"),
                                 class = "btn-default btn-lg", width = "100%"))),
                  column(8,
                    conditionalPanel(
                      "output.step < 5",
                      actionButton("w_next", "Next", icon = icon("arrow-right"),
                                   class = "btn-primary btn-lg", width = "100%")),
                    conditionalPanel(
                      "output.step == 5",
                      actionButton("clin_go", "Work out my dose",
                                   icon = icon("play"),
                                   class = "btn-success btn-lg", width = "100%"))))
            )
          )
        )),

        conditionalPanel(
          "output.showres == 'yes'",
          fluidRow(
            column(width = 7,
              uiOutput("clin_verdict"),
              fluidRow(valueBoxOutput("clin_vb1", 6), valueBoxOutput("clin_vb2", 6)),
              box(title = "Where the two definitions put the optimum",
                  width = 12, status = "primary", solidHeader = TRUE,
                  plotOutput("clin_plot", height = "360px"))),
            column(width = 5,
              box(title = "Take this to your statistician", width = 12,
                  status = "success", solidHeader = TRUE,
                  p(class = "lead-note",
                    "Written from your own numbers: the estimand, the analysis
                     that matches it, the dose it implies, and why the trial
                     will not settle the question on its own."),
                  downloadButton("clin_dl", "Download the summary",
                                 class = "btn-success btn-lg"),
                  br(), br(),
                  verbatimTextOutput("clin_protocol")),
              box(title = "How well the model matched your numbers", width = 12,
                  status = "info", solidHeader = TRUE, collapsible = TRUE,
                  collapsed = TRUE, uiOutput("clin_fitnote")))
          ))
      ),

      ## ------------------------------------------------------------ EXPLORE --
      tabItem(
        "explore",
        fluidRow(
          valueBoxOutput("vb_dB", width = 3),
          valueBoxOutput("vb_dA", width = 3),
          valueBoxOutput("vb_gap", width = 3),
          valueBoxOutput("vb_level", width = 3)
        ),
        fluidRow(
          column(
            width = 8,
            box(title = "The two population objectives", status = "primary",
                solidHeader = TRUE, width = 12,
                p(class = "lead-note",
                  "Each curve is normalised to its own maximum, so what you are",
                  "comparing is", strong("where they peak"), "and not how high",
                  "they rise. Antigen potency is recalibrated on every change so",
                  "that the breadth optimum stays at 396 µg."),
                plotOutput("p_obj", height = "420px")),
            box(title = "Three patients", status = "info", solidHeader = TRUE,
                width = 12,
                p(class = "lead-note",
                  "Expected breadth for patients with different panel sizes,",
                  "each normalised to its own maximum. Larger panels peak at",
                  "lower doses — that is what makes the two objectives",
                  "separate."),
                plotOutput("p_pat", height = "340px")),
            box(title = "Separation against τ", status = "warning",
                solidHeader = TRUE, width = 12,
                p(class = "lead-note",
                  "Computed with competition switched off. The open circle is",
                  "your current setting. At τ = 0 the two estimands coincide",
                  "exactly, which is Proposition 1."),
                plotOutput("p_tau", height = "320px"))
          ),
          column(width = 4, slider_box)
        )
      ),

      ## -------------------------------------------------------- SENSITIVITY --
      tabItem(
        "sens",
        fluidRow(box(
          title = "One factor at a time", status = "primary", solidHeader = TRUE,
          width = 12,
          p(class = "lead-note",
            "Each parameter is varied about your current setting on the Explore",
            "tab, with antigen potency recalibrated in every row so the breadth",
            "optimum is held at 396 µg. This is Table 3 of the paper, computed",
            "live."),
          actionButton("run_sens", "Compute sensitivity table",
                       icon = icon("play"), class = "btn-success btn-lg"),
          br(), br(),
          tableOutput("t_sens"),
          div(class = "callout",
              tags$b("What to look at:"),
              " σ", tags$sub("b"), " barely moves the separation, while τ, the",
              " panel-size contrast and the quality coupling move it a great",
              " deal. Between-patient spread in immune capacity determines how",
              " much individualised dosing would buy; it is close to irrelevant",
              " to the gap between the two population estimands.")
        ))
      ),

      ## --------------------------------------------------------------- GUIDE --
      tabItem(
        "guide",
        fluidRow(box(
          width = 12, status = "primary", solidHeader = TRUE,
          title = "How to use this dashboard", div(class = "guide",

          h3("What problem this shows"),
          p("A multi-epitope vaccine gives each patient several peptides at",
            "once, and the number given is fixed by the patient's HLA type.",
            "Immunogenicity is scored epitope by epitope, so the efficacy",
            "outcome arrives as a count out of a denominator that differs",
            "between patients."),
          p("There are two natural ways to summarise that across a trial:"),
          tags$ul(
            tags$li(strong("Breadth, Q_B."), " Total the responses across all",
                    "epitopes. A patient with twenty epitopes counts twenty",
                    "times as heavily as one with a single epitope. This asks:",
                    tags$em("which dose covers the most epitopes across the",
                            "population?")),
            tags$li(strong("Cluster-averaged, Q_A."), " Average the per-epitope",
                    "response rate within each patient, then average across",
                    "patients equally. This asks:", tags$em("which dose works",
                    "best for a typical epitope in a typical patient?"))),
          p("Both are legitimate. They are maximised at", strong("different"),
            "doses, and no amount of data will reconcile them, because they",
            "are the maximisers of different functions."),

          h3("How to read the main plot"),
          p("The solid blue curve is the breadth objective, the dashed red",
            "curve is the cluster-averaged objective. Vertical lines mark each",
            "maximum. When the lines sit apart, a trial reporting",
            "\"the optimal biological dose\" is reporting a number that depends",
            "on which analysis was run."),

          h3("What each slider does"),
          tags$ul(
            tags$li(strong("τ"), " — how much more responsive the common HLA",
                    "supertype is. This is the master switch. At τ = 0 with",
                    "κ = 0 the two estimands coincide exactly."),
            tags$li(strong("κ"), " — how much better the epitopes are in the",
                    "common supertype, on top of the capacity gain."),
            tags$li(strong("σ_b"), " — spread of immune capacity between",
                    "patients. Try moving this across its whole range: the",
                    "separation hardly changes. That is a real finding, not a",
                    "bug."),
            tags$li(strong("Panel means"), " — average number of epitopes in",
                    "each supertype. Widening the gap between them widens the",
                    "separation."),
            tags$li(strong("φ and c"), " — antigenic competition, where",
                    "co-administered epitopes compete for the same immune",
                    "resource. Switching it on is slower to compute.")),

          h3("Using this at the protocol stage"),
          tags$ol(
            tags$li(strong("Decide which estimand you mean"), " before the",
                    "trial, and write it into the protocol. \"Optimal",
                    "biological dose\" on its own is not a specification."),
            tags$li(strong("Match the analysis to it."), " A binomial",
                    "regression on patient-level counts targets breadth. A",
                    "marginal regression on the per-patient rate, or a",
                    "cluster-weighted analysis, targets the cluster-averaged",
                    "summary. An exchangeable GEE lands somewhere between the",
                    "two, at a position the estimated intracluster correlation",
                    "chooses for you — which means it conceals the decision",
                    "rather than making it."),
            tags$li(strong("Check the dose grid."), " The",
                    tags$em("Dose level"), " box tells you whether the two",
                    "estimands would actually send different doses forward from",
                    "100 / 250 / 600 / 1500 µg. Often they land on the same",
                    "level, and then the choice is academic for that grid. Move",
                    "the sliders and it stops being academic."),
            tags$li(strong("Do not expect the trial to settle it."), " In the",
                    "paper's simulation the confidence interval for either dose",
                    "spans a 2.7-fold range at 60 patients, and the interval",
                    "around one estimand contains the other in more than four",
                    "replicates in five. The data cannot adjudicate; the",
                    "protocol has to.")),

          h3("Where to start"),
          p("If you are running the trial rather than analysing it, start on",
            "the", strong("Clinic"), "tab. It walks through five questions, one",
            "screen at a time, each about a cohort you have already run — your",
            "dose levels, which one covered the most epitopes, average panel",
            "sizes in each HLA group, and the response rate you saw in each. It",
            "returns the dose to carry forward, the analysis that matches your",
            "question, and a written summary."),
          p("That summary is written to be", tags$em("checked"), "— take it to",
            "your trial statistician and settle the estimand together before",
            "the protocol is finalised. The dashboard makes the choice visible",
            "and quantifies what it costs; it does not make it for you. This",
            "tab is for seeing how the answer moves once you want to look",
            "behind it."),

          h3("What this dashboard is not"),
          p("Every number here follows from the generative model in Section 4.1",
            "with parameters you have set. It is a tool for thinking about",
            "magnitude and direction before a trial, not an analysis of any",
            "real dataset. The parameters that matter most — the association",
            "between HLA supertype and panel size, and between supertype and",
            "response capacity — should be estimated from existing phase I",
            "immunogenicity series before this is used to justify a design.")
        ))))
      ,

      ## --------------------------------------------------------------- ABOUT --
      tabItem(
        "about",
        fluidRow(box(
          width = 12, status = "primary", solidHeader = TRUE, title = "About",
          div(class = "guide",
            p(strong("Which optimal biological dose? Estimand ambiguity under",
                     "informative cluster size, with application to",
                     "multi-epitope cancer vaccines.")),
            p("Atanu Bhattacharjee, Division of Population Health and Genomics,",
              "School of Medicine, University of Dundee."),
            h3("Method"),
            p("The objectives are evaluated exactly rather than simulated.",
              "Conditional on HLA supertype the latent shift is normal and the",
              "panel size is independent of it, so each objective is a finite",
              "sum over two strata; integration over epitope quality and over",
              "immune capacity uses Gauss–Hermite quadrature with",
              "Golub–Welsch nodes, matching the epiOBD package. Because there",
              "is no Monte Carlo, τ = 0 returns a separation of exactly zero."),
            h3("Reproducibility"),
            p("The same model, the six analyses, the manuscript figures and",
              "tables and the interval-estimation results are in the R scripts",
              "accompanying the paper: epiOBD_core.R, epiOBD_analyses.R,",
              "epiOBD_figures.R, epiOBD_tables.R, epiOBD_sensitivity.R and",
              "epiOBD_inference.R.")
          )))
      )
    )
  )
)

## ============================ SERVER =========================================
server <- function(input, output, session) {

  params <- reactive({
    list(tau = input$tau, kappa = input$kappa, sigma_b = input$sigma_b,
         lam_common = input$lam_common, lam_rare = input$lam_rare,
         p_common = input$p_common, phi = input$phi, cval = input$cval)
  })

  res <- reactive({
    p <- params()
    if (p$phi > 0)
      withProgress(message = "Computing with antigenic competition",
                   value = 0.5, solve_model(p))
    else solve_model(p)
  })

  observeEvent(input$reset, {
    updateSliderInput(session, "tau", value = 0.5)
    updateSliderInput(session, "kappa", value = 0.67)
    updateSliderInput(session, "sigma_b", value = 0.8)
    updateSliderInput(session, "lam_common", value = 20)
    updateSliderInput(session, "lam_rare", value = 11)
    updateSliderInput(session, "p_common", value = 0.6)
    updateSliderInput(session, "phi", value = 0)
    updateSliderInput(session, "cval", value = 0.35)
  })

  ## ---- value boxes ----------------------------------------------------------
  output$vb_dB <- renderValueBox({
    valueBox(sprintf("%.0f µg", res()$dB), HTML("d*<sub>B</sub> &mdash; breadth"),
             icon = icon("layer-group"), color = "blue")
  })
  output$vb_dA <- renderValueBox({
    valueBox(sprintf("%.0f µg", res()$dA),
             HTML("d*<sub>A</sub> &mdash; cluster-averaged"),
             icon = icon("user"), color = "red")
  })
  output$vb_gap <- renderValueBox({
    g <- res()$gap
    if (abs(g) < 0.05) g <- 0
    valueBox(sprintf("%.1f%%", g),
             if (abs(g) < 0.05) "the estimands coincide" else "separation",
             icon = icon("arrows-left-right"),
             color = if (abs(g) < 0.05) "green" else if (g < 10) "yellow" else "orange")
  })
  output$vb_level <- renderValueBox({
    r <- res()
    lB <- nearest_level(r$dB); lA <- nearest_level(r$dA)
    same <- lB == lA
    valueBox(if (same) sprintf("%g µg", lB) else sprintf("%g / %g", lB, lA),
             if (same) "both estimands pick this dose level"
             else "the estimands pick different levels",
             icon = icon("flask"), color = if (same) "aqua" else "purple")
  })

  ## ---- plots ---------------------------------------------------------------
  base_par <- function()
    par(mar = c(4.6, 4.8, 1.0, 1.2), mgp = c(3.0, 0.8, 0), las = 1,
        cex.axis = 1.15, cex.lab = 1.3, font.lab = 2, bty = "n", col.axis = COL_INK)

  output$p_obj <- renderPlot({
    r <- res()
    keep <- r$u >= log(140) & r$u <= log(1250)
    nB <- max(r$QB); nA <- max(r$QA)
    yb <- range(c(r$QB[keep] / nB, r$QA[keep] / nA))
    base_par()
    plot(exp(r$u[keep]), r$QB[keep] / nB, type = "n", log = "x", xaxt = "n",
         ylim = c(max(yb[1], 0.80), 1.005),
         xlab = "Dose per peptide (µg)", ylab = "Objective (own maximum = 1)")
    at <- c(150, 250, 400, 600, 1000)
    axis(1, at = at, labels = at)
    abline(h = seq(0.8, 1, by = 0.02), col = COL_GRID, lwd = 1)
    abline(v = at, col = COL_GRID, lwd = 1)
    abline(v = r$dB, col = COL_B, lwd = 2, lty = 3)
    abline(v = r$dA, col = COL_A, lwd = 2, lty = 3)
    lines(exp(r$u[keep]), r$QB[keep] / nB, col = COL_B, lwd = 4)
    lines(exp(r$u[keep]), r$QA[keep] / nA, col = COL_A, lwd = 4, lty = 2)
    legend("bottom", horiz = TRUE, bty = "n", cex = 1.15, seg.len = 2.4,
           text.font = 2, lwd = 4, lty = c(1, 2), col = c(COL_B, COL_A),
           legend = c(expression(bold(Q[B]~"breadth")),
                      expression(bold(Q[A]~"cluster-averaged"))))
    mtext(sprintf("d*B = %.0f", r$dB), side = 3, adj = 0, col = COL_B,
          font = 2, cex = 1.15)
    mtext(sprintf("d*A = %.0f  (%+.1f%%)", r$dA, r$gap), side = 3, adj = 1,
          col = COL_A, font = 2, cex = 1.15)
  })

  output$p_pat <- renderPlot({
    r <- res()
    keep <- r$u >= log(90) & r$u <= log(1900)
    nB <- max(r$QB)
    base_par()
    plot(exp(r$u[keep]), r$QB[keep] / nB, type = "n", log = "x", xaxt = "n",
         ylim = c(0.70, 1.02), xlab = "Dose per peptide (µg)",
         ylab = "Expected breadth (own maximum = 1)")
    at <- c(100, 200, 400, 800, 1600); axis(1, at = at, labels = at)
    abline(v = at, col = COL_GRID); abline(h = seq(0.7, 1, 0.1), col = COL_GRID)
    lines(exp(r$u[keep]), r$QB[keep] / nB, col = COL_POP, lwd = 8)
    cases <- list(
      list(m = round(r$st[[1]]$mbar * 1.3), s = r$st[[1]], col = "#0b4f8a"),
      list(m = round(r$st[[1]]$mbar * 0.85), s = r$st[[1]], col = "#f0a202"),
      list(m = max(4, round(r$st[[2]]$mbar * 0.7)), s = r$st[[2]], col = "#2e8b57"))
    labs <- character(0)
    for (i in seq_along(cases)) {
      cc <- cases[[i]]
      br <- cc$m * Ffun(r$u[keep] + r$logA + cc$s$shift)
      lines(exp(r$u[keep]), br / max(br), col = cc$col, lwd = 3.5)
      points(exp(r$u[keep])[which.max(br)], 1, pch = 19, col = cc$col, cex = 1.4)
      labs[i] <- sprintf("m = %d (%s)", cc$m,
                         if (cc$s$S == 1) "common" else "rare")
    }
    legend("bottom", horiz = TRUE, bty = "n", cex = 1.1, lwd = 3.5,
           text.font = 2, col = vapply(cases, function(z) z$col, ""),
           legend = labs)
  })

  output$p_tau <- renderPlot({
    p <- params(); p$phi <- 0
    tg <- seq(0, 1, length.out = 13)
    gaps <- vapply(tg, function(t) { q <- p; q$tau <- t; solve_model(q)$gap }, 0)
    cur <- { q <- p; solve_model(q)$gap }
    base_par()
    plot(tg, gaps, type = "n", xlab = expression(bold(tau)),
         ylab = "Separation (%)", ylim = range(c(0, gaps)) * c(1, 1.08))
    abline(h = pretty(c(0, gaps)), col = COL_GRID)
    lines(tg, gaps, col = COL_A, lwd = 4)
    points(tg, gaps, pch = 19, col = COL_A, cex = 1.1)
    points(p$tau, cur, pch = 21, bg = "white", col = COL_INK, cex = 2.4, lwd = 3)
  })

  ## ---- clinic wizard --------------------------------------------------------
  wstep <- reactiveVal(1L)
  observeEvent(input$w_next, wstep(min(5L, wstep() + 1L)))
  observeEvent(input$w_back, wstep(max(1L, wstep() - 1L)))

  output$step <- renderText(as.character(wstep()))
  outputOptions(output, "step", suspendWhenHidden = FALSE)

  shown <- reactiveVal("no")
  observeEvent(input$clin_go, shown("yes"))
  output$showres <- renderText(shown())
  outputOptions(output, "showres", suspendWhenHidden = FALSE)

  STEP_TITLES <- c("Step 1 of 5 — what you are optimising",
                   "Step 2 of 5 — your dose levels",
                   "Step 3 of 5 — your best dose so far",
                   "Step 4 of 5 — peptides per patient",
                   "Step 5 of 5 — response rates")
  output$w_title <- renderUI(HTML(STEP_TITLES[wstep()]))
  output$w_progress <- renderUI({
    pct <- 100 * wstep() / 5
    div(style = "margin:2px 0 14px 0;",
        div(style = "height:10px;background:#e4e8ee;border-radius:5px;",
            div(style = sprintf(
              "height:10px;width:%.0f%%;background:#00a65a;border-radius:5px;", pct))))
  })

  ## ---- clinic ---------------------------------------------------------------
  ## the dose menu follows whatever the investigator typed in Step 2
  clin_levels <- reactive({
    req(input$clin_levels)
    lv <- suppressWarnings(as.numeric(trimws(strsplit(input$clin_levels, ",")[[1]])))
    sort(unique(lv[is.finite(lv) & lv > 0]))
  })
  observeEvent(clin_levels(), {
    lv <- clin_levels()
    if (length(lv) >= 2) {
      keep <- if (!is.null(input$clin_anchor_sel) &&
                  as.numeric(input$clin_anchor_sel) %in% lv)
        input$clin_anchor_sel else lv[max(1, round(length(lv) / 2))]
      updateSelectInput(session, "clin_anchor_sel", choices = lv, selected = keep)
    }
  })

  ## three starting points, so nobody faces an empty form
  preset <- function(levels, anchor, lc, lr, pc, rate, gap, n) {
    updateTextInput(session, "clin_levels", value = paste(levels, collapse = ", "))
    updateSelectInput(session, "clin_anchor_sel", choices = levels,
                      selected = anchor)
    updateNumericInput(session, "clin_lam_c", value = lc)
    updateNumericInput(session, "clin_lam_r", value = lr)
    updateSliderInput(session, "clin_pc", value = pc)
    updateSliderInput(session, "clin_n", value = n)
    updateRadioButtons(session, "clin_mode", selected = "rough")
    updateSliderInput(session, "clin_rate_all", value = rate)
    updateRadioButtons(session, "clin_gap_desc", selected = gap)
    wstep(1L)
  }
  observeEvent(input$preset_neo,
    preset(c(100, 250, 600, 1500), 600, 20, 11, 0.6, 0.48, "moderate", 60))
  observeEvent(input$preset_panel,
    preset(c(150, 300, 600, 1200), 600, 14, 8, 0.55, 0.40, "small", 80))
  observeEvent(input$preset_small,
    preset(c(400, 600, 800), 600, 22, 9, 0.6, 0.55, "large", 30))

  ## response rates, however the investigator chose to supply them
  clin_rates <- reactive({
    req(input$clin_mode, input$clin_pc)
    if (identical(input$clin_mode, "rates"))
      return(c(input$clin_rate_c, input$clin_rate_r))
    g <- switch(input$clin_gap_desc, none = 0, small = 0.05,
                moderate = 0.10, large = 0.18)
    p <- input$clin_rate_all; w <- input$clin_pc
    c(min(0.95, p + (1 - w) * g), max(0.02, p - w * g))
  })

  clin <- eventReactive(input$clin_go, ignoreNULL = FALSE, {
    req(input$clin_lam_c, input$clin_lam_r, input$clin_pc)
    lv <- clin_levels()
    validate(need(length(lv) >= 2, "Give at least two dose levels in Step 2."))
    anchor <- suppressWarnings(as.numeric(input$clin_anchor_sel))
    if (!is.finite(anchor)) anchor <- lv[max(1, round(length(lv) / 2))]
    r <- clin_rates()
    validate(need(all(r > 0 & r < 1),
                  "Response rates must be between 0 and 1."),
             need(r[1] >= r[2],
                  "The common HLA group should not respond worse than the other
                   group; if it does, swap the two groups over."))
    withProgress(message = "Fitting to your data", value = 0.5,
      clinic_solve(list(levels = lv, anchor = anchor,
                        lam_common = input$clin_lam_c,
                        lam_rare = input$clin_lam_r,
                        p_common = input$clin_pc,
                        rate_common = r[1], rate_rare = r[2])))
  })

  output$clin_dl <- downloadHandler(
    filename = function()
      paste0("dose-estimand-summary-", Sys.Date(), ".txt"),
    content = function(file) {
      cl <- clin()
      writeLines(c(
        "WHICH OPTIMAL BIOLOGICAL DOSE? -- design-stage summary",
        paste("Generated", format(Sys.time(), "%Y-%m-%d %H:%M")),
        strrep("-", 70), "",
        "YOUR INPUTS",
        paste("  Dose levels (ug/peptide):", paste(cl$levels, collapse = ", ")),
        paste("  Dose with the most responding epitopes in total:",
              input$clin_anchor_sel),
        paste("  Patients:", input$clin_n),
        sprintf("  Mean peptides: %g (common HLA group), %g (other group)",
                input$clin_lam_c, input$clin_lam_r),
        sprintf("  Share in the common group: %.0f%%", 100 * input$clin_pc),
        sprintf("  Per-epitope response: %.2f (common), %.2f (other)",
                cl$fit$obs_common, cl$fit$obs_rare),
        "",
        "RESULT",
        sprintf("  Best for coverage              : %.0f ug  -> dose level %g",
                cl$dB, cl$levB),
        sprintf("  Best for per-epitope response  : %.0f ug  -> dose level %g",
                cl$dA, cl$levA),
        sprintf("  Separation between them        : %.1f%%", cl$gap),
        if (cl$levA == cl$levB)
          "  Both definitions select the same administered level on this grid."
        else
          "  THE TWO DEFINITIONS SELECT DIFFERENT LEVELS ON THIS GRID.",
        "", strrep("-", 70), "",
        protocol_text(cl, input$clin_choice, input$clin_n),
        strrep("-", 70),
        "Model fitted to two quantities only: the gap in per-epitope response",
        "between HLA groups, and the overall response level. Projection for",
        "design; not an analysis of trial data."), file)
    })

  output$clin_verdict <- renderUI({
    cl <- clin()
    same <- cl$levB == cl$levA
    chosen <- if (input$clin_choice == "breadth") cl$levB else cl$levA
    other  <- if (input$clin_choice == "breadth") cl$levA else cl$levB
    box(width = 12, status = if (same) "success" else "danger",
        solidHeader = TRUE,
        title = if (same) "Your dose grid is not sensitive to the choice"
                else "The choice changes which dose you take forward",
        div(class = "guide",
          h3(sprintf("Take %g µg per peptide forward.", chosen)),
          if (same)
            p("Both definitions land on the same administered level here, so on",
              tags$b("this"), "grid the decision is not at risk. That is a",
              "property of where your dose levels happen to sit, not a general",
              "result — add or move a level and it can change. Pre-specify the",
              "estimand anyway, because the next trial's grid may differ.")
          else
            p("The other definition would have selected", tags$b(sprintf("%g µg", other)),
              "instead. Two competent statisticians analysing your data could",
              "hand you different doses, and neither would be wrong. This is why",
              "the estimand has to be written into the protocol before the",
              "analysis, not chosen afterwards."),
          p(class = "callout",
            tags$b("Use this analysis: "),
            if (input$clin_choice == "breadth")
              "binomial regression on patient-level counts (k_i out of m_i), quadratic in log dose."
            else
              "regression on the per-patient mean response rate, or a cluster-weighted analysis with weight 1/m_i.")))
  })

  output$clin_vb1 <- renderValueBox({
    cl <- clin()
    valueBox(sprintf("%.0f µg", cl$dB),
             HTML("Best for <b>coverage</b> (unrounded)"),
             icon = icon("layer-group"), color = "blue")
  })
  output$clin_vb2 <- renderValueBox({
    cl <- clin()
    valueBox(sprintf("%.0f µg", cl$dA),
             HTML("Best for <b>per-epitope response</b> (unrounded)"),
             icon = icon("user"), color = "red")
  })

  output$clin_plot <- renderPlot({
    cl <- clin()
    keep <- cl$u >= log(min(cl$levels) * 0.7) & cl$u <= log(max(cl$levels) * 1.4)
    nB <- max(cl$QB); nA <- max(cl$QA)
    base_par()
    plot(exp(cl$u[keep]), cl$QB[keep] / nB, type = "n", log = "x", xaxt = "n",
         ylim = c(min(c(cl$QB[keep] / nB, cl$QA[keep] / nA)), 1.02),
         xlab = "Dose per peptide (µg)", ylab = "Objective (own maximum = 1)")
    axis(1, at = cl$levels, labels = cl$levels)
    abline(v = cl$levels, col = COL_GRID, lwd = 2)
    rug(cl$levels, side = 1, lwd = 4, ticksize = 0.04, col = COL_INK)
    lines(exp(cl$u[keep]), cl$QB[keep] / nB, col = COL_B, lwd = 4)
    lines(exp(cl$u[keep]), cl$QA[keep] / nA, col = COL_A, lwd = 4, lty = 2)
    abline(v = cl$dB, col = COL_B, lwd = 2, lty = 3)
    abline(v = cl$dA, col = COL_A, lwd = 2, lty = 3)
    legend("bottom", horiz = TRUE, bty = "n", cex = 1.15, seg.len = 2.4,
           text.font = 2, lwd = 4, lty = c(1, 2), col = c(COL_B, COL_A),
           legend = c("coverage", "per-epitope response"))
    mtext("| = your dose levels", side = 3, adj = 0, cex = 1.05, font = 2,
          col = COL_INK)
  })

  output$clin_protocol <- renderText({
    protocol_text(clin(), input$clin_choice, input$clin_n)
  })

  output$clin_fitnote <- renderUI({
    cl <- clin(); f <- cl$fit
    div(class = "guide",
      p("Two quantities were fitted to your numbers: the latent difference",
        "between HLA groups, set by the", tags$b("gap"), "in response rate, and",
        "the spread of immune capacity, set by the overall", tags$b("level"),
        "of response. Nothing else was estimated."),
      tags$ul(
        tags$li(sprintf("Your figures: %.2f (common) and %.2f (less common).",
                        f$obs_common, f$obs_rare)),
        tags$li(sprintf("Model reproduces: %.2f and %.2f.",
                        f$rate_common, f$rate_rare)),
        tags$li(sprintf("Fitted separation parameter %.2f; capacity spread %.2f.",
                        f$delta, f$sigma_b))),
      if (!is.null(f$clamp))
        p(class = "callout", tags$b("Note: "),
          "your observed response level sits ", f$clamp,
          "er than this model can reach with a realistic spread of immune",
          " capacity, so the spread was held at its limit. The gap between",
          " groups is matched, and the gap is what drives the separation",
          " between the two optimal doses — but treat the absolute doses as",
          " indicative rather than exact."),
      p("The separation between the two optima is almost insensitive to the",
        "capacity spread, which is why fitting it to the response level does",
        "not distort the dose recommendation."))
  })

  ## ---- sensitivity table ----------------------------------------------------
  sens <- eventReactive(input$run_sens, {
    p <- params()
    grid <- list(
      list(nm = "sigma_b (capacity spread)", key = "sigma_b",
           vals = c(0.4, 0.6, 0.8, 1.0, 1.2)),
      list(nm = "tau (informativeness)", key = "tau",
           vals = c(0, 0.25, 0.5, 0.75, 1.0)),
      list(nm = "kappa (quality coupling)", key = "kappa",
           vals = c(0, 0.35, 0.67, 1.0, 1.35)),
      list(nm = "prevalence, common supertype", key = "p_common",
           vals = c(0.4, 0.5, 0.6, 0.7, 0.8)),
      list(nm = "mean panel, common supertype", key = "lam_common",
           vals = c(16, 18, 20, 22, 24))
    )
    out <- NULL
    withProgress(message = "Computing sensitivity", value = 0, {
      n <- sum(vapply(grid, function(g) length(g$vals), 0L)); i <- 0
      for (g in grid) for (v in g$vals) {
        q <- p; q[[g$key]] <- v; q$phi <- p$phi
        r <- solve_model(q)
        out <- rbind(out, data.frame(
          Parameter = if (v == g$vals[1]) g$nm else "",
          Value = format(v),
          `d*_A` = sprintf("%.0f", r$dA),
          `Separation (%)` = sprintf("%.1f", if (abs(r$gap) < 0.05) 0 else r$gap),
          `|obs - pred|` = sprintf("%.4f", abs(r$uA - r$uB - r$pred)),
          check.names = FALSE))
        i <- i + 1; incProgress(1 / n)
      }
    })
    out
  })

  output$t_sens <- renderTable(sens(), striped = TRUE, hover = TRUE,
                               bordered = TRUE, width = "100%", align = "llrrr")
}

shinyApp(ui, server)
