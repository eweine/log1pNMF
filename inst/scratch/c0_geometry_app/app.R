# Interactive companion to note_c0_geometry.tex: the achievable set of
# lambda-directions for p = 3 features, comparing
#
#   * log1p NMF, c -> 0+, with an INTERCEPT plus two structural factors:
#     the closure of the softmax family sigma(l1 f1 + l2 f2), l >= 0.
#     The shaded cloud shows the region reachable with structural
#     loadings bounded by T (the slider); as T grows the region fills
#     out toward its boundary at infinity -- the red frontier, which is
#     exactly the sparse segment set of the appendix theorem (segments
#     between basis vectors of adjacent upper-right-hull vertices).
#
#   * additive Poisson NMF (c = infinity) with the SAME two factors and
#     no intercept: normalized rates are convex combinations of the two
#     normalized factor columns -- the dark blue SEGMENT. Optionally add
#     a constant factor to the additive model too (toggle): the set
#     becomes the triangle spanned by the uniform vector and the two
#     normalized factors.
#
# Run with:  Rscript -e "shiny::runApp('inst/scratch/c0_geometry_app')"
# then open the printed http://127.0.0.1:PORT in a browser.

library(shiny)

## ---- geometry helpers -------------------------------------------------------
V <- rbind(c(0, 0), c(1, 0), c(1/2, sqrt(3)/2))   # e1, e2, e3 in the plane
proj <- function(lam) as.numeric(t(V) %*% (lam / sum(lam)))
softmax <- function(e) { e <- e - max(e); exp(e) / sum(exp(e)) }

## frontier of the theorem: sweep loading directions mu >= 0 and record which
## feature attains argmax_j mu . f_j; consecutive distinct argmaxes are
## adjacent upper-right-hull vertices, and the frontier is the union of
## simplex edges between them
frontier_pairs <- function(F) {
  th  <- seq(0, pi/2, length.out = 2001)
  amx <- sapply(th, function(a) which.max(cos(a) * F[, 1] + sin(a) * F[, 2]))
  seqv <- rle(amx)$values
  if (length(seqv) < 2) return(list(vertices = seqv, pairs = NULL))
  list(vertices = seqv,
       pairs = cbind(seqv[-length(seqv)], seqv[-1]))
}

## ---- ui ---------------------------------------------------------------------
num <- function(id, val) numericInput(id, NULL, val, min = 0, step = 0.1,
                                      width = "70px")
ui <- fluidPage(
  titlePanel("Achievable directions of λ (p = 3): log1p c → 0⁺ with intercept vs. additive NMF"),
  sidebarLayout(
    sidebarPanel(width = 3,
      h5(strong("Structural factors (rows = features)")),
      tags$table(
        tags$tr(tags$th(""), tags$th("f·1"), tags$th("f·2")),
        tags$tr(tags$td("j = 1"), tags$td(num("f11", 3)),   tags$td(num("f12", 0.5))),
        tags$tr(tags$td("j = 2"), tags$td(num("f21", 2)),   tags$td(num("f22", 2))),
        tags$tr(tags$td("j = 3"), tags$td(num("f31", 0.5)), tags$td(num("f32", 3)))
      ),
      checkboxInput("Tinf", "T = ∞: draw the FULL achievable set (exact boundary)", TRUE),
      sliderInput("Tmax", "Structural loading budget T (log1p cloud shows ||l|| ≤ T)",
                  min = 0.5, max = 40, value = 8, step = 0.5),
      checkboxInput("show_frontier", "Show c → 0⁺ frontier (theorem's sparse set)", TRUE),
      checkboxInput("show_additive", "Show additive model, two factors (segment)", TRUE),
      checkboxInput("add_const", "Give the additive model an intercept too (triangle; rank-fair comparison)", TRUE),
      helpText("The shaded cloud is the log1p (c → 0⁺) achievable set with an",
               "intercept: softmax(l₁f·1 + l₂f·2) over 0 ≤ l ≤ T.",
               "Increase T to watch it fill toward the red frontier -- the",
               "appendix theorem's segment set, i.e. the boundary at infinity.",
               "The dark blue segment is everything the 2-factor additive model",
               "can reach.")
    ),
    mainPanel(width = 9,
      fluidRow(
        column(7, plotOutput("simplex", height = "540px")),
        column(5, plotOutput("factorspace", height = "540px"))
      )
    )
  )
)

## ---- server -----------------------------------------------------------------
server <- function(input, output, session) {
  Fmat <- reactive({
    F <- rbind(c(input$f11, input$f12),
               c(input$f21, input$f22),
               c(input$f31, input$f32))
    validate(need(all(is.finite(F)) && all(F >= 0), "entries must be nonnegative"))
    F
  })

  output$simplex <- renderPlot({
    F <- Fmat(); T <- input$Tmax
    op <- par(mar = c(1, 1, 2.5, 1)); on.exit(par(op))
    plot(NA, xlim = c(-0.09, 1.09), ylim = c(-0.09, 0.96), asp = 1,
         axes = FALSE, xlab = "", ylab = "",
         main = "Achievable directions on the simplex")
    polygon(V[c(1, 2, 3), 1], V[c(1, 2, 3), 2], border = "grey60")
    text(V[1, 1] - 0.04, V[1, 2] - 0.03, expression(e[1]))
    text(V[2, 1] + 0.04, V[2, 2] - 0.03, expression(e[2]))
    text(V[3, 1],        V[3, 2] + 0.04, expression(e[3]))

    ## T = infinity: the full achievable set, drawn exactly. Its boundary is
    ## the single-factor curve sigma(t f1) (uniform -> e_argmax f1), the
    ## frontier edges through the upper-right-hull vertices, and the
    ## single-factor curve sigma(t f2) reversed (e_argmax f2 -> uniform).
    if (input$Tinf) {
      tt <- seq(0, 1, length.out = 300)^2 * 80
      c1 <- t(sapply(tt, function(t) proj(softmax(t * F[, 1]))))
      c2 <- t(sapply(tt, function(t) proj(softmax(t * F[, 2]))))
      fp <- frontier_pairs(F)
      mid <- V[fp$vertices, , drop = FALSE]
      poly <- rbind(c1, mid, c2[rev(seq_len(nrow(c2))), ])
      polygon(poly[, 1], poly[, 2], border = NA,
              col = adjustcolor("#9BB7BD", 0.45))
      lines(c1, col = "#0E5C6B", lwd = 2)
      lines(c2, col = "#0E5C6B", lwd = 2)
    }

    ## finite-T cloud: loading box [0, T]^2 (grid warped toward the corners
    ## so the frontier approach is visible)
    g  <- seq(0, 1, length.out = 90)^2 * T
    ll <- expand.grid(l1 = g, l2 = g)
    pts <- t(apply(ll, 1, function(l) proj(softmax(F %*% as.numeric(l)))))
    points(pts, pch = 16, cex = 0.35, col = adjustcolor("#5F8790", 0.35))
    u <- proj(rep(1, 3))
    points(u[1], u[2], pch = 16, cex = 1.1)
    text(u[1] - 0.09, u[2], "uniform", cex = 0.9)

    ## additive model(s), same structural factors
    f1 <- proj(F[, 1]); f2 <- proj(F[, 2])
    if (input$add_const) {
      polygon(c(u[1], f1[1], f2[1]), c(u[2], f1[2], f2[2]),
              border = "#0072B2", col = adjustcolor("#0072B2", 0.12), lwd = 2)
    }
    if (input$show_additive) {
      segments(f1[1], f1[2], f2[1], f2[2], col = "#0072B2", lwd = 3)
      points(rbind(f1, f2), pch = 17, col = "#0072B2", cex = 1.2)
      text(f1[1], f1[2] - 0.045, expression(hat(f)[1]), col = "#0072B2", cex = 0.9)
      text(f2[1], f2[2] - 0.045, expression(hat(f)[2]), col = "#0072B2", cex = 0.9)
    }

    ## the theorem's frontier
    if (input$show_frontier) {
      fp <- frontier_pairs(F)
      if (!is.null(fp$pairs))
        for (r in seq_len(nrow(fp$pairs)))
          segments(V[fp$pairs[r, 1], 1], V[fp$pairs[r, 1], 2],
                   V[fp$pairs[r, 2], 1], V[fp$pairs[r, 2], 2],
                   col = "#C2410C", lwd = 4)
      points(V[unique(fp$vertices), , drop = FALSE], pch = 15,
             col = "#C2410C", cex = 1.2)
    }

    legend("topleft", bty = "n", cex = 0.85,
           legend = c("log1p, c → 0⁺: cloud = ||l|| ≤ T; light fill = T = ∞",
                      "c → 0⁺ frontier (theorem)",
                      "additive, 2 factors",
                      if (input$add_const) "additive + intercept"),
           col = c(adjustcolor("#9BB7BD", 0.9), "#C2410C", "#0072B2",
                   if (input$add_const) adjustcolor("#0072B2", 0.4)),
           pch = c(16, 15, NA, if (input$add_const) 15),
           lty = c(NA, 1, 1, if (input$add_const) NA),
           lwd = c(NA, 4, 3, if (input$add_const) NA))
  })

  output$factorspace <- renderPlot({
    F <- Fmat()
    op <- par(mar = c(4, 4, 2.5, 1)); on.exit(par(op))
    lim <- c(0, max(F) * 1.15 + 0.1)
    plot(F[, 1], F[, 2], pch = 16, cex = 1.4, xlim = lim, ylim = lim,
         asp = 1, xlab = expression(f[j1]), ylab = expression(f[j2]),
         main = "Structural rows & upper right hull")
    text(F[, 1], F[, 2], labels = paste0("j=", 1:3), pos = 4, offset = 0.5)
    fp <- frontier_pairs(F)
    vs <- fp$vertices
    if (length(vs) >= 2)
      lines(F[vs, 1], F[vs, 2], col = "#C2410C", lwd = 2)
    points(F[unique(vs), 1], F[unique(vs), 2], pch = 15, col = "#C2410C",
           cex = 1.2)
    legend("topright", bty = "n", cex = 0.85, lwd = 2, col = "#C2410C",
           legend = "upper right hull (vertices = frontier features)")
  })
}

shinyApp(ui, server)
