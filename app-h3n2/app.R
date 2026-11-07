# ==========================================================================
# BIO463 - Molecular Epidemiology of H3N2
# Interactive companion to Molecular-Epi-H3N2.Rmd
#
# All phylogenetics (NJ tree, treedater clock, lineages-through-time) is
# PRECOMPUTED by prep/02_prepare_h3n2.R, so this app needs only
# shiny, bslib, ggplot2 and dplyr and can run in the browser via shinylive.
# ==========================================================================

library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)

CASES <- readRDS("data/h3n2_cases.rds")
TREE  <- readRDS("data/h3n2_tree.rds")
SUM   <- readRDS("data/summary.rds")

# ---- matte green palette (shared with the norovirus app) ----------------
GREEN_DARK  <- "#3F5A47"
GREEN       <- "#5A7A62"
GREEN_LIGHT <- "#8CA393"
CREAM       <- "#F7F8F4"
INK         <- "#2E3A32"
MUTED       <- "#7C8A80"
GRID        <- "#E2E6DF"
CLAY        <- "#B4756B"   # southern hemisphere / contrast series
BLUE        <- "#6C8CA8"   # lineages

HEMI_COLS <- c(Northern = GREEN, Southern = CLAY)

theme_matte <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      plot.background  = element_rect(fill = CREAM, colour = NA),
      panel.background = element_rect(fill = CREAM, colour = NA),
      panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      text             = element_text(colour = INK),
      axis.text        = element_text(colour = MUTED, size = base - 2),
      plot.title       = element_text(face = "bold", colour = GREEN_DARK,
                                      size = base + 1),
      plot.subtitle    = element_text(colour = MUTED, size = base - 1),
      legend.position  = "bottom",
      legend.key.size  = unit(10, "pt")
    )
}

sci <- function(x, d = 2) formatC(x, format = "e", digits = d)

DATE_MIN <- min(CASES$week)
DATE_MAX <- max(CASES$week)

# ==========================================================================
ui <- page_navbar(
  title = "Molecular Epidemiology of H3N2",
  id = "nav",
  theme = bs_theme(
    version = 5, bg = "#FBFCFA", fg = INK, primary = GREEN,
    base_font = font_google("Inter", local = FALSE),
    heading_font = font_google("Inter", local = FALSE)
  ),
  header = tags$head(tags$style(HTML(sprintf("
    .navbar { background-color: %s !important; border:0; }
    .navbar .navbar-brand, .navbar .nav-link { color: #EDF1EA !important; }
    .navbar .nav-link.active { color:#FFF !important; font-weight:600;
      border-bottom:2px solid #C9D6C5; }
    body { background:#FBFCFA; }
    .card { border:1px solid %s; border-radius:6px; box-shadow:none; }
    .lead-note { color:%s; font-size:.95rem; line-height:1.6; max-width:60rem; }
    h4,h5 { color:%s; font-weight:620; }
    .q-box { background:%s; border-left:3px solid %s; padding:.85rem 1.1rem;
      border-radius:4px; margin-bottom:.75rem; color:%s; font-size:.93rem;
      line-height:1.55; }
    .form-label, label { font-size:.85rem; color:%s; font-weight:550; }
    .irs--shiny .irs-bar, .irs--shiny .irs-handle>i:first-child { background:%s; }
    .irs--shiny .irs-single, .irs--shiny .irs-from, .irs--shiny .irs-to {
      background:%s; }
  ", GREEN_DARK, GRID, MUTED, CREAM, GRID, GREEN_DARK, MUTED, MUTED,
     GREEN_DARK, CREAM, GREEN_LIGHT, INK, GREEN_DARK, GREEN, GREEN_DARK)))),

  # ----------------------------------------------------------- overview
  nav_panel(
    "Overview",
    div(class = "container-fluid", style = "max-width:1100px;",
      h4("Reading an epidemic off a phylogeny"),
      p(class = "lead-note",
        "Seasonal influenza H3N2 circulates globally and reinvents its surface
         proteins fast enough that last year's immunity is only partly useful
         this year. That leaves two very different traces: case counts
         reported week by week through surveillance networks, and the genetic
         record carried in the virus itself."),
      p(class = "lead-note",
        sprintf("This app puts the two side by side. The genetic side is %d
                 hemagglutinin (HA) sequences from the northern hemisphere
                 sampled between %s and %s; the epidemiological side is %s
                 weeks of WHO sentinel surveillance data. The question to hold
                 on to: how much of the epidemic can you recover from the
                 sequences alone?",
                SUM$n_tip, substr(SUM$date_min, 1, 4), substr(SUM$date_max, 1, 4),
                format(SUM$n_weeks, big.mark = ","))),
      p(class = "lead-note",
        "Each step is built from the one before it, so the order matters. ",
        "Start with the surveillance record, because that is the epidemic as ",
        "the clinic sees it, and it is what the genetics will eventually be ",
        "judged against. Then turn to the sequences: the tree shows how these ",
        "viruses are related, and its shape alone says something about how ",
        "influenza persists. Putting that tree on a time axis is what makes a ",
        "rate measurable, so the molecular clock comes next. A dated tree can ",
        "then be read as a count of how many lineages were alive at each ",
        "moment — a measure of diversity that owes nothing to case ",
        "reports. Only at the end do the two independent records go on the ",
        "same axes."),
      br(),
      div(class = "q-box",
          strong("Start with the surveillance data. "),
          "Each tab has controls on the left and a short note under the ",
          "figure that interprets what you are looking at. Work through them ",
          "in order; the questions on the last tab refer back to things you ",
          "can change along the way."),
      div(style = "margin: 1.2rem 0 .4rem 0;",
          actionButton("go_cases", "Begin → the surveillance data",
                       class = "btn btn-primary")),
      p(class = "lead-note", style = "margin-top:1.4rem; font-size:.85rem;",
        sprintf("Sequences are HA segment 4 from GenBank; case data is WHO
                 sentinel surveillance. Built %s.", SUM$built))
    )
  ),

  # -------------------------------------------------------- surveillance
  nav_panel(
    "Surveillance data",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        checkboxGroupInput("c_hemi", "Hemisphere",
                           choices = c("Northern", "Southern"),
                           selected = c("Northern", "Southern")),
        radioButtons("c_metric", "Measure",
                     c("Influenza-positive specimens" = "positive",
                       "Percent of specimens positive" = "pct_positive",
                       "A(H3) subtype positives" = "h3",
                       "Specimens tested" = "tested"), "positive"),
        sliderInput("c_range", "Date range",
                    min = DATE_MIN, max = DATE_MAX,
                    value = c(DATE_MIN, DATE_MAX),
                    timeFormat = "%Y-%m"),
        checkboxInput("c_log", "Log scale (y)", FALSE),
        checkboxInput("c_points", "Show weekly points", FALSE),
        hr(),
        div(style = sprintf("font-size:.8rem; color:%s; line-height:1.5;", MUTED),
            "Sentinel surveillance: a fixed network of reporting sites that
             sample systematically, so counts track real activity rather than
             changes in who happened to get tested.")
      ),
      card(full_screen = TRUE,
           card_header(textOutput("c_title", inline = TRUE)),
           plotOutput("c_plot", height = "420px")),
      card(card_body(div(class = "lead-note", uiOutput("c_note")))),
      div(style = "margin:.2rem 0 1.2rem 0;",
          actionButton("go_tree", "Next \u2192 build the phylogeny",
                       class = "btn btn-primary"))
    )
  ),

  # ----------------------------------------------------------- phylogeny
  nav_panel(
    "Phylogeny",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        checkboxInput("t_colour", "Color tips by sampling year", TRUE),
        checkboxInput("t_points", "Show tips", TRUE),
        sliderInput("t_years", "Restrict tips shown to years",
                    min = min(TREE$years), max = max(TREE$years),
                    value = c(min(TREE$years), max(TREE$years)),
                    step = 1, sep = ""),
        hr(),
        div(style = sprintf("font-size:.8rem; color:%s; line-height:1.5;", MUTED),
            "A neighbor-joining tree scaled to time with a strict molecular
             clock. The x-axis is calendar time, so a tip's horizontal
             position is when that virus was sampled.")
      ),
      card(full_screen = TRUE,
           card_header(sprintf("Northern hemisphere H3N2 HA — dated phylogeny (%d sequences)",
                               SUM$n_tip)),
           plotOutput("t_plot", height = "520px")),
      card(card_body(div(class = "lead-note", uiOutput("t_note")))),
      div(style = "margin:.2rem 0 1.2rem 0;",
          actionButton("go_clock", "Next \u2192 fit the molecular clock",
                       class = "btn btn-primary"))
    )
  ),

  # --------------------------------------------------------------- clock
  nav_panel(
    "Molecular clock",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        checkboxInput("r_line", "Show regression line", TRUE),
        checkboxInput("r_colour", "Color by sampling year", FALSE),
        hr(),
        div(style = sprintf("font-size:.8rem; color:%s; line-height:1.5;", MUTED),
            "Root-to-tip regression: how far each tip sits from the root of
             the tree, against when it was sampled. The slope is the
             evolutionary rate.")
      ),
      card(full_screen = TRUE,
           card_header("Root-to-tip regression"),
           plotOutput("r_plot", height = "420px")),
      card(card_body(div(class = "lead-note", uiOutput("r_note")))),
      div(style = "margin:.2rem 0 1.2rem 0;",
          actionButton("go_ltt", "Next \u2192 count lineages over time",
                       class = "btn btn-primary"))
    )
  ),

  # ------------------------------------------------------------ lineages
  nav_panel(
    "Lineages over time",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        checkboxInput("l_log", "Log scale (y)", TRUE),
        sliderInput("l_range", "Date range",
                    min = DATE_MIN, max = DATE_MAX,
                    value = c(DATE_MIN, DATE_MAX),
                    timeFormat = "%Y-%m"),
        hr(),
        div(style = sprintf("font-size:.8rem; color:%s; line-height:1.5;", MUTED),
            "Lineages through time counts how many branches of the tree are
             alive at each moment. It is a measure of standing genetic
             diversity, and it rises when many distinct viruses are
             transmitting at once.")
      ),
      card(full_screen = TRUE,
           card_header("Lineages through time"),
           plotOutput("l_plot", height = "420px")),
      card(card_body(div(class = "lead-note", uiOutput("l_note")))),
      div(style = "margin:.2rem 0 1.2rem 0;",
          actionButton("go_both", "Next \u2192 compare with the cases",
                       class = "btn btn-primary"))
    )
  ),

  # ------------------------------------------------------ cases+lineages
  nav_panel(
    "Cases vs lineages",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        radioButtons("b_mode", "How to overlay",
                     c("Dual axis (as in the tutorial)" = "dual",
                       "Both rescaled to 0–1" = "norm"), "dual"),
        checkboxGroupInput("b_hemi", "Case series",
                           choices = c("Northern", "Southern"),
                           selected = c("Northern", "Southern")),
        sliderInput("b_range", "Date range",
                    min = DATE_MIN, max = DATE_MAX,
                    value = c(as.Date("2014-01-01"), DATE_MAX),
                    timeFormat = "%Y-%m"),
        hr(),
        div(style = sprintf("font-size:.8rem; color:%s; line-height:1.5;", MUTED),
            "The sequences and the case counts are independent data. If the
             genetic diversity curve tracks the epidemic curve, then sequence
             data alone carries epidemiological information.")
      ),
      card(full_screen = TRUE,
           card_header("H3N2 cases and viral lineages"),
           plotOutput("b_plot", height = "440px")),
      card(card_body(div(class = "lead-note", uiOutput("b_note")))),
      div(style = "margin:.2rem 0 1.2rem 0;",
          actionButton("go_q", "Next \u2192 the questions",
                       class = "btn btn-primary"))
    )
  ),

  # ----------------------------------------------------------- questions
  nav_panel(
    "Questions",
    div(class = "container-fluid", style = "max-width:900px;",
      h4("Work through these"),
      p(class = "lead-note",
        "Each question can be answered by changing something in the tabs."),
      br(),
      h5("1. Surveillance data"),
      div(class = "q-box",
          "On the Surveillance tab, show both hemispheres with ",
          strong("Influenza-positive specimens"), ". Describe the timing of ",
          "the two curves relative to each other. What causes that pattern, ",
          "and what does it imply for which hemisphere's season informs ",
          "vaccine strain selection for the other?"),
      div(class = "q-box",
          "Now switch to ", strong("Percent of specimens positive"),
          ". The shape changes. Why can raw positive counts mislead when ",
          "comparing across years or between hemispheres, and what does ",
          "dividing by specimens tested correct for?"),
      h5("2. The shape of the tree"),
      div(class = "q-box",
          "On the Phylogeny tab, look at the overall shape. Influenza HA ",
          "trees are described as ", strong("ladder-like"), ": a single trunk ",
          "with short side branches that die out. Contrast this with a tree ",
          "where many deep lineages persist for years. What does the ladder ",
          "tell you about immune escape and lineage turnover?"),
      div(class = "q-box",
          "Use the year slider to show only tips from a single season, then ",
          "widen it. Do sequences from one season cluster together on the ",
          "tree, or scatter across it? What would each pattern mean?"),
      h5("3. Rate of evolution"),
      div(class = "q-box",
          sprintf("The fitted clock rate is %s substitutions/site/year,
                   close to published estimates for influenza A HA of around
                   4 × 10⁻³. Compare it with the within-host norovirus rate
                   (2.1 × 10⁻²) from the other app in this course. Which is
                   faster, and does a faster per-site rate necessarily mean
                   faster antigenic change?", sci(SUM$rate))),
      div(class = "q-box",
          sprintf("The root-to-tip regression uses %d sequences over 16 years.
                   The trend is clearly upward, but R² is only %.2f, and the
                   regression slope (%s) and treedater's clock rate (%s) do not
                   agree. Give two reasons why individual tips scatter around
                   the clock line even when the virus really does evolve at a
                   constant average rate. Does that scatter undermine the dated
                   tree, or only the precision of the slope?",
                  SUM$n_tip, SUM$r2, sci(SUM$slope), sci(SUM$rate))),
      h5("4. Diversity and epidemiology"),
      div(class = "q-box",
          "On the Lineages tab, turn the log scale off and on. The curve rises ",
          "steeply at the right-hand edge. How much of that is real growth in ",
          "diversity, and how much is an artefact of recent sequences being ",
          "overrepresented in the sample?"),
      div(class = "q-box",
          "On the Cases vs lineages tab, compare the two overlay modes. The ",
          "dual-axis version requires choosing a scaling factor, which can be ",
          "set to make almost any two curves look aligned. Does the ",
          "rescaled version support the same conclusion? What would you need ",
          "in order to make this a quantitative rather than visual claim?"),
      div(class = "q-box",
          "If lineage diversity tracks the epidemic curve, sequence data can ",
          "substitute for case data in places where surveillance is weak. ",
          "Name one strength and one serious limitation of that substitution."),
      br(), br()
    )
  )
)

# ==========================================================================
server <- function(input, output, session) {

  # step-by-step navigation: each tab hands off to the next
  observeEvent(input$go_cases, nav_select("nav", "Surveillance data", session))
  observeEvent(input$go_tree,  nav_select("nav", "Phylogeny", session))
  observeEvent(input$go_clock, nav_select("nav", "Molecular clock", session))
  observeEvent(input$go_ltt,   nav_select("nav", "Lineages over time", session))
  observeEvent(input$go_both,  nav_select("nav", "Cases vs lineages", session))
  observeEvent(input$go_q,     nav_select("nav", "Questions", session))

  # ------------------------------------------------------ surveillance
  METRIC_LAB <- c(positive = "Influenza-positive specimens",
                  pct_positive = "Specimens positive (%)",
                  h3 = "A(H3) positive specimens",
                  tested = "Specimens tested")

  cases_f <- reactive({
    d <- CASES[CASES$week >= input$c_range[1] & CASES$week <= input$c_range[2], ]
    if (length(input$c_hemi)) d <- d[d$hemisphere %in% input$c_hemi, ]
    else d <- d[0, ]
    d$value <- d[[input$c_metric]]
    d <- d[!is.na(d$value), ]
    # a log axis cannot show zeros; drop them explicitly and report how many,
    # rather than letting ggplot discard them silently
    dropped <- 0L
    if (isTRUE(input$c_log)) {
      dropped <- sum(d$value <= 0)
      d <- d[d$value > 0, ]
    }
    attr(d, "dropped") <- dropped
    d
  })

  output$c_title <- renderText(METRIC_LAB[[input$c_metric]])

  output$c_plot <- renderPlot({
    d <- cases_f()
    if (!nrow(d)) {
      return(ggplot() + annotate("text", 0, 0, label = "Select a hemisphere",
                                 colour = MUTED, size = 5) +
               theme_void() +
               theme(plot.background = element_rect(fill = CREAM, colour = NA)))
    }
    p <- ggplot(d, aes(week, value, colour = hemisphere)) +
      geom_line(linewidth = 0.6)
    if (isTRUE(input$c_points)) p <- p + geom_point(size = 0.7, alpha = 0.55)
    if (isTRUE(input$c_log)) {
      p <- p + scale_y_log10()
    }
    p + scale_colour_manual(values = HEMI_COLS, name = NULL) +
      labs(x = NULL, y = METRIC_LAB[[input$c_metric]]) +
      theme_matte()
  })

  output$c_note <- renderUI({
    d <- cases_f()
    if (!nrow(d)) return(HTML("Nothing selected."))
    pk <- d[which.max(d$value), ]
    dropped <- attr(d, "dropped")
    HTML(paste0(sprintf(
      "Showing <b>%s</b> weekly records from %s to %s. The largest value in
       this window is <b>%s</b> (%s, week of %s). The northern and southern
       series peak about six months apart, which is why the WHO holds two
       vaccine composition meetings each year — one hemisphere's season is
       evidence about the strains the other will face.",
      format(nrow(d), big.mark = ","), min(d$week), max(d$week),
      format(round(pk$value, 1), big.mark = ","), pk$hemisphere, pk$week),
      if (!is.null(dropped) && dropped > 0) sprintf(
        "<br><br><span style='font-size:.88rem;'><b>%d</b> week%s in this
         window recorded zero, and a logarithmic axis cannot show zero, so
         those points are not plotted — which is why the lines break up.
         Worth noticing: the gaps are the inter-season troughs, and a log
         axis quietly removes exactly the weeks when nothing was
         circulating.</span>",
        dropped, if (dropped == 1) "" else "s") else ""))
  })

  # --------------------------------------------------------- phylogeny
  output$t_plot <- renderPlot({
    nd <- TREE$layout$nodes
    sg <- TREE$layout$segs
    tips <- nd[nd$is_tip, ]
    tips <- tips[!is.na(tips$year) &
                 tips$year >= input$t_years[1] &
                 tips$year <= input$t_years[2], ]

    p <- ggplot() +
      geom_segment(data = sg, aes(x = x, xend = xend, y = y, yend = yend),
                   colour = GREEN_DARK, linewidth = 0.22, lineend = "round")

    if (isTRUE(input$t_points) && nrow(tips)) {
      if (isTRUE(input$t_colour)) {
        p <- p + geom_point(data = tips, aes(x, y, colour = year), size = 1.1) +
          scale_colour_gradient(low = GREEN_LIGHT, high = CLAY, name = "Year")
      } else {
        p <- p + geom_point(data = tips, aes(x, y), colour = GREEN, size = 1.1)
      }
    }

    p + labs(x = "Year", y = NULL) +
      theme_matte() +
      theme(panel.grid.major.y = element_blank(),
            axis.text.y = element_blank(),
            legend.key.width = unit(34, "pt"))
  })

  output$t_note <- renderUI({
    HTML(sprintf(
      "The tree is strikingly <b>ladder-like</b>: one long trunk running left
       to right through time, with short side branches that appear and then
       stop. That is the signature of rapid lineage turnover under immune
       selection — each season's dominant viruses descend from the
       previous trunk, and the lineages that fail to escape existing immunity
       go extinct rather than persisting. The inferred common ancestor of
       this sample sits at about <b>%d</b>, and the %d sequences run to %s.",
      floor(TREE$tmrca), SUM$n_tip, substr(SUM$date_max, 1, 4)))
  })

  # ------------------------------------------------------------- clock
  output$r_plot <- renderPlot({
    d <- TREE$rtt
    p <- if (isTRUE(input$r_colour)) {
      ggplot(d, aes(decimal_date, distance, colour = year)) +
        geom_point(size = 1.4, alpha = 0.85) +
        scale_colour_gradient(low = GREEN_LIGHT, high = CLAY, name = "Year")
    } else {
      ggplot(d, aes(decimal_date, distance)) +
        geom_point(size = 1.4, alpha = 0.8, colour = GREEN)
    }
    if (isTRUE(input$r_line)) {
      p <- p + geom_abline(slope = TREE$slope, intercept = TREE$intercept,
                           colour = GREEN_DARK, linewidth = 0.6,
                           linetype = "22")
    }
    p + labs(x = "Sampling date (decimal years)",
             y = "Root-to-tip distance (subs/site)",
             subtitle = sprintf("slope = %s subs/site/year   ·   R² = %.3f",
                                sci(TREE$slope), TREE$r2)) +
      theme_matte() +
      theme(legend.key.width = unit(34, "pt"))
  })

  output$r_note <- renderUI({
    HTML(sprintf(
      "Each point is one sequence. The upward trend is clear — later samples
       sit further from the root, which is what licenses putting dates on the
       tree at all — but the scatter is substantial: R² = <b>%.2f</b>, so
       sampling time explains only about %.0f%% of the variation in
       root-to-tip distance. The regression slope is <b>%s</b>
       subs/site/year, while treedater's fitted clock rate is <b>%s</b> —
       the latter very close to the canonical influenza A HA figure of
       roughly 4 × 10⁻³.<br><br>
       <span style='font-size:.88rem;'>Why so much scatter, and why do the two
       estimates differ? Neighbour-joining distances are noisy; a ladder-like
       tree puts most tips on short side branches of uneven length; and
       sampling is uneven across years. The regression fits a line through the
       tips only, whereas treedater fits the whole tree, so they are not
       obliged to agree.</span>",
      TREE$r2, 100 * TREE$r2, sci(TREE$slope), sci(TREE$rate)))
  })

  # ---------------------------------------------------------- lineages
  ltt_f <- reactive({
    d <- TREE$ltt
    d[d$week >= input$l_range[1] & d$week <= input$l_range[2], ]
  })

  output$l_plot <- renderPlot({
    d <- ltt_f()
    if (!nrow(d)) return(NULL)
    p <- ggplot(d, aes(week, N)) +
      geom_step(colour = BLUE, linewidth = 0.9)
    if (isTRUE(input$l_log)) p <- p + scale_y_log10()
    p + labs(x = NULL, y = "Number of lineages") + theme_matte()
  })

  output$l_note <- renderUI({
    d <- ltt_f()
    HTML(sprintf(
      "At the right-hand edge <b>%d</b> lineages are counted as co-circulating.
       Read this curve carefully: a lineages-through-time plot built from a
       sample always rises steeply toward the present, partly because
       recent viruses are better represented and partly because recent
       branches have not yet had time to go extinct. The informative signal is
       in the <i>shape and timing of the wobbles</i>, not the absolute height.
       %s",
      max(d$N),
      if (isTRUE(input$l_log))
        "The log scale makes the earlier, smaller fluctuations visible; turn it
         off to see how much the recent rise dominates."
      else
        "Turn the log scale on to see the structure in the earlier years, which
         the recent rise otherwise flattens."))
  })

  # ----------------------------------------------------- cases+lineages
  output$b_plot <- renderPlot({
    cs <- CASES[CASES$week >= input$b_range[1] & CASES$week <= input$b_range[2], ]
    if (length(input$b_hemi)) cs <- cs[cs$hemisphere %in% input$b_hemi, ]
    else cs <- cs[0, ]
    lt <- TREE$ltt[TREE$ltt$week >= input$b_range[1] &
                   TREE$ltt$week <= input$b_range[2], ]
    if (!nrow(lt) || !nrow(cs)) {
      return(ggplot() + annotate("text", 0, 0, colour = MUTED, size = 5,
               label = "Nothing to show in this window") + theme_void() +
               theme(plot.background = element_rect(fill = CREAM, colour = NA)))
    }

    if (input$b_mode == "norm") {
      cs$v <- cs$positive / max(cs$positive, na.rm = TRUE)
      lt$v <- lt$N / max(lt$N, na.rm = TRUE)
      ggplot() +
        geom_line(data = cs, aes(week, v, colour = hemisphere),
                  linewidth = 0.6) +
        geom_step(data = lt, aes(week, v, colour = "Lineages"),
                  linewidth = 0.9) +
        scale_colour_manual(values = c(HEMI_COLS, Lineages = BLUE),
                            name = NULL) +
        labs(x = NULL, y = "Each series scaled to its own maximum") +
        theme_matte()
    } else {
      sf <- max(cs$positive, na.rm = TRUE) / max(lt$N, na.rm = TRUE)
      ggplot() +
        geom_line(data = cs, aes(week, positive, colour = hemisphere),
                  linewidth = 0.6) +
        geom_step(data = lt, aes(week, N * sf, colour = "Lineages"),
                  linewidth = 0.9) +
        scale_colour_manual(values = c(HEMI_COLS, Lineages = BLUE),
                            name = NULL) +
        scale_y_continuous(
          name = "Influenza-positive specimens",
          sec.axis = sec_axis(~ . / sf, name = "Number of lineages")) +
        labs(x = NULL) +
        theme_matte() +
        theme(axis.title.y.right = element_text(colour = BLUE))
    }
  })

  output$b_note <- renderUI({
    if (input$b_mode == "dual") {
      HTML("The lineage curve rises and falls roughly in step with reported
            cases: when many distinct viruses are transmitting, surveillance
            is finding many positives. That correspondence is the point of the
            tutorial — genomic data carries epidemiological signal, which
            is what structured coalescent and birth–death models exploit
            to estimate population size through time.<br><br>
            <span style='font-size:.88rem;'>But be skeptical of this particular
            figure. The right-hand axis was scaled so the two maxima coincide;
            a different scaling factor would tell a different visual story.
            Switch to the rescaled view and decide whether the agreement
            survives.</span>")
    } else {
      HTML("With both series divided by their own maximum, no scaling factor
            is doing any persuasive work. Judge the timing of the peaks rather
            than their heights. Where the curves genuinely agree, the genetic
            and epidemiological records are describing the same waves; where
            they diverge, ask whether the sequences for that period were
            sampled evenly — sampling effort, not biology, is the usual
            culprit.")
    }
  })
}

shinyApp(ui, server)
