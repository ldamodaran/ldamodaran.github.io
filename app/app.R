# ==========================================================================
# BIO463 - Within-Host Norovirus Evolution
# Interactive companion to Norovirus-Prolonged-Infection.Rmd
#
# All phylogenetics (NJ trees, treedater molecular clock, seqinr dN/dS) is
# PRECOMPUTED by prep/01_prepare_data.R. This app only reads small data
# files and plots them, so it runs entirely in the browser via shinylive.
#
# Packages here must stay limited to: shiny, bslib, ggplot2, dplyr
# ==========================================================================

library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)

# ---- data ---------------------------------------------------------------
ALN  <- readRDS("data/p33_alignment.rds")
TREE <- readRDS("data/p33_tree.rds")
DNDS <- readRDS("data/p33_dnds.rds")
SYD  <- readRDS("data/sydney_tree.rds")
SUM  <- readRDS("data/summary.rds")

META   <- ALN$meta
LABELS <- META$label
N_SEQ  <- nrow(META)

# ---- matte green palette ------------------------------------------------
GREEN_DARK  <- "#3F5A47"
GREEN       <- "#5A7A62"
GREEN_LIGHT <- "#8CA393"
CREAM       <- "#F7F8F4"
INK         <- "#2E3A32"
MUTED       <- "#7C8A80"
GRID        <- "#E2E6DF"

NT_COLS <- c(A = "#7E9B76", C = "#6C8CA8", G = "#C8A05E", T = "#B4756B",
             `-` = "#D8DAD4", N = "#BFC3BC")

# amino acids grouped by property, all desaturated
AA_COLS <- c(
  A = "#9FB497", V = "#9FB497", I = "#9FB497", L = "#9FB497", M = "#8AA384",
  F = "#C2A878", W = "#A98E63", P = "#CFD6C6",
  S = "#8FAFA2", T = "#A6C0B5", N = "#BFD2C9", Q = "#D5E0DA",
  R = "#7D96B2", K = "#95ABC2", H = "#AFC0D1",
  D = "#C08A84", E = "#D3A8A2",
  C = "#A892B0", G = "#D2D4CE", Y = "#7FA8B8",
  `*` = "#2E3A32", `-` = "#D8DAD4", X = "#E8EAE5"
)

theme_matte <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      plot.background   = element_rect(fill = CREAM, colour = NA),
      panel.background  = element_rect(fill = CREAM, colour = NA),
      panel.grid.major  = element_line(colour = GRID, linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      text              = element_text(colour = INK),
      axis.text         = element_text(colour = MUTED, size = base - 2),
      plot.title        = element_text(face = "bold", colour = GREEN_DARK,
                                       size = base + 1),
      plot.subtitle     = element_text(colour = MUTED, size = base - 1),
      legend.position   = "bottom",
      legend.key.size   = unit(10, "pt")
    )
}

sci <- function(x, d = 2) formatC(x, format = "e", digits = d)

# small stat card
stat_card <- function(value, label, note = NULL) {
  div(class = "stat-card",
      div(class = "stat-value", value),
      div(class = "stat-label", label),
      if (!is.null(note)) div(class = "stat-note", note))
}

# ---- long-format helpers ------------------------------------------------
melt_aln <- function(mat, cols_keep = NULL, value_name = "base") {
  if (!is.null(cols_keep)) mat <- mat[, cols_keep, drop = FALSE]
  pos <- if (is.null(cols_keep)) seq_len(ncol(mat)) else cols_keep
  df <- data.frame(
    taxon    = rep(rownames(mat), times = ncol(mat)),
    position = rep(pos, each = nrow(mat)),
    value    = as.vector(mat),
    stringsAsFactors = FALSE
  )
  df$taxon <- factor(df$taxon, levels = rev(rownames(mat)))
  names(df)[3] <- value_name
  df
}

# ==========================================================================
# UI
# ==========================================================================
ui <- page_navbar(
  title = "Within-Host Norovirus Evolution",
  id    = "nav",
  theme = bs_theme(
    version    = 5,
    bg         = "#FBFCFA",
    fg         = INK,
    primary    = GREEN,
    base_font  = font_google("Inter", local = FALSE),
    heading_font = font_google("Inter", local = FALSE)
  ),
  header = tags$head(tags$style(HTML(sprintf("
    .navbar { background-color: %s !important; border:0; }
    .navbar .navbar-brand, .navbar .nav-link { color: #EDF1EA !important; }
    .navbar .nav-link.active { color: #FFFFFF !important; font-weight: 600;
      border-bottom: 2px solid #C9D6C5; }
    body { background: #FBFCFA; }
    .card { border: 1px solid %s; border-radius: 6px; box-shadow: none; }
    .lead-note { color: %s; font-size: .95rem; line-height: 1.6;
      max-width: 60rem; }
    .stat-card { background: %s; border: 1px solid %s; border-radius: 6px;
      padding: .9rem 1rem; height: 100%%; }
    .stat-value { font-size: 1.45rem; font-weight: 650; color: %s;
      line-height: 1.15; }
    .stat-label { font-size: .78rem; text-transform: uppercase;
      letter-spacing: .06em; color: %s; margin-top: .25rem; }
    .stat-note { font-size: .78rem; color: %s; margin-top: .4rem; }
    h4, h5 { color: %s; font-weight: 620; }
    .q-box { background: %s; border-left: 3px solid %s;
      padding: .85rem 1.1rem; border-radius: 4px; margin-bottom: .75rem;
      color: %s; font-size: .93rem; line-height: 1.55; }
    table.dataframe, .table { font-size: .85rem; }
    .form-label, label { font-size: .85rem; color: %s; font-weight: 550; }
    .irs--shiny .irs-bar, .irs--shiny .irs-handle>i:first-child { background: %s; }
    .irs--shiny .irs-single, .irs--shiny .irs-from, .irs--shiny .irs-to {
      background: %s; }
  ", GREEN_DARK, GRID, MUTED, CREAM, GRID, GREEN_DARK, MUTED, MUTED,
     GREEN_DARK, CREAM, GREEN_LIGHT, INK, GREEN_DARK, GREEN, GREEN_DARK)))),

  # ---------------------------------------------------------------- overview
  nav_panel(
    "Overview",
    div(class = "container-fluid", style = "max-width:1100px;",
      h4("Tracking a virus as it evolves inside one person"),
      p(class = "lead-note",
        "Most of what we know about viral evolution comes from comparing ",
        "sequences collected from ", em("different"), " people. But a virus ",
        "also evolves inside a single host, and in a patient who cannot ",
        "clear the infection we can watch that happen directly."),
      p(class = "lead-note",
        sprintf("This app follows patient P33, who had a prolonged norovirus
                 GII.4 infection. %d capsid sequences were collected over
                 %d days. Nothing else about the outcome is given away here:
                 you will build it up one step at a time.",
                SUM$p33_n, SUM$p33_span_days)),
      p(class = "lead-note",
        "The order matters, because each step depends on the one before it. ",
        "First look at the alignment and find the handful of positions that ",
        "changed at all. Those differences are what a phylogeny is built ",
        "from, so the tree comes next. Putting the tree on a time axis is ",
        "what makes a rate measurable, so the molecular clock follows. And ",
        "once you know how fast change accumulates, you can ask the harder ",
        "question of whether that change is being favored or removed — ",
        "which is where dN/dS comes in."),
      br(),
      div(class = "q-box",
          strong("Start with the alignment. "),
          "Each tab has controls on the left and a short note under the ",
          "figure that interprets what you are looking at. Work through them ",
          "in order; the questions on the last tab refer back to things you ",
          "can change along the way."),
      div(style = "margin: 1.2rem 0 .4rem 0;",
          actionButton("go_align", "Begin → the alignment",
                       class = "btn btn-primary")),
      p(class = "lead-note", style = "margin-top:1.5rem; font-size:.85rem;",
        sprintf("Dates for P33 are shifted to a placeholder year to protect
                 patient privacy; only the intervals between samples are real.
                 Built %s.", SUM$built))
    )
  ),

  # --------------------------------------------------------------- alignment
  nav_panel(
    "Alignment",
    layout_sidebar(
      sidebar = sidebar(
        width = 290,
        radioButtons("aln_level", "Sequence type",
                     c("Nucleotide" = "nt", "Amino acid" = "aa"), "nt"),
        radioButtons("aln_view", "Show",
                     c("Whole gene" = "all",
                       "Variable sites only" = "var",
                       "Zoom to a window" = "win"), "all"),
        conditionalPanel(
          "input.aln_view == 'win'",
          uiOutput("win_slider"),
          helpText("Tip: the window 1020–1170 nt (codons 340–390) ",
                   "contains five variable sites.")
        ),
        checkboxInput("aln_outline", "Outline variable sites", TRUE),
        conditionalPanel(
          "input.aln_view != 'all'",
          checkboxInput("aln_letters", "Show letters in cells", TRUE)
        ),
        hr(),
        div(style = sprintf("font-size:.8rem; color:%s; line-height:1.5;", MUTED),
            "Rows are ordered earliest at top. Each column is one position; ",
            "a column of a single color means every time point agrees.")
      ),
      card(
        full_screen = TRUE,
        card_header(textOutput("aln_title", inline = TRUE)),
        plotOutput("aln_plot", height = "430px")
      ),
      card(
        card_body(
          div(class = "lead-note", uiOutput("aln_note"))
        )
      ),
      div(style = "margin:.2rem 0 1.2rem 0;",
          actionButton("go_tree", "Next \u2192 build the phylogeny",
                       class = "btn btn-primary"))
    )
  ),

  # --------------------------------------------------------------- phylogeny
  nav_panel(
    "Phylogeny",
    layout_sidebar(
      sidebar = sidebar(
        width = 290,
        radioButtons("tree_which", "Dataset",
                     c("P33 — within one patient" = "p33",
                       "Sydney 2012 — global" = "syd"), "p33"),
        conditionalPanel("input.tree_which == 'p33'",
                         checkboxInput("tree_tips", "Show tip labels", TRUE)),
        hr(),
        div(style = sprintf("font-size:.8rem; color:%s; line-height:1.5;", MUTED),
            "Both trees are neighbor-joining trees scaled to time with a ",
            "strict molecular clock (treedater). The x-axis is calendar time, ",
            "so the left-hand root is the inferred common ancestor.")
      ),
      card(
        full_screen = TRUE,
        card_header(textOutput("tree_title", inline = TRUE)),
        plotOutput("tree_plot", height = "480px")
      ),
      card(card_body(div(class = "lead-note", uiOutput("tree_note")))),
      div(style = "margin:.2rem 0 1.2rem 0;",
          actionButton("go_clock", "Next \u2192 fit the molecular clock",
                       class = "btn btn-primary"))
    )
  ),

  # ------------------------------------------------------------------- clock
  nav_panel(
    "Molecular clock",
    layout_sidebar(
      sidebar = sidebar(
        width = 290,
        radioButtons("rtt_which", "Dataset",
                     c("P33 — within one patient" = "p33",
                       "Sydney 2012 — global" = "syd",
                       "Rate comparison" = "both"), "p33"),
        conditionalPanel("input.rtt_which != 'both'",
                         checkboxInput("rtt_line", "Show regression line", TRUE),
                         checkboxInput("rtt_lab", "Label points", FALSE)),
        hr(),
        div(style = sprintf("font-size:.8rem; color:%s; line-height:1.5;", MUTED),
            "Root-to-tip regression: distance from the root of the tree ",
            "against sampling date. A straight line means clock-like ",
            "evolution, and the slope is the rate.")
      ),
      card(
        full_screen = TRUE,
        card_header(textOutput("rtt_title", inline = TRUE)),
        plotOutput("rtt_plot", height = "430px")
      ),
      card(card_body(div(class = "lead-note", uiOutput("rtt_note")))),
      div(style = "margin:.2rem 0 1.2rem 0;",
          actionButton("go_sel", "Next \u2192 measure selection",
                       class = "btn btn-primary"))
    )
  ),

  # --------------------------------------------------------------- selection
  nav_panel(
    "Selection (dN/dS)",
    layout_sidebar(
      sidebar = sidebar(
        width = 290,
        selectInput("sel_a", "Earlier sequence", LABELS, LABELS[1]),
        selectInput("sel_b", "Later sequence", LABELS, LABELS[N_SEQ]),
        radioButtons("sel_show", "Codon table",
                     c("Non-synonymous only" = "non",
                       "All differences" = "all"), "non"),
        hr(),
        div(style = sprintf("font-size:.8rem; color:%s; line-height:1.5;", MUTED),
            HTML("<b>&omega; &lt; 1</b> purifying selection &mdash; amino-acid
                  changes are being removed.<br>
                  <b>&omega; &asymp; 1</b> neutral.<br>
                  <b>&omega; &gt; 1</b> positive selection &mdash; change is
                  being favored.")),
        hr(),
        actionLink("sel_reset", "Reset to first vs last",
                   style = sprintf("font-size:.82rem; color:%s;", GREEN))
      ),
      uiOutput("sel_cards"),
      layout_columns(
        col_widths = c(5, 7),
        card(card_header("Synonymous vs non-synonymous"),
             plotOutput("sel_bar", height = "300px")),
        card(card_header(textOutput("sel_tab_title", inline = TRUE)),
             div(style = "max-height:300px; overflow-y:auto;",
                 tableOutput("sel_table")))
      ),
      card(card_body(div(class = "lead-note", uiOutput("sel_note")))),
      div(style = "margin:.2rem 0 1.2rem 0;",
          actionButton("go_q", "Next \u2192 the questions",
                       class = "btn btn-primary"))
    )
  ),

  # --------------------------------------------------------------- questions
  nav_panel(
    "Questions",
    div(class = "container-fluid", style = "max-width:900px;",
      h4("Work through these"),
      p(class = "lead-note",
        "Each question can be answered by changing something in the tabs."),
      br(),
      h5("1. Reading an alignment"),
      div(class = "q-box",
          "On the Alignment tab, switch between ", strong("Whole gene"),
          " and ", strong("Variable sites only"), ". Most of the capsid is ",
          "identical across all nine time points. Why is that expected over ",
          "three months inside one person, when the same gene differs ",
          "substantially between GII.4 variants sampled years apart?"),
      div(class = "q-box",
          "Switch the sequence type from nucleotide to amino acid with ",
          strong("Variable sites only"), " selected. ",
          sprintf("Of the %d variable nucleotide sites, %d change the amino
                   acid. ", SUM$n_var_nt, SUM$n_var_aa),
          "If mutations were accumulating without regard to their effect on ",
          "the protein, roughly a quarter of random changes would be silent. ",
          "What does this ratio already hint about selection, before you ",
          "compute dN/dS on the next tab?"),
      h5("2. Trees and time"),
      div(class = "q-box",
          "On the Phylogeny tab, compare the P33 tree with the Sydney 2012 ",
          "tree. The P33 tree spans months and the Sydney tree spans more ",
          "than a decade. Describe the shape difference. In which tree do the ",
          "sampling dates line up with the tip positions more tightly, and why?"),
      div(class = "q-box",
          "The root of the P33 tree is dated before the first sample. What ",
          "biological event does that inferred root correspond to?"),
      h5("3. Rates"),
      div(class = "q-box",
          sprintf("On the Molecular clock tab, choose Rate comparison. The
                   within-host rate is about %.1f times the population-level
                   rate for the same variant. Explain this using the idea of a
                   transmission bottleneck — which mutations get counted
                   in each case?", SUM$fold)),
      div(class = "q-box",
          "The P33 regression uses only nine points over three months. What ",
          "would make you cautious about the precision of that slope? Look at ",
          "the R-squared values for both datasets."),
      h5("4. Selection"),
      div(class = "q-box",
          "On the Selection tab, keep the earliest sequence fixed and step the ",
          "later sequence through each time point in turn. Do the counts of ",
          "synonymous and non-synonymous differences grow at the same pace?"),
      div(class = "q-box",
          "Compare the earliest sequence with one of the early time points ",
          "(day 18, 22 or 29). dN/dS comes back undefined. Read the note under ",
          "the figure and explain, in terms of how dS is calculated, why that ",
          "happens — and why it is biologically interesting rather than ",
          "merely a numerical failure."),
      div(class = "q-box",
          sprintf("For the first versus last comparison, ω = %.2f, with %d
                   non-synonymous changes against only %d synonymous ones.
                   State in one sentence what that means for selection on the
                   capsid. Why is this the <i>opposite</i> of what is usually
                   found for a structural gene at the population level?",
                  SUM$fl_omega, SUM$fl_non, SUM$fl_syn) |> HTML()),
      div(class = "q-box",
          "Most published GII.4 capsid analyses report dN/dS well below 1 ",
          "across the gene. This patient gives a value above 1 over three ",
          "months. What is different about the selective environment inside a ",
          "single immunocompromised host, and which part of the capsid would ",
          "you expect those amino-acid changes to fall in?"),
      div(class = "q-box",
          "A single ratio averaged across 540 codons is a blunt instrument: it ",
          "can hide a few strongly selected sites in a sea of conserved ones, ",
          "and it can also be pulled above 1 by a small number of changes when ",
          "the total count is low. Given that only 12 codons differ here, how ",
          "much confidence do you place in ω = ",
          sprintf("%.2f", SUM$fl_omega),
          "? What analysis would you run to test selection site by site?"),
      br(), br()
    )
  ),

  nav_spacer(),
  nav_item(tags$a("BIO463", href = "../", class = "nav-link",
                  style = "opacity:.85;"))
)

# ==========================================================================
# Server
# ==========================================================================
server <- function(input, output, session) {

  # step-by-step navigation: each tab hands off to the next
  observeEvent(input$go_align, nav_select("nav", "Alignment", session))
  observeEvent(input$go_tree,  nav_select("nav", "Phylogeny", session))
  observeEvent(input$go_clock, nav_select("nav", "Molecular clock", session))
  observeEvent(input$go_sel,   nav_select("nav", "Selection (dN/dS)", session))
  observeEvent(input$go_q,     nav_select("nav", "Questions", session))

  # ------------------------------------------------------- alignment
  is_nt <- reactive(input$aln_level == "nt")
  cur_mat  <- reactive(if (is_nt()) ALN$nt else ALN$aa)
  cur_var  <- reactive(if (is_nt()) ALN$var_nt else ALN$var_aa)
  cur_unit <- reactive(if (is_nt()) "Nucleotide" else "Amino acid")

  output$win_slider <- renderUI({
    n <- ncol(cur_mat())
    if (is_nt()) {
      sliderInput("aln_win", "Nucleotide positions", 1, n,
                  value = c(1020, 1170), step = 10, width = "100%")
    } else {
      sliderInput("aln_win", "Codon positions", 1, n,
                  value = c(340, 390), step = 5, width = "100%")
    }
  })

  cols_shown <- reactive({
    n <- ncol(cur_mat())
    switch(input$aln_view,
      all = seq_len(n),
      var = cur_var(),
      win = {
        w <- input$aln_win
        if (is.null(w)) seq_len(n) else seq(max(1, w[1]), min(n, w[2]))
      })
  })

  output$aln_title <- renderText({
    k <- length(cols_shown())
    sprintf("P33 capsid — %s alignment, %d position%s shown",
            tolower(cur_unit()), k, if (k == 1) "" else "s")
  })

  output$aln_plot <- renderPlot({
    cols <- cols_shown()
    mat  <- cur_mat()
    df   <- melt_aln(mat, cols, "value")
    pal  <- if (is_nt()) NT_COLS else AA_COLS
    # any unexpected character -> grey, so scale never errors
    extra <- setdiff(unique(df$value), names(pal))
    if (length(extra)) pal <- c(pal, setNames(rep("#C9CCC6", length(extra)), extra))

    discrete_x <- length(cols) <= 220
    df$xpos <- if (discrete_x) factor(df$position, levels = cols) else df$position

    p <- ggplot(df, aes(x = xpos, y = taxon, fill = value))

    if (discrete_x) {
      p <- p + geom_tile(colour = CREAM, linewidth = 0.25)
      brk <- if (input$aln_view == "var") as.character(cols) else
        as.character(cols[seq(1, length(cols),
                              length.out = min(12, length(cols)))])
      p <- p + scale_x_discrete(breaks = brk)
    } else {
      p <- p + geom_tile()
      p <- p + scale_x_continuous(expand = c(0, 0))
    }

    if (isTRUE(input$aln_outline)) {
      vdf <- df[df$position %in% cur_var(), ]
      if (nrow(vdf) && discrete_x) {
        p <- p + geom_tile(data = vdf, aes(x = xpos, y = taxon),
                           fill = NA, colour = INK, linewidth = 0.7,
                           inherit.aes = FALSE)
      } else if (nrow(vdf)) {
        p <- p + geom_vline(xintercept = unique(vdf$position),
                            colour = INK, linewidth = 0.25, alpha = 0.55)
      }
    }

    if (isTRUE(input$aln_letters) && input$aln_view != "all" &&
        length(cols) <= 120) {
      p <- p + geom_text(aes(label = value), size = 2.6,
                         colour = INK, fontface = "bold")
    }

    p + scale_fill_manual(values = pal, name = NULL,
                          guide = guide_legend(nrow = 1)) +
      labs(x = sprintf("%s position", cur_unit()), y = NULL) +
      theme_matte() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(size = 8, colour = INK),
            axis.text.x = element_text(
              angle = if (discrete_x && length(cols) > 40) 90 else 0,
              vjust = 0.5, hjust = if (discrete_x && length(cols) > 40) 1 else 0.5))
  })

  output$aln_note <- renderUI({
    k <- length(cols_shown())
    nv <- sum(cols_shown() %in% cur_var())
    HTML(sprintf(
      "Showing <b>%d</b> of %d %s positions, of which <b>%d</b> vary between
       time points. Across the whole gene only %d of %d nucleotide positions
       change during the infection, and just %d of those alter the protein —
       the other changes are synonymous.",
      k, ncol(cur_mat()), tolower(cur_unit()), nv,
      SUM$n_var_nt, SUM$p33_nt, SUM$n_var_aa))
  })

  # ------------------------------------------------------- phylogeny
  tre <- reactive(if (input$tree_which == "p33") TREE else SYD)

  output$tree_title <- renderText({
    if (input$tree_which == "p33")
      sprintf("P33 capsid — dated phylogeny (%d sequences)", SUM$p33_n)
    else
      sprintf("Sydney 2012 GII.4 capsid — dated phylogeny (n = %d)", SYD$n)
  })

  output$tree_plot <- renderPlot({
    t  <- tre()
    nd <- t$layout$nodes
    sg <- t$layout$segs
    tips <- nd[nd$is_tip, ]
    p33 <- input$tree_which == "p33"

    rng  <- range(nd$x)
    padr <- if (p33 && isTRUE(input$tree_tips)) diff(rng) * 0.42
            else diff(rng) * 0.04

    p <- ggplot() +
      geom_segment(data = sg, aes(x = x, xend = xend, y = y, yend = yend),
                   colour = GREEN_DARK,
                   linewidth = if (p33) 0.45 else 0.22,
                   lineend = "round") +
      geom_point(data = tips, aes(x = x, y = y),
                 colour = GREEN, size = if (p33) 2 else 0.5) +
      scale_x_continuous(limits = c(rng[1], rng[2] + padr)) +
      labs(x = if (p33) "Time (decimal years within the infection)" else "Year",
           y = NULL) +
      theme_matte() +
      theme(panel.grid.major.y = element_blank(),
            axis.text.y = element_blank(),
            legend.position = "none")

    if (p33 && isTRUE(input$tree_tips)) {
      p <- p + geom_text(data = tips,
                         aes(x = x, y = y, label = label),
                         hjust = 0, nudge_x = diff(rng) * 0.02,
                         size = 3.1, colour = INK)
    }
    p
  })

  output$tree_note <- renderUI({
    if (input$tree_which == "p33") {
      HTML(sprintf(
        "Nine sequences diverging progressively from a common ancestor dated to
         roughly the start of the infection. The tree is time-scaled, so
         horizontal position is calendar time rather than genetic distance.
         Inferred clock rate: <b>%s</b> substitutions/site/year.",
        sci(TREE$rate)))
    } else {
      HTML(sprintf(
        "A random subsample of <b>%d</b> of the %d available Sydney 2012 capsid
         sequences, collected worldwide between %s and %s. This tree covers
         more than a decade of global transmission, and its clock rate
         (<b>%s</b> subs/site/year) reflects only those mutations that
         successfully passed from host to host.<br><br>
         <span style='font-size:.86rem;'>A caveat worth raising in class: this
         is a neighbor-joining tree under a strict clock, chosen because it is
         fast and transparent, not because it is the best available model. It
         puts the common ancestor around %d, earlier than published estimates
         for this variant, and treedater reports a poor clock fit for part of
         the tree. A maximum-likelihood tree with a relaxed clock would date it
         better.</span>",
        SYD$n, SYD$n_available, substr(SYD$date_min, 1, 4),
        substr(SYD$date_max, 1, 4), sci(SYD$rate), floor(SYD$tmrca)))
    }
  })

  # ------------------------------------------------------- molecular clock
  output$rtt_title <- renderText({
    switch(input$rtt_which,
      p33  = "P33 — root-to-tip regression",
      syd  = "Sydney 2012 — root-to-tip regression",
      both = "Within-host vs population-level evolutionary rate")
  })

  output$rtt_plot <- renderPlot({
    if (input$rtt_which == "both") {
      df <- data.frame(
        Dataset = factor(c("P33\n(within-host)", "Sydney 2012\n(population)"),
                         levels = c("P33\n(within-host)",
                                    "Sydney 2012\n(population)")),
        Rate = c(SUM$p33_rate, SUM$syd_rate))
      return(
        ggplot(df, aes(Dataset, Rate, fill = Dataset)) +
          geom_col(width = 0.5) +
          geom_text(aes(label = sci(Rate)), vjust = -0.6,
                    size = 4, colour = INK) +
          scale_fill_manual(values = c(GREEN, GREEN_LIGHT)) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                             labels = function(z) sci(z, 1)) +
          labs(x = NULL, y = "Substitutions / site / year") +
          theme_matte(13) + theme(legend.position = "none")
      )
    }

    p33 <- input$rtt_which == "p33"
    t   <- if (p33) TREE else SYD
    df  <- t$rtt

    p <- ggplot(df, aes(decimal_date, distance)) +
      geom_point(colour = GREEN, size = if (p33) 2.6 else 1.1,
                 alpha = if (p33) 1 else 0.6)

    if (isTRUE(input$rtt_line)) {
      p <- p + geom_abline(slope = t$slope, intercept = t$intercept,
                           colour = GREEN_DARK, linewidth = 0.6,
                           linetype = "22")
    }
    if (isTRUE(input$rtt_lab) && p33) {
      p <- p + geom_text(aes(label = label), size = 2.9, colour = MUTED,
                         vjust = -0.9)
    }

    p + labs(x = if (p33) "Sampling time (decimal years)" else "Sampling year",
             y = "Root-to-tip distance (subs/site)",
             subtitle = sprintf("slope = %s subs/site/year   ·   R² = %.3f",
                                sci(t$slope), t$r2)) +
      scale_y_continuous(labels = function(z) sci(z, 1)) +
      theme_matte()
  })

  output$rtt_note <- renderUI({
    switch(input$rtt_which,
      p33 = HTML(sprintf(
        "The upward trend confirms clock-like evolution within the host: later
         samples sit further from the root. The fitted slope is <b>%s</b>
         subs/site/year (R² = %.3f). With only nine points over %d days,
         treat the precision cautiously — the direction is more reliable
         than the exact value.",
        sci(TREE$slope), TREE$r2, SUM$p33_span_days)),
      syd = HTML(sprintf(
        "With %d sequences over more than a decade the regression is far better
         constrained (R² = %.3f). Slope: <b>%s</b> subs/site/year.",
        SYD$n, SYD$r2, sci(SYD$slope))),
      both = HTML(sprintf(
        "The within-host rate is about <b>%.1f×</b> the population-level
         rate for the same GII.4 variant. Inside a host, mutations accumulate
         with little filtering, including many that are mildly deleterious. At
         transmission, only a tiny fraction of the viral population founds the
         next infection, and that bottleneck removes most of those variants.
         The population-level rate therefore counts only the changes that
         survived repeated bottlenecks over a decade of global spread.",
        SUM$fold)))
  })

  # ------------------------------------------------------- selection
  observeEvent(input$sel_reset, {
    updateSelectInput(session, "sel_a", selected = LABELS[1])
    updateSelectInput(session, "sel_b", selected = LABELS[N_SEQ])
  })

  pair <- reactive({
    ia <- match(input$sel_a, LABELS)
    ib <- match(input$sel_b, LABELS)
    if (is.na(ia) || is.na(ib) || ia == ib) return(NULL)
    i <- min(ia, ib); j <- max(ia, ib)
    key <- paste(i, j, sep = "-")
    list(row = DNDS$pairs[DNDS$pairs$key == key, ],
         changes = DNDS$changes[[key]],
         flipped = ia > ib)
  })

  output$sel_cards <- renderUI({
    pr <- pair()
    if (is.null(pr)) {
      return(div(class = "q-box",
                 "Choose two different sequences to compare."))
    }
    r <- pr$row
    om <- r$omega
    verdict <- if (!is.finite(om)) "dS = 0, undefined" else
      if (om < 1) "purifying selection" else
      if (om > 1) "positive selection" else "neutral"
    div(
      layout_columns(
        col_widths = c(2, 2, 2, 2, 2, 2),
        stat_card(r$days_apart, "Days apart"),
        stat_card(r$total, "Codon differences"),
        stat_card(r$syn, "Synonymous", "silent"),
        stat_card(r$nonsyn, "Non-synonymous", "protein changing"),
        stat_card(sprintf("%.4f", r$dS), "dS",
                  sprintf("dN = %.4f", r$dN)),
        stat_card(if (is.na(om)) "—" else sprintf("%.3f", om),
                  "dN/dS (ω)", verdict)
      ),
      br()
    )
  })

  output$sel_bar <- renderPlot({
    pr <- pair(); if (is.null(pr)) return(NULL)
    df <- data.frame(
      Type = factor(c("Synonymous", "Non-synonymous"),
                    levels = c("Synonymous", "Non-synonymous")),
      Count = c(pr$row$syn, pr$row$nonsyn))
    ggplot(df, aes(Type, Count, fill = Type)) +
      geom_col(width = 0.5) +
      geom_text(aes(label = Count), vjust = -0.5, size = 4.6, colour = INK) +
      scale_fill_manual(values = c(GREEN_LIGHT, GREEN)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(x = NULL, y = "Number of codon differences") +
      theme_matte(12) + theme(legend.position = "none")
  })

  output$sel_tab_title <- renderText({
    pr <- pair(); if (is.null(pr)) return("Codon differences")
    if (input$sel_show == "non") "Non-synonymous codon changes"
    else "All codon differences"
  })

  output$sel_table <- renderTable({
    pr <- pair(); if (is.null(pr)) return(NULL)
    ch <- pr$changes
    if (!nrow(ch)) return(data.frame(` ` = "No differences between these two sequences",
                                     check.names = FALSE))
    if (input$sel_show == "non") ch <- ch[ch$type == "Non-synonymous", ]
    if (!nrow(ch)) return(data.frame(` ` = "No non-synonymous changes between these two",
                                     check.names = FALSE))
    if (pr$flipped) {
      ch <- transform(ch, codon_from = ch$codon_to, codon_to = ch$codon_from,
                      aa_from = ch$aa_to, aa_to = ch$aa_from)
      ch$change <- ifelse(ch$type == "Synonymous",
                          sprintf("%s%d (silent)", ch$aa_from, ch$codon),
                          sprintf("%s%d%s", ch$aa_from, ch$codon, ch$aa_to))
    }
    out <- ch[order(ch$codon), c("codon", "nt_pos", "codon_from", "codon_to",
                                 "aa_from", "aa_to", "change", "type")]
    names(out) <- c("Codon", "nt pos", "Codon A", "Codon B",
                    "AA A", "AA B", "Change", "Type")
    out
  }, striped = TRUE, spacing = "xs", width = "100%", align = "c")

  output$sel_note <- renderUI({
    pr <- pair(); if (is.null(pr)) return(NULL)
    r <- pr$row
    om <- r$omega
    if (!is.finite(om)) {
      return(HTML(sprintf(
        "These two sequences differ at <b>%d</b> codon%s, and <b>every one of
         those changes alters the protein</b> — there are no synonymous
         differences at all, so dS = 0 and the ratio dN/dS cannot be computed
         (dividing by zero). That is not a failure of the method so much as a
         signal: over only %d days there has not yet been enough time to
         accumulate silent changes, yet amino-acid changes are already present.
         Choose time points further apart to get a defined ω.",
        r$total, if (r$total == 1) "" else "s", r$days_apart)))
    }
    HTML(sprintf(
      "Over %d days these two sequences differ at <b>%d</b> codons: <b>%d</b>
       synonymous and <b>%d</b> non-synonymous. Normalising by the number of
       sites at which each kind of change could have occurred gives
       dN = %.4f, dS = %.4f, so <b>ω = %.3f</b>. %s Note that a
       non-zero dS is itself evidence of continuing replication and mutation
       throughout the infection — the virus was never static, even where
       the protein stayed the same.",
      r$days_apart, r$total, r$syn, r$nonsyn, r$dN, r$dS, om,
      if (om < 1)
        "Because ω is below 1, most amino-acid-changing mutations are being
         removed by purifying selection, consistent with structural and
         functional constraints on the capsid."
      else
        "Because ω exceeds 1, amino-acid change is being <i>favored</i>
         over this interval. In a prolonged infection this is the expected
         signature of antigenic escape: the capsid is the target of the
         patient's antibody response, so mutations that change the epitopes
         on its protruding domain are advantageous to the virus."))
  })
}

shinyApp(ui, server)
