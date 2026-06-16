#' Default content for the intro module
#'
#' Returns the list of app-specific content the intro page renders.
#' Override any field by passing a `content` list to `mod_intro_ui()`;
#' missing fields fall back to these placeholders.
#'
#' @noRd
default_intro_content <- function() {
  list(
    # TODO: study title + one-line plain-language summary
    title = "Study Title Goes Here",
    summary = "One-line plain-language summary of the study and what this portal lets you explore.",

    # Markdown files for the long-form prose (NULL -> a placeholder line)
    abstract_file = NULL,
    methods_file = NULL,

    # "How to use" body (any tag/tagList); NULL -> a placeholder line
    howto = NULL,

    # Static label shown in the "Groups compared" stat box
    groups_label = "Groups compared",

    # Each link entry may have: text (prefix), href, link_text.
    # Empty by default so a calling app fully controls which links show
    # (modifyList merges the links sublist, so defaults would otherwise
    # leak through).
    links = list(),

    # Provenance / independence disclaimer (NULL -> omit the card).
    # Accepts a tagList / HTML / character.
    disclaimer = NULL,

    # Small muted attribution line under the disclaimer (NULL -> omit).
    built_by = NULL
  )
}

# Render one citation/resource list item from a link entry.
intro_link_item <- function(label, entry) {
  tags$li(
    tags$b(paste0(label, ": ")),
    if (!is.null(entry$text)) paste0(entry$text, " "),
    if (!is.null(entry$href)) tags$a(href = entry$href, entry$link_text)
  )
}

# Render markdown from a file path, or a placeholder if none supplied.
# Uses shiny::markdown() (commonmark) to avoid the extra `markdown` pkg.
intro_markdown <- function(path, placeholder) {
  if (!is.null(path) && file.exists(path)) {
    shiny::markdown(paste(readLines(path, warn = FALSE), collapse = "\n"))
  } else {
    p(placeholder)
  }
}

#' intro UI Function
#'
#' @description A reusable landing/about page: study summary, dataset
#'   stat boxes, abstract, how-to guide, methods, and citation/links.
#'   All app-specific content is supplied via `content`.
#'
#' @param id Internal parameter for {shiny}.
#' @param content A named list of content; see `default_intro_content()`.
#'   Any omitted field falls back to the default placeholder.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_intro_ui <- function(id, content = list()) {
  ns <- NS(id)
  content <- utils::modifyList(default_intro_content(), content)

  tagList(

    # ----------------------------
    # HERO HEADER
    # ----------------------------
    div(
      class = "mb-4",
      if (!is.null(content$title)) h2(content$title),
      p(class = "text-muted", content$summary)
    ),

    # ----------------------------
    # DATASET-AT-A-GLANCE STAT BOXES
    # ----------------------------
    bslib::layout_columns(
      fill = FALSE,
      bslib::value_box(
        title = "Samples",
        value = textOutput(ns("n_samples")),
        showcase = icon("vials")
      ),
      bslib::value_box(
        title = "Proteins quantified",
        value = textOutput(ns("n_proteins")),
        showcase = icon("layer-group")
      ),
      bslib::value_box(
        title = "Significant (default cutoff)",
        value = textOutput(ns("n_significant")),
        showcase = icon("chart-simple")
      ),
      bslib::value_box(
        title = content$groups_label,
        value = textOutput(ns("groups")),
        showcase = icon("code-compare")
      )
    ),

    # ----------------------------
    # ABSTRACT (collapsible)
    # ----------------------------
    bslib::accordion(
      class = "my-4",
      open = TRUE,
      bslib::accordion_panel(
        "Abstract",
        intro_markdown(content$abstract_file, "TODO: abstract not yet provided.")
      )
    ),

    # ----------------------------
    # HOW TO USE
    # ----------------------------
    bslib::card(
      bslib::card_header("How to use this app"),
      if (is.null(content$howto)) {
        p(class = "text-muted", "TODO: usage notes not yet provided.")
      } else {
        content$howto
      }
    ),

    # ----------------------------
    # CITATION + LINKS
    # ----------------------------
    bslib::card(
      class = "mt-2",
      bslib::card_header("Citation & resources"),
      tags$ul(
        lapply(names(content$links), function(label) {
          intro_link_item(label, content$links[[label]])
        })
      )
    ),

    # ----------------------------
    # DISCLAIMER / PROVENANCE
    # ----------------------------
    if (!is.null(content$disclaimer)) {
      bslib::card(
        class = "mt-2 border-warning",
        bslib::card_header("Disclaimer"),
        content$disclaimer,
        if (!is.null(content$built_by)) {
          p(class = "text-muted fst-italic mt-2 mb-0", content$built_by)
        }
      )
    }
  )
}

#' intro Server Functions
#'
#' @param id Module id.
#' @param reactive_de A reactive returning the differential-expression
#'   data.frame.
#' @param fc_col,p_col Column names used for the default significance count.
#' @param groups_value Text shown in the "Groups compared" stat box.
#'
#' @noRd
mod_intro_server <- function(id,
                             reactive_de,
                             fc_col = "logFC",
                             p_col = "P.Value",
                             groups_value = "KRAS mutant vs wildtype") {
  moduleServer(id, function(input, output, session) {

    # Intensity columns follow the M#####a sample-id pattern. Fall back
    # to counting numeric columns if none match.
    n_samples_val <- reactive({
      df <- reactive_de()
      n <- sum(grepl("^M[0-9]+a$", colnames(df)))
      if (n == 0) n <- sum(vapply(df, is.numeric, logical(1)))
      n
    })

    output$n_samples <- renderText({
      format(n_samples_val(), big.mark = ",")
    })

    output$n_proteins <- renderText({
      format(nrow(reactive_de()), big.mark = ",")
    })

    output$n_significant <- renderText({
      df <- reactive_de()
      n <- sum(
        df[[p_col]] < 0.05 & abs(df[[fc_col]]) >= 1,
        na.rm = TRUE
      )
      format(n, big.mark = ",")
    })

    output$groups <- renderText({
      groups_value
    })

  })
}
