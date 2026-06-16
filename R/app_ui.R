#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    bslib::page_fillable(

      # ----------------------------
      # THEME (optional but recommended)
      # ----------------------------
      theme = bslib::bs_theme(
        version = 5,
        bootswatch = "flatly"
      ),

      # ----------------------------
      # HEADER (title + icon)
      # ----------------------------
      div(
        style = "display: flex; align-items: center; gap: 10px; margin-bottom: 10px;",

        div(
          style = "font-size: 20px; font-weight: 600;",
          "Proteomics of wild-type and mutant KRAS human colorectal cancer samples"
        )
      ),

      # ----------------------------
      # MAIN AREA WITH TABS
      # ----------------------------
      bslib::navset_tab(

        bslib::nav_panel(
          "About",
          mod_intro_ui(
            "intro",
            content = list(
              title = NULL,
              summary = "Explore differential protein expression in KRAS-mutant vs wildtype colorectal cancer tumors.",
              abstract_file = app_sys("app/about/abstract.md"),
              methods_file = app_sys("app/about/methods.md"),
              howto = tags$ul(
                tags$li(tags$b("Volcano"), " — differential expression; adjust the cutoffs to change what's significant."),
                tags$li(tags$b("Heatmap"), " — enter genes and pick a grouping to compare expression."),
                tags$li("Click any point or tile to open a detail panel for that protein.")
              ),
              links = list(
                Data = list(
                  text = "Wu H. et al. Quantitative Proteomics of wild-type and mutant KRAS human colorectal cancer samples. Soochow University. PRIDE",
                  href = "https://www.ebi.ac.uk/pride/archive/projects/PXD073413",
                  link_text = "PXD073413"
                ),
                License = list(
                  text = "Dataset reused under",
                  href = "https://creativecommons.org/publicdomain/zero/1.0/",
                  link_text = "CC0 1.0 (Public Domain)"
                )
              ),
              disclaimer = tagList(
                shiny::p(
                  "This is an independent demonstration application built for ",
                  "portfolio purposes. It reanalyzes a publicly available dataset ",
                  "from PRIDE/ProteomeXchange (accession ",
                  tags$a(
                    href = "https://www.ebi.ac.uk/pride/archive/projects/PXD073413",
                    "PXD073413"
                  ),
                  ", Wu H. et al., Soochow University). It is ",
                  tags$b("not affiliated with, reviewed by, or endorsed by"),
                  " the original authors. Figures shown here come from an ",
                  "independent reanalysis pipeline and may differ from any ",
                  "results in the original study. Data are reused under the ",
                  tags$a(
                    href = "https://creativecommons.org/publicdomain/zero/1.0/",
                    "CC0 1.0 Public Domain Dedication"
                  ),
                  "."
                )
              ),
              built_by = "Built by Jacob Epstein as a demonstration of interactive bioinformatics app development (R / golem / bslib / plotly / limma)."
            )
          )
        ),

        bslib::nav_panel(
          "Volcano",
          h3("Differential Analysis"),
          mod_volcano_ui("volcano")
        ),

        bslib::nav_panel(
          "Heatmap",
          h3("Expression Heatmap"),
          mod_heatmap_ui("heatmap")
        )
      ),

      # Shared slide-in detail panel (opened by volcano OR heatmap clicks)
      mod_protein_detail_ui("detail")
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "KRAS CRC Proteomics"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
