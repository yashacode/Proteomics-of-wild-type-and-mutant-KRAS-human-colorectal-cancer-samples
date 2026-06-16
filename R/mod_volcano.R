#' volcano UI Function
#'
#' @description A shiny Module for an interactive volcano plot with
#'   module-owned threshold controls, a significant-gene table, and a
#'   TSV download of the original table.
#'
#' @param id Internal parameter for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_volcano_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Volcano controls",
      sliderInput(
        ns("neglog10p_cutoff"),
        "Significance cutoff: -log10(p)",
        min = 0,
        max = 10,
        value = round(-log10(0.05), 2),
        step = 0.1
      ),
      sliderInput(
        ns("fc_cutoff"),
        "Fold change cutoff: |log2 fc|",
        min = 0,
        max = 5,
        value = 1,
        step = 0.1
      ),
      downloadButton(ns("download_tsv"), "Download table (TSV)")
    ),
    plotly::plotlyOutput(ns("plot"), height = "550px")
  )
}

#' volcano Server Functions
#'
#' @param id Module id.
#' @param reactive_de A reactive returning the differential-expression
#'   data.frame (the original, unfiltered table).
#' @param fc_col,p_col,feature_col Column names in the data.frame.
#' @param up_label,down_label Legend labels for the positive / negative
#'   fold-change groups (e.g. group-A vs group-B enriched).
#'
#' @return A list of reactives: `clicked_gene` (character or NULL) and
#'   `sig_table` (the significant-feature data.frame).
#'
#' @noRd
mod_volcano_server <- function(id,
                               reactive_de,
                               fc_col = "logFC",
                               p_col = "P.Value",
                               feature_col = "Genes",
                               up_label = "Up",
                               down_label = "Down") {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Convert the -log10(p) slider back to a raw p-value cutoff
    p_cutoff <- reactive({
      10^(-input$neglog10p_cutoff)
    })

    # Single source of truth for "what is significant"
    sig_data <- reactive({
      df <- reactive_de()
      df[
        df[[p_col]] < p_cutoff() &
          abs(df[[fc_col]]) >= input$fc_cutoff,
        ,
        drop = FALSE
      ]
    })

    output$plot <- plotly::renderPlotly({
      plot_volcano(
        reactive_de(),
        fc_col = fc_col,
        p_col = p_col,
        feature_col = feature_col,
        fc_cutoff = input$fc_cutoff,
        p_cutoff = p_cutoff(),
        up_label = up_label,
        down_label = down_label,
        interactive = TRUE,
        plotly_source = ns("plot"),
        click_input_id = ns("click")
      )
    })

    output$download_tsv <- downloadHandler(
      filename = function() {
        paste0("differential_expression_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        utils::write.table(
          reactive_de(),
          file,
          sep = "\t",
          row.names = FALSE,
          quote = FALSE
        )
      }
    )

    # Gene of the clicked point. Driven by a custom click handler (see
    # plot_volcano's click_input_id) that fires on every click, including
    # repeats of the same point.
    clicked_gene <- reactive({
      ev <- input$click
      if (is.null(ev) || is.null(ev$value)) {
        return(NULL)
      }
      ev$value
    })

    list(
      clicked_gene = clicked_gene,
      sig_table = sig_data
    )
  })
}
