#' heatmap UI Function
#'
#' @description Heatmap of group-averaged, row z-scored expression. Users
#'   enter genes (one per line) and choose colData column(s) to group
#'   samples by.
#'
#' @param id Internal parameter for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_heatmap_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Heatmap controls",
      width = 280,
      textAreaInput(
        ns("genes"),
        "Genes (one per line)",
        rows = 10,
        placeholder = "BHMT\nCPS1\nSERPINB3"
      ),
      selectInput(
        ns("group_by"),
        "Group samples by",
        choices = NULL,
        multiple = FALSE
      ),
      checkboxInput(ns("cluster_rows"), "Cluster rows", value = TRUE)
    ),
    plotly::plotlyOutput(ns("plot"), height = "600px")
  )
}

#' heatmap Server Functions
#'
#' @param id Module id.
#' @param se A SummarizedExperiment.
#' @param default_genes Character vector pre-filled in the textbox.
#' @param default_group_by colData column selected by default.
#' @param gene_col rowData column to match genes against.
#'
#' @return A list with reactive `clicked_gene` (gene of the clicked tile,
#'   or NULL).
#'
#' @noRd
mod_heatmap_server <- function(id,
                               se,
                               default_genes = character(0),
                               default_group_by = NULL,
                               gene_col = "Genes") {
  moduleServer(id, function(input, output, session) {

    coldata_cols <- colnames(SummarizedExperiment::colData(se))

    # Populate grouping choices + defaults from the SE.
    updateSelectInput(
      session,
      "group_by",
      choices = coldata_cols,
      selected = if (is.null(default_group_by)) coldata_cols[1] else default_group_by
    )
    if (length(default_genes) > 0) {
      updateTextAreaInput(
        session,
        "genes",
        value = paste(default_genes, collapse = "\n")
      )
    }

    genes <- reactive({
      raw <- strsplit(input$genes, "\n", fixed = TRUE)[[1]]
      trimws(raw)
    })

    # The heatmap object (also carries attr "gene_order" for click mapping)
    hm <- reactive({
      gl <- genes()
      gl <- gl[nzchar(gl)]
      validate(
        need(length(gl) > 0, "Enter at least one gene (one per line)."),
        need(length(input$group_by) > 0, "Select a grouping column.")
      )
      tryCatch(
        plot_heatmap(
          se,
          genes = gl,
          group_by = input$group_by,
          gene_col = gene_col,
          cluster_rows = input$cluster_rows,
          interactive = TRUE,
          plotly_source = session$ns("plot"),
          click_input_id = session$ns("click")
        ),
        error = function(e) validate(need(FALSE, conditionMessage(e)))
      )
    })

    output$plot <- plotly::renderPlotly({
      hm()
    })

    # Heatmap traces don't carry a `key`; a custom click handler (see
    # plot_heatmap's click_input_id) pushes the clicked row position (y,
    # bottom = 1) on every click. Resolve it to a gene via the stored order.
    clicked_gene <- reactive({
      ev <- input$click
      if (is.null(ev) || is.null(ev$value)) {
        return(NULL)
      }
      order <- attr(hm(), "gene_order")
      idx <- round(as.numeric(ev$value))
      if (is.null(order) || is.na(idx) || idx < 1 || idx > length(order)) {
        return(NULL)
      }
      order[idx]
    })

    list(clicked_gene = clicked_gene)

  })
}
