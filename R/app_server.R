#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  de_data <- readRDS(app_sys("app/data/cancer_de.rds"))
  sample_meta <- readRDS(app_sys("app/data/sample_meta.rds"))
  cancer_se <- readRDS(app_sys("app/data/cancer_se.rds"))

  mod_intro_server(
    "intro",
    reactive_de = reactive(de_data),
    groups_value = "KRAS mutant vs wildtype"
  )

  # Returns list(clicked_gene, sig_table) for wiring other tabs later
  volcano_out <- mod_volcano_server(
    "volcano",
    reactive_de = reactive(de_data),
    up_label = "KRAS mutant enriched",
    down_label = "KRAS wild-type enriched"
  )

  heatmap_out <- mod_heatmap_server(
    "heatmap",
    se = cancer_se,
    default_genes = c("BHMT", "CPS1", "SERPINB3", "MAT1A", "STARD4"),
    default_group_by = "mutation"
  )

  # Merge clicks from the volcano and the heatmap into one selection that
  # drives the shared detail panel (most recent click wins). A bump counter
  # is bundled in so re-clicking the *same* gene still invalidates and
  # re-opens the panel (a plain reactiveVal would dedupe identical values).
  selection <- reactiveVal(NULL)
  register_click <- function(gene) {
    if (is.null(gene)) return()
    n <- if (is.null(isolate(selection()))) 1L else isolate(selection())$n + 1L
    selection(list(gene = gene, n = n))
  }
  observeEvent(volcano_out$clicked_gene(), register_click(volcano_out$clicked_gene()))
  observeEvent(heatmap_out$clicked_gene(), register_click(heatmap_out$clicked_gene()))

  # Slide-in detail panel driven by volcano OR heatmap clicks
  mod_protein_detail_server(
    "detail",
    clicked_gene = reactive(selection()$gene),
    reactive_de = reactive(de_data),
    sample_meta = sample_meta,
    group_colors = c(Mutant = "#D55E00", Wildtype = "#0072B2"),
    pubmed_context = "KRAS"
  )

}
