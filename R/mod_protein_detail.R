#' protein_detail UI Function
#'
#' @description A slide-in (Bootstrap offcanvas) panel that shows details
#'   for a single protein/gene: metadata, external database links, and a
#'   per-group intensity boxplot. Opened programmatically when a point is
#'   clicked elsewhere (e.g. the volcano plot).
#'
#' @param id Internal parameter for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_protein_detail_ui <- function(id) {
  ns <- NS(id)
  tagList(

    # The offcanvas panel itself (hidden until opened)
    div(
      class = "offcanvas offcanvas-end",
      tabindex = "-1",
      id = ns("panel"),
      style = "width: 440px;",

      div(
        class = "offcanvas-header",
        h5(class = "offcanvas-title", textOutput(ns("title"), inline = TRUE)),
        tags$button(
          type = "button",
          class = "btn-close",
          `data-bs-dismiss` = "offcanvas",
          `aria-label` = "Close"
        )
      ),

      div(
        class = "offcanvas-body",
        uiOutput(ns("info")),
        hr(),
        h6("Abundance by group"),
        plotOutput(ns("boxplot"), height = "300px"),
        hr(),
        h6("External resources"),
        uiOutput(ns("links"))
      )
    ),

    # Tiny handler so the server can open the offcanvas on click.
    tags$script(shiny::HTML(sprintf(
      "Shiny.addCustomMessageHandler('%s', function(id){
         var el = document.getElementById(id);
         if (el) { bootstrap.Offcanvas.getOrCreateInstance(el).show(); }
       });",
      ns("open_offcanvas")
    )))
  )
}

# Build the external-resource links for a gene / UniProt accession.
# `pubmed_context` (optional) adds an "AND <context>" term to the PubMed
# search, e.g. "KRAS" -> a study-specific co-search.
protein_detail_links <- function(gene, accession, pubmed_context = NULL) {
  gene_enc <- utils::URLencode(gene, reserved = TRUE)

  if (is.null(pubmed_context) || !nzchar(pubmed_context)) {
    pubmed_label <- "PubMed"
    pubmed_term <- gene_enc
  } else {
    ctx_enc <- utils::URLencode(pubmed_context, reserved = TRUE)
    pubmed_label <- sprintf("PubMed (gene + %s)", pubmed_context)
    pubmed_term <- paste0(gene_enc, "+AND+", ctx_enc)
  }

  items <- list(
    UniProt = if (!is.na(accession) && nzchar(accession)) {
      paste0("https://www.uniprot.org/uniprotkb/", accession, "/entry")
    },
    `Human Protein Atlas` = paste0("https://www.proteinatlas.org/search/", gene_enc),
    GeneCards = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene_enc)
  )
  items[[pubmed_label]] <- paste0("https://pubmed.ncbi.nlm.nih.gov/?term=", pubmed_term)
  items <- items[!vapply(items, is.null, logical(1))]

  tags$ul(
    class = "list-unstyled",
    lapply(names(items), function(label) {
      tags$li(
        tags$a(href = items[[label]], target = "_blank", rel = "noopener", label),
        " ↗"
      )
    })
  )
}

#' protein_detail Server Functions
#'
#' @param id Module id.
#' @param clicked_gene A reactive returning the clicked gene symbol (or NULL).
#' @param reactive_de A reactive returning the DE data.frame.
#' @param sample_meta A data.frame with columns SampleID and group.
#' @param fc_col,p_col,feature_col,accession_col Column names.
#' @param group_colors Optional named character vector mapping group
#'   values to fill colors. If NULL, ggplot's default palette is used.
#' @param pubmed_context Optional term added to the PubMed search
#'   (e.g. the study's gene/condition); NULL searches the gene alone.
#'
#' @noRd
mod_protein_detail_server <- function(id,
                                      clicked_gene,
                                      reactive_de,
                                      sample_meta,
                                      fc_col = "logFC",
                                      p_col = "P.Value",
                                      feature_col = "Genes",
                                      accession_col = "Protein.Group",
                                      group_colors = NULL,
                                      pubmed_context = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # The DE row for the clicked gene (first match).
    selected <- reactive({
      g <- clicked_gene()
      if (is.null(g) || !nzchar(g)) {
        return(NULL)
      }
      df <- reactive_de()
      hit <- df[df[[feature_col]] == g, , drop = FALSE]
      if (nrow(hit) == 0) {
        return(NULL)
      }
      hit[1, , drop = FALSE]
    })

    # Open the offcanvas whenever a new gene is clicked.
    observeEvent(clicked_gene(), {
      if (!is.null(selected())) {
        session$sendCustomMessage(ns("open_offcanvas"), ns("panel"))
      }
    })

    output$title <- renderText({
      row <- selected()
      if (is.null(row)) "Protein detail" else row[[feature_col]]
    })

    output$info <- renderUI({
      row <- selected()
      if (is.null(row)) {
        return(p(class = "text-muted", "Click a significant point to inspect it."))
      }
      desc <- if ("First.Protein.Description" %in% colnames(row)) {
        row[["First.Protein.Description"]]
      } else {
        NULL
      }
      tagList(
        if (!is.null(desc)) p(class = "text-muted", desc),
        tags$ul(
          class = "list-unstyled mb-0",
          tags$li(tags$b("Accession: "), row[[accession_col]]),
          tags$li(tags$b("log2 FC: "), round(row[[fc_col]], 3)),
          tags$li(tags$b("P-value: "), signif(row[[p_col]], 3)),
          if ("adj.P.Val" %in% colnames(row)) {
            tags$li(tags$b("Adj. P: "), signif(row[["adj.P.Val"]], 3))
          }
        )
      )
    })

    output$links <- renderUI({
      row <- selected()
      if (is.null(row)) {
        return(NULL)
      }
      protein_detail_links(
        row[[feature_col]],
        row[[accession_col]],
        pubmed_context = pubmed_context
      )
    })

    output$boxplot <- renderPlot({
      row <- selected()
      if (is.null(row)) {
        return(NULL)
      }
      sample_cols <- intersect(sample_meta$SampleID, colnames(row))
      vals <- as.numeric(row[1, sample_cols])
      # Treat 0 (missing/undetected) as NA so it doesn't drag the plot
      # down. The boxplot shows detected values only -- detection counts
      # are surfaced on the axis so missingness stays visible (the DE
      # logFC is computed on imputed data and can differ for proteins
      # with high missingness).
      vals[vals == 0] <- NA

      plot_df <- data.frame(
        intensity = vals,
        group = sample_meta$group[match(sample_cols, sample_meta$SampleID)]
      )

      # Detection counts per group (detected = non-missing)
      detect <- lapply(split(plot_df$intensity, plot_df$group), function(v) {
        c(detected = sum(!is.na(v)), total = length(v))
      })
      group_labels <- vapply(names(detect), function(g) {
        sprintf("%s\n(%d/%d detected)", g, detect[[g]]["detected"], detect[[g]]["total"])
      }, character(1))
      names(group_labels) <- names(detect)

      plot_df <- plot_df[!is.na(plot_df$intensity), , drop = FALSE]

      ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = group, y = intensity, fill = group)
      ) +
        ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.5) +
        ggplot2::geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
        (if (!is.null(group_colors)) {
          ggplot2::scale_fill_manual(values = group_colors)
        }) +
        ggplot2::scale_x_discrete(labels = group_labels) +
        ggplot2::labs(
          x = NULL,
          y = "Intensity",
          caption = "Detected (non-zero) values only; logFC is computed on imputed data."
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
    })

  })
}
