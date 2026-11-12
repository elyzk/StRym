# about.R -- the landing page.

about_ui <- function() {
  div(class = "container-fluid", style = "max-width:1100px",
    div(class = "hero", tags$img(src = "leaf.svg", class = "leafbg", alt = ""),
      h1(class = "herologo", tags$img(src = "strym-logo-white.png", alt = "StRym")),
      div(class = "tag", "Streamlined transcriptomic Regulatory mapping"),
      p("From RNA-seq data to regulatory networks, in one place. Tell StRym about your experiment and it plans the analysis, ",
        "runs it, and saves everything in a tidy, numbered folder you can refer back to."),
      actionButton("go_design", "Describe your experiment  →", class = "btn")),

    h3("How it works"),
    div(class = "how",
      div(class = "step", div(class = "num", "1"), h4("Describe your experiment"),
          p("Whole organ/tissue RNA-seq or single cell/single nucleus RNA-seq, how many conditions and time points, which analyses and plots you would like.")),
      div(class = "step", div(class = "num", "2"), h4("Review the plan"),
          p("StRym lists the steps for your choices and checks that your data fit. Tick or untick anything.")),
      div(class = "step", div(class = "num", "3"), h4("Get an organized folder"),
          p("Results land in numbered sub-folders with a summary page. Run everything at once, or one step at a time."))),

    br(),
    div(class = "card", h4("Prefer to do it step by step?"),
      p("Every analysis module is also available on its own under ", strong("Analysis modules"), " in the top bar, with all of its parameters:"),
      div(class = "modgrid",
        div(class = "mod", strong("TUX"), br(), "RNA-seq reads → counts → differentially expressed genes (HISAT2, Cufflinks, Cuffdiff)."),
        div(class = "mod", strong("GENIST"), br(), "Gene regulatory networks from time-course data (Bayesian scoring)."),
        div(class = "mod", strong("RTP-STAR"), br(), "Directed regulatory networks with random forests (GENIE3), clustering and time-course direction."))),
    div(class = "card", h4("Good to know"),
      tags$ul(
        tags$li(HTML("Your input data are only ever <em>read</em>; results are written to the folder you choose.")),
        tags$li(HTML("Hard drive plugged in after the app started? Open the <strong>Drives</strong> tab and press <em>Mount drive</em>.")),
        tags$li("Long analyses run in the background: you can close this page and come back to see the progress."))),
    div(class = "card", h4("Credits"),
      p(class = "hint", "StRym is an R re-implementation of the TuxNet modules. GENIST: M. Angels de Luis Balaguer. RTP-STAR: Natalie M. Clark. ",
        "GENIE3: Huynh-Thu, Irrthum, Wehenkel & Geurts (PLoS ONE 2010). TuxNet: https://github.com/rspurney/TuxNet. ",
        "Example data and the Arabidopsis TF / gene-name lists come from the TuxNet repository.")))
}
