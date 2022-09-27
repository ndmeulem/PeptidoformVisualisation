#' Shiny app server function
#'
#' @param input provided by shiny
#' @param output provided by shiny
#' @param session provided by shiny
#' @import tidyverse shiny shinyhelper rmarkdown knitr msqrob2 QFeatures limma plotly ggplot2 DT wesanderson BiocManager utils ExploreModelMatrix

library(shiny)
library(shinyhelper)
library(QFeatures)
library(tidyverse)
library(DT)
library(ggplot2)
library(plotly)
library(limma)

server <- (function(input, output, session) {

    #Activate the helper files
    shiny::addResourcePath("PeptidoformViz", system.file("helpfiles", package="PeptidoformViz"))
    shinyhelper::observe_helpers(help_dir = system.file("helpfiles", package="PeptidoformViz"),
                    withMathJax = TRUE)

    #add variables to work with
    variables <- reactiveValues(pe = NULL)

    #When the user clicks read data or example data, it will trigger a number of events:
    #get ecols, read in files, get coldata, do filtering steps

    observeEvent(input$example_data, {


          intensity_file = paste0(system.file("example_data", package="PeptidoformViz"),
                                  "/example_intensitiesfile.csv")

          metadata_file = paste0(system.file("example_data", package="PeptidoformViz"),
                                  "/example_metadatafile.csv")

          #Get intensity columns
          ecols <- grep("2021",
                        names(read.delim(intensity_file,
                                         skip = 2,
                                         sep = ",", header = TRUE)))

          #Read in peptides intensity file into QFeatures object
          pe <- readQFeatures(intensity_file,
                              ecol = ecols,
                              name = "peptideRaw", sep = ",",
                              skip = 2)

          #Read in metadatafile
          metadataFile <- read.delim(metadata_file,
                                     sep = ",")

          variables$pe <- pe
          variables$metadataFile <- metadataFile

          output$read_in_example_data <- renderText({
            "Data read in complete. You may now continue to the preprocessing tab."
          })
    })

    observeEvent(input$go, {
        #make sure files are actually uploaded
        req(input$data$name, input$metadata$name, input$intensityIdentifier)

        #Get intensity columns
        ecols <- grep(input$intensityIdentifier,
                      names(read.delim(input$data$datapath, skip = input$skip,
                                       sep = input$separator, header = input$header)))

        #Read in peptides intensity file into QFeatures object
        pe <- readQFeatures(input$data$datapath,
                            ecol = ecols,
                            name = "peptideRaw", sep = input$separator,
                            skip = input$skip)

        #Read in metadatafile
        metadataFile <- read.delim(input$metadata$datapath,
                                   sep = input$separatorMetadata)

        variables$pe <- pe
        variables$metadataFile <- metadataFile

        output$read_in_data <- renderText({
          "Data read in complete. You may now continue to the preprocessing tab."
        })
    })

    observeEvent({input$go
                  input$example_data}, {

        #Make coldata for the pe object based on metadatafile
        idcols <- colnames(variables$metadataFile)[2:length(colnames(variables$metadataFile))]
        for (col in idcols){
            colData(variables$pe)[as.character(col)] <- as.factor(variables$metadataFile[[as.character(col)]])
        }

        #Update selectInput for arranging of x-axis
        updateSelectInput(session, "x_axis",
                          label = "select variable",
                          choices = variables$idcols,
                          selected = variables$idcols[1])

        #Filtering steps: already calculate nNonZero, zero -> NA is necessary
        rowData(variables$pe[["peptideRaw"]])$nNonZero <- rowSums(assay(variables$pe[["peptideRaw"]]) > 0)
        pe <- zeroIsNA(variables$pe, i = "peptideRaw")

        #Calculate some quick stats about the data
        features = paste(rowData(pe[["peptideRaw"]])[,2], rowData(pe[["peptideRaw"]])[,3],sep="_")
        stats = tibble(stats = c("number of unique proteins",
                                 "number of unique peptides",
                                 "number of unique peptidoforms"))
        stats$values = c(length(unique(rowData(pe[["peptideRaw"]])$Protein)),
                         length(unique(rowData(pe[["peptideRaw"]])$Sequence)),
                         length(unique(features)))

        output$statstable <- renderTable(stats)

        #Plot the before preprocessing plots
        #for ggplot: pivot the intensity assay to long format
        assay_pe_long <- assay(pe[["peptideRaw"]]) %>%
          as_tibble() %>%
          pivot_longer(cols = everything(), names_to = "sample", values_to = "intensity")
        #set the colour palette
        pal <- wesanderson::wes_palette("Darjeeling1",
                                        length(unique(assay_pe_long$sample)),
                                        type = "continuous")
        #Density plot
        output$densityBefore <- renderPlotly({

          p1 <- ggplot(assay_pe_long, aes(x=intensity, color=sample)) +
            geom_density(show.legend = F) +
            scale_colour_manual(values = pal)
          p1 <- ggplotly(p1) %>% layout(showlegend=F) %>%
            config(toImageButtonOptions = list(
              format = "png", filename = "densityplot_before", width = 1920, height = 1080
            ))
          p1
        })
        #Boxplot
        output$boxplotBefore <- renderPlot({
          p1 = ggplot(assay_pe_long, aes(x = sample, y=intensity, color=sample)) +
            geom_boxplot(show.legend = F) +
            scale_color_manual(values=pal)+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size =10))
          p1
        })

        variables$pe <- pe
        variables$idcols <- idcols
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    #---------------
    #Preprocessing
    #--------------
    #User will click on preprocessing button, then do preprocessing as indicated by user
    observeEvent(input$preprocess, {
      #Each time the preprocessing changes, this will be triggered,
      #so data needs to start from scratch every time
      pe2 <- variables$pe

      #logtransformation
      if (input$logTransform != "none"){
        pe2 <- logTransform(pe2, base = as.numeric(input$logTransform),
                           i = "peptideRaw", name = "peptideLog")
      }
      else if (input$logTransform == "none"){
        #This probably means the data are already logtransformed, or it is not necessary
        #So I will manually add the peptideRaw as peptideLog
        pe2 <- addAssay(pe2, pe2[["peptideRaw"]], "peptideLog")
      }
      colData(pe2[["peptideLog"]]) <- colData(pe2)

      #filter on number of nonzero columns
      pe2 <- pe2[rowData(pe2[["peptideLog"]])$nNonZero>=input$nnonzero,,]

      #Global normalisation
      #Either it is a robust normalisation, a centering or no normalisation
      if (input$normalisationMethodGlobal == "robust"){
        #list to store summarised values per protein
        pe_robustval <- list()

        protein_column <- colnames(rowData(pe2[["peptideLog"]]))[input$proteinColumn]

        #Do not change pe2 yet, gives problems further down the line
        pe3 <- aggregateFeatures(pe2,
                                i = "peptideLog",
                                fcol = protein_column,
                                na.rm = TRUE,
                                name = "proteinRobust",
                                fun = MsCoreUtils::robustSummary)

        #Fill out list with the summarised values
        for (prot in unique(rowData(pe2[["peptideLog"]])[,protein_column])){
          pe_robustval[[prot]] = assay(pe3[["proteinRobust"]])[prot,]
        }

        #Use the summarised values to get the normalised assay
        y <- as_tibble(assay(pe2[["peptideLog"]]))
        y_new <- tibble()
        for (prot in unique(rowData(pe2[["peptideLog"]])[,protein_column])){
          #pe_sub <- filterFeatures(pe2,~grepl(protein_column,pattern=prot))
          #pe_sub <- pe2[["peptideLog"]][rowData(pe2[["peptideLog"]])[,protein_column]==prot,]
          pe_sub <- pe2[grepl(prot,rowData(pe2[["peptideLog"]])[[protein_column]])]
          #as_tibble is belangrijk hier, anders worden er random rijen verschoven
          y_ <- as_tibble(assay(pe_sub[["peptideLog"]]))
          #center assay based on the corresponding protein value
          y_scale <- base::scale(y_, center = pe_robustval[[prot]], scale = FALSE)
          y_new <- rbind(y_new, y_scale)
        }
        #Add the normalised assay as a new assay to the existing pe
        y_assay <- SummarizedExperiment(assays=as.matrix(y_new), rowData=rowData(pe2[["peptideLog"]]), colData=colData(pe2[["peptideLog"]]))
        pe2 <- addAssay(pe2, y_assay, name = "peptideLogNorm")
        #Add assay link -> peptideLogNorm is linked to peptideLog
        pe2 <- addAssayLinkOneToOne(pe2, from = "peptideLog", to = "peptideLogNorm")
      }

      else if (input$normalisationMethodGlobal != "none"){
        pe2 <- normalize(pe2,
                        i = "peptideLog",
                        name = "peptideLogNorm",
                        method = input$normalisationMethodGlobal)
      }
      else if (input$normalisationMethodGlobal == "none"){
        #I do continue with the name peptideLognorm, so I will have to add this
        #in this case as well, it is then just unnormalised
        pe2 <- addAssay(pe2, pe2[["peptideLog"]], "peptideLogNorm")
      }

      #Update selectInput for data table as some of proteins may have been filtered out
      updateSelectInput(session, "protein",
                        label = "Protein",
                        choices = rowData(pe2[["peptideLogNorm"]])[[input$proteinColumn]],
                        selected = rowData(pe2[["peptideLogNorm"]])[[input$proteinColumn]][1])

      #Update selectInput for arranging of x-axis
      updateSelectInput(session, "x_axis",
                        label = "select variable",
                        choices = variables$idcols,
                        selected = variables$idcols[1])

      #Calculate some quick stats about the data after preprocessing
      features = paste(rowData(pe2[["peptideLogNorm"]])[,2], rowData(pe2[["peptideLogNorm"]])[,3],sep="_")
      stats = tibble(stats = c("number of unique proteins",
                               "number of unique peptides",
                               "number of unique peptidoforms"))
      stats$values = c(length(unique(rowData(pe2[["peptideLogNorm"]])$Protein)),
                       length(unique(rowData(pe2[["peptideLogNorm"]])$Sequence)),
                       length(unique(features)))


      output$statstableafter <- renderTable(stats)

      #Plot the after preprocessing plots
      #for ggplot: pivot the intensity assay to long format
      assay_pe_long <- assay(pe2[["peptideLogNorm"]]) %>%
        as_tibble() %>%
        pivot_longer(cols = everything(), names_to = "sample", values_to = "intensity")
      #set the colour palette
      pal <- wesanderson::wes_palette("Darjeeling1",
                                      length(unique(assay_pe_long$sample)),
                                      type = "continuous")
      #Density plot
      output$densityAfter <- renderPlotly({

        p1 <- ggplot(assay_pe_long, aes(x=intensity, color=sample)) +
          geom_density(show.legend = F) +
          scale_colour_manual(values = pal)
        p1 <- ggplotly(p1) %>% layout(showlegend=F) %>%
          config(toImageButtonOptions = list(
            format = "png", filename = "densityplot_after", width = 1920, height = 1080
          ))
        p1
      })
      #Boxplot
      output$boxplotAfter <- renderPlot({
        p1 = ggplot(assay_pe_long, aes(x = sample, y=intensity, color=sample)) +
          geom_boxplot(show.legend = F) +
          scale_color_manual(values=pal)+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size =10))
        p1
      })

      variables$pe2 <- pe2
    })

    #If user clicks on radiobuttons for normalisation or selects different protein,
    #or user changes preprocessing settings, or chooses to arrange x-axis differently on plots,
    #variables$proteindf needs to update to the correct version
    observeEvent({input$protein
                 input$normalisationMethod
                 input$preprocess
                 input$x_axis}, {

        req(input$x_axis)
        #Get data for particular protein
        print(variables$pe2)
        proteinpe <- variables$pe2[grepl(input$protein,
                                rowData(variables$pe2[["peptideLogNorm"]])[[input$proteinColumn]])]
        print(proteinpe)

        #Normalisation + get dataset instead of QFeatures object
        #get all metadata + intensity values together in df
        if( input$normalisationMethod != "none"){
            proteinpe <- normalize(proteinpe,
                            i = "peptideLogNorm",
                            name = "peptideLogNormProtein",
                            method = input$normalisationMethod)
            df <- longFormat(proteinpe[["peptideLogNormProtein"]])
        } else if (input$normalisationMethod == "none"){
            df <- longFormat(proteinpe[["peptideLogNorm"]])
        }

        colH <- colData(proteinpe)
        for (col in variables$idcols){
            df[as.character(col)] = as.factor(colH[df$colname, col])
        }
        #add id column
        df <- df %>% as_tibble() %>% tidyr::unite("id", all_of(variables$idcols),sep="_", remove = F)
        #add biorepeat column
        df$biorepeat <- df$id %>% as.factor %>% as.double
        #add features column
        df$features <-
            rep(paste(rowData(proteinpe[["peptideLogNorm"]])[,input$sequenceColumn],
                      rowData(proteinpe[["peptideLogNorm"]])[,input$modificationsColumn],sep="_"),
                      length(unique(df$biorepeat)))
        #arrange according to input$x_axis
        df <- as.data.frame(df)
        df <- df[order(df[,input$x_axis]),]
        df$id = factor(df$id, unique(df$id))
        #Save dataset
        variables$proteindf <- df
        }, ignoreInit = TRUE)


    #Visualisation
    #Plot data table: wide format so that users can easily see the features
    output$proteinDataTable <- DT::renderDataTable(server = FALSE, {
        #Transform dataset into wide format
        proteindf_wide <- variables$proteindf %>% tibble::as_tibble() %>%
            tidyr::pivot_wider(id_cols = c("rowname", "features"), names_from = "id", values_from = "value")
        variables$proteindf_wide = proteindf_wide
        DT::datatable(proteindf_wide,
            extensions = "Buttons",
            options = list(
            paging = TRUE,
            pageLength = 4,
            searching = TRUE,
            dom = "Bfrtip",
            buttons = list("copy", list(
              extend = "collection",
              buttons = c("csv", "excel"),
              text = "Download",
              exportOptions = list(
                modifier = list(page="all")
              ))
            ))
            ) %>% DT::formatStyle(names(proteindf_wide), lineHeight="80%")
    })

    #If user clicks deselect button -> all selected rows are deselected
    proxy = dataTableProxy('proteinDataTable')
    observeEvent(input$deselect, {
      proxy %>% selectRows(NULL)
    })


    #delay input$proteinDataTable_rows_selected so that it does not
    #recalculate every time a row is clicked
    rows_selected <- reactive(input$proteinDataTable_rows_selected)
    rows_selected_d <- debounce(rows_selected, 1000)

    #Plot lineplot
    #Use plotly to make the plot interactive
    output$lineplot <- renderPlotly({
      #set colour palette
      pal <- wesanderson::wes_palette("Darjeeling1",
                                      length(variables$proteindf_wide$features),
                                      type = "continuous")

      #base lineplot (ggplot)
      p1 <- ggplot(
        data = variables$proteindf %>% as.data.frame,
        aes(x = as.factor(id), y = value, group = rowname,colour=rowname)) +
        geom_line(show.legend = F)  +
        scale_colour_manual(values = pal) +
        ggtitle("log2-normalised data") +
        xlab("sample id") +
        ylab("Intensity (log2)") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


      #Plot base lineplot already when nothing is selected
      if (is.null(rows_selected_d())){
        #turn ggplot into plotly object
        p1 <- ggplotly(p1) %>% layout(showlegend=F) %>%
          config(toImageButtonOptions = list(
            format = "png", filename = "lineplot", width = 1920, height = 1080
          ))
        print(p1)
        }

      #If a row is selected, highlight that row
        else if (!is.null(rows_selected_d())){
          #Get dataset with selected rows
          features_selected <- variables$proteindf_wide[rows_selected_d(),]$features
          filter_proteindf <- variables$proteindf %>% as_tibble %>%
              filter(features %in% features_selected)

          hlp1 <- p1 +
              geom_line(
                  aes(x=as.factor(id),
                      y=value),
                  size=2.3,
                  color="palevioletred4",
                  show.legend=FALSE,
                  data = as_tibble(filter_proteindf)) +
              geom_line(
                  aes(x=as.factor(id),
                      y=value),
                  size=0.7,
                  color="grey85",
                  show.legend=FALSE,
                  data = as_tibble(filter_proteindf))

          hlp1 <- ggplotly(hlp1) %>% layout(showlegend=F) %>%
            config(toImageButtonOptions = list(
              format = "png", filename = "lineplot", width = 1920, height = 1080
            ))
          print(hlp1)
        }
    })

    #Plot boxplot
    output$boxplot <- renderPlotly({
      #set colour palette
      pal <- wesanderson::wes_palette("Darjeeling1",
                                      length(unique(variables$proteindf$id)),
                                      type = "continuous")
      #base boxplot
      boxplot <- variables$proteindf %>% as_tibble() %>%
        ggplot(aes(x=as.factor(id),y=value,col=as.factor(id))) +
        geom_boxplot() +
        #geom_jitter(aes(col=features),size=0.7) +
        xlab("sample id") +
        ylab("Intensity (log2)") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              legend.title = element_blank(), legend.position = "None") +
        scale_colour_manual(values = pal)

      #If nothing is selected, already print base boxplot
      if(is.null(rows_selected_d())){
        boxplot <- ggplotly(boxplot) %>% layout(showlegend=F) %>%
          config(toImageButtonOptions = list(
            format = "png", filename = "boxplot", width = 1920, height = 1080
          ))
        print(boxplot)
      }
      #If something is selected, highlight selected dots
       else if (!is.null(rows_selected_d())){
         #Get dataset with selected rows
            features_selected <- variables$proteindf_wide[rows_selected_d(),]$features
            filter_proteindf <- variables$proteindf %>% as_tibble %>%
                filter(features %in% features_selected)

            boxplot_select <- boxplot +
              geom_point(data=filter_proteindf,
                         aes(x=as.factor(id),y=value), color = "red", size = 1.2)

            boxplot_select <- ggplotly(boxplot_select) %>% layout(showlegend=F) %>%
              config(toImageButtonOptions = list(
                format = "png", filename = "boxplot", width = 1920, height = 1080
              ))
            print(boxplot_select)
            }
    })

    #Modeling tab
    ##output possible modeling variables
    output$available_modelvariables <- renderText({
      colnames(variables$metadataFile)},
      sep = ", ")

    ##output metadata
    output$designVariables <- renderDataTable(variables$metadataFile,
                                              options = list(
                                                pageLength = 5
                                              ))

    ##output designmatrix
    design <- reactive(ExploreModelMatrix::VisualizeDesign(variables$metadataFile,input$designformula))
    ##Dit is niet wat het moet zijn, ik weet niet hoe dat gewoon de exploremodelmatrix moet worden
    output$designmatrix <- renderUI({
      lapply(1:length(design()[[2]]),
             function(j){
               renderPlot({design()[[2]][[j]]}, height = 800, width = 1400)
             })
    })

    ##model the data when user clicks button
    observeEvent(input$fitModel,{
      pe2 <- variables$pe2
      #check whether peptideLogNorm actually exists
      if (!"peptideLogNorm" %in% names(pe2)){
        output$model_fitted <- renderText({"Please go through the preprocessing step first"})
      }
      req(pe2[["peptideLogNorm"]])
      pe2 <- msqrob2::msqrob(object = pe2, i = "peptideLogNorm",
                             formula = stats::as.formula(input$designformula),
                             robust=FALSE)
      output$model_fitted <- renderText({"Model fitting complete"})
      variables$pe2 <- pe2
    })

    #Inference tab
    ##output possible parameter names
    output$available_parameters <- renderText({
      colnames(design()$designmatrix)
      }, sep = ", ")

    ##hypothesistest
})
