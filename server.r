library(shiny)
library(shinyjs)
library(tidyverse)
library(DT)
library(ggplot2)

setwd("Y:/IRB_00096551/egfr-table-files")

#shiny::runExample("10_download")
#master.file <- read_tsv("/Users/u6013142/Documents/Projects/eGFR/tables/vitals/bp_bmi_slope_cohort_summary.tsv", na = "NA")
# THIS DATA COULD BE A PULL FROM AROUND THE TIME SCOTT WAS WORKING ON HIS PAPER
#master.file <- read_tsv("slope_cohort_summary_right_censored_vitals_labs_meds_insulin.tsv")
master.file <- master.file.temp
bmi.file <- read_tsv("bmi_study_pop.tsv")
bp.file <- read_tsv("bp_study_pop.tsv")
outpatient.egfr.3om <- read_tsv("slope_cohort_egfr.tsv")
#master.file <- read_tsv("/Users/u6013142/Documents/Projects/eGFR/tables/masterFile18AndOlder.tsv", na = ".")
#outpatient.egfr.3om <- read_tsv("/Users/u6013142/Documents/Projects/eGFR/tables/outpatient_3OrMoreEgfr_18AndOlder.tsv")
#master.file <- read_tsv("/Users/u6013142/Documents/Projects/eGFR/tables/not_ckd3_3months_summary.tsv")
#outpatient.egfr.3om <- read_tsv("/Users/u6013142/Documents/Projects/eGFR/tables/not_ckd3_3months_egfr.tsv")
#master.file <- master.file %>% filter(followup_yrs > 0.08333)
#idstoexclude <- unique(outpatient.egfr.3om$dist_id) %in% master.file$dist_id
#outpatient.egfr.3om <- outpatient.egfr.3om[idstoexclude,]
outpatient.egfr.3om$dist_id <- as.character(outpatient.egfr.3om$dist_id)
bmi.file$DIST_ID <- as.character(bmi.file$DIST_ID)
bp.file$DIST_ID <- as.character(bp.file$DIST_ID)

#master.file <- master.file %>% select(-c(uncorrected_adj_r_2, corrected_adj_r_2))

eGFR_outpatient_data <- data.frame(id = as.integer(master.file$dist_id),
                                   sex = as.character(master.file$sex),
                                   ethnicity = as.character(master.file$ethn),
                                   first_age = master.file$age_first_egfr,
                                   first_egfr = master.file$first_egfr,
                                   first_ckd_stage = as.character(master.file$first_ckd_stage),
                                   last_age = master.file$age_last_egfr,
                                   last_egfr = master.file$last_egfr,
                                   last_ckd_stage = as.character(master.file$last_ckd_stage),
                                   number_egfrs = master.file$n_egfr,
                                   followup_years = master.file$followup_yrs,
                                   #slope = master.file$corrected_slope,
                                   slope = master.file$uncorrected_slope,
                                   r2 = master.file$corrected_r2,
                                   obsv_esrd_followup = master.file$obsv_esrd_fu_yrs)
                                   #esrd_followup = master.file$icd_esrd_fu_yrs,
                                   #dialysis_followup = master.file$icd_dialysis_fu_yrs,
                                   #kidney_followup = master.file$icd_kidney_tx_fu_yrs)

eGFR_outpatient_data$sex <- as.character(eGFR_outpatient_data$sex)
eGFR_outpatient_data$ethnicity <- as.character(eGFR_outpatient_data$ethnicity)
totalPatients <- nrow(eGFR_outpatient_data)
GFR_outpatient_data <- eGFR_outpatient_data[1:100,]

shinyServer(function(input, output, session) {
  
  
  #output$masterTable <- renderDataTable(eGFR_outpatient_data, 
  #                                      options = list(paging = TRUE, columnDefs = list(list(className = 'dt-center', targets = 0:15)), scrollX = TRUE), 
  #                                      filter = 'top', 
  #                                      rownames = FALSE)
  output$masterTable <- DT::renderDataTable(eGFR_outpatient_data, 
                                        options = list(paging = TRUE, columnDefs = list(list(className = 'dt-center', targets = 0:15)), scrollX = TRUE))
  
  proxyMasterTable = dataTableProxy("masterTable")
  
  shinyjs::hide("sbsOptions")
  fullList <- list("Slope" = "slope", "Number of eGFR Measurements" = "number_egfrs","Number of Follow Up Years" = "followup_years","Sex (0 = Male, 1 = Female)" = "sex","Ethnicity (0 = White, 1 = Black)" = "ethnicity","First Age" = "first_age","First eGFR Measurement" = "first_egfr","First CKD Stage" = "first_ckd_stage","Last Age" = "last_age","Last eGFR Measurement" = "last_egfr","Last CKD Stage" = "last_ckd_stage","Renal ICD" = "renal_icd","Diabetes Type" = "dm_type")
  sexList <- list("Slope" = "slope", "eGFR" = "ckd_epi_egfr", "Number of eGFR Measurements" = "number_egfrs", "Number of Follow Up Years" = "followup_years", "Ethnicity (0 = White, 1 = Black)" = "ethnicity", "First Age" = "first_age", "First eGFR Measurement" = "first_egfr", "First CKD Stage" = "first_ckd_stage", "Last Age" = "last_age", "Last eGFR Measurement" = "last_egfr", "Last CKD Stage" = "last_ckd_stage", "Renal ICD" = "renal_icd","Diabetes Type" = "dm_type")
  ethnList <- list("Slope" = "slope", "eGFR" = "ckd_epi_egfr", "Number of eGFR Measurements" = "number_egfrs", "Number of Follow Up Years" = "followup_years", "Sex (0 = Male, 1 = Female)" = "sex", "First Age" = "first_age", "First eGFR Measurement" = "first_egfr", "First CKD Stage" = "first_ckd_stage", "Last Age" = "last_age", "Last eGFR Measurement" = "last_egfr", "Last CKD Stage" = "last_ckd_stage", "Renal ICD" = "renal_icd","Diabetes Type" = "dm_type")
                           
  observe({
    if (input$histType == "sbs") {
      shinyjs::show("sbsOptions")
      if (input$sbsOptions == "sex") {
        updateRadioButtons(session, "histOptions", "Histogram Options", choices = sexList)
      } else {
        updateRadioButtons(session, "histOptions", "Histogram Options", choices = ethnList)
      }
    } else {
      shinyjs::hide("sbsOptions")
      updateRadioButtons(session, "histOptions", "Histogram Options", choices = fullList)
    }
  })
  
  observeEvent(input$clear, {
    proxyMasterTable %>% selectRows(NULL)
  })
  observeEvent(input$all, {
    proxyMasterTable %>% selectRows(input$masterTable_rows_all)
  })
  observeEvent(input$page, {
    proxyMasterTable %>% selectRows(input$masterTable_rows_current)
  })

  
  output$debug <- renderPrint({
    selectedPatients = input$masterTable_rows_selected
    pagePatients = input$masterTable_rows_all
    totalFiltered <- length(pagePatients)
    cat(paste0("Showing ", totalFiltered, " out of ", totalPatients))
    cat("\n\n")
    cat('Selected Patients:\n\n')
    if (length(selectedPatients) < 100) {
      cat(master.file$dist_id[selectedPatients], sep = ", ")
    } else {
      cat(paste0(master.file$dist_id[selectedPatients[1:100]]), "and more....")
    }
    cat('\n\nNumber of selected patients:\n\n')
    cat(length(selectedPatients), sep = ', ')
  })
  
  currentPlot <- reactive({
    pagePatients <- input$masterTable_rows_all
    selectedPatients <- input$masterTable_rows_selected
    radioChoice <- input$histOptions
    radioFactor <- input$sbsOptions
    binWidth <- input$binwidth
    if (length(selectedPatients) == 0) {
      if (input$histType == "single") {
        if (radioChoice != "ckd_epi_egfr" & radioChoice != "slope" & radioChoice != "sex" & radioChoice != "ethnicity" & radioChoice != "first_ckd_stage" & radioChoice != "last_ckd_stage" & radioChoice != "renal_icd" & radioChoice != "dm_type") {
          print(ggplot() + aes(eGFR_outpatient_data[pagePatients,radioChoice]) + geom_histogram(binwidth = binWidth))
        } else if (radioChoice == "ckd_epi_egfr") {
          print(ggplot() + aes(outpatient.egfr.3om[outpatient.egfr.3om$dist_id %in% eGFR_outpatient_data[pagePatients,]$ID, radioChoice]) + geom_histogram(binwidth = binWidth))
        } else if (radioChoice == "slope") {
          downPlot = function() {print(ggplot(eGFR_outpatient_data[pagePatients,], aes(x = slope)) + geom_histogram(breaks = c(min(master.file$uncorrected_slope), seq(-100, 100, by=binWidth))) + xlim(c(-100,100)) + theme_bw())}
          downPlot()
        } else if (radioChoice == "sex") {
          print(ggplot(eGFR_outpatient_data[pagePatients,], aes(sex, ..count..)) + geom_bar(aes(fill = sex)) + scale_fill_discrete(name = radioChoice, breaks = c("0","1"), labels = c("Male","Female")))
        } else if (radioChoice == "ethnicity") {
          print(ggplot(eGFR_outpatient_data[pagePatients,], aes(ethnicity, ..count..)) + geom_bar(aes(fill = ethnicity)) + scale_fill_discrete(name = radioChoice, breaks = c("0","1"), labels = c("White","Black")))
        } else if (radioChoice == "first_ckd_stage") {
          downPlot = function() {print(ggplot(eGFR_outpatient_data[pagePatients,], aes(first_ckd_stage, ..count..)) + geom_bar(aes(fill = first_ckd_stage)) + scale_fill_discrete(name = "CKD Stage"))}#, breaks = c("0","1"), labels = c("","")))
          downPlot()
        } else if (radioChoice == "last_ckd_stage") {
          print(ggplot(eGFR_outpatient_data[pagePatients,], aes(last_ckd_stage, ..count..)) + geom_bar(aes(fill = last_ckd_stage)) + scale_fill_discrete(name = "CKD Stage"))#, breaks = c("0","1"), labels = c("","")))
        } else if (radioChoice == "renal_icd") {
          print(ggplot(eGFR_outpatient_data[pagePatients,], aes(renal_icd, ..count..)) + geom_bar(aes(fill = renal_icd)) + scale_fill_discrete(name = radioChoice, breaks = c("0","1"), labels = c("No","Yes")))
        } else if (radioChoice == "dm_type") {
          print(ggplot(eGFR_outpatient_data[pagePatients,], aes(dm_type, ..count..)) + geom_bar(aes(fill = dm_type)) + scale_fill_discrete(name = "DM Type"))#, breaks = c("0","1"), labels = c("White","Black")))
        }
      } else {
        if (radioChoice != "ckd_epi_egfr" & radioChoice != "slope") {
          #print(ggplot() + aes(eGFR_outpatient_data[pagePatients,c(radioChoice,"sex")], x = radioChoice) + geom_histogram(binwidth = binWidth) + facet_grid(~sex)) + theme_bw()
          print(ggplot(eGFR_outpatient_data[pagePatients,], aes(x = radioChoice)) + stat_count(width = binWidth) + facet_grid(~sex) + theme_bw())
        } else if (radioChoice == "ckd_epi_egfr") {
          print(ggplot() + aes(outpatient.egfr.3om[outpatient.egfr.3om$dist_id %in% eGFR_outpatient_data[pagePatients,]$id, c(radioChoice,"sex")]) + geom_histogram(binwidth = binWidth) + facet_grid(~sex))
        } else {
          print(ggplot(eGFR_outpatient_data[pagePatients,], aes(x = slope)) + geom_histogram(breaks = c(min(master.file$uncorrected_slope), seq(-100, 100, binWidth), max(master.file$uncorrected_slope))) + facet_grid(as.formula(paste0("~",radioFactor))) + theme_bw())
        }
      }
    } else if (length(selectedPatients) == 1) {
      print(ggplot() + aes(outpatient.egfr.3om[outpatient.egfr.3om$dist_id %in% eGFR_outpatient_data[selectedPatients,]$id, "ckd_epi_egfr"]) + geom_histogram(binwidth = binWidth))
    } else if (length(selectedPatients) > 1) {
      print(ggplot() + aes(outpatient.egfr.3om[outpatient.egfr.3om$dist_id %in% eGFR_outpatient_data[selectedPatients,]$id, "ckd_epi_egfr"]) + geom_histogram(binwidth = binWidth))
    }
  })
  
  output$hist <- renderPlot({
    currentPlot()
  })
   
  output$download <- downloadHandler(
    filename = function() {
      paste("test.png")
    },
    content = function(file) {
      pagePatients <- input$masterTable_rows_all
      binWidth <- input$binwidth
      ggsave(file, currentPlot(), device = "png", dpi = 500)
      #print(ggplot(eGFR_outpatient_data[pagePatients,], aes(x = slope)) + geom_histogram(breaks = c(min(master.file$uncorrected_slope), seq(-100, 100, by=binWidth))) + xlim(c(-100,100)) + theme_bw())
      #ggsave(file, plotInput(), device = "png")
      #dev.off()
    }
  )
  
  output$histSum <- renderPrint({
    selectedPatients <- input$masterTable_rows_selected
    if (length(selectedPatients) == 0) {
      cat(paste0("Here's a histogram of everyones ", input$histOptions))
    } else if (length(selectedPatients) == 1) {
      print(paste0("Here's a histogram of the eGFR measurements of patient ", eGFR_outpatient_data[selectedPatients,]$id))
    } else if (length(selectedPatients) > 1) {
      if (length(selectedPatients) < 100) {
        print(eGFR_outpatient_data[selectedPatients,]$id, sep = ", ")
      } else {
        print(paste0(eGFR_outpatient_data[selectedPatients[1:100],]$id, sep = ", ", " and more...."))
      }
    }
  })
  
  output$line <- renderPlot({
    selectedPatients <- input$masterTable_rows_selected
    xAxis <- input$lineOptions
    lab <- input$labOptions
    ids <- eGFR_outpatient_data[selectedPatients, ]$id
    if (lab == "egfr") {
      if (length(selectedPatients) == 1) {
        print(ggplot(outpatient.egfr.3om[outpatient.egfr.3om$dist_id %in% ids,], aes(x=futime_yrs, y=ckd_epi_egfr)) + 
          geom_line() + geom_point() +
          ylim(c(0, 200)) +
          labs(x = "Years of Follow-up", y = "Estimated eGFR") +
          ggtitle(paste0("Patient ", ids)))
      } else if (length(selectedPatients) > 1) {
        print(ggplot(outpatient.egfr.3om[outpatient.egfr.3om$dist_id %in% ids,], aes(x=futime_yrs, y=ckd_epi_egfr, color=dist_id)) +
          geom_line() + geom_point() +
          ylim(c(0, 200)) +
          labs(x = "Years", y = "eGFR") +
          ggtitle(paste0("eGFR Decline Starting at Each CKD Stage, Patients ",paste(ids, collapse = ", "))))
      }
    } else if (lab == "bmi") {
      if (length(selectedPatients) == 1) {
        print(ggplot(bmi.file[bmi.file$DIST_ID %in% ids, ], aes(x = BMI_N, y = BMI_CALC_BMI)) +
                geom_line() + geom_point() +
                ylim(c(0, 100)) +
                labs(x = "Record N", y = "BMI (kg/m^2)") +
                ggtitle(paste0("Patient ", ids)))
      } else if (length(selectedPatients) > 1) {
        print(ggplot(bmi.file[bmi.file$DIST_ID %in% ids, ], aes(x = BMI_N, y = BMI_CALC_BMI, color = DIST_ID)) +
                geom_line() + geom_point() +
                ylim(c(0, 100)) +
                labs(x = "Record Number", y = "BMI (kg/m^2)") +
                ggtitle(paste0("BMI of Patients ", paste(ids, collapse = ", "))))
      }
    }
  })
  
  output$lineExplain <- renderPrint({
    selectedPatients <- input$masterTable_rows_selected
    lab <- input$labOptions
    if (length(selectedPatients) == 0) {
      cat("Cannot generate line plot until you select at least one patient from the table")
    } else if (lab == "egfr") {
      outpatient.egfr.3om[outpatient.egfr.3om$dist_id %in% eGFR_outpatient_data$id[selectedPatients],] %>% print(n=nrow(.))
    } else if (lab == "bmi") {
      bmi.file[bmi.file$DIST_ID %in% eGFR_outpatient_data$id[selectedPatients], ] %>% print(n=nrow(.))
    } else if (lab == "bp") {
      bp.file[bp.file$DIST_ID %in% eGFR_outpatient_data$id[selectedPatients], ] %>% print(n=nrow(.))
    }
  })
  
  output$regression <- renderPlot({
    selectedPatients <- input$masterTable_rows_selected
    ids <- eGFR_outpatient_data[selectedPatients,]$id
    #ids <- eGFR_outpatient_data[40,]$id
    #ids <- eGFR_outpatient_data[121:124,]$id
    selected.sum <- eGFR_outpatient_data[selectedPatients,]
    #selected.sum <- eGFR_outpatient_data[40,]
    #selected.sum <- eGFR_outpatient_data[121:124,]
    ##### START HERE GETTING THE REGRESSIONS TO WORK AND THEN ADDING THE CONDITIONALS INTO ALL THE PLOTS!!!!
    pat.egfr <- outpatient.egfr.3om[outpatient.egfr.3om$dist_id %in% ids,]
    
    #browser()
    timeFrameChoice <- input$regressionOptions
    time = pat.egfr$futime_yrs
    if(timeFrameChoice == "days") {
      time = pat.egfr$futime_days
    }
    print(timeFrameChoice)
    if (length(selectedPatients) == 1) {
      #if (master.file[selectedPatients,]$esrd_icd == 1 | master.file[selectedPatients,]$dialysis_icd == 1 | master.file[selectedPatients,]$kidney_tx_icd == 1) {
      if (!is.na(selected.sum$obsv_esrd_followup) | !is.na(selected.sum$esrd_followup) | !is.na(selected.sum$dialysis_followup) | !is.na(selected.sum$kidney_followup)) {
        earliest.event <- min(c(selected.sum$obsv_esrd_followup, selected.sum$esrd_followup, selected.sum$dialysis_followup, selected.sum$kidney_followup), na.rm = TRUE)
        
        ggplot(pat.egfr, aes(futime_yrs, ckd_epi_egfr)) +
        #ggplot(pat.egfr, aes(time, ckd_epi_egfr)) +
          geom_point(data = subset(pat.egfr, futime_yrs <= earliest.event)) +
          geom_point(data = subset(pat.egfr, futime_yrs > earliest.event), alpha = 0.35) +
          #geom_vline(aes(xintercept = as.numeric(earliest.event)), linetype=4, col="red") +
          #geom_smooth(data = subset(pat.egfr, futime_yrs <= earliest.event), size = 1, color = "royalblue3", method = "lm", se = FALSE) +
          geom_smooth(data = subset(pat.egfr, futime_yrs <= earliest.event), size = 1, color = "royalblue3", method = "lm", se = FALSE) +
          labs(x = "Years", y = "eGFR") +
          ylim(c(0, 200)) +
          ggtitle(paste0("Estimated Linear eGFR Decline of Patient ", ids)) + {
            if (!is.na(selected.sum$obsv_esrd_followup)) {
              geom_vline(aes(xintercept = as.numeric(selected.sum$obsv_esrd_followup), col = "Observed ESRD"), size = 0.75, linetype = 4)
            }
          } + {
            if (!is.na(selected.sum$esrd_followup)) {
              geom_vline(aes(xintercept = as.numeric(selected.sum$esrd_followup), col = "ICD ESRD"), size = 0.75, linetype = 4)
            }
          } + {
            if (!is.na(selected.sum$dialysis_followup)) {
              geom_vline(aes(xintercept = as.numeric(selected.sum$dialysis_followup), col = "ICD Dialysis"), size = 0.75, linetype = 4)
            }
          } + {
            if (!is.na(selected.sum$kidney_followup)) {
              geom_vline(aes(xintercept = as.numeric(selected.sum$kidney_followup), col = "ICD Kidney Tx"), size = 0.75, linetype = 4)
            }
          } + { 
            scale_color_manual(name = "Event", values = c(`Observed ESRD` = "mediumorchid4", `ICD ESRD` = "tomato2", `ICD Dialysis` = "springgreen4", `ICD Kidney Tx` = "tan2"))
          }
        
      } else {
        #ggplot(outpatient.egfr.3om[outpatient.egfr.3om$dist_id %in% ids,], aes(futime_yrs, ckd_epi_egfr)) +
        #  geom_point() +
        #  geom_smooth(method = "lm",se = FALSE, color = "royalblue3") +
        #  ylim(c(0, 200)) +
        #  labs(x = "Years", y = "eGFR") +
        #  ggtitle(paste0("Linear Regression Plot of Patient ", ids))
        ggplot(outpatient.egfr.3om[outpatient.egfr.3om$dist_id %in% ids,], aes(time, ckd_epi_egfr)) +
          geom_point() +
          geom_smooth(method = "lm",se = FALSE, color = "royalblue3") +
          ylim(c(0, 200)) +
          labs(x = str_to_title(timeFrameChoice), y = "eGFR") +
          ggtitle(paste0("Linear Regression Plot of Patient ", ids))
      }
    } else if (length(selectedPatients) > 1) {
      selected.sum$earliest_event <- apply(selected.sum, 1, function(x) min(c(x[[14]], x[[15]], x[[16]]), na.rm = TRUE))
      selected.sum$earliest_event <- sapply(selected.sum$earliest_event, function(x) if(is.infinite(x)){ x <- NA}else{x <- x})
      selected.sum$earliest_event <- as.numeric(selected.sum$earliest_event)
      distid_earliestevent <- as.data.frame(cbind(selected.sum$id, selected.sum$earliest_event))
      colnames(distid_earliestevent) <- c("dist_id","earliest_event")
      regression.egfrs <- merge(pat.egfr, distid_earliestevent, all = TRUE)
      regression.egfrs <- regression.egfrs %>% filter(!is.na(earliest_event) & futime_yrs <= earliest_event | is.na(earliest_event))

      ggplot(regression.egfrs) +
        geom_jitter(aes(x=futime_yrs,y=ckd_epi_egfr,color=dist_id)) + 
        geom_smooth(aes(x=futime_yrs,y=ckd_epi_egfr,color=dist_id), method = lm, se = FALSE) +
        facet_wrap(~dist_id, scales = "free_x") +
        labs(x = "Years", y = "eGFR") +
        ylim(c(0, 200)) +
        ggtitle(paste0("Non-Linear eGFR Decline in Patients ", paste(ids, collapse = ", ")))
    }
  })
  
  output$regressionExplain<- renderPrint({
    selectedPatients <- input$masterTable_rows_selected
    ids <- eGFR_outpatient_data[selectedPatients,]$ID
    selected.sum <- eGFR_outpatient_data[selectedPatients,]
    if (length(selectedPatients) == 0) {
      cat("Cannot generate line plot until you select at least one patient from the table")
    } else if (length(selectedPatients == 1)) {
      outpatient.egfr.3om[outpatient.egfr.3om$dist_id %in% eGFR_outpatient_data$id[selectedPatients],] %>% print(n=nrow(.))
    } else {
      cat("Can only do one regression at a time at this point. Please select only 1 patient.") 
    }
  })
})
