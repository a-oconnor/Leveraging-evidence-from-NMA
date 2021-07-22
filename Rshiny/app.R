library("netmeta")
library("NlcOptim")
library("kableExtra")
library("dplyr")
library(shiny)
library("DEoptimR")
source(file = "SampleSize.R")

# contrast-level dat input only currently 
ui <-
  fluidPage(titlePanel("Sample Size Calculation to reach a predefined power for a future three-arm study with a new treatment, 
                         based on the existing network"),
            tabsetPanel(
              tabPanel("sample data",
                       fluidRow(
                         
                         column(4,style="background-color:lightgoldenrodyellow",
                                h4("Sample network data"),
                                "This is a sample data consisting 20 studies and 6 treatments for investigators to refer, which is accessible on",
                                tags$a(href="https://github.com/dapengh/Leveraging-evidence-from-NMA/blob/main/Rshiny/sampledat.csv", 
                                       "Github"),
                                ".",
                                tags$br(),
                                "This dataset contains 5 columns that are:",
                                tags$ul(
                                  tags$li("studlab: study label/id"), 
                                  tags$li("treat1: label/number for first treatment"), 
                                  tags$li("treat2: label/number for second treatment"),
                                  tags$li("TE: estimate of treatment effect measured in log odds ratio"), 
                                  tags$li("seTE: standard error of treatment estimate")
                                ),
                                helpText("Note that for trials with more than two arms, all pairwise comparisons are required in the dataset. For example, a three-arm trial should have 3 rows of comparisons, a four-arm trial shold have 6 rows, etc."),
                                # Horizontal line ----
                                tags$hr(),
                                h4("Design the new three arm trial"),
                                tags$b("Please select two treatments used in your future trial"),
                                splitLayout(
                                  uiOutput("treatment1_s"),
                                  uiOutput("treatment2_s")
                                ),
                                tags$b("Please select the type of events"),
                                helpText("If higher probability of event means better, select \"Good events\"; if the opposite, select \"Bad events\"."),
                                radioButtons("eventType_s",
                                             label = NULL,
                                             c("Good events" = "good", "Bad events" = "bad"),
                                             selected = "bad"),
                                tags$b("Please select the type of testing"),
                                helpText("Whether you want to show superiority or non-inferiority of the new treatment to the existing one."),
                                radioButtons("testType_s",
                                             label = NULL,
                                             c("Superiority" = "sup", "Non-Inferiority" = "noninf"),
                                             selected = "sup"),
                                uiOutput("Margin_s")
                         ),
                         
                         column(4,style="background-color:seashell",
                                h4("Confirm the parameters"),
                                radioButtons("trialType_s",
                                             "Do you want to achieve a desired power or to use a fixed sample size?",
                                             c("Desired power" = "power","Fixed sample size" = "size"),
                                             selected = "size"),
                                uiOutput("Restriction_s"),
                                uiOutput("Risk_s"),
                                radioButtons("newtrt_s", "How would you input the effect size of the new treatment?",
                                             c("Risk" = "arisk", "Log odds ratio" = "LOR"), selected = "arisk"),
                                uiOutput("newtrt_size_s"),
                                tags$hr(),
                                br(),
                                h4("Estimate the cost of the treatments"),
                                radioButtons("cost_s",
                                             "Do you want to calculate the cost for each allocation plan?",
                                             c("Yes",
                                               "No"),
                                             selected = "No"),
                                uiOutput("cost_number_s")
                         ),
                         
                         column(4,style="background-color:aliceblue",
                                h4("Information from the previous network"),
                                htmlOutput("txtOutput_s"),
                                tags$hr(),
                                conditionalPanel(condition = "output.powertest_s == 1",
                                                 h4(textOutput("Text_s")),
                                                 tableOutput("tabOutput_s")),
                                conditionalPanel(condition = "output.powertest_s == 0",
                                                 h4("The power you set cannot be reached."))
                         )
                       )),
              tabPanel("your data",
                       fluidRow(
                         
                         column(4,style="background-color:lightgoldenrodyellow",
                                h4("Step1: Upload your network data"),
                                helpText("When using the shiny app to calculate the sample size, please upload your dataset using the same format and column names as the sample data"),
                                fileInput("file1", "Choose CSV File",
                                          multiple = F,
                                          accept = ".csv"),
                                # Horizontal line ----
                                tags$hr(),
                                h4("Step2: Design the new three arm trial"),
                                helpText("This part would appear after the data with correct format is uploaded"),
                                tags$b("Please select two treatments used in your future trial"),
                                splitLayout(
                                  uiOutput("treatment1"),
                                  uiOutput("treatment2")
                                ),
                                tags$b("Please select the type of events"),
                                helpText("If higher probability of event means better, select \"Good events\"; if the opposite, select \"Bad events\"."),
                                radioButtons("eventType",
                                             label = NULL,
                                             c("Good events" = "good", "Bad events" = "bad"),
                                             selected = "bad"),
                                tags$b("Please select the type of testing"),
                                helpText("Whether you want to show superiority or non-inferiority of the new treatment to the existing one."),
                                radioButtons("testType",
                                             label = NULL,
                                             c("Superiority" = "sup", "Non-Inferiority" = "noninf"),
                                             selected = "sup"),
                                uiOutput("Margin")
                         ),
                         
                         column(4,style="background-color:seashell",
                                h4("Step3: Confirm the parameters"),
                                radioButtons("trialType",
                                             "Do you want to achieve a desired power or to use a fixed sample size?",
                                             c("Desired power" = "power","Fixed sample size" = "size"),
                                             selected = "size"),
                                uiOutput("Restriction"),
                                uiOutput("Risk"),
                                radioButtons("newtrt", "How would you input the effect size of the new treatment?",
                                             c("Risk" = "arisk", "Log odds ratio" = "LOR"), selected = "arisk"),
                                uiOutput("newtrt_size"),
                                tags$hr(),
                                br(),
                                h4("Step4: Estimate the cost of the treatments"),
                                radioButtons("cost",
                                             "Do you want to calculate the cost for each allocation plan?",
                                             c("Yes",
                                               "No"),
                                             selected = "No"),
                                uiOutput("cost_number")
                         ),
                         
                         column(4,style="background-color:aliceblue",
                                h4("Information from the previous network"),
                                htmlOutput("txtOutput"),
                                tags$hr(),
                                conditionalPanel(condition = "output.powertest == 1",
                                                 h4(textOutput("Text")),
                                                 tableOutput("tabOutput")),
                                conditionalPanel(condition = "output.powertest == 0",
                                                 h4("The power you set cannot be reached."))
                         )
                       ))
            )
  )

server <- function(input, output,session) {
  
  
  filedata <- reactive({
    infile <- input$file1
    if (!is.null(infile)) {
      # User has not uploaded a file yet
      read.csv(infile$datapath)
    }
  })
  
  arm <- reactive({
    dat <- filedata()
    if (!is.null(dat)){
      nma_old <- netmeta(TE,seTE,treat1,treat2,studlab,data=dat,
                         sm="OR",comb.fixed = T,comb.random = F)
      arms <- unique(c(dat$treat1,dat$treat2))
      list(nma_old = nma_old, arms = arms)
    }
  })
  
  baseline <- reactive({
    dat <- filedata()
    if (!is.null(dat)){
      nma_old <- arm()$nma_old
      trt1 <- input$choice1
      trt2 <- input$choice2
      
      # get baseline risk
      base_trt <- input$choice2
      lor_1 <- nma_old$TE.fixed[base_trt,trt1]
      lor_2 <- nma_old$TE.fixed[base_trt,trt2]
      p1 <- lor2prob(input$baseline_risk,lor_2)
      p2 <- lor2prob(input$baseline_risk,lor_1)
      list(p1 = p1, p2 = p2) 
    }
    
  })
  
  sigma_nma_old <- reactive({
    dat <- filedata()
    if (!is.null(dat)){
      nma_old <- arm()$nma_old
      trt1 <- input$choice1
      trt2 <- input$choice2
      nma_old$seTE.fixed[trt1,trt2]
    }
    
  })
  
  
  output$treatment1 <- renderUI({
    
    dat <- filedata()
    
    if(!is.null(dat)){
      radioButtons(inputId = "choice1",
                   label = "comparator",
                   choices = arm()$arms,
                   selected = arm()$arms[2])
      
    }
  })
  
  arm2 <- reactive({
    dat <- filedata()
    if (!is.null(dat)){
      setdiff(arm()$arms,input$choice1)
    }
  })
  
  output$treatment2 <- renderUI({
    
    dat <- filedata()
    
    if(!is.null(dat)){
      radioButtons(inputId="choice2", 
                   label="baseline",
                   choices= arm2(),
                   selected = character(0))
      
    }
  })
  
  arm_comp <- reactive({
    dat <- filedata()
    if (!is.null(dat)){
      c(input$choice1, input$choice2)
    }
  })
  
  # output$comptrt <- renderUI({
  #   
  #   dat <- filedata()
  #   
  #   if(!is.null(dat)){
  #     radioButtons(inputId="trt_comp", label="Please select the treatment to be compared with the new treatment",
  #                  choices= arm_comp(),
  #                  selected = character(0))
  #     
  #   }
  # })
  
  name_trt1 <- reactive({
    dat <- filedata()
    if (!is.null(dat)){
      input$choice1
    }
  })
  
  name_trt2 <- reactive({
    dat <- filedata()
    if (!is.null(dat)){
      input$choice2
    }
  })
  
  output$Margin <- renderUI({
    if(input$testType == 'noninf'){
      numericInput("margin","Margin",step=0.000001,value = 0.2,max = 0.99999,min = 0.01)
    }
  })
  
  output$Restriction <- renderUI({
    if(input$trialType == 'power'){
      numericInput("power_level","Predefind Power",step=0.000001,value = 0.8,max = 0.99999,min = 0.1)
    }else{
      numericInput("sample_size","Total sample size",step=1,value = 300,max = Inf,min = 1)
    }
  })
  
  
  output$Risk <- renderUI({
    numericInput("baseline_risk","Please input the risk of the baseline treatment",
                 step=0.000001,value = 0.80,max = 0.99999,min = 0.00001)
  })
  
  output$newtrt_size <- renderUI({
    if(input$newtrt == "arisk"){
      tagList(
        numericInput("risk3",paste0("Please input the risk of event of the new treatment"),
                     step=0.000001,value = 0.10,max = 0.99999,min = 0.00001)
      )
    }else{
      tagList(
        numericInput("lor_new",paste0("The log odds ratio of ",name_trt1(), " to the new treatment"),
                     step=0.000001,value = 0,max = Inf,min = -Inf)
      )
    }
  })
  
  output$cost_number <- renderUI({
    if(input$cost=="Yes"){
      tagList(
        numericInput("cost1",paste0("Cost($) per treatment(",name_trt2(),")"),
                     step=0.000001,value = 2,min = 0.00001),
        numericInput("cost2",paste0("Cost($) per treatment(",name_trt1(),")"),
                     step=0.000001,value = 2,min = 0.00001),
        numericInput("cost3",paste0("Cost($) per new treatment"),
                     step=0.000001,value = 2,min = 0.00001),
        numericInput("cost4","Cost($) per animal",
                     step=0.000001,value = 2,min = 0.00001)
        
      ) 
    }
  })
  
  output$txtOutput  = renderUI({
    dat <- filedata()
    if(!is.null(dat)){
      str <- paste0("Standard error of the estimated effect size 
                          between selected two treatments by the previous network is ",
                    round(sigma_nma_old(),4))
      str1 <- paste0("The risk of ",name_trt1()," estimated by the previous network is ",
                     round(baseline()$p2,4))
      HTML(paste(str, str1, sep = '<br/>'))
    }
  })
  
  output$Text = renderText({
    dat <- filedata()
    if(!is.null(dat)){
      if(input$trialType == "power"){
        str <- ("The optimal sample size for each treatment in the future trial")
      }
      if(input$trialType == "size"){
        str <- (paste0("The power in the future trial given total sample size = ", input$sample_size))
      }
      HTML(str)
    }
  })
  
  output$tabOutput <- function() {
    dat <- filedata()
    if(!is.null(dat) & !is.null(input$choice1) & !is.null(input$choice2)){
      eventtype <- input$eventType
      if(input$trialType == "power"){
        sigma <- sigma_nma_old()
        power_level <- input$power_level
        
        risk1 <- baseline()$p1
        risk2 <- baseline()$p2
        risk3 <- ifelse(input$newtrt == 'arisk', input$risk3, lor2prob(baseline()$p2,input$lor_new))
        
        if(input$testType == "sup"){
          samplesize_even = rep(SolveSampleSize_Withprev_equal_sup(risk1,risk2,risk3,sigma,power_level,eventtype)/3,3)
          samplesize = SolveSampleSize_Withprev_sup(risk1,risk2,risk3,sigma,power_level,eventtype)
          #samplesize_single = SolveSampleSize_Single(risk1,risk2,power_level)
          samplesize_single_even = rep(SolveSampleSize_Single_equal_sup(risk1,risk2,risk3,power_level,eventtype)/3,3)
        }else{
          margin <- input$margin
          testtype <- input$testType
          samplesize_even = rep(SolveSampleSize_Withprev_equal(risk1,risk2,risk3,sigma,power_level, margin = margin, testtype = testtype,eventtype)/3,3)
          samplesize = SolveSampleSize_Withprev(risk1,risk2,risk3,sigma,power_level, margin = margin, testtype = testtype,eventtype)
          samplesize_single_even = rep(SolveSampleSize_Single_equal(risk1,risk2,risk3,power_level, margin = margin, testtype = testtype,eventtype)/3,3)
        }
        
        output_dat <- data.frame(n1 = c(samplesize_even[1],samplesize[1],
                                        samplesize_single_even[1]),
                                 n2 = c(samplesize_even[2],samplesize[2],
                                        samplesize_single_even[2]),
                                 n3 = c(samplesize_even[3],samplesize[3],
                                        samplesize_single_even[3]))
        
        output_dat$total <- output_dat$n1+output_dat$n2+output_dat$n3
        
        collapse_rows_dt <- cbind(C1 = c(rep("with the existing network", 2), "without the existing network"),
                                  C2 = c("even","uneven","even"),
                                  output_dat)
        colnames(collapse_rows_dt) <- c("","",name_trt2(),name_trt1(),"New treatment","Total")
        if(input$cost=="Yes"){
          cost <- c(input$cost1,input$cost2,input$cost3,input$cost4)
          costs <- output_dat$n1*cost[1] + output_dat$n2*cost[2] + output_dat$n3*cost[3] + output_dat$total*cost[4]
          collapse_rows_dt[,7] <- costs
          colnames(collapse_rows_dt)[1:2] <- c("","")
          colnames(collapse_rows_dt)[7] <- "Costs ($)"
          collapse_rows_dt %>%
            kbl() %>%
            kable_styling() %>%
            collapse_rows(columns = 1:2, valign = "middle")%>%
            add_header_above(c("","","Sample Size" = 4,"")) 
        }else{
          collapse_rows_dt %>%
            kbl() %>%
            kable_styling() %>%
            collapse_rows(columns = 1:2, valign = "middle") %>%
            add_header_above(c("","","Sample Size" = 4)) 
          
        }
        
      }else{
        sigma <- sigma_nma_old()
        samplesize <- input$sample_size
        
        risk1 <- baseline()$p1
        risk2 <- baseline()$p2
        risk3 <- ifelse(input$newtrt == 'arisk', input$risk3, lor2prob(baseline()$p2,input$lor_new))
        
        if(input$testType == "sup"){
          power_even = SolvePower_Withprev_equal_sup(risk1,risk2,risk3,sigma,samplesize,eventtype)
          power = SolvePower_Withprev_sup(risk1,risk2,risk3,sigma,samplesize,eventtype)
          #samplesize_single = SolveSampleSize_Single(risk1,risk2,power_level)
          power_single_even = SolvePower_Single_equal_sup(risk1,risk2,risk3,samplesize,eventtype)
        }else{
          margin <- input$margin
          testtype <- input$testType
          power_even = SolvePower_Withprev_equal(risk1,risk2,risk3,sigma,samplesize, margin = margin, testtype = testtype,eventtype)
          power = SolvePower_Withprev(risk1,risk2,risk3,sigma,samplesize, margin = margin, testtype = testtype,eventtype)
          power_single_even = SolvePower_Single_equal(risk1,risk2,risk3,samplesize, margin = margin, testtype = testtype,eventtype)
        }
        
        output_dat <- data.frame(n1 = c(power_even[1],power[1],
                                        power_single_even[1]),
                                 n2 = c(power_even[2],power[2],
                                        power_single_even[2]),
                                 n3 = c(power_even[3],power[3],
                                        power_single_even[3]),
                                 power = c(power_even[4],power[4],
                                           power_single_even[4]))
        
        
        collapse_rows_dt <- cbind(C1 = c(rep("with the existing network", 2), "without the existing network"),
                                  C2 = c("even","uneven","even"),
                                  output_dat)
        colnames(collapse_rows_dt) <- c("","",name_trt2(),name_trt1(),"New treatment","Power")
        
        
        if(input$cost=="Yes"){
          cost <- c(input$cost1,input$cost2,input$cost3,input$cost4)
          costs <- output_dat$n1*cost[1] + output_dat$n2*cost[2] + output_dat$n3*cost[3] + input$sample_size*cost[4]
          collapse_rows_dt[,7] <- costs
          colnames(collapse_rows_dt)[1:2] <- c("","")
          colnames(collapse_rows_dt)[7] <- "Costs ($)"
          collapse_rows_dt %>%
            kbl() %>%
            kable_styling() %>%
            collapse_rows(columns = 1:2, valign = "middle")%>%
            add_header_above(c("","","Sample Allocation" = 3,"","")) 
        }else{
          collapse_rows_dt %>%
            kbl() %>%
            kable_styling() %>%
            collapse_rows(columns = 1:2, valign = "middle") %>%
            add_header_above(c("","","Sample Allocation" = 3,"")) 
          
        }
      }
    }
  }
  
  output$powertest <- reactive({
    if(input$trialType == "power"){
      sigma <- sigma_nma_old()
      power_level <- input$power_level
      
      risk1 <- baseline()$p1
      risk2 <- baseline()$p2
      risk3 <- ifelse(input$newtrt == 'arisk', input$risk3, lor2prob(baseline()$p2,input$lor_new))
      if(input$testType == "sup"){
        tmp <- try(SolveSampleSize_Single_equal_sup(risk1,risk2,risk3,power_level,event_type = input$eventType))
        x <- !inherits(tmp, "try-error")
      }else{
        margin <- input$margin
        testtype <- input$testType
        tmp <- try(SolveSampleSize_Single_equal(risk1,risk2,risk3,power_level, margin = margin, testtype = testtype,event_type = input$eventType))
        x <- !inherits(tmp, "try-error")
      }
    }else{
      x <- 1
    }
    return(x)
  })  
  outputOptions(output, 'powertest', suspendWhenHidden=FALSE)
  
  
  
  
  
  
  
  
  
  ##### functions for sample data
  filedata_s <- reactive({
    infile_s <- input$file1_s
    if(is.null(infile_s)){
      read.csv(file = "sampledat.csv")
    }
  })
  
  arm_s <- reactive({
    dat_s <- filedata_s()
    if (!is.null(dat_s)){
      nma_old <- netmeta(TE,seTE,treat1,treat2,studlab,data=dat_s,
                         sm="OR",comb.fixed = T,comb.random = F)
      arms <- unique(c(dat_s$treat1,dat_s$treat2))
      list(nma_old = nma_old, arms = arms)
    }
  })
  
  baseline_s <- reactive({
    dat <- filedata_s()
    if (!is.null(dat)){
      nma_old <- arm_s()$nma_old
      trt1 <- input$choice1_s
      trt2 <- input$choice2_s
      
      # get baseline risk
      base_trt <- input$choice2_s
      lor_1 <- nma_old$TE.fixed[base_trt,trt1]
      lor_2 <- nma_old$TE.fixed[base_trt,trt2]
      p1 <- lor2prob(input$baseline_risk_s,lor_2)
      p2 <- lor2prob(input$baseline_risk_s,lor_1)
      list(p1 = p1, p2 = p2) 
    }
    
  })
  
  sigma_nma_old_s <- reactive({
    dat <- filedata_s()
    if (!is.null(dat)){
      nma_old <- arm_s()$nma_old
      trt1 <- input$choice1_s
      trt2 <- input$choice2_s
      nma_old$seTE.fixed[trt1,trt2]
    }
    
  })
  
  
  output$treatment1_s <- renderUI({
    
    dat <- filedata_s()
    
    if(!is.null(dat)){
      radioButtons(inputId = "choice1_s",
                   label = "comparator",
                   choices = arm_s()$arms,
                   selected = arm_s()$arms[2])
      
    }
  })
  
  arm2_s <- reactive({
    dat <- filedata_s()
    if (!is.null(dat)){
      setdiff(arm_s()$arms,input$choice1_s)
    }
  })
  
  output$treatment2_s <- renderUI({
    
    dat <- filedata_s()
    
    if(!is.null(dat)){
      radioButtons(inputId="choice2_s", 
                   label="baseline",
                   choices= arm2_s(),
                   selected = character(0))
      
    }
  })
  
  arm_comp_s <- reactive({
    dat <- filedata_s()
    if (!is.null(dat)){
      c(input$choice1_s, input$choice2_s)
    }
  })
  
  # output$comptrt_s <- renderUI({
  #   
  #   dat <- filedata_s()
  #   
  #   if(!is.null(dat)){
  #     radioButtons(inputId="trt_comp_s", label=NULL,
  #                  choices= arm_comp_s(),
  #                  selected = character(0))
  #     
  #   }
  # })
  
  name_trt1_s <- reactive({
    dat <- filedata_s()
    if (!is.null(dat)){
      input$choice1_s
    }
  })
  
  name_trt2_s <- reactive({
    dat <- filedata_s()
    if (!is.null(dat)){
      input$choice2_s
    }
  })
  
  output$Margin_s <- renderUI({
    if(input$testType_s == 'noninf'){
      numericInput("margin_s","Margin",step=0.000001,value = 0.2,max = 0.99999,min = 0.01)
    }
  })
  
  output$Restriction_s <- renderUI({
    if(input$trialType_s == 'power'){
      numericInput("power_level_s","Predefind Power",step=0.000001,value = 0.8,max = 0.99999,min = 0.1)
    }else{
      numericInput("sample_size_s","Total sample size",step=1,value = 300,max = Inf,min = 1)
    }
  })
  
  
  output$Risk_s <- renderUI({
    numericInput("baseline_risk_s","Please input the risk of the baseline treatment",
                 step=0.000001,value = 0.80,max = 0.99999,min = 0.00001)
  })
  
  output$newtrt_size_s <- renderUI({
    if(input$newtrt_s == "arisk"){
      tagList(
        numericInput("risk3_s",paste0("Please input the risk of event of the new treatment"),
                     step=0.000001,value = 0.10,max = 0.99999,min = 0.00001)
      )
    }else{
      tagList(
        numericInput("lor_new_s",paste0("The log odds ratio of ",name_trt1_s(), " to the new treatment"),
                     step=0.000001,value = 0,max = Inf,min = -Inf)
      )
    }
  })
  
  output$cost_number_s <- renderUI({
    if(input$cost_s=="Yes"){
      tagList(
        numericInput("cost1_s",paste0("Cost($) per treatment(",name_trt2_s(),")"),
                     step=0.000001,value = 2,min = 0.00001),
        numericInput("cost2_s",paste0("Cost($) per treatment(",name_trt1_s(),")"),
                     step=0.000001,value = 2,min = 0.00001),
        numericInput("cost3_s",paste0("Cost($) per new treatment"),
                     step=0.000001,value = 2,min = 0.00001),
        numericInput("cost4_s","Cost($) per animal",
                     step=0.000001,value = 2,min = 0.00001)
        
      ) 
    }
  })
  
  output$txtOutput_s  = renderUI({
    dat <- filedata_s()
    if(!is.null(dat)){
      str <- paste0("Standard error of the estimated effect size 
                          between selected two treatments by the previous network is ",
                    round(sigma_nma_old_s(),4))
      str1 <- paste0("The risk of ",name_trt1_s()," estimated by the previous network is ",
                     round(baseline_s()$p2,4))
      
      
      HTML(paste(str, str1, sep = '<br/>'))
    }
  })
  
  output$Text_s = renderText({
    dat <- filedata_s()
    if(!is.null(dat)){
      if(input$trialType_s == "power"){
        str <- ("The optimal sample size for each treatment in the future trial")
      }
      if(input$trialType_s == "size"){
        str <- (paste0("The power in the future trial given total sample size = ", input$sample_size_s))
      }
      HTML(str)
    }
  })
  
  output$tabOutput_s <- function() {
    dat <- filedata_s()
    if(!is.null(dat) & !is.null(input$choice1_s) & !is.null(input$choice2_s)){
      eventtype <- input$eventType_s
      if(input$trialType_s == "power"){
        sigma <- sigma_nma_old_s()
        power_level <- input$power_level_s
        
        risk1 <-baseline_s()$p1
        risk2 <- baseline_s()$p2
        risk3 <- ifelse(input$newtrt_s == 'arisk', input$risk3_s, lor2prob(baseline_s()$p2,input$lor_new_s))
        
        if(input$testType_s == "sup"){
          samplesize_even = rep(SolveSampleSize_Withprev_equal_sup(risk1,risk2,risk3,sigma,power_level,eventtype)/3,3)
          samplesize = SolveSampleSize_Withprev_sup(risk1,risk2,risk3,sigma,power_level,eventtype)
          #samplesize_single = SolveSampleSize_Single(risk1,risk2,power_level)
          samplesize_single_even = rep(SolveSampleSize_Single_equal_sup(risk1,risk2,risk3,power_level,eventtype)/3,3)
        }else{
          margin <- input$margin_s
          testtype <- input$testType_s
          samplesize_even = rep(SolveSampleSize_Withprev_equal(risk1,risk2,risk3,sigma,power_level, margin = margin, testtype = testtype,eventtype)/3,3)
          samplesize = SolveSampleSize_Withprev(risk1,risk2,risk3,sigma,power_level, margin = margin, testtype = testtype,eventtype)
          samplesize_single_even = rep(SolveSampleSize_Single_equal(risk1,risk2,risk3,power_level, margin = margin, testtype = testtype,eventtype)/3,3)
        }
        
        output_dat <- data.frame(n1 = c(samplesize_even[1],samplesize[1],
                                        samplesize_single_even[1]),
                                 n2 = c(samplesize_even[2],samplesize[2],
                                        samplesize_single_even[2]),
                                 n3 = c(samplesize_even[3],samplesize[3],
                                        samplesize_single_even[3]))
        
        output_dat$total <- output_dat$n1+output_dat$n2+output_dat$n3
        
        collapse_rows_dt <- cbind(C1 = c(rep("with the existing network", 2), "without the existing network"),
                                  C2 = c("even","uneven","even"),
                                  output_dat)
        colnames(collapse_rows_dt) <- c("","",name_trt2_s(),name_trt1_s(),"New treatment","Total")
        if(input$cost_s=="Yes"){
          cost <- c(input$cost1_s,input$cost2_s,input$cost3_s,input$cost4_s)
          costs <- output_dat$n1*cost[1] + output_dat$n2*cost[2] + output_dat$n3*cost[3] + output_dat$total*cost[4]
          collapse_rows_dt[,7] <- costs
          colnames(collapse_rows_dt)[1:2] <- c("","")
          colnames(collapse_rows_dt)[7] <- "Costs ($)"
          collapse_rows_dt %>%
            kbl() %>%
            kable_styling() %>%
            collapse_rows(columns = 1:2, valign = "middle")%>%
            add_header_above(c("","","Sample Size" = 4,"")) 
        }else{
          collapse_rows_dt %>%
            kbl() %>%
            kable_styling() %>%
            collapse_rows(columns = 1:2, valign = "middle") %>%
            add_header_above(c("","","Sample Size" = 4)) 
          
        }
        
      }else{
        sigma <- sigma_nma_old_s()
        samplesize <- input$sample_size_s
        
        risk1 <-baseline_s()$p1
        risk2 <- baseline_s()$p2
        risk3 <- ifelse(input$newtrt_s == 'arisk', input$risk3_s, lor2prob(baseline_s()$p2,input$lor_new_s))
        
        if(input$testType_s == "sup"){
          power_even = SolvePower_Withprev_equal_sup(risk1,risk2,risk3,sigma,samplesize,eventtype)
          power = SolvePower_Withprev_sup(risk1,risk2,risk3,sigma,samplesize,eventtype)
          #samplesize_single = SolveSampleSize_Single(risk1,risk2,power_level)
          power_single_even = SolvePower_Single_equal_sup(risk1,risk2,risk3,samplesize,eventtype)
        }else{
          margin <- input$margin_s
          testtype <- input$testType_s
          power_even = SolvePower_Withprev_equal(risk1,risk2,risk3,sigma,samplesize, margin = margin, testtype = testtype,eventtype)
          power = SolvePower_Withprev(risk1,risk2,risk3,sigma,samplesize, margin = margin, testtype = testtype,eventtype)
          power_single_even = SolvePower_Single_equal(risk1,risk2,risk3,samplesize, margin = margin, testtype = testtype,eventtype)
        }
        
        output_dat <- data.frame(n1 = c(power_even[1],power[1],
                                        power_single_even[1]),
                                 n2 = c(power_even[2],power[2],
                                        power_single_even[2]),
                                 n3 = c(power_even[3],power[3],
                                        power_single_even[3]),
                                 power = c(power_even[4],power[4],
                                           power_single_even[4]))
        
        
        collapse_rows_dt <- cbind(C1 = c(rep("with the existing network", 2), "without the existing network"),
                                  C2 = c("even","uneven","even"),
                                  output_dat)
        colnames(collapse_rows_dt) <- c("","",name_trt2_s(),name_trt1_s(),"New treatment","Power")
        
        
        if(input$cost_s=="Yes"){
          cost <- c(input$cost1_s,input$cost2_s,input$cost3_s,input$cost4_s)
          costs <- output_dat$n1*cost[1] + output_dat$n2*cost[2] + output_dat$n3*cost[3] + input$sample_size_s*cost[4]
          collapse_rows_dt[,7] <- costs
          colnames(collapse_rows_dt)[1:2] <- c("","")
          colnames(collapse_rows_dt)[7] <- "Costs ($)"
          collapse_rows_dt %>%
            kbl() %>%
            kable_styling() %>%
            collapse_rows(columns = 1:2, valign = "middle")%>%
            add_header_above(c("","","Sample Allocation" = 3,"","")) 
        }else{
          collapse_rows_dt %>%
            kbl() %>%
            kable_styling() %>%
            collapse_rows(columns = 1:2, valign = "middle") %>%
            add_header_above(c("","","Sample Allocation" = 3,"")) 
          
        }
      }
    }
  }
  
  output$powertest_s <- reactive({
    if(input$trialType_s == "power"){
      sigma <- sigma_nma_old_s()
      power_level <- input$power_level_s
      
      risk1 <-baseline_s()$p1
      risk2 <- baseline_s()$p2
      risk3 <- ifelse(input$newtrt_s == 'arisk', input$risk3_s, lor2prob(baseline_s()$p2,input$lor_new_s))
      if(input$testType_s == "sup"){
        tmp <- try(SolveSampleSize_Single_equal_sup(risk1,risk2,risk3,power_level,event_type = input$eventType_s))
        x <- !inherits(tmp, "try-error")
      }else{
        margin <- input$margin_s
        testtype <- input$testType_s
        tmp <- try(SolveSampleSize_Single_equal(risk1,risk2,risk3,power_level, margin = margin, testtype = testtype,event_type = input$eventType_s))
        x <- !inherits(tmp, "try-error")
      }
    }else{
      x <- 1
    }
    return(x)
  })  
  outputOptions(output, 'powertest_s', suspendWhenHidden=FALSE)
  
  
  
}

shinyApp(ui = ui, server = server)

