#   TIDD: peptide rescoring s/w
#
#   Written by H. Li <hllee@hanyang.ac.kr>
#
#   Copyright (C) 2022 BIS Labs, Hanyang Univ. Korea
#
#   TIDD is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   TIDD is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with TIDD.  If not, see <http://www.gnu.org/licenses/>.


library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinyFiles)


# Shiny UI -------

feature_list<-
  list(
    h5("Feature selction:"),
    tags$div(align = 'left', 
             class = 'multicol', 
             checkboxGroupInput(inputId  = 'check1', 
                                label    = "Select the features:", 
                                choices  = c(""),
                                selected = ,
                                inline   = FALSE))
  )
ui <- fluidPage(
  
  navbarPage("TIDD v1.1",
             
             theme=shinytheme("flatly"), 
             id = "main_navbar",
             
             tabPanel(
               "About",
               tabBox(width = 12,
                      tabPanel("About TIDD",
                               p("Tool-Independent and Data-Dependent PSM rescoring method."),
                               br(),
                               p("TIDD is a universal post-processing tool which supports confident peptideidentifications regardless of the search engine adopted. TIDD can work for any(including newly developed) search engines because it calculates universalfeatures that assess peptide-spectrum match quality while it allows additional featuresprovided by search engines (or users) as well."),
                               p("Here, we support two types of TIDD version -- a simple GUI based TIDD and command version of TIDD.
"),
                               
                               hr(style = "border-top: 1px solid #000000;"),
                               
                               p("GUI"),
                               p("1. Data Management: extract features and check feature densities"),
                               p("2. Model fitting: fit SVM models"),
                               
                               br(),
                               h6("*Computational resources supported by free version of \"R studio cloud\" is limited, so it takes a little time for feature extraction."),
                               h6("*Sample file description:"),
                               h6("1.test_required.tsv: contain 6 required features for feature extraction"),
                               h6("2.*_feature.stv: contain common features except calculated Xcorr"),
                               h6("3.*_xcorrColAdded.tsv: add calculated Xcorr feature"),
                               h6("4.*_ModelFitting_Results.tsv: add rescoring score feature"),
                               br(),
                               
                               hr(style = "border-top: 1px solid #000000;"),
                               p("Required:"),
                               p("1. R (v4.1 or above) packages: shiny, shinythemes, shinydashboard,e1071,ROCR"),
                               p("2. java jre/jdk v1.8 or above"),
                               
                               hr(style = "border-top: 1px solid #000000;"),
                               
                               p("Tutorial:"),
                               tags$a(href="https://docs.google.com/document/d/168oGySS15xrobeqTY54_QLMgWQofnPIApx6ma5ntGa8/edit?usp=sharing", "How to run TIDD-GUI"),
                               br(),
                               p("Download: "),
                               tags$a(href="https://github.com/HonglanLi/TIDD.git", "1.Download code from github"),
                               br(),
                               tags$a(href="https://rstudio.cloud/spaces/178915/project/2994889", "2.Download code from R studio cloud"),
                               br(),
                               tags$a(href="https://honglan-li.shinyapps.io/project/", "3.Test GUI from shinyapps.io"),
                               
                      ),
                      # tabPanel("Install & Download",
                      #          
                      # )
                      
               )
             ),
             tabPanel(
               "File upload & download",
               
               column(width=6,
                      
                      h5("Upload your own files to web"),
                      hr(style = "border-top: 1px solid #000000;"),
                      
                      h6("*.Tab-deliminated PSM file must contains:"),
                      h6("1.file name col.: must same as the corresponding mgf file name "),
                      h6("2.scan number col: scan number in the *mgf file"),
                      h6("3.charge"),
                      h6("4.precursor mass: not precursor mass over charge"),
                      h6("5.peptide: must contain modification mass"),
                      h6("6.protein: list proteins with separator \";\" and must contain the pre&next amino acid info at least at one protein" ),
                      
                      h6("E.g."), 
                      h6("file_name /scan_num charge / precursor_mass    /   peptide   /              protein"),
                      h6("test.mgf / 5169   /   2  / 944.481541  /  AAC+57.021464NLLQR  /  sp|Q01813|PFKAP_HUMAN(pre=K,next=G);sp|Q01813-2|PFKAP_HUMAN   "),
                      
                      br(), 
                      
                      h6("*.MGF format spectra file must contains:"),
                      h6("1.TITLE=$filename+\".\"+$scan+\".\"+$scan+\".\"+$charge"),
                      
                      h6("E.g. at \"test.mgf\" file"),
                      h6("BEGIN IONS"),
                      h6("..."),
                      h6("TITLE=test.5169.5169.2"),
                      
                      br(),
                      radioButtons("filetype","Select one file type to upload (Max size 2G per one time)",choices=c("PSM file (tab-deliminated file)"="1","MS/MS spectra file (MGF)"="2")),
                      fileInput("upload","Upload", multiple=TRUE),
                      h6("Uploaded PSM or MGF files will be stored in \"PSM or Local\" directory respectively."),
                      
                      tableOutput("files")
               ),
               
               column(width=6,
                      
                      h5("Download files from web"),
                      hr(style = "border-top: 1px solid #000000;"),
                      
                      h5("Select files to download:"),
                      selectInput('downloadpsm',"PSM file: ",choices = c( list.files('data/PSM/.'))),
                      downloadButton('downloadData', 'Download *.tsv'),
                      
                       h5("Select mgf files to remove:"),
                       selectInput('downloadmgf','MGF files:', choices = c("",list.files('data/MGF/test/.')),selected=c(""), multiple = TRUE),
                       downloadButton('downloadmgfData', 'Download *.mgf'),
                      
                      
                      br(),
                      br(),
                      textOutput("downloadinfo"),
                      
                      
               )
               
               
               
             ),
             tabPanel(
               "Data Management",
               fluidRow(column(
                 3, uiOutput("dart_currentTime", container = span)
               ),
               column(
                 3, uiOutput("dart_last_update", container = span)
               )), 
               sidebarLayout(
                 sidebarPanel(
                   width = 3,
                   
                   h5("File Loading"),
                   hr(style = "border-top: 1px solid #000000;"),
                   h6("Peptide-spectrum matching key info.:"),
                   h6("1. filename"),
                   h6("2. scan number"),
                   h6("3. charge"),
                   
                   br(),
                   selectInput('psmfile',"Select one PSM file (tab-deliminated): ",choices = c("", list.files('data/PSM/.')),selected=c("test_required.tsv")),
                   
                   
                   selectInput('msfile','Select one directory contained spectrum files(*.mgf):', choices = c("",list.files('data/MGF/.')),selected="test"),
                   hr(style = "border-top: 1px solid #000000;"),
                   
                   br(),
                   
                   #h5("Additional info"),
                   
                   # tags$style(type='text/css', "text-input { font-size: 8px; line-height: 8px;} "),
                   # textInput("dec_pre","Decoy prefix: ", "XXX_"),
                   # textInput("prot_id","Protein column index: ", "6"),
                   # textInput("frag_tol","Fragment tol.(Da): ", "0.025"),
                   # 
                   br(),
                   actionButton("loadButton", "First. Load Files", class="btn-success")
                 ),
                 mainPanel
                 (
                   fluidRow(
                     dashboardBody
                     (
                       tags$style(HTML("


                .box.box-solid.box-primary>.box-header {

                }

                .box.box-solid.box-primary{

                border-left-color:#f7f8f9;
                background:#f7f8f9
                }

                ")),
                       
                       column(width=8,
                              fluidRow(
                                tabBox
                                (   width=12,
                                  id="psms", height="550",
                                  tabPanel("PSMs",
                                           {
                                             box(width=12, height="100",
                                             )
                                             box(width=12, height="450",
                                                 
                                                 p("Tutorial:"),
                                                 tags$a(href="https://docs.google.com/document/d/168oGySS15xrobeqTY54_QLMgWQofnPIApx6ma5ntGa8/edit?usp=sharing", "How to run TIDD-GUI: data management"),
                                                 
                                                 textOutput('psm_filename'),
                                                 
                                                 h6("Features in the psm result file."),
                                                 br(),
                                                 tags$div(
                                                   align='left',
                                                   class='multicol',
                                                   style="font-size: 10px; padding:0px",
                                                   radioButtons("psmf","Column info (index.feature name):", choices=c("-"), inline=TRUE)
                                                   
                                                 )
                                                 
                                             )
                                             
                                             
                                             
                                           }),
                                  tabPanel("Preview (table)",
                                           
                                           dataTableOutput('psmlist_preview'),
                                           
                                  ),
                                  
                                ),),
                              fluidRow(
                                tabBox
                                (width=12,
                                  id="ms2", height="450px",
                                  
                                  tabPanel("Extraction log",
                                           box(width=12, height="50",
                                               textOutput('tidd0'),
                                           ),
                                           box(width=12,height="50",
                                               verbatimTextOutput('tidd1') ),
                                           
                                  ),
                                  
                                  tabPanel("Summary&Plot",
                                           
                                           fluidRow(width=12,
                                                    verbatimTextOutput("summary_feature"),
                                                    
                                           ),
                                           fluidRow(width=12,
                                                    plotOutput('feature_plot'))
                                  ),
                                  
                                ),),
                       ),
                     ),
                     
                     column(width = 4,
                            box( 
                              width=12,
                              theme="flatly",
                              h5("6 required features for PSM feature extraction:"),
                              hr(style = "border-top: 1px solid #000000;"),
                              h6("Write the corresponding column indices in the same order as below:"),
                              h6("(separator between indices: \",\")"),
                              h6("1. filename"),
                              h6("2. scan number"),
                              h6("3. charge"),
                              h6("4. precursor mass"),
                              h6("5. peptide "),
                              h6("6. protein (separator \";\")"),
                              
                              h6(""),
                              
                              textInput("required_feature","required column indices:","1,2,3,4,5,6"),
                              textInput("dec_pre","Decoy prefix: ", "XXX_"),
                              textInput("frag_tol","Fragment tol.(Da): ", "0.025"),
                              
                              
                              actionButton("featureButton", " Second.Extract!", class="btn-success"),
                              h6("Check the \"Extraction log\" tab whether it is finished or not")
                              
                              
                              #uiOutput('feature_extract'),
                              
                              
                            )
                            
                     )
                   )
                 )
               )
             ),
             tabPanel(
               "Model Fitting",
               
               sidebarLayout(
                 sidebarPanel(
                   width = 3,
                   
                   h5("PSMs for model fitting"),
                   hr(style = "border-top: 1px solid #000000;"),
                   h6("Tab-deliminated file:"),
                   h6("1. head: do not start with '#' symbol"),
                   h6("2. one column for one feature"),
                   
                   br(),
                   selectInput('psm_feature_file',"Select one PSM file : ",choices = c(list.files('data/PSM/.')),
                   ),
                   actionButton("load2Button", " LOAD <GO> !", class="btn-success")
                 ),
                 mainPanel
                 (
                   tags$head(tags$style(
                     
                     
                     HTML("
                           
                            .multicol .shiny-options-group{
                            -webkit-column-count: 3;* Chrome, Safari, Opera */
                            -moz-column-count: 3;    /* Firefox */
                            column-count: 3;
                            -moz-column-fill: balanced;
                            -column-fill: balanced;
                            }
                            .checkbox{margin: 0;}
                            .checkbox p{margin:0px;}
                            .shiny-input-cantainer{margin-bottom:0px;}
                            ")
                   )),
                   
                   fluidRow(
                     
                     dashboardBody
                     (
                       column(width=12,
                              fluidRow(
                                br(),
                                tabBox(width=12,height="450px",
                                       tabPanel("Feature selection",
                                                br(),
                                                p("Tutorial:"),
                                                tags$a(href="https://docs.google.com/document/d/168oGySS15xrobeqTY54_QLMgWQofnPIApx6ma5ntGa8/edit?usp=sharing", "How to run TIDD-GUI: model fitting"),
                                                br(),
                                                textOutput('feature_filename'),
                                                h5("Please select x, y and init score correctly."),
                                                # h5("y: selected factor type of feature, such as \"TorD\""),
                                                # h5("init score: select numeric feature, such as 
                                                #    \"CalXcorr\""),
                                                # 
                                                h5("x:"),
                                                
                                                column(width=9, height = "450px",
                                                       tags$div(
                                                         align='left',
                                                         class='multicol',
                                                         
                                                         checkboxGroupInput("check1","",choices = c("")),
                                                         style="font-size:11px ! important; text-align; left"
                                                       )),
                                                column(
                                                  width=3, height="450px",
                                                  textOutput("x"),
                                                  tags$style(type='text/css', ".selectize-input { font-size: 11px; line-height: 11px;} .selectize-dropdown { font-size: 11px; line-height: 11px; }"),
                                                  selectInput("y"," y (tar/dec): (Select TorD feature)", multiple = FALSE, choices = list("====")),
                                                  selectInput("model"," model:", multiple = FALSE, choices = c("iterative SVM")),
                                                  selectInput("init_score"," init score: (Select xcorr or calXcorr feature)", multiple = FALSE, choices = c("====")),
                                                  selectInput("k_cross"," k-CV:", multiple = FALSE, choices = c("3","5","10")),
                                                  tags$style(type='text/css', "text-input { font-size: 11px; line-height: 11px;} "),
                                                  textInput("tain_size","Max train size: ", "1000"),
                                                  textInput("iteration","Max iteration (SVM): ", "1"),
                                                  actionButton("modelbutton", " Fitting <GO> !", class="btn-success")
                                                  
                                                )
                                       ),),),
                              fluidRow(
                                tabBox
                                (width=12,
                                  id="model_fit", height="450px",
                                  tabPanel("Model fitting",
                                           
                                           fluidRow(
                                             column(width=6,
                                                    h5("log:"),
                                                    verbatimTextOutput("confirm"),
                                                    textOutput("static"),
                                                    textOutput("initID"),
                                                    verbatimTextOutput("log"),
                                                    
                                             ),
                                             column(width=6,
                                                    h5("plot:"),
                                                    plotOutput("IDplot")
                                             )
                                           )
                                  ),
                                  
                                ),),
                       ),
                     ),
                     
                   )
                 )
               )
             )
             
  )
)
