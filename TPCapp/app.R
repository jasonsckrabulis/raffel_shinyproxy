# Libraries
library(shiny)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)
library(rmarkdown)
library(knitr)

# Define UI
ui <- fluidPage(

    titlePanel(
        h1("Fitting Thermal Performance Curves")
    ),
    h6("Developed by Jason Sckrabulis"),
    h3("This Shiny-R app was developed to enable users to fit generic thermal performance curves for any response variable."),
    hr(),
    
    # Tabs
    tabsetPanel(id="tabset",
        
        # Analysis
        tabPanel("Analysis",
            h4("Upload a .csv file and receive a thermal performance curve output by following the directions in the sidebar below. The analysis defaults to a log-normal error distribution without mass-scaling."),
            hr(),
            
            sidebarLayout(
                
                # Sidebar options panel
                sidebarPanel(
                    
                    helpText("(1) -- Upload your .csv file for your data."),
                    helpText("Your .csv needs to be formatted with the following column headers as written: 'temperature', 'mass', 'response' (all lowercase and no single quotations). If your data does not contain organism mass values (i.e., chemical reaction rate data), please insert mass=1 for all datapoints (the mass-correction option will not change your model output if this is the case)."),
                    fileInput(inputId="file",label="Choose .csv file:",multiple=FALSE,
                              accept=c("text/csv","text/comma-separated-values,text/plain",".csv")),
                    hr(),
                    
                    helpText("(2) -- Options: Check the respective box to enable specific options. Otherwise, skip this step."),
                    checkboxInput(inputId="massCheck",label="Correct for organism mass"),
                    checkboxInput(inputId="normalCheck",label="Normal error distribution"),
                    hr(),
                    
                    helpText("(3) -- Click button below to start the analysis"),
                    actionButton(inputId="confirm",label="Start Analysis"),
                    hr(),
                    
                    downloadButton(outputId="downloadReport",label="Download HTML Report")
                    
                ),
                
                # Main output panel
                mainPanel(
                    
                    plotOutput("arrGraphC"),
                    plotOutput("arrGraphInv"),
                    tableOutput("arrTable")
                    
                )
                
            )
            
        ),
        
        # Notes and Help
        tabPanel("Notes and Help",
            h4("This tab contains notes about using the current version of this Shiny app and commonly asked questions about results."),
            hr(),
            
            mainPanel(
                
                p("The default model assumes a log-normal error distribution and does not correct for mass."),
                p("This app currently does not support datasets that have any values of zero for the response. Data are log-transformed as is, and do not have a constant added prior to log-transformation (i.e., is log(response), NOT log(response+1)."),
                p("In this app, mass is incorporated into the model as mass^(-1/4). Parameters are estimated as normal, but the response is plotted as the 'mass-corrected response' by dividing by mass^(-1/4). This is done prior to log-transformation.")
            
            ),
            
        ),
        
        # Development notes
        tabPanel("Development",
            hr(),
            
            mainPanel(
                
                # Known Issues
                h3("Known Issues:"),
                h4("None"),
                hr(),
                
                # Planned Updates
                h3("Planned Updates (Listed by priority):"),
                h4("Help text description tab dynamic UI"),
                h4("Arrhenius-type thermal performance curve - Boltzmann-Arrhenius model (direct comparison with future Sharpe-Schoolfield model)"),
                h4("Dynamically exclude mass UI elements if mass is not present in uploaded data"),
                h4("Unimodal thermal performance curve option - Sharpe-Schoolfield model"),
                hr(),
                
                # Version History
                h3("Version History:"),
                h4("v0.0.4 (Jul. 30, 2020): Basic HTML report generator and AICc calculations"),
                h4("v0.0.3 (Jul. 22, 2020): Mass-corrected and normal error distribution for Arrhenius"),
                h4("v0.0.2 (Jul. 14, 2020): Default Arrhenius model"),
                h4("v0.0.1 (Jul. 10, 2020): Basic UI")
                
            )
            
        )
        
    )
    
)

# Define server logic
server <- function(input, output) {
    
    # Reactives
    v<-reactiveValues(doPlot=FALSE)
    
    # Observe 'confirm' button to display plots
    observeEvent(input$confirm,{
        
        v$doPlot<-TRUE
    
    })
    
    # Events
    tpc<-eventReactive(input$confirm,{
        
        # Load data
        inFile<-input$file
        if(is.null(inFile)){
            
            return(NULL)
        
        }
        
        # Read into R and log-transform
        dat<-read.csv(inFile$datapath,header=TRUE)
        dat<-subset(dat,response!='NA')
        dat<-subset(dat,response>=0)
        dat$lnresponse<-log(dat$response)
        dat$lnmassresponse<-log(dat$response/(dat$mass^(-.25)))
        
        # Constants
        K<-273.15
        k<-8.62*10^-5
        
        # Starting estimates for Arrhenius
        arrStart<-lm(lnresponse~I(1/(k*(temperature+K))),data=dat)
        starting<-c(pto=exp(summary(arrStart)$coef[1]),Ea=abs(summary(arrStart)$coef[2]))
        
        # Equations
        # Arrhenius
        arrEq<-function(x,m,pto,Ea){
            
            if(input$massCheck==TRUE){
                
                mass<-m
                
            } else if(input$massCheck==FALSE){
                
                mass<-1
                
            }
            
            if(input$normalCheck==TRUE){
                    
                y<-(mass^(-.25))*pto*exp(-Ea/k*(1/(x+K)))
                
            } else if(input$normalCheck==FALSE){
                    
                y<-log((mass^(-.25))*pto*exp(-Ea/k*(1/(x+K))))
                
            }

        }
        
        # Arrhenius predictions from model
        arrPreds<-function(x,m,model){
            
            pto<-summary(model)$coef[1]
            Ea<-summary(model)$coef[2]
            
            if(input$massCheck==TRUE){
                
                mass<-m
            
            } else if(input$massCheck==FALSE){
            
                mass<-1
                    
            }
            
            if(input$normalCheck==TRUE){
                
                y<-mass^(-.25)*pto*exp(-Ea/k*(1/(x+K)))
                
            } else if(input$normalCheck==FALSE){
                
                y<-log(mass^(-.25)*pto*exp(-Ea/k*(1/(x+K))))
            }
            
        }
        
        # AICc calculation based on Xiao et al. 2011
        calcAICc<-function(numParams,loglik,n){
            
            numParams<-numParams+1
            2*numParams-2*loglik+2*numParams*(numParams+1)/(n-numParams-1)
            
        }
        
        if(input$normalCheck==TRUE){
            
            if(input$massCheck==TRUE){
                
                arr<-nls(response~arrEq(x=temperature,m=mass,pto,Ea),start=starting,data=dat,control=c(warnOnly=TRUE))
                
                arrSD<-sd(dat$response-arrPreds(x=dat$temperature,m=dat$mass,model=arr))
                arrLL<-sum(log(dnorm(dat$response,arrPreds(x=dat$temperature,m=dat$mass,model=arr),arrSD)))
                arrAICc<-calcAICc(numParams=2,loglik=arrLL,n=length(dat$temperature))
                    
                # Plot temp as C
                arrPlotC<-ggplot(data=dat,aes(x=temperature,y=(response/mass^(-.25))))+
                    geom_point(size=2)+
                    stat_function(fun=function(x)coef(arr)[1]*exp(-coef(arr)[2]/k*(1/(x+K))),geom="line",color="black",size=1,linetype="solid")+
                    labs(x="Temperature (C)",y="Mass-corrected Response")+
                    theme_cowplot(12)
                
                # Plot temp as 1/kT
                arrPlotInv<-ggplot(data=dat,aes(x=1/(k*(temperature+K)),y=lnmassresponse))+
                    geom_point(size=2)+
                    stat_function(fun=function(x)log(coef(arr)[1]*exp(-coef(arr)[2]*x)),geom="line",color="black",size=1,linetype="solid")+
                    labs(x="Temperature (1/kT)",y="Ln Mass-corrected Response")+
                    theme_cowplot(12)
                
                # Parameter estimates
                arrTable<-data.frame("Parameter"=c("A","Ea","AICc"),"Estimate"=c(summary(arr)$coef[1,1],summary(arr)$coef[2,1],arrAICc),"Standard Error"=c(summary(arr)$coef[1,2],summary(arr)$coef[2,2],"--"))
                
            } else if(input$massCheck==FALSE){
                
                arr<-nls(response~arrEq(x=temperature,m=mass,pto,Ea),start=starting,data=dat,control=c(warnOnly=TRUE))
                
                arrSD<-sd(dat$response-arrPreds(x=dat$temperature,m=dat$mass,model=arr))
                arrLL<-sum(log(dnorm(dat$response,arrPreds(x=dat$temperature,m=dat$mass,model=arr),arrSD)))
                arrAICc<-calcAICc(numParams=2,loglik=arrLL,n=length(dat$response))
                
                # Plot temp as C
                arrPlotC<-ggplot(data=dat,aes(x=temperature,y=response))+
                    geom_point(size=2)+
                    stat_function(fun=function(x)coef(arr)[1]*exp(-coef(arr)[2]/k*(1/(x+K))),geom="line",color="black",size=1,linetype="solid")+
                    labs(x="Temperature (C)",y="Response")+
                    theme_cowplot(12)
                
                # Plot temp as 1/kT
                arrPlotInv<-ggplot(data=dat,aes(x=1/(k*(temperature+K)),y=lnresponse))+
                    geom_point(size=2)+
                    stat_function(fun=function(x)log(coef(arr)[1]*exp(-coef(arr)[2]*x)),geom="line",color="black",size=1,linetype="solid")+
                    labs(x="Temperature (1/kT)",y="Ln Response")+
                    theme_cowplot(12)
                
                # Parameter estimates
                arrTable<-data.frame("Parameter"=c("A","Ea","AICc"),"Estimate"=c(summary(arr)$coef[1,1],summary(arr)$coef[2,1],arrAICc),"Standard Error"=c(summary(arr)$coef[1,2],summary(arr)$coef[2,2],"--"))
                
            }
            
        } else if(input$normalCheck==FALSE){
            
            if(input$massCheck==TRUE){
                
                arr<-nls(lnresponse~arrEq(x=temperature,m=mass,pto,Ea),start=starting,data=dat,control=c(warnOnly=TRUE))
                arrSD<-sd(dat$lnmassresponse-arrPreds(x=dat$temperature,m=dat$mass,model=arr))
                arrLL<-sum(log(dlnorm(dat$response,arrPreds(x=dat$temperature,m=dat$mass,model=arr),arrSD)))
                
                arrAICc<-calcAICc(numParams=2,loglik=arrLL,n=length(dat$temperature))
                
                # Plot temp as C
                arrPlotC<-ggplot(data=dat,aes(x=temperature,y=lnmassresponse))+
                    geom_point(size=2)+
                    stat_function(fun=function(x)log(coef(arr)[1]*exp(-coef(arr)[2]/k*(1/(x+K)))),geom="line",color="black",size=1,linetype="solid")+
                    labs(x="Temperature (C)",y="Ln Mass-corrected Response")+
                    theme_cowplot(12)
                
                # Plot temp as 1/kT
                arrPlotInv<-ggplot(data=dat,aes(x=1/(k*(temperature+K)),y=lnmassresponse))+
                    geom_point(size=2)+
                    stat_function(fun=function(x)log(coef(arr)[1]*exp(-coef(arr)[2]*x)),geom="line",color="black",size=1,linetype="solid")+
                    labs(x="Temperature (1/kT)",y="Ln Mass-corrected Response")+
                    theme_cowplot(12)
                
                # Parameter estimates
                arrTable<-data.frame("Parameter"=c("A","Ea","AICc"),"Estimate"=c(summary(arr)$coef[1,1],summary(arr)$coef[2,1],arrAICc),"Standard Error"=c(summary(arr)$coef[1,2],summary(arr)$coef[2,2],"--"))
                
            } else if(input$massCheck==FALSE){

                arr<-nls(lnresponse~arrEq(x=temperature,m=mass,pto,Ea),start=starting,data=dat,control=c(warnOnly=TRUE))
                arrSD<-sd(dat$lnresponse-arrPreds(x=dat$temperature,m=dat$mass,model=arr))
                arrLL<-sum(log(dlnorm(dat$response,arrPreds(x=dat$temperature,m=dat$mass,model=arr),arrSD)))
                
                arrAICc<-calcAICc(numParams=2,loglik=arrLL,n=length(dat$temperature))
                
                # Plot temp as C
                arrPlotC<-ggplot(data=dat,aes(x=temperature,y=lnresponse))+
                    geom_point(size=2)+
                    stat_function(fun=function(x)log(coef(arr)[1]*exp(-coef(arr)[2]/k*(1/(x+K)))),geom="line",color="black",size=1,linetype="solid")+
                    labs(x="Temperature (C)",y="Ln Response")+
                    theme_cowplot(12)
                
                # Plot temp as 1/kT
                arrPlotInv<-ggplot(data=dat,aes(x=1/(k*(temperature+K)),y=lnresponse))+
                    geom_point(size=2)+
                    stat_function(fun=function(x)log(coef(arr)[1]*exp(-coef(arr)[2]*x)),geom="line",color="black",size=1,linetype="solid")+
                    labs(x="Temperature (1/kT)",y="Ln Response")+
                    theme_cowplot(12)
                
                # Parameter estimates
                arrTable<-data.frame("Parameter"=c("A","Ea","AICc"),"Estimate"=c(summary(arr)$coef[1,1],summary(arr)$coef[2,1],arrAICc),"Standard Error"=c(summary(arr)$coef[1,2],summary(arr)$coef[2,2],"--"))
                
            }
            
        }
        
        # List of outputs
        list(arrPlotC=arrPlotC,arrPlotInv=arrPlotInv,arrTable=arrTable)
    })
    
    # Outputs
    # Arrhenius models
    output$arrGraphC<-renderPlot({
        if(v$doPlot==FALSE){
            return()
        } else{
            tpc()$arrPlotC
        }
    })
    
    output$arrGraphInv<-renderPlot({
        if(v$doPlot==FALSE){
            return()
        } else{
            tpc()$arrPlotInv
        }
    })
    
    output$arrTable<-renderTable({
        if(v$doPlot==FALSE){
            return()
        } else{
            tpc()$arrTable
        }
    })
    
    # Download PDF
    output$downloadReport<-downloadHandler(
        filename="TPC Ouput.html",
        content=function(file){
            tempReport<-file.path(tempdir(),"report.Rmd")
            file.copy("report.Rmd",tempReport,overwrite=TRUE)
            params<-list(arrPlotC=tpc()$arrPlotC,
                arrPlotInv=tpc()$arrPlotInv,
                arrTable=tpc()$arrTable)
            rmarkdown::render(tempReport,output_file=file,
                params=params,
                envir=new.env(parent=globalenv())
            )
        }
    )
}

# Run the application 
shinyApp(ui=ui,server=server)