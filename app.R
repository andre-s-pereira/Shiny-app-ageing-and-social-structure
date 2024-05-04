library(shiny)
library(igraph)
library(visreg)
library(plyr)
library(tidyverse)
library(magrittr)
library(lubridate)
library(readr)
library(cowplot)
library(ggpubr)
library(ggeffects)
library(ggplot2)
library(tidygraph)
library(ggraph)

ui <- fluidPage(
  
  titlePanel("Siracusa et al. (2023). Ageing in a collective: The impact of ageing individuals on social network structure."),
  sidebarLayout(
    sidebarPanel(
      numericInput(
        inputId = "n_simulations",
        label = "Number of simulations (likely to crash above 1000)",
        min = 1, max = 1000,
        value=100
      ),
        numericInput(
          inputId = "n_kin",
          label = "Number of matrilines/kin groups",
          min = 1, max = 50,
          value= 5
        ),
      sliderInput(
        inputId = "percentage_of_old_individuals",
        label = "Percentage of old individuals",
        min = 0, max = 100,
        value=25
      ),
      selectInput(
        inputId= "net_measure", 
        label="Plot", 
        choices=c("Mean degree", "Transitivity", "Diameter", "Social network (only plots the last simulated network, similar to setting number of simulations=1)")
      ),
      strong("ShinyApp description"),
      p("This ShinyApp allows readers to use the agent-based model presented in Siracusa et al., 2023, “Ageing in a collective: The impact of ageing individuals on social network structure”. This agent-based model simulates how changing proportions of old individuals in a network of female rhesus macaques affects network structure. Users can control the number of simulations and the percentage of old individuals in the simulated rhesus macaque group."),
      p("All simulated groups features 50 adult females, which mirrors the mean number of adult females in real groups on Cayo Santiago. Each simulated group has 10 kin clusters, with 5 individuals in each cluster. Individuals in the same kin cluster are close kin and individuals in different kin clusters are non-kin. The kinship structure of the simulated networks mirrors that of real Cayo Santiago groups."),
      p("The linking probability of the simulated individuals (i.e., the probability that two simulated individuals have an edge, representing a grooming link) depends on the age and kinship of the individuals and is based on the empirical Cayo Santiago data. The linking probability for dyads that consist of old individuals who are related to each other (i.e., old/old kin) is 0.33; for old/old non-kin dyads it is 0.02; for old/young kin dyads it is 0.37; for old/young non-kin dyads it is 0.05; for young/young kin dyads it is 0.27 and for young/young non-kin dyads it is 0.08."),
    ),
    mainPanel(
      plotOutput("myPlot")
    )
  )
)

server <- function(input, output, session) {
  
  output$myPlot <- renderPlot({
    
    nsims <- input$n_simulations
    nkin <- input$n_kin
    pold <- input$percentage_of_old_individuals
    netm <- input$net_measure
    
    #################### 
    #### Producing the graphs
    ####################
    
    source("age_model_app.R")
    age_model(pold, nkin, nsims)
    load(file="age_simulation_dataset.Rdata")
    
    #############
    
    if(netm == "Mean degree") {
      if(length(unique(outputdata$mean.degree))==1) {
        md<-round(outputdata$mean.degree[1],digits=2)
        text = paste("Mean degree was =", md, "in all simulations")
        ggplot() + 
          annotate("text", x = 4, y = 25, size=8, label = text) + 
          theme_void()
      } else{ 
          ggplot(outputdata, aes(x = mean.degree)) + 
            geom_density(alpha = 0.5, show.legend = FALSE, size=1.1, fill="#009E73") +
            ylim (0,2) +
            xlim(1,7) +
            labs(x="Network mean degree", y="Density") +
            theme_classic(base_size = 33, base_family = "", base_line_size=1.1, base_rect_size=1.1) 
          
      }  
     
    } else if(netm == "Transitivity") {
      if(length(unique(outputdata$transitivity))==1) {
        tr<-round(outputdata$transitivity[1],digits=2)
        text = paste("Transitivity was =", tr, "in all simulations")
        ggplot() + 
          annotate("text", x = 4, y = 25, size=8, label = text) + 
          theme_void()
      } else{ 
      ggplot(outputdata, aes(x = transitivity)) + 
        geom_density(alpha = 0.5, show.legend = FALSE, size=1.1, fill="#009E73") +
        ylim(0,25) +
        xlim(0,0.3) +
        labs(x="Network transitivity", y="Density") +
        theme_classic(base_size = 33, base_family = "", base_line_size=1.1, base_rect_size=1.1)
          
      }
      
    } else if(netm == "Diameter") { 
      if(length(unique(outputdata$diameter))==1) {
        dm<-round(outputdata$diameter[1],digits=2)
        text = paste("Diameter was =", dm, "in all simulations")
        ggplot() + 
          annotate("text", x = 4, y = 25, size=8, label = text) + 
          theme_void()
      } else{ 
      ggplot(outputdata, aes(x = diameter)) + 
        geom_density(alpha = 0.5, show.legend = FALSE, size=1.1, fill="#009E73") +
        ylim(0,2) +
        xlim(3,20) +
        labs(x="Network diameter", y="Density") +
        theme_classic(base_size = 33, base_family = "", base_line_size=1.1, base_rect_size=1.1)
      }
      
    } else if(netm == "Social network (only plots the last simulated network, similar to setting number of simulations=1)") {
      load(file="age_network_plot.Rdata")
      load(file="nodedata.Rdata")
      if(0 %in% nodedata$agegroup & 1 %in% nodedata$agegroup) {
        inet$agegroup <- as.factor(inet$agegroup)
        ggraph(inet, layout = "fr") +
        geom_edge_link0(edge_colour = "dark gray", show.legend = F) +
        scale_edge_width(range = c(0.5, 2)) + #control size of edge
        scale_color_viridis(discrete=TRUE, option = "magma", begin = 0, end = 0.9) +
        geom_node_point(aes(colour = kingroup, shape = agegroup, size = 2), show.legend = T) +
        theme_graph() +
        coord_fixed() +
        scale_shape_manual(values = c(17, 19),
                           labels = c("Young", "Old")) +
        scale_size(guide = "none") +
        guides(color= "none") +
        labs(shape = "Age", caption = "Colours = kingroups                                               ") +
        theme_graph(base_family="sans") +
        theme(legend.position = c(0.9, 0.05),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 9),
              legend.box.background = element_rect(color="black", size=0.5))
        } else if(0 %in% nodedata$agegroup){
          inet$agegroup <- as.factor(inet$agegroup)
          ggraph(inet, layout = "fr") +
            geom_edge_link0(edge_colour = "dark gray", show.legend = F) +
            scale_edge_width(range = c(0.5, 2)) + #control size of edge
            scale_color_viridis(discrete=TRUE, option = "magma", begin = 0, end = 0.9) +
            geom_node_point(aes(colour = kingroup, shape = agegroup, size = 2), show.legend = T) +
            theme_graph() +
            coord_fixed() +
            scale_shape_manual(values = c(17),labels = c("Young")) +
            scale_size(guide = "none") +
            guides(color= "none") +
            labs(shape = "Age", caption = "Colours = kingroups                                               ") +
            theme_graph(base_family="sans") +
            theme(legend.position = c(0.9, 0.05),
                  legend.title = element_text(size = 10),
                  legend.text = element_text(size = 9),
                  legend.box.background = element_rect(color="black", size=0.5))
        } else if(1 %in% nodedata$agegroup){
          inet$agegroup <- as.factor(inet$agegroup)
          ggraph(inet, layout = "fr") +
            geom_edge_link0(edge_colour = "dark gray", show.legend = F) +
            scale_edge_width(range = c(0.5, 2)) + #control size of edge
            scale_color_viridis(discrete=TRUE, option = "magma", begin = 0, end = 0.9) +
            geom_node_point(aes(colour = kingroup, shape = agegroup, size = 2), show.legend = T) +
            theme_graph() +
            coord_fixed() +
            scale_shape_manual(values = c(19),labels = c("Old")) +
            scale_size(guide = "none") +
            guides(color= "none") +
            labs(shape = "Age", caption = "Colours = kingroups                                               ") +
            theme_graph(base_family="sans") +
            theme(legend.position = c(0.9, 0.05),
                  legend.title = element_text(size = 10),
                  legend.text = element_text(size = 9),
                  legend.box.background = element_rect(color="black", size=0.5))
        }
    }
  })
}

shinyApp(ui = ui, server = server)