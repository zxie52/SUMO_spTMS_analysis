library(tidyverse, quietly = TRUE);
df <- data.frame(
                 condtion = c("Probe1 Accuracy", 
                              "Probe1 Accuracy", 
                              "Probe2 Stay Accuracy", 
                              "Probe2 Stay Accuracy", 
                              "Probe2 Switch Accuracy",
                              "Probe2 Switch Accuracy"),
                 accuracy = c(68.42, 70.97, 
                              70.31, 71.11, 
                              59.68, 58.22),
                 error = c(2.54, 2.69,
                           3, 3.09,
                           2.8, 2.61),
                 items= c("left items", "right items",
                          "left items", "right items",
                          "left items", "right items"))

library(scales)
#df %>%
  #group_by(items) %>%
ggplot(data = df, aes(x=condtion, y=accuracy, fill=items)) +
geom_bar(stat="identity", position = position_dodge()) +
coord_cartesian(ylim = c(50, 80))+
geom_text(aes(label=accuracy), vjust=-4.5, color="black",
          position = position_dodge(0.9), size=3.5) +
scale_fill_manual(values = c("Red", "Blue")) +
geom_errorbar(aes(ymin = accuracy - error, ymax = accuracy + error), 
              width = .2, colour = "purple",
              position=position_dodge(.9))+
theme_bw()+
labs(title = "Memory Performance")+
ggeasy::easy_center_title()
