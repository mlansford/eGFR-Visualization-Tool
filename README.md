# eGFR-Visualization-Tool

This is an RShiny app originally developed by Scott Frodsham that gives users a few different options for plotting patient eGFR data. The data as a whole can be plotted on a histogram, eGFR measurements can be plotted per patient on a line plot, they can be plotted with a linear regression indicating the slope of rapid or slow decliners, and the app also displays a table of all or selected patient information. Some functions of the app remain limited or are missing completely and will be subject to future development.

### USAGE:
This app should and can only be run on the University of Utah CHPC redwood node. Electronic Health records required by the app to run should only be stored in the protected environment. Thus the current method for running this app is that it will only work within the protected environment by those that have access to the IRB directory where the EHR files are stored (user must be on the IRB). This ensures that the patient data is secure and is used within the parameters set by the IRB.

#### A Couple of Ways to Run the App within the PE:
```
# Install the following R packages:
# • shiny
# • shinyjs
# • DT
# • tidyverse
# • ggplot2

Then use the following command:
runGitHub( "<this repository name>", "<your user name>")
```

```
# First clone the repository with git. If you have cloned it into
# ~/shiny_example, first go to that directory, then use runApp():
setwd("~/shiny_example")
runApp()
```
