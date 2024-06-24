####
#### R script for Ohigashi et al (2024)
#### 
#### 2024.06.24 Ohigashi
#### R 4.3.3

dir.create("00_SessionInfo")
dir.create("Data")
dir.create("FigCode")
dir.create("Function")

writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/01_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
