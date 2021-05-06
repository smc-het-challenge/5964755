FROM ubuntu

RUN apt-get update 
RUN apt-get -y install build-essential 
RUN apt-get -y install g++ 

RUN apt-get -y install r-base

#setup R configs
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('combinat')"
RUN Rscript -e "install.packages('DPpackage')"

RUN mkdir -p /opt/BayClone3

ADD parseInputData.cpp /opt/BayClone3

ADD BayClone3-CC-v2_dream.cpp /opt/BayClone3
ADD BayClone3-CC-v2_dream.R /opt/BayClone3
ADD fn_BayClone3-CC-v2_dream.R /opt/BayClone3
ADD RobustFunctionsV5_ForICGC_v4_dream.R /opt/BayClone3

WORKDIR /opt/BayClone3

RUN g++ -o parseInputData parseInputData.cpp 

RUN R CMD SHLIB BayClone3-CC-v2_dream.cpp

