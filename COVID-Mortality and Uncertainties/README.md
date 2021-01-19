# COVID-19 
### R13B-Group 4

## Introduction 
In a year we all thought would be the best yet, the Coronavirus pandemic as put the whole world on hold. With nearly a year since the first case, this virus, making its way all around the globe has affected the way we live, businesses and most of all taken loved ones from us. As this virus spreads in transmission through respiratory droplets and medical researchers aim to find a vaccine, data analysts are using tools to understand who this virus affects the most. We have heard similar news around the world but also different information that makes us ask questions like "why did it affect them and not me"? 

In this repository, the team aims to tackle the **driving question: How deadly is COVID-19 and how can we present data about this question so that the uncertainty is made very clear to the user of the visualisation of the results?** By answering this driving question, we aim to improve the misconception illustrated by poor analytics and miscommunication and therefore provide a concise conclusion to the question. To address this driving question we explored different datasets and conducted analysis on different countries to come up with methods to deal with multiple uncertainties that present in both calculations and presentations.


## Table of Contents

* [Product Notebooks](https://github.sydney.edu.au/swan9801/R13B-Group4-COVID/blob/master/Product_Notebook_Outline.ipynb):This is our final Product Notebook

* [Process Notebooks](https://github.sydney.edu.au/swan9801/R13B-Group4-COVID/tree/master/Process%20Notebooks): cleaned process notebooks (enhanced version based on checkpoint 2) into different phases of the full data analysis process, with preecisely one notebook for each part of the process including Data Engineering, Exploration Data Analysis (EDA), Case Fatality Rate (CFR), Infection Fatality Rate (IFR) and Uncertainties Visualisation.

It includes: 
  * [Code Review, Think alouds and Retrospective](https://github.sydney.edu.au/swan9801/R13B-Group4-COVID/tree/master/Process%20Notebooks/Code%20Review) 
  * [Checkpoint 1](https://github.sydney.edu.au/swan9801/R13B-Group4-COVID/tree/master/Process%20Notebooks/Checkpoint%201): Exploratory Data Analysis for COVID-19 in Australia

  * [Checkpoint 2](https://github.sydney.edu.au/swan9801/R13B-Group4-COVID/tree/master/Process%20Notebooks/Checkpoint%202): Exploratory Data Analysis for COVID-19 in India, United States, Russia, South Africa and Brazil
   * [Datasets](https://github.sydney.edu.au/swan9801/R13B-Group4-COVID/tree/master/Process%20Notebooks/Checkpoint%202/Model%20Datasets) used for Checkpoint 2


* [Group Report](https://github.sydney.edu.au/swan9801/R13B-Group4-COVID/blob/master/GROUP%20REPORT.pdf): This is the team's group work process report

* [Colour Simulator](https://github.sydney.edu.au/swan9801/R13B-Group4-COVID/tree/master/ColourSimulator): This folder contains the images we generated for our visualisations through the colour blind generator.

## Understanding the Notebooks

* **Key files**

The main COVID-19 data file we used is obtained from [Our World in Data](https://ourworldindata.org/coronavirus). The data is a collection of the COVID-19 data maintained and updated daily, which includes key data that we used in our analysis, including country codes, confirmed cases, deaths, tests and many other factors that are of our interest.

The second file we used for analysis is the [John Hopkins University](https://coronavirus.jhu.edu/map.html) COVID-19 Data Repository. This file was mainly used to source the missing number of recovered cases for each country that *Our World in Data* does not have.

The following [model datasets](https://github.sydney.edu.au/swan9801/R13B-Group4-COVID/tree/master/Process%20Notebooks/Checkpoint%202/Model%20Datasets) were obtained from [Our World in Data](https://ourworldindata.org/covid-models) in order to estimate the true number of infection cases for each country. This is the fundamental datasets we used for calculation of the Infection Fatality Rates (IFR). 5 files were imported from this site, with 4 different models and one combined file. 
Short description of the models:

- Imperial College London (ICL):

Age-structured SEIR model focuses on low-and middle-income countries
Model uses age and country-specific data on demographics patterns of social contact, hospital availability and risk of hospitalisation and death
Assumes sufficient access to healthcare
Assumes change in transmission over time is a function of average mobility trends

- The Institute for Health Metrics and Evaluation (IHME)

Model uses different data to simulate transmission and disease progression: mobility, social distancing policies, population density, pneumonia seasonality and death rate, air pollution, altitude, smoking rates and self-reported contacts and mask use.
Death model assumes relationship between confirmed deaths, confirmed cases and testing levels. 

- Youyang Gu (YYG)

Model created and optimised for the US
assumptions on how reopening will affect social distancing and ultimately transmission.


- The London School of Hydiene & Tropicala Medicine (LSHTM)

Assumes delay adjusted CFR of 1.4% baseline that does not account for different age distributions outside China.  Overestimation for younger populations and underestimation in countries with older populations. 


* **Key variables**
Some key variables used and explored in the notebooks are defined below:

Variables in the 'Our World in Data' data set:

- `new_cases_smoothed` New confirmed cases of COVID-19 (7-day smoothed)

- `new_deaths_smoothed` New deaths attributed to COVID-19 (7-day smoothed)

- `new_tests_smoothed` New tests for COVID-19 (7-day smoothed). For countries that don't report testing data on a daily basis, we assume that testing changed equally on a daily basis over any periods in which no data was reported. This produces a complete series of daily figures, which is then averaged over a rolling 7-day window	

- `total_deaths` New deaths attributed to COVID-19 (cumulative value)

- `total_cases` Total confirmed cases of COVID-19 (cumulative value)

- `location` Country name

- `date` Reported date of observation

- `aged_65_older` Share of the population that is 65 years and older, most recent year available

- `human_development_index` Summary measure of average achievement in key dimensions of human development: a long and healthy life, being knowledgeable and have a decent standard of living

Variables in the 'John Hopkins' data set:

- `total_recovered` New recovered attributed to COVID-19 (cumulative value)

## Getting Started

Below are the steps to run the project:
1. Download the repository or simply clone the repository to your GitHub Desktop containing all the files and folders
2. Find Product_Notebook_Outline.ipynb file (this is the Product Notebook)
3. Make sure you have the up to date libraries installed to run the file
- Use pandas `1.1.3`
- Use numpy `1.15.4`
- Use plotly.express `4.11.0`
- Use sklearn.linear_model `0.20.1`
4. Run the Product Notebook


## Authors

| Team Member | Role | Individual Country | Contact details |
| --- | --- | --- | --- |
| Cynthia Mather | Manager |  India  | cmat7799@uni.sydney.edu.au |
| Sherry Wang | Tracker |  United States  | swan9801@uni.sydney.edu.au |
| Hendrix Wang| Co-ordinator |  South Africa  | zwan6697@uni.sydney.edu.au |
| Johnson Yun| Editor |  Russia  | zyun3912@uni.sydney.edu.au |
| Charles Zhang | Formatter |  Brazil  | zzha3884@uni.sydney.edu.au |

## Times to meet outside class

Saturday 3:00 - 5:00 pm
