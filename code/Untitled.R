rm(list=ls())
options(stringsAsFactors = F)
gc()
library(ggthemes)
library(dplyr)
library(readxl)
library(writexl)
library(tableone)
library(tidyr)
library(lubridate)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(survminer)
library(survival)
GSE134520=read_excel("USE-GSE134520.xlsx",sheet = 1)
GSE150290=read_excel("USE-GSE150290.xlsx" )
GSE183904=read_excel("USE-GSE183904.xlsx",sheet = 1)
HRA000704=read_excel("USE-HRA000704.xlsx",sheet = 1)
colnames(GSE134520)=stringr::str_split(colnames(GSE134520),"_",simplify = T)[,1]
colnames(GSE150290)=stringr::str_split(colnames(GSE150290),"_",simplify = T)[,1]
colnames(GSE183904)=stringr::str_split(colnames(GSE183904),"_",simplify = T)[,1]
colnames(HRA000704)=stringr::str_split(colnames(HRA000704),"_",simplify = T)[,1]






