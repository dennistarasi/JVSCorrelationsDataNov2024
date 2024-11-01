#We performed post-sampling calculations to quantify several topographic and biotic variables of interest in predicting ESS. We calculated each plot’s straight-line distance to the nearest stream in QGIS using the EPA National Hydrography Dataset Plus stream data (regions 02 and 03N). We calculated each plot’s straight-line distance to the nearest road centerline in the study region using state DOT maps of primary and secondary roads in Virginia, North Carolina and South Carolina. These maps depicted all public roads with at least one lane of traffic in each direction, thereby serving as a common corridor for vehicle movement. 
#Note: Because many plots included federally listed species or were collected on private lands, calculations of proximity measures utilized GPS coordinates with more accuracy than what is publicly available and were completed prior to the start of this code.

library(plyr)
setwd("")
Abio<-read.csv("PiedmontPredictors.csv",head=T)
Abio$Slope<-ifelse((Abio$Slope==-2),NA,Abio$Slope)
Abio$Aspect<-ifelse((Abio$Aspect==-2),NA,Abio$Aspect)
Abio$Aspect<-ifelse((Abio$Aspect==-1),NA,Abio$Aspect)

#Some plots include soil data for multiple individual subplots; for these cases we used the average value of each soil characteristic for the whole plot. Negative aspect and slope values were assigned during sampling as unable to be determined, so these were converted to NA.

#Using CVS species records, we calculated the total native species richness of each plot. 
Plot<-read.csv("PlotSpecies.csv", head=T)
Origins<-read.csv("SpeciesOrigins.csv", head=T)
Species<-merge(Plot,Origins, all.x=T)
Native<-Species[Species[,"Origin"]=="Native",]
NRich<-ddply(Native, "Plot", summarize, NRich=length(Origin))
dat1<-merge(Abio,NRich, all.x=T, all.y=F)

#Where multiscale species records were available, we calculated the 400 m2 native species richness. A total of 118 plots >400 m2 did not include sub-plot species records. Using the Arrhenius log-log species area model (Fridley et al. 2005), we modeled the species accumulation from 400 to 1000 m2 in the 382 plots in the dataset with species data at each spatial grain. We then extrapolated that training data to the 118 plots lacking in 400 m2 records to scale down from the plot level to a standardized 400 m2 (Supplementary File 1).

Intensives<-read.csv("IntensiveSpeciesOnly.csv",head=T)
Sp2<-merge(Intensives,Origins, all.x=T)
Nat<-Sp2[Sp2[,"Origin"]=="Native",]
NatRich<-ddply(Nat, "Plot", summarize, NatRich=length(Origin))
dat1a<-merge(dat1,NatRich, all.x=T, all.y=F)

#Arrhenius log–log models(log S = Z log A + log c);See Supp File 1 for Calculations. 600 Scalar is 0.9292, 800 is 0.882, 1000 is 0.847 (See Supplemental File 1)
dat1b<-subset(dat1a,Area==400)
dat1b$NatRich<-dat1b$NRich
dat1c<-subset(dat1a,Area==600)
dat1c$ProxyRich<-dat1c$NRich*0.9292
dat1c$NatRich<-ifelse(is.na(dat1c$NatRich),dat1c$ProxyRich,dat1c$NatRich)
dat1d<-subset(dat1a,Area==800)
dat1d$ProxyRich<-dat1d$NRich*0.882
dat1d$NatRich<-ifelse(is.na(dat1d$NatRich),dat1d$ProxyRich,dat1d$NatRich)
dat1e<-subset(dat1a,Area==1000)
dat1e$ProxyRich<-dat1e$NRich*0.847
dat1e$NatRich<-ifelse(is.na(dat1e$NatRich),dat1e$ProxyRich,dat1e$NatRich)
dat1f<-subset(dat1a,Area==700)
dat1g<-subset(dat1a,Area==500)
dat1h<-rbind(dat1b [,c(1:27,30)],dat1c[,c(1:27,30)],dat1d[,c(1:27,30)],dat1e[,c(1:27,30)],dat1f[,c(1:27,30)],dat1g[,c(1:27,30)])

#Because soil moisture may be a suitable predictor of exotic species success but is not measured in CVS protocols, we used species records to calculate a proxy measure for plot moisture and wetness. For this purpose, we calculated the total visual cover of species in each plot that are listed as obligate wetland or facultative wetland species according to the 2020 National Wetland Plant List for the Eastern Mountains and Piedmont Region. Cover values were assigned as the geometric means of the percentage cover for that cover class (Table 2), as species occurrences tend to occur disproportionately in the lower portion of each cover class (Peet and Roberts 2013).
Wetland<-read.csv("EMPWetland.csv",head=T)
tempa<-merge(Species,Wetland)
Geo<-read.csv("GeoMean.csv",head=T)
temp<-merge(tempa,Geo)
Wet<-ddply(temp, "Plot", summarize, Wet=sum(GeoMean))
dat2<-merge(dat1h,Wet, all.x=T)
dat2$Wet[is.na(dat2$Wet)]<-0

#write.csv(dat2, "PiedmontRawData.csv")
#dat2<-read.csv("PiedmontRawData.csv",header=T)
#We tested each continuous predictor variable for normality of distribution, and transformed non-normal distributions as needed, such as aspect and soil nutrient content (see table 1 for transformations).  
dat2$cosAspect<-cos(dat2$Aspect)
dat2$logOrganic<-log(dat2$Organic)
dat2$logS<-log(dat2$S)
dat2$logP<-log(dat2$P)
dat2$logCa<-log(dat2$Ca)
dat2$logMg<-log(dat2$Mg)
dat2$logK<-log(dat2$K)
dat2$logNa<-log(dat2$Na)
dat2$logBasalArea<-log(dat2$BasalArea)
dat2$logStemDensity<-log(dat2$StemDensity)
dat2$logB<-log(dat2$B)
dat2$logFe<-log(dat2$Fe)
dat2$logMn<-log(dat2$Mn)
dat2$logCu<-log(dat2$Cu)
dat2$logZn<-log(dat2$Zn)
dat2$logAl<-log(dat2$Al)


#We then standardized each variable individually, such that the mean measured value for each variable was converted to zero. All values for that variable were likewise converted to their corresponding z score to demonstrate the deviation of variables with ESS, for more uniform comparison of impacts between variables.
dat3<-dat2[,c(1,4,6:8,10:15,28:45)]

dat3[,30:57]<-sapply(dat3[,2:29], function(dat3)(dat3-mean(dat3,na.rm=TRUE))/sd(dat3,na.rm=TRUE))

colnames(dat3)[30:57]<-c("zElevation","zSlope","zWater","zRoad","zSand","zSilt","zClay","zPH","zExchangeCapacity","zBaseSaturation","zNatRich","zWet","zcosAspect","zlogOrganic","zlogS","zlogP","zlogCa","zlogMg","zlogK","zlogNa","zlogBasal_Area","zlogStemDensity", "zlogB", "zlogFe","zlogMn","zlogCu","zlogZn","zlogAl")

dat4<-dat3[,c(1,30:57)]

#We calculated several continuous and binary metrics of exotic species success for each plot.  First, we summed the total cover of species considered “exotic” to the southeastern United States, according to Weakley (2022). Cover values were assigned as the geometric means of the percentage cover for that cover class (Table 2), as species occurrences tend to occur disproportionately in the lower portion of each cover class (Peet and Roberts 2013). All exotic species cover means were summed for each plot. Because different species may occupy different zones of the canopy and understory, some geometric means exceed 100% cover. Plots were also binned according to their cover for binomial regression.
	
Geo<-read.csv("GeoMean.csv",head=T)

temp<-merge(Species,Origins, all.x=T)
temp2<-merge(temp, Geo)
temp3<-temp2[temp2$Origin=="Alien",]
Exotics<-ddply(temp3, "Plot", summarize, Exotics=sum(GeoMean))

dat5<-merge(dat4, Exotics, all.x=T)
dat5$Exotics[is.na(dat5$Exotics)]<-0
dat5$Ex2[dat5$Exotics>2]<-1
dat5$Ex2[is.na(dat5$Ex2)]<-0
dat5$Ex5[dat5$Exotics>5]<-1
dat5$Ex5[is.na(dat5$Ex5)]<-0
dat5$Ex10[dat5$Exotics>10]<-1
dat5$Ex10[is.na(dat5$Ex10)]<-0
	
#write.csv(dat5, "PiedmontAnalysis.csv")
	
#To determine which abiotic and biotic variables correlate with ESS, we performed logistic regression of each continuous variable with the percent exotic cover, using the MASS package in R version 2.15.1 (R Core Team 2012). Because 601 plots in the dataset of 1425 plots did not have any exotic species cover, we utilized a generalized linear model with a negative binomial regression.

dat6<-dat5[,2:30]

library(MASS)
dat7<-matrix(nrow=28,ncol=12)
for(i in 1:28){
	z<-glm.nb(Exotics~dat6[,i], dat6)
	dat7[i,1]<-coefficients(z)[1]
	dat7[i,2]<-coefficients(z)[2]
	dat7[i,3]<-summary(z)$coefficients[8]
	}


#We also performed binomial regressions of each continuous variable with each exotic cover threshold  for a total of 28 logistic and 84 binary models. Plots meeting the threshold for exotic cover were designated with a “1” at the given threshold, while all other plots were designated with a “0”. We calculated the slope of the binomial regression line and the Wald test statistic to determine significance of each model (Legendre and Legendre 2012). 

datEx2<-dat5[,c(2:29,31)]

for(i in 1:28){
	z<-glm(Ex2~datEx2[,i], datEx2, family=binomial)
	dat7[i,4]<-coefficients(z)[1]
	dat7[i,5]<-coefficients(z)[2]
	dat7[i,6]<-summary(z)$coefficients[8]
	}

datEx5<-dat5[,c(2:29,32)]

for(i in 1:28){
	z<-glm(Ex5~datEx5[,i], datEx5, family=binomial)
	dat7[i,7]<-coefficients(z)[1]
	dat7[i,8]<-coefficients(z)[2]
	dat7[i,9]<-summary(z)$coefficients[8]
	}
	
	datEx10<-dat5[,c(2:29,33)]

for(i in 1:28){
	z<-glm(Ex10~datEx10[,i], datEx10, family=binomial)
	dat7[i,10]<-coefficients(z)[1]
	dat7[i,11]<-coefficients(z)[2]
	dat7[i,12]<-summary(z)$coefficients[8]
	}

dat8<-as.data.frame(dat7)
colnames(dat8)<-c("NBInt","NBSlope","NBSig","Bin2Int","Bin2Slope","Bin2Sig","Bin5Int","Bin5Slope","Bin5Sig","Bin10Int","Bin10Slope","Bin10Sig")
rownames(dat8)<-c("Elevation","Slope","Water","Road","Sand","Silt","Clay","PH","ExchangeCapacity","BaseSaturation","NativeRichness","Wetland","Aspect","OrganicMatter","Sulfur","Phosphorus","Calcium","Magnesium","Potassium","Sodium","BasalArea","StemDensity","Boron", "Iron","Manganese","Copper","Zinc","Aluminum")

#write.csv(dat8, "Regressions.csv")

#We then developed a classification and regression tree. Plots with at least 5% cover of exotic species were designated with a “1” and all others were designated with a “0”.

library(rpart)
library(rpart.plot)
temp<-dat5[,c(1,32)]
dat9<-merge(dat2,temp)
ExTree<-rpart(Ex5~Elevation+cosAspect+Slope+Water+Road+Organic+Sand+Silt+Clay+pH+ExchangeCapacity+BaseSaturation+NatRich+Wet+BasalArea+StemDensity+S+P+Ca+Mg+K+Na+B+Fe+Mn+Cu+Zn+Al, minsplit=30, dat9)
rpart.plot(ExTree, extra=1, uniform=T)

#Because soil texture data was most sparse in the dataset, we also developed a classification tree on a partial dataset, constrained to partition only the 944 plots that contained all 25 non-soil texture predictor variables, 79 of which expressed 5% exotic cover. However, this tree did not meaningfully differ from the tree using all plots (Supplementary File 2), so we maintained robust geographic representation with all plots included in the tree.
#Tree using only plots with complete data (soil texture omitted), representing 1071 total plots, 94 of which are invaded. No dramatic differences from tree using all plots.
#dat10<-dat9[,c(1:4,6:9,13:30,46)]
#dat11<-na.omit(dat10)
#write.csv(dat11,"CompletePlots.csv")
#NewTree<-rpart(Ex5~Elevation+cosAspect+Slope+Water+Road+Organic+pH+ExchangeCapacity+BaseSaturation+NatRich+Wet+BasalArea+StemDensity+S+P+Ca+Mg+K+Na+B+Fe+Mn+Cu+Zn+Al, minsplit=30, dat11)
#rpart.plot(NewTree, extra=1, uniform=T)

