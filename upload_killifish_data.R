
upload_killifish_data <- function(population){
  
  # population options: "naranjo", "guanapo", "both"

data <- read.csv("all_meso_data.csv", header=T,as.is = "embryo.dry.mass")

#Rearrange and make some calculations
#Select Experiment

if(population=="both"){
data <- data[which(data$spp=='killifish' &
                   #data$drainage=="Naranjo" &
                   (data$experiment=='KilliPhenoxDens' | data$experiment== "KilliPhenoxSize")
                   ),]
}

if(population=="naranjo"){
  data <- data[which(data$spp=='killifish' &
                       data$drainage=="Naranjo" &
                       (data$experiment=='KilliPhenoxDens' | data$experiment== "KilliPhenoxSize")
  ),]
}

if(population=="guanapo"){
  data <- data[which(data$spp=='killifish' &
                       data$drainage=="Guanapo" &
                       (data$experiment=='KilliPhenoxDens' | data$experiment== "KilliPhenoxSize")
  ),]
}

if(population!="naranjo"){
# exclude Guanapo data for KilliPhenoxSize for channels 15 & 16 (a tree fell on them)
data <- data[-which(data$experiment=="KilliPhenoxSize" & (data$channel==15 | data$channel==16) & data$drainage=="Guanapo"),]
}

data$e1 <- 0
data$e2 <- 0

#Number the experiment
data$exper <- 0
data[which(data$experiment=='KilliPhenoxDens'),'exper'] <- 1
data[which(data$experiment=='KilliPhenoxSize'),'exper'] <- 2


#Add in immature eggs
data$num.embryos2 <- data$num.embryos
data[which(is.na(data$num.embryos2)==TRUE),'num.embryos2'] <- 0

data$num.imm.eggs2 <- data$num.imm.eggs
data[which(is.na(data$num.imm.eggs2)==TRUE),'num.imm.eggs2'] <- 0

data$num.embryos <- data$num.embryos2 + data$num.imm.eggs2

data[which(data$final.sex=='M'),'num.embryos'] <- NA

#Calculate growth
#only survivors are measured
data$growth <- NA
data[which(data$survival==1),'growth'] <- data[which(data$survival==1),'recap.total.length'] - data[which(data$survival==1),'initial.total.length']
data[which(data$survival==0),'num.embryos'] <- NA
data[which(data$survival==0),'embryo.dry.mass'] <- NA





#Add in blocks
data[which(data$channel>=1 & data$channel<=4),'block'] <- 1
data[which(data$channel>=5 & data$channel<=8),'block'] <- 2
data[which(data$channel>=9 & data$channel<=12),'block'] <- 3
data[which(data$channel>=13 & data$channel<=16),'block'] <- 4
data$block <- as.factor(data$block)


#Create Channel.id which is unique to the experiment, drainage, replicate, and channel
data$channel.id <- do.call(paste, c(data[c("experiment", "drainage","replicate","channel")], sep = ".")) 

#Order the data
data <- data[order(data$channel.id,data$initial.total.length),]



data <- data[which(data$initial.fish==1),c('experiment','drainage','block','channel','channel.id','density','initial.fish','analysis.use','died.early',
                                           'final.sex','color.combo','riv.pheno','growth','survival','recap.total.length','pregnant','num.embryos',
                                           'stage','off.length','survival','initial.total.length','e1','e2','exper')]

#Create KO and LP columns. These will get used as dummy variables
data$KO <- 0
data$LP <- 0
data[which(data$riv.pheno=='KO'),'KO'] <- 1
data[which(data$riv.pheno=='LP'),'LP'] <- 1

data$drainage <- droplevels(data$drainage)

unique(data$drainage)

data[which(data$drainage==unique(data$drainage)[1]),'drain.num'] <- 1
data[which(data$drainage==unique(data$drainage)[2]), 'drain.num'] <- 2

#Create gdata, which will be the data set that ultimately gets analyzed
touse <- which( is.na(data$died.early)==TRUE  )
gdata <- data[touse,c('exper','drain.num','channel.id','block','final.sex','growth','num.embryos','pregnant','off.length','stage','survival','initial.total.length','KO','LP')]

#Create a unique number for each channel
for (i in unique(gdata$channel.id)){
  gdata[which((gdata$channel.id)==i),'channel.num'] <- which(unique(gdata$channel.id)==i)
}

#Calculate the number of fish in each channel
gdata$count <- 1
densities <- ddply(gdata,c('channel.num','KO'),summarise,density=sum(count))
densities <- reshape2::dcast(densities,channel.num~KO,value.var='density')
names(densities) <- c('channel.num','density.LP','density.KO')
densities[which(is.na(densities$density.LP)==T),'density.LP'] <- 0
densities[which(is.na(densities$density.KO)==T),'density.KO'] <- 0

gdata <- merge(x=gdata,y=densities,by=c('channel.num'),all.x=TRUE)


#To be able to include experiments where there is only one phenotype, need to add the density of the same phenotype.
#This will not affect estimate are these are marked with 'both.pheno' in the model statement. It is just to keep it from throwing errors.
#by including NA's
gdata$both.pheno <- 0
gdata[which(gdata$density.KO>0 & gdata$density.LP>0),'both.pheno'] <- 1


gdata <- gdata[,c('exper','drain.num','block','final.sex','growth','num.embryos','pregnant','off.length','stage','survival','density.LP','density.KO','initial.total.length','channel.num','KO','LP','both.pheno')]


#Not used for killifish
gdata[which(is.na(gdata$off.length)==TRUE),'stage'] <- 0

return(gdata)
}
