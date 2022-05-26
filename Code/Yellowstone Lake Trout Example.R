#### Fuzzy Cognitve Maps - 8/21/19
#### T.J. Clark

#### AIM: run FCM on yellowstone lake scenario, round2

###############################################################################

#### STEP 1: Build Matrix for null-model
#### based off of "yellowstone matrix per cap"
# flux lakers to simulate introduction

w <- matrix(0,nrow=14,ncol=14)
species <- c("sun","phyto","amphipod","zooSM","zooLG","cuts","sucker",
             "otter","bear","osprey","eagle","gullother",
             "lakers","elk")
dimnames(w) <- list(species,species)

# Sun
# added in b/c we want to monitor changes in phytoplankton
w["sun","sun"] <- 0
w["phyto","sun"] <- 1

# Phytoplankton
w["phyto","phyto"] <- 0.5
w["amphipod","phyto"] <- 1 
w["zooSM","phyto"] <- 1 
w["zooLG","phyto"] <- 1 
w["sucker","phyto"] <- 0.5

# Amphipods
w["phyto","amphipod"] <- -0.25
w["cuts","amphipod"] <- 0.9 
w["sucker","amphipod"] <- 0.1
w["lakers","amphipod"] <- 0.25

# Small Zooplankton
w["phyto","zooSM"] <- -0.5
w["cuts","zooSM"] <- 0.1
w["sucker","zooSM"] <- 0.1
w["lakers","zooSM"] <- 0.1

# Large Zooplankton
w["phyto","zooLG"] <- -0.5
w["cuts","zooLG"] <- 1 
w["sucker","zooLG"] <- 0.1
w["lakers","zooLG"] <- 0.1

# Cutthroat Trout
w["amphipod","cuts"] <- -0.9
w["zooSM","cuts"] <- -0.1
w["zooLG","cuts"] <- -0.75
w["otter","cuts"] <- 0.75
w["bear","cuts"] <- 0.1
w["osprey","cuts"] <- 0.9
w["eagle","cuts"] <- 0.25
w["gullother","cuts"] <- 0.1
w["lakers","cuts"] <- 0.75

# Longnose Suckers
w["phyto","sucker"] <- -0.1
w["amphipod","sucker"] <- -0.1
w["zooSM","sucker"] <- -0.1
w["zooLG","sucker"] <- -0.1
w["otter","sucker"] <- 0.75
w["bear","sucker"] <- 0.1
w["osprey","sucker"] <- 0.9
w["eagle","sucker"] <- 0.25
w["gullother","sucker"] <- 0.1
w["lakers","sucker"] <- 0.1

# Otter
w["cuts","otter"] <- -0.1
w["sucker","otter"] <- -0.1

# Bears
w["cuts","bear"] <- -0.1
w["sucker","bear"] <- -0.1
w["elk","bear"] <- -0.1
w["bear","bear"] <- 1 

# Osprey
w["cuts","osprey"] <- -0.1
w["sucker","osprey"] <- -0.1

# Eagle
w["cuts","eagle"] <- -0.1
w["sucker","eagle"] <- -0.1
w["eagle","eagle"] <- 1 

# Gulls, Terns, others...
w["cuts","gullother"] <- -0.1
w["sucker","gullother"] <- -0.1
w["gullother","gullother"] <- 1 

# Lakers
w["amphipod","lakers"] <- -0.25
w["zooSM","lakers"] <- -0.1
w["zooLG","lakers"] <- -0.1
w["cuts","lakers"] <- -0.9
w["sucker","lakers"] <- -0.1

# Elk
w["bear","elk"] <- 0.1
w["elk","elk"] <- 1.5 # make it survive

# null model run - do things work?
# holding lakers so that they can't affect things
s<- c(1,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0,0.5)
fix<- fix<- c(1,0,0,0,0,0,0,0,0,0,0,0,1,0)
af <- rep(c(1,rep(2,13)))

# note, c changed to 1, which follows rules
test<- fcm(s,w,tol=0.001,lambda=0.5,c=1,fix=fix,af=af)
colnames(test$state) <- c("Sun","Phyto","Amphipod","ZooSM","ZooLG","Cuts",
                          "Sucker","Otter","Bear","Osprey",
                          "Eagle","Gull-other","Laker","Elk")

n<- test$iter

windows();
plot(c(1,n),c(0,1),type='n',xlab="Time (Iteration #)",
     ylab="Relative abundance (Activation level)",
     bty='n',tck=0.02)
for(i in 2:ncol(test$state)) {
  lines(1:n,test$state[1:n,i],lty=1,col=i,lwd=2)
}
legend("topright",c("Phyto","Amphipod","ZooSM","ZooLG","Cuts",
                    "Sucker","Otter","Bear","Osprey",
                    "Eagle","Gull-other","Laker","Elk"),
       lty=1,col=1:13,bty='n',lwd=2)

View(test$state)

#########################################################################

#### Step 2) Build Membership Functions
# based off of the manuscript
# we will stick to these triangular functions...
# RULE: start is 50% in most cases
# high = 50%, 80%, 1
# middle = 30%, 50%, 70%
# low = 0%, 20%, 50%

# 1) Phyto abundance: beg = 2.2; end = 0.5
# based off of chlorophyll a concentration (ug/L)
high<- matrix(cbind(c(2.2,3.52,4.4),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(1.76,2.2,2.64),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,0.88,2.2),c(1,1,0)),ncol=2)

phyto.mf<- list(high=high,mod=mod,low=low)

# 2) Amphipod abundance: NO INFO
# using basic membership functions, starting at 0.5
high<- matrix(cbind(c(0.5,0.8,1),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(0.3,0.5,0.7),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,0.2,0.5),c(1,1,0)),ncol=2)

amphipod.mf<- list(high=high,mod=mod,low=low)

# 3) Small Zoo: beg = 67.8, end = 30.3
# based off of small zoo biomass in mg/L
high<- matrix(cbind(c(67.8,108.48,135.6),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(54.24,67.8,81.36),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,27.12,67.8),c(1,1,0)),ncol=2)

zooSM.mf<- list(high=high,mod=mod,low=low)

# 4) Large Zoo: beg = 8.5, end = 102.9
high<- matrix(cbind(c(21.25,34,42.5),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(12.75,21.25,29.75),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,8.5,21.25),c(1,1,0)),ncol=2)

zooLG.mf<- list(high=high,mod=mod,low=low)

# 5) Cutthroat Trout: beg = 44.6, end = 22.1
high<- matrix(cbind(c(44.6,71.36,89.2),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(35.68,44.6,53.52),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,17.84,44.6),c(1,1,0)),ncol=2)

cuts.mf<- list(high=high,mod=mod,low=low)

#6) Longnose Sucker: beg = 30, end = 4.6
high<- matrix(cbind(c(30,48,60),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(24,30,36),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,12,30),c(1,1,0)),ncol=2)

sucker.mf<- list(high=high,mod=mod,low=low)

#7) Otters: beg = 0.73, end = 0.14
# prevalence of scat (so we see a behavioral shift)
# this is pretty close to our normal proportional fuzz...
high<- matrix(cbind(c(0.73,0.8,1),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(0.6,0.73,0.8),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,0.2,0.73),c(1,1,0)),ncol=2)

otter.mf<- list(high=high,mod=mod,low=low)

#8) Bears: beg = 0.5, end = 0.2
# proportion of visits
high<- matrix(cbind(c(0.5,0.8,1),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(0.4,0.5,0.6),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,0.2,0.5),c(1,1,0)),ncol=2)

bear.mf<- list(high=high,mod=mod,low=low)

#9) Osprey: beg = 37.6, end = 3.2
# total nest count
high<- matrix(cbind(c(37.6,60.16,75.2),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(30.08,37.6,45.12),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,15.04,37.6),c(1,1,0)),ncol=2)

osprey.mf<- list(high=high,mod=mod,low=low)

#10) Bald Eagle: beg = 6.4, end = 7.6
# total nest count
# this is obscured by recovery from DDT
# range is 0 to 16 over the time period, made 18 the max

high<- matrix(cbind(c(9,14.4,18),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(4.8,9,12.6),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,3.6,8),c(1,1,0)),ncol=2)

eagle.mf<- list(high=high,mod=mod,low=low)

#11) Gull Other???
# same as amphipods...we don't have info, starting at 0.8
high<- matrix(cbind(c(0.5,0.8,1),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(0.4,0.5,0.6),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,0.2,0.5),c(1,1,0)),ncol=2)

gullother.mf<- list(high=high,mod=mod,low=low)

# 12) Lake Trout: beg = 0, end = 953 (in abundance)
# we don't need fuzzy for this...but let's do it anyway
# let's make 1000 = 1?

high<- matrix(cbind(c(500,800,1000),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(300,500,700),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,200,500),c(1,1,0)),ncol=2)

lakers.mf<- list(high=high,mod=mod,low=low)

# 13) Elk Calves: beg = .1, end = 0.5
# based off of proportion of elk calf in bear diet
# because high elk calf mortality doesn't necessarily translate into
# reduced lambda...

high<- matrix(cbind(c(0.25,0.4,.5),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(0.15,0.25,0.35),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,0.1,0.25),c(1,1,0)),ncol=2)

elk.mf<- list(high=high,mod=mod,low=low)

####
# BASE ONE IGNORE
high<- matrix(cbind(c(0.5,0.8,1),c(0,1,1)),ncol=2)
mod<- matrix(cbind(c(0.4,0.5,0.6),c(0,1,0)),ncol=2)
low<- matrix(cbind(c(0,0.2,0.5),c(1,1,0)),ncol=2)

crossbill.mf<- list(high=high,mod=mod,low=low)

# Test membership functions with graph
mf<- sucker.mf
max<- mf[[1]][3,1]
mean<- mf[[2]][2,1]
min<- 0
tmp<- rescale.mf(mf,mean,min,max)
n<- length(mf)

windows();
par(mfrow=c(1,2))
plot(c(0,max),c(0,1),type='n',tck=0.02,las=1,
     xlab="Lake Trout Abundance",ylab="Membership Value")
for(i in 1:n){
  lines(mf[[i]][,1],mf[[i]][,2])
}
title("Fuzzy membership functions (original scale)")
plot(c(0,1),c(0,1),type='n',tck=0.02,las=1,
     xlab="Relative Lake Trout Abundance",ylab="Membership Value")
for(i in 1:n){
  lines(tmp[[i]][,1],tmp[[i]][,2])
}
title("Fuzzy membership functions (rescaled)")

########################################################################

#### Step 3: RUN FCM under two scenarios
#### A) Pre Lake Trout (0 lake trout in sim)
#### B) Post Lake Trout (~950 lake trout in sim)

# set up data.frame!
result <- data.frame(matrix(NA,2,14))
result[,1] <- c("Pre","Post")
result[,13] <- c(0,1)
colnames(result) <- c("Treatment","Phyto","Amphipod","ZooSM","ZooLG",
                      "Cuts","Sucker","Otter","Bear","Osprey",
                      "Eagle","Gullother","Lakers","Elk")

# determine starting values via fuzzify, then defuzzify
# create a starting value fxn
start_fxn <- function(mf,start){
  start.fuzz <- fuzzify(mf, start, rescale=F, # rescale MUST BE FALSE
                        min=0,max=mf[[1]][3,1],
                        mean=mf[[2]][2,1])
  start.val <- defuzzify(mf,start.fuzz,rescale=T,
                       min=0,max=mf[[1]][3,1],
                       mean=mf[[2]][2,1])
  start.val[5]
}

phyto.s <- start_fxn(phyto.mf,2.2)
amphipod.s <- start_fxn(amphipod.mf,0.5)
zooSM.s <- start_fxn(zooSM.mf,67.8)
zooLG.s <- start_fxn(zooLG.mf,8.5)
cuts.s <- start_fxn(cuts.mf,44.6)
sucker.s <- start_fxn(sucker.mf,30)
otter.s <- start_fxn(otter.mf,0.73)
bear.s <- start_fxn(bear.mf, 0.5)
osprey.s <- start_fxn(osprey.mf, 37.6)
eagle.s <- start_fxn(eagle.mf,6.4)
gullother.s <- start_fxn(gullother.mf,0.5)
elk.s <- start_fxn(elk.mf,.1)

# Fuzzy Model Runs
# what happens when we introduce lake trout into Yellowstone Lake?
for(m in 1:nrow(result)){

  s <- c(5,phyto.s,amphipod.s,zooSM.s,zooLG.s,cuts.s,sucker.s,
         otter.s,bear.s,osprey.s,eagle.s,gullother.s,
         result[m,13],elk.s)
  fix<- c(1,0,0,0,0,0,0,0,0,0,0,0,1,0)
  af <- rep(c(1,rep(2,13)))

  test<- fcm(s,w,tol=0.001,lambda=0.5,c=1,fix=fix,af=af)

  result[m,2] <- test$state[test$iter,2] #phyto
  result[m,3] <- test$state[test$iter,3] #amphi
  result[m,4] <- test$state[test$iter,4] #zooSM
  result[m,5] <- test$state[test$iter,5] #zooLG
  result[m,6] <- test$state[test$iter,6] #cuts
  result[m,7] <- test$state[test$iter,7] #sucker
  result[m,8] <- test$state[test$iter,8] #otter
  result[m,9] <- test$state[test$iter,9] #bear
  result[m,10] <- test$state[test$iter,10] #osprey
  result[m,11] <- test$state[test$iter,11] #eagle
  result[m,12] <- test$state[test$iter,12] #gullother
  result[m,14] <- test$state[test$iter,14] #elk

}
View(result)

#########################################################################

#### Step 4: Bring Together Results

ling <- data.frame(t(result))
ling <- ling[-1,]
colnames(ling) <- c("Pre","Post")

# does direction of change match real results?
ling$direction <- c("Yes","No","Yes","Yes","Yes",
                    "Yes","Yes","Yes","Yes","Yes",
                    "Yes","NA","Yes")

# does fuzzyness match real results?
ling$fuzzysim <- c("High","NA","Medium","Medium","Low",
                   "Low","Low","Low","Low","Low",
                   "NA","NA","Low")

ling$fuzzyreal <- c("Low","NA","Low","High","Low",
                    "Low","Low","Low","Low","Medium",
                    "NA","NA","Low")

fuzzify(bear.mf, 0.1, rescale=F)

# does fuzziness match real results?
ling$final <- c("No","NA","No","No","Yes",
                "Yes","Yes","Yes","Yes","No",
                "NA","NA","Yes")

write.csv(ling,"Yellowstone Lake Results_round2.csv")
ling <- read.csv("Yellowstone Lake Results_round2.csv")

#########################################################################

#### Step 5: Let's dumb down everything: no tweaking/abund info, or sun

w <- matrix(0,nrow=13,ncol=13)
species <- c("phyto","amphipod","zooSM","zooLG","cuts","sucker",
             "otter","bear","osprey","eagle","gullother",
             "lakers","elk")
dimnames(w) <- list(species,species)

# Phytoplankton
w["phyto","phyto"] <- 0.5
w["amphipod","phyto"] <- 0.25
w["zooSM","phyto"] <- 0.25
w["zooLG","phyto"] <- 0.25
w["sucker","phyto"] <- 0.1

# Amphipods
w["phyto","amphipod"] <- -0.25
w["cuts","amphipod"] <- 0.25
w["sucker","amphipod"] <- 0.1
w["lakers","amphipod"] <- 0.25

# Small Zooplankton
w["phyto","zooSM"] <- -0.5
w["cuts","zooSM"] <- 0.1
w["sucker","zooSM"] <- 0.1
w["lakers","zooSM"] <- 0.1

# Large Zooplankton
w["phyto","zooLG"] <- -0.5
w["cuts","zooLG"] <- 0.25
w["sucker","zooLG"] <- 0.1
w["lakers","zooLG"] <- 0.1

# Cutthroat Trout
w["amphipod","cuts"] <- -0.9
w["zooSM","cuts"] <- -0.1
w["zooLG","cuts"] <- -0.75
w["otter","cuts"] <- 0.75
w["bear","cuts"] <- 0.1
w["osprey","cuts"] <- 0.9
w["eagle","cuts"] <- 0.25
w["gullother","cuts"] <- 0.1
w["lakers","cuts"] <- 0.75

# Longnose Suckers
w["phyto","sucker"] <- -0.1
w["amphipod","sucker"] <- -0.1
w["zooSM","sucker"] <- -0.1
w["zooLG","sucker"] <- -0.1
w["otter","sucker"] <- 0.75
w["bear","sucker"] <- 0.1
w["osprey","sucker"] <- 0.9
w["eagle","sucker"] <- 0.25
w["gullother","sucker"] <- 0.1
w["lakers","sucker"] <- 0.1

# Otter
w["cuts","otter"] <- -0.1
w["sucker","otter"] <- -0.1

# Bears
w["cuts","bear"] <- -0.1
w["sucker","bear"] <- -0.1
w["elk","bear"] <- -0.1

# Osprey
w["cuts","osprey"] <- -0.1
w["sucker","osprey"] <- -0.1

# Eagle
w["cuts","eagle"] <- -0.1
w["sucker","eagle"] <- -0.1

# Gulls, Terns, others...
w["cuts","gullother"] <- -0.1
w["sucker","gullother"] <- -0.1

# Lakers
w["amphipod","lakers"] <- -0.25
w["zooSM","lakers"] <- -0.1
w["zooLG","lakers"] <- -0.1
w["cuts","lakers"] <- -0.9
w["sucker","lakers"] <- -0.1

# Elk
w["bear","elk"] <- 0.1

# null model run - do things work?
# holding lakers so that they can't affect things
s<- c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
fix<- c(0,0,0,0,0,0,0,0,0,0,0,0,0)
af <- rep(c(1,rep(2,12)))

# note, c changed to 1, which follows rules
test<- fcm(s,w,tol=0.001,lambda=0.5,c=1,fix=fix,af=af)
colnames(test$state) <- c("Phyto","Amphipod","ZooSM","ZooLG","Cuts",
                          "Sucker","Otter","Bear","Osprey",
                          "Eagle","Gull-other","Laker","Elk")

n<- test$iter

windows();
plot(c(1,n),c(0,1),type='n',xlab="Time (Iteration #)",
     ylab="Relative abundance (Activation level)",
     bty='n',tck=0.02)
for(i in 2:ncol(test$state)) {
  lines(1:n,test$state[1:n,i],lty=1,col=i,lwd=2)
}
legend("topright",c("Phyto","Amphipod","ZooSM","ZooLG","Cuts",
                    "Sucker","Otter","Bear","Osprey",
                    "Eagle","Gull-other","Laker","Elk"),
       lty=1,col=1:13,bty='n',lwd=2)

##

result <- data.frame(matrix(NA,2,14))
result[,1] <- c("Pre","Post")
result[,13] <- c(0,1)
colnames(result) <- c("Treatment","Phyto","Amphipod","ZooSM","ZooLG",
                      "Cuts","Sucker","Otter","Bear","Osprey",
                      "Eagle","Gullother","Lakers","Elk")

for(m in 1:nrow(result)){

  s <- c(0.5,0.5,0.5,0.5,0.5,0.5,
         0.5,0.5,0.5,0.5,0.5,
         result[m,13],0.5)
  fix<- c(0,0,0,0,0,0,0,0,0,0,0,1,0)
  af <- rep(c(1,rep(2,12)))

  test<- fcm(s,w,tol=0.001,lambda=0.5,c=1,fix=fix,af=af)

  result[m,2] <- test$state[test$iter,1] #phyto
  result[m,3] <- test$state[test$iter,2] #amphi
  result[m,4] <- test$state[test$iter,3] #zooSM
  result[m,5] <- test$state[test$iter,4] #zooLG
  result[m,6] <- test$state[test$iter,5] #cuts
  result[m,7] <- test$state[test$iter,6] #sucker
  result[m,8] <- test$state[test$iter,7] #otter
  result[m,9] <- test$state[test$iter,8] #bear
  result[m,10] <- test$state[test$iter,9] #osprey
  result[m,11] <- test$state[test$iter,10] #eagle
  result[m,12] <- test$state[test$iter,11] #gullother
  result[m,14] <- test$state[test$iter,13] #elk

}
View(result)

##

ling$dir_noabund <- c("No","Yes","Yes","No","Yes","Yes",
                      "Yes","Yes","Yes","Yes","Yes",
                      "NA","Yes") # 10/12
ling$fuzzysim_noabund <- c("Medium","NA","Low","Low","Low",
                           "Low","Low","Low","Low","Low",
                           "NA","NA","Low")
ling$final_noabund <- c("No","NA","Yes","No","Yes","Yes",
                        "Yes","Yes","Yes","No","NA","NA","Yes")

##########################################################################


#### Step 6: no info on interaction coefficients

w <- matrix(0,nrow=13,ncol=13)
species <- c("phyto","amphipod","zooSM","zooLG","cuts","sucker",
             "otter","bear","osprey","eagle","gullother",
             "lakers","elk")
dimnames(w) <- list(species,species)

# Phytoplankton
w["phyto","phyto"] <- 1
w["amphipod","phyto"] <- 1
w["zooSM","phyto"] <- 1
w["zooLG","phyto"] <- 1
w["sucker","phyto"] <- 1

# Amphipods
w["phyto","amphipod"] <- -1
w["cuts","amphipod"] <- 1
w["sucker","amphipod"] <- 1
w["lakers","amphipod"] <- 1

# Small Zooplankton
w["phyto","zooSM"] <- -1
w["cuts","zooSM"] <- 1
w["sucker","zooSM"] <- 1
w["lakers","zooSM"] <- 1

# Large Zooplankton
w["phyto","zooLG"] <- -1
w["cuts","zooLG"] <- 1
w["sucker","zooLG"] <- 1
w["lakers","zooLG"] <- 1

# Cutthroat Trout
w["amphipod","cuts"] <- -1
w["zooSM","cuts"] <- -1
w["zooLG","cuts"] <- -1
w["otter","cuts"] <- 1
w["bear","cuts"] <- 1
w["osprey","cuts"] <- 1
w["eagle","cuts"] <- 1
w["gullother","cuts"] <- 1
w["lakers","cuts"] <- 1

# Longnose Suckers
w["phyto","sucker"] <- -1
w["amphipod","sucker"] <- -1
w["zooSM","sucker"] <- -1
w["zooLG","sucker"] <- -1
w["otter","sucker"] <- 1
w["bear","sucker"] <- 1
w["osprey","sucker"] <- 1
w["eagle","sucker"] <- 1
w["gullother","sucker"] <- 1
w["lakers","sucker"] <- 1

# Otter
w["cuts","otter"] <- -1
w["sucker","otter"] <- -1

# Bears
w["cuts","bear"] <- -1
w["sucker","bear"] <- -1
w["elk","bear"] <- -1

# Osprey
w["cuts","osprey"] <- -1
w["sucker","osprey"] <- -1

# Eagle
w["cuts","eagle"] <- -1
w["sucker","eagle"] <- -1

# Gulls, Terns, others...
w["cuts","gullother"] <- -1
w["sucker","gullother"] <- -1

# Lakers
w["amphipod","lakers"] <- -1
w["zooSM","lakers"] <- -1
w["zooLG","lakers"] <- -1
w["cuts","lakers"] <- -1
w["sucker","lakers"] <- -1

# Elk
w["bear","elk"] <- 1

# null model run - do things work?
# holding lakers so that they can't affect things
s<- c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
fix<- c(0,0,0,0,0,0,0,0,0,0,0,0,0)
af <- rep(c(1,rep(2,12)))

# note, c changed to 1, which follows rules
test<- fcm(s,w,tol=0.001,lambda=.1,c=1,fix=fix,af=af)
# NOTE: set lambda to 0.1 to speed convergence
colnames(test$state) <- c("Phyto","Amphipod","ZooSM","ZooLG","Cuts",
                          "Sucker","Otter","Bear","Osprey",
                          "Eagle","Gull-other","Laker","Elk")

n<- test$iter

windows();
plot(c(1,n),c(0,1),type='n',xlab="Time (Iteration #)",
     ylab="Relative abundance (Activation level)",
     bty='n',tck=0.02)
for(i in 2:ncol(test$state)) {
  lines(1:n,test$state[1:n,i],lty=1,col=i,lwd=2)
}
legend("topright",c("Phyto","Amphipod","ZooSM","ZooLG","Cuts",
                    "Sucker","Otter","Bear","Osprey",
                    "Eagle","Gull-other","Laker","Elk"),
       lty=1,col=1:13,bty='n',lwd=2)

##

result <- data.frame(matrix(NA,2,14))
result[,1] <- c("Pre","Post")
result[,13] <- c(0,1)
colnames(result) <- c("Treatment","Phyto","Amphipod","ZooSM","ZooLG",
                      "Cuts","Sucker","Otter","Bear","Osprey",
                      "Eagle","Gullother","Lakers","Elk")

for(m in 1:nrow(result)){

  s <- c(0.5,0.5,0.5,0.5,0.5,0.5,
         0.5,0.5,0.5,0.5,0.5,
         result[m,13],0.5)
  fix<- c(0,0,0,0,0,0,0,0,0,0,0,1,0)
  af <- rep(c(1,rep(2,12)))

  test<- fcm(s,w,tol=0.001,lambda=0.1,c=1,fix=fix,af=af)

  result[m,2] <- test$state[test$iter,1] #phyto
  result[m,3] <- test$state[test$iter,2] #amphi
  result[m,4] <- test$state[test$iter,3] #zooSM
  result[m,5] <- test$state[test$iter,4] #zooLG
  result[m,6] <- test$state[test$iter,5] #cuts
  result[m,7] <- test$state[test$iter,6] #sucker
  result[m,8] <- test$state[test$iter,7] #otter
  result[m,9] <- test$state[test$iter,8] #bear
  result[m,10] <- test$state[test$iter,9] #osprey
  result[m,11] <- test$state[test$iter,10] #eagle
  result[m,12] <- test$state[test$iter,11] #gullother
  result[m,14] <- test$state[test$iter,13] #elk

}
View(result)

##

ling$dir_noint <- c("No","Yes","Yes","No","Yes","Yes",
                      "Yes","Yes","Yes","Yes","Yes",
                      "NA","Yes") # still 10/12
ling$fuzzysim_noint <- c("High","NA","Low","Low","Low",
                           "Low","Low","Low","Low","Low",
                           "NA","NA","Low")
ling$final_noint <- c("No","NA","Yes","No","Yes","Yes",
                      "Yes","Yes","Yes","No","NA","NA","Yes")

#### Caution: still pretty bad at predicting, Low is just Low...
#### be cautious in how you interpret these results!

new_FCMresults <- write.csv(ling,"YellowstoneLake_results3.csv")

