#### T.J. Clark - 2/14/19
#### Newfoundland Red squirrel/crossbill example

#### Two objectives:
# a) build FCMs with squirrel-crossbill-marten-spruce web

#################################################################################
# FULL MATRIX

w <- matrix(0,nrow=12,ncol=12)
species <- c("spruce","other conifer","new crossbill",
             "other crossbill", "seed eating bird",
             "squirrel", "marten","goshawk", "rodent",
             "grouse","hare","lynx")
dimnames(w) <- list(species,species)

## PRODUCERS
# black spruce
w["spruce","spruce"] <- 0.5
w["other conifer","spruce"] <- -0.1
w["new crossbill","spruce"] <- 0.9
w["other crossbill","spruce"] <- 0.5
w["squirrel","spruce"] <- 0.75
w["seed eating bird","spruce"] <- 0.5

# other conifers
w["other conifer","other conifer"] <- 0.5
w["spruce","other conifer"] <- -0.1
w["other crossbill","other conifer"] <- 0.5
w["seed eating bird","other conifer"] <- 0.5
w["squirrel","other conifer"] <- 0.75

## HERBIVORES
# Newfoundland crossbill
w["spruce","new crossbill"] <- -0.25
w["marten","new crossbill"] <- 0.1
w["goshawk","new crossbill"] <- 0.1

# Other crossbills (white-winged and red)
w["spruce","other crossbill"] <- -0.25
w["other conifer","other crossbill"] <- -0.25
w["marten","other crossbill"] <- 0.1
w["goshawk","other crossbill"] <- 0.1

# Other Seed-eating Birds (grosbeaks, finches, etc.)
w["spruce","seed eating bird"] <- -0.25
w["other conifer","seed eating bird"] <- -0.25
w["marten","seed eating bird"] <- 0.1
w["goshawk","seed eating bird"] <- 0.1

# Red Squirrels
w["spruce","squirrel"] <-  -.5
w["other conifer","squirrel"] <- -0.5
w["new crossbill","squirrel"] <- -0.25 # these are trait-based indirect
w["goshawk","squirrel"] <- 0.5
w["marten","squirrel"] <- 0.25

# Rodents
w["rodent","rodent"] <- 1
w["marten","rodent"] <- 0.75

## PREDATORS
# martens
w["squirrel","marten"] <- -0.25
w["rodent","marten"] <- -0.25
w["new crossbill","marten"] <- -0.1
w["other crossbill","marten"] <- -0.1
w["seed eating bird","marten"] <- -0.1

# goshawks
w["squirrel","goshawk"] <- -0.5
w["new crossbill","goshawk"] <- -0.1
w["other crossbill","goshawk"] <- -0.1
w["seed eating bird","goshawk"] <- -0.1
w["grouse", "goshawk"] <- -0.25
w["hare","goshawk"] <- -0.25

# grouse
w["goshawk","grouse"] <- 0.25
w["lynx","grouse"] <- 0.1
w["grouse","grouse"] <- 1

# hare
w["goshawk","hare"] <- 0.25
w["lynx","hare"] <- 0.75
w["hare","hare"] <- 1

#lynx
w["grouse","lynx"] <- -0.1
w["hare","lynx"] <- -0.35
w["lynx","lynx"] <- 0

### Run null-model
s<- rep(0.5,12)
fix<- rep(0,12)
#af<- c(1,1,2,2,2,2,2,2,2,2,2,2)
af<- c(1,1,2,2,2,2,2,2,1,1,1,2) # try fixing the populations...

# note, c changed to 1, which follows rules
test<- fcm(s,w,tol=0.001,lambda=0.5,c=1,fix=fix,af=af)

n<- test$iter

windows();
plot(c(1,n),c(0,1),type='n',xlab="Time (Iteration #)",
     ylab="Relative abundance (Activation level)",
     main="Null Model Runs for Expanded Food Web",
     bty='n',tck=0.02)
for(i in 1:ncol(test$state)) {
  lines(1:n,test$state[1:n,i],lty=1,col=i,lwd=2)
}
legend("topright",c("Black Spruce","Other Conifers",
                    "Newfoundland Crossbill","Other Crossbills",
                    "Other Seed-eating Birds",
                    "Red Squirrel","Pine Marten","Goshawk","Rodents",
                    "Grouse","Hare","Lynx"),
       lty=1,col=1:12,bty='n',lwd=2)

View(test$state)

## Introduce the Squirrels
# iterate number of squirrels and watch crossbill abundance
squirrel <- seq(0,1,0.1)

# other variable of interest
result <- data.frame(squirrel)
result[,2] <- rep(NA,nrow(result)) # spruce
result[,3] <- rep(NA,nrow(result)) # other conifer
result[,4] <- rep(NA,nrow(result)) # new crossbill
result[,5] <- rep(NA,nrow(result)) # other crossbills
result[,6] <- rep(NA,nrow(result)) # seed-eating birds
result[,7] <- rep(NA,nrow(result)) # marten
result[,8] <- rep(NA,nrow(result)) # goshawks
result[,9] <- rep(NA,nrow(result)) # rodents
result[,10] <- rep(NA,nrow(result)) # grouse
result[,11] <- rep(NA,nrow(result)) # hares
result[,12] <- rep(NA,nrow(result)) # lynx

# what happens when we introduce more squirrels?
for(m in 1:nrow(result)){

  s<- c(0.5,0.5,0.5,0.5,0.5,result[m,1],0.5,0.5,0.5,0.5,0.5,0.5)
  fix<- c(0,0,0,0,0,1,0,0,0,0,0,0)
  af<- c(1,1,2,2,2,2,2,2,2,2,2,2)
  #af<- c(1,1,2,2,2,2,2,2,1,1,1,2) # try fixing the populations...

  test<- fcm(s,w,tol=0.001,lambda=0.5,c=1,fix=fix,af=af)

  result[m,1] <- test$state[test$iter,6] # squirrels
  result[m,2] <- test$state[test$iter,1]
  result[m,3] <- test$state[test$iter,2]
  result[m,4] <- test$state[test$iter,3]
  result[m,5] <- test$state[test$iter,4]
  result[m,6] <- test$state[test$iter,5]
  result[m,7] <- test$state[test$iter,7]
  result[m,8] <- test$state[test$iter,8]
  result[m,9] <- test$state[test$iter,9]
  result[m,10] <- test$state[test$iter,10]
  result[m,11] <- test$state[test$iter,11]
  result[m,12] <- test$state[test$iter,12]

}
# results = end state of each converged sample...no fuzziness yet
names(result) <- c("squirrel_intro", "spruce", "other_conifer",
                   "new_crossbill","other_crossbill","seed_eaters",
                   "marten","goshawks","rodents",
                   "grouse","hare","lynx")

# Lets make the crossbills fuzzy!
# Fuzzify and defuzzify
crossbill.f<- fuzzify(crossbill.mf,result$lynx,rescale=T,
                      min=0,max=1,mean=0.5)
crossbill.df<- defuzzify(crossbill.mf,crossbill.f,rescale=T,
                         min=0,max=1,mean=0.5)
result$lynx.df <- crossbill.df[,5]

# plot
windows();
plot(result$squirrel_intro,result$spruce.df,
     type="b",lwd=2,col="black",ylim=c(0,1),
     xlab="Relative Red Squirrel Abundance (at Introduction)",
     ylab="Relative Species Abundance (at Equilibrium)",
     main="Effect of the Introduction of Red Squirrels on
     the Newfoundland Ecosystem")
lines(result$squirrel_intro,result$otherconifer.df,lwd=2,type="b",col="red")
lines(result$squirrel_intro,result$newcrossbill.df,lwd=2,type="b",col="blue")
lines(result$squirrel_intro,result$othercrossbill.df,lwd=2,type="b",col="green")
lines(result$squirrel_intro,result$seedeaters.df,lwd=2,type="b",col="purple")
lines(result$squirrel_intro,result$marten.df,lwd=2,type="b",col="pink")
lines(result$squirrel_intro,result$goshawk.df,lwd=2,type="b",col="yellow")
lines(result$squirrel_intro,result$rodent.df, lwd=2, type="b",col="coral")
lines(result$squirrel_intro,result$grouse.df, lwd=2, type="b",col="grey")
lines(result$squirrel_intro,result$hare.df, lwd=2, type="b",col="orange")
lines(result$squirrel_intro,result$lynx.df, lwd=2, type="b",col="brown")
legend("topright",c("Black Spruce","Other Conifers",
                    "Newfoundland Crossbill","Other Crossbills and Birds",
                    "Pine Marten","Goshawk","Rodent","Grouse","Hare","Lynx"),
       lty=1,col=c("black","red","blue","purple","pink","yellow","coral","grey",
                   "orange","brown"),
       bty='n',lwd=2)

#################################################################################

#### SENSITIVITY/UNCERTAINTY ANALYSIS - 3/28/22
# two ideas: 1) percent change in abundance following full intro of squirrels; 2) % of sims in which crossbills decrease...

# uncertainty analysis for loop, following Baker et al. 2018
# going to change interactions that affected red squirrels or newfoundland crossbills between -1 to 0 or 0 to 1
# then show responses of species...or maybe just stick to red squirrels...

sims <- 10000 # 10000

uncert <- rep(NA,sims)

for (i in 1:sims){

w <- matrix(0,nrow=12,ncol=12)
species <- c("spruce","other conifer","new crossbill",
             "other crossbill", "seed eating bird",
             "squirrel", "marten","goshawk", "rodent",
             "grouse","hare","lynx")
dimnames(w) <- list(species,species)

## PRODUCERS
# black spruce
w["spruce","spruce"] <- 0.5
w["other conifer","spruce"] <- -0.1
w["new crossbill","spruce"] <- runif(1,0.25,1)#0.9
w["other crossbill","spruce"] <- 0.5
w["squirrel","spruce"] <- runif(1,0.25,1)#0.75
w["seed eating bird","spruce"] <- 0.5
# other conifers
w["other conifer","other conifer"] <- 0.5
w["spruce","other conifer"] <- -0.1
w["other crossbill","other conifer"] <- 0.5
w["seed eating bird","other conifer"] <- 0.5
w["squirrel","other conifer"] <- 0.75
# Newfoundland crossbill
w["spruce","new crossbill"] <- runif(1,-0.75,-0.01)#-0.25
w["marten","new crossbill"] <- 0.1
w["goshawk","new crossbill"] <- 0.1
# Other crossbills (white-winged and red)
w["spruce","other crossbill"] <- -0.25
w["other conifer","other crossbill"] <- -0.25
w["marten","other crossbill"] <- 0.1
w["goshawk","other crossbill"] <- 0.1
# Other Seed-eating Birds (grosbeaks, finches, etc.)
w["spruce","seed eating bird"] <- -0.25
w["other conifer","seed eating bird"] <- -0.25
w["marten","seed eating bird"] <- 0.1
w["goshawk","seed eating bird"] <- 0.1
# Red Squirrels
w["spruce","squirrel"] <-  runif(1,-1,-0.25)#-.5
w["other conifer","squirrel"] <- -0.5
w["new crossbill","squirrel"] <- runif(1,-0.75,-0.25)#-0.25 # these are trait-based indirect
w["goshawk","squirrel"] <- 0.5
w["marten","squirrel"] <- 0.25
# Rodents
w["rodent","rodent"] <- 1
w["marten","rodent"] <- 0.75
# martens
w["squirrel","marten"] <- -0.25
w["rodent","marten"] <- -0.25
w["new crossbill","marten"] <- -0.1
w["other crossbill","marten"] <- -0.1
w["seed eating bird","marten"] <- -0.1
# goshawks
w["squirrel","goshawk"] <- -0.5
w["new crossbill","goshawk"] <- -0.1
w["other crossbill","goshawk"] <- -0.1
w["seed eating bird","goshawk"] <- -0.1
w["grouse", "goshawk"] <- -0.25
w["hare","goshawk"] <- -0.25
# grouse
w["goshawk","grouse"] <- 0.25
w["lynx","grouse"] <- 0.1
w["grouse","grouse"] <- 1
# hare
w["goshawk","hare"] <- 0.25
w["lynx","hare"] <- 0.75
w["hare","hare"] <- 1
#lynx
w["grouse","lynx"] <- -0.1
w["hare","lynx"] <- -0.35
w["lynx","lynx"] <- 0

## RUN CODE
# iterate number of squirrels and watch crossbill abundance
squirrel <- seq(0,1,0.1)

# other variable of interest
result <- data.frame(squirrel)
result[,2] <- rep(NA,nrow(result)) # spruce
result[,3] <- rep(NA,nrow(result)) # other conifer
result[,4] <- rep(NA,nrow(result)) # new crossbill
result[,5] <- rep(NA,nrow(result)) # other crossbills
result[,6] <- rep(NA,nrow(result)) # seed-eating birds
result[,7] <- rep(NA,nrow(result)) # marten
result[,8] <- rep(NA,nrow(result)) # goshawks
result[,9] <- rep(NA,nrow(result)) # rodents
result[,10] <- rep(NA,nrow(result)) # grouse
result[,11] <- rep(NA,nrow(result)) # hares
result[,12] <- rep(NA,nrow(result)) # lynx

# what happens when we introduce more squirrels?
for(m in 1:nrow(result)){
  
  s<- c(0.5,0.5,0.5,0.5,0.5,result[m,1],0.5,0.5,0.5,0.5,0.5,0.5)
  fix<- c(0,0,0,0,0,1,0,0,0,0,0,0)
  af<- c(1,1,2,2,2,2,2,2,2,2,2,2)
  #af<- c(1,1,2,2,2,2,2,2,1,1,1,2) # try fixing the populations...
  
  test<- fcm(s,w,tol=0.001,lambda=0.5,c=1,fix=fix,af=af)
  
  result[m,1] <- test$state[test$iter,6] # squirrels
  result[m,2] <- test$state[test$iter,1]
  result[m,3] <- test$state[test$iter,2]
  result[m,4] <- test$state[test$iter,3]
  result[m,5] <- test$state[test$iter,4]
  result[m,6] <- test$state[test$iter,5]
  result[m,7] <- test$state[test$iter,7]
  result[m,8] <- test$state[test$iter,8]
  result[m,9] <- test$state[test$iter,9]
  result[m,10] <- test$state[test$iter,10]
  result[m,11] <- test$state[test$iter,11]
  result[m,12] <- test$state[test$iter,12]
  
}
# results = end state of each converged sample...no fuzziness yet
names(result) <- c("squirrel_intro", "spruce", "other_conifer",
                   "new_crossbill","other_crossbill","seed_eaters",
                   "marten","goshawks","rodents",
                   "grouse","hare","lynx")

# only record results from red squirrels
uncert[i] <- (result$new_crossbill[11]-result$new_crossbill[1])/result$new_crossbill[1]

}

# red squirrel results (don't need a graph)
mean(uncert) # -0.9526 (normal analysis was a -0.99% change)
quantile(uncert, c(0.05,0.95)) # -0.9968687, -0.6588288




