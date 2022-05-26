#### Simple fuzzy cognitive map
#### This sets up FCM and the null argument that nothing changes

fcm<- function(s,w,fix=rep(0,length(s)),
               lambda=0.5, max.iter=500,
               af=rep(1,length(s)), tol=0.001,c=5) {
#
# Fuzzy cognitive map algorithm
# s is the initial state vector variable [0,1]
# w is the matrix of edge weights [-1,1]
# fix is the vector indicating which element of s are fixed boundary conditions
# lambda is the relaxation parameter to speed convergence
# activation functions are (1)logistic, (2) threshold exponential or (3) threshold linear
# af - vector specifying which activation function to be used on each element of s
#
# D. Ramsey 30/4/04
#

state<- matrix(nrow=max.iter,ncol=length(s))
raw<- matrix(nrow=max.iter,ncol=length(s))
state[1,]<- s
raw[1,]<- s

iter<- 1

for(z in 2:max.iter) {

                tmp<- w %*% state[z-1,]
    for(i in 1:nrow(w)) {


        if(af[i]==1){
            # logistic
            state[z,i] <- lambda* (1/(1+exp(-c*tmp[i]))) + (1-lambda)*state[z-1,i]
            raw[z,i]<- tmp[i]
        }
        else if(af[i]==2) {
            # exponential
            state[z,i]<- lambda * max(0,1-exp(-c*tmp[i])) + (1-lambda)*state[z-1,i]
            raw[z,i]<- tmp[i]
        }
        else if(af[i]==3){
            # linear
            state[z,i] <- lambda* (min(1,max(0,c*tmp[i]))) + (1-lambda)*state[z-1,i]
            raw[z,i]<- tmp[i]
        }

    }
        state[z,fix==1]<- state[1,fix==1] # reset boundary conditions (control)


    ind<- abs(state[z,] - state[z-1,]) < tol

    if(length(state[z,ind])==length(s)) break
        iter<- iter + 1
}

    if(iter == max.iter) cat(" WARNING ! Convergence not reached in ",max.iter," iterations",'\n')
        else cat("convergence reached after ",iter," iterations",'\n')

list(state=state[1:iter,],raw=raw[1:iter,],iter=iter)

}

#============================================================
#make adjacency matrix (w)
# NULL MODEL: interaction matrix
# also includes results from the Hebbian learning algorithm

w<- matrix(0,nrow=5,ncol=5)
species<- c("trees","kokako","rat","possum","stoat")
dimnames(w)<- list(species,species)

w['possum','trees']<- 0.5
w['rat','trees']<- 0.5
w['kokako','trees']<- 0.5
w['stoat','rat']<- 0.5

w['kokako','possum']<- -0.17
w['kokako','rat']<- -0.67
w['kokako','stoat']<- -0.2

w['trees','possum']<- -0.20
w['rat','stoat']<- -0.20
w['trees','kokako']<- -0.20
w['trees','rat']<- -0.20
w['stoat','kokako']<- 0.20
w['rat','kokako']<- 0.20
w['possum','kokako']<- 0.20

diag(w)<- -0.2
w['kokako','kokako']<- 0
w['trees','trees']<- 0.4

#=======================================
# Run NULL model (no management)

s<- c(0.5,0.5,0.5,0.5,0.5)
fix<- c(0,0,0,0,0)
af<- c(1,2,2,2,2)

test<- fcm(s,w,tol=0.001,lambda=0.5,c=5,fix=fix,af=af)

n<- test$iter
windows();
plot(c(1,n),c(0,1),type='n',xlab="Iteration number",ylab="Activation level",
bty='n',tck=0.02)
for(i in 1:ncol(test$state)) {
lines(1:n,test$state[1:n,i],lty=1,col=i)
}
legend(n/2,0.9,c("fruit/foliage","kokako","rat",'possum','stoat'),lty=1,col=1:5,bty='n')


