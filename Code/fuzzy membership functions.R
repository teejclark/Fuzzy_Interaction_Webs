#
# Fuzzy membership functions


fuzzify<- function(mf, x, rescale=F, min=NULL,mean=NULL,max=NULL){
	
# find fuzzy membership of x in each fuzzy set mf
# mf is a list of fuzzy set values (membership functions)
# triangular, trapezoid and shelf functions only
# output is x and the membership of x in each fuzzy set in mf 
# rescale first rescales mf on the unit interval using rescale.mf()
#
# D. Ramsey 25/4/04

if(rescale) {
	mf<- rescale.mf(mf,mean,min,max)
}
n<- length(mf) # we have n fuzzy sets
numx<- length(x)
result<- matrix(NA,nrow=numx,ncol=(n+1))
result[,1]<- x

mu.x<- rep(NA,n) # membership value of x in each set 

for(i in 1:numx) {
	xi<- x[i]
	
	for(j in 1:n){

	# 1st and last elements of mf are shelf functions

		if(j==1) {
			left<- mf[[j]][1,1]
			right<- mf[[j]][2,1]

			if(xi >= right) mu.x[j]<- 1
				else
			if(xi > left & xi < right) mu.x[j]<- 1-(right - xi)/(right-left)
				else mu.x[j]<- 0
				}
		else if(j==n)
			{# last is a shelf
	
			left<- mf[[j]][2,1]
			right<- mf[[j]][3,1]

			if(xi <= left) mu.x[j]<- 1
				else
			if(xi > left & xi < right) mu.x[j]<- (right - xi)/(right-left)
				else mu.x[j]<- 0
			}

		else if(nrow(mf[[j]])==3)
			{# triangular
			left<- mf[[j]][1,1]
			centre<- mf[[j]][2,1]
			right<- mf[[j]][3,1]

			if(xi > left & xi <= centre) mu.x[j]<- 1-(centre - xi)/(centre-left)
				else
			if(xi > centre & xi < right) mu.x[j]<- (right - xi)/(right-centre)
				else 	mu.x[j]<- 0
				}
		else if(nrow(mf[[j]])==4)
			{# trapezoid
			left<- mf[[j]][1,1]
			centre1<- mf[[j]][2,1]
			centre2<- mf[[j]][3,1]
			right<- mf[[j]][4,1]

			if(xi > left & xi < centre1) mu.x[j]<- 1-(centre1 - xi)/(centre1-left)
				else
			if(xi > centre2 & xi < right) mu.x[j]<- (right - xi)/(right-centre2)
				else
			if(xi >=centre1 & xi <= centre2) mu.x[j]<- 1
				else 	mu.x[j]<- 0
			}# end trapezoid
					
	
		} #j

result[i,c(2:(n+1))]<- mu.x

	}#i
	return(result)
}

#------------------------------------------------------------


defuzzify<- function(mf,fs, rescale=F,min=NULL,mean=NULL,max=NULL,nint=100){
#
# Defuzzify using "centre of gravity" method (numerical integration of piecewise fuzzy sets)
# fs is output of fuzzify()
# rescale - scale fs to the unit interval and defuzzify 
# or else map unit interval back to real scale				
# nint = delta x - number of intervals used for integration

tmp<- mf
n<- length(tmp)

if(rescale) {
	tmp<- 	rescale.mf(mf,mean,min,max)
	dx<- seq(0,1,1/nint)
	}
else {
	tmp<- mf
	dx<- seq(0,max,max/nint)
	} 
		
out<- rep(NA,nrow(fs))

for(i in 1:nrow(fs)){
	
	mu.x<- fs[i,2:(n+1)]
	sumfx<- rep(NA,length(dx))
	sumfx.j<- rep(NA,n)
	
	
	for(j in 1:length(dx)){
			sumfx.j<- rep(0,n)
		for(k in 1:n){
			# 1st j is  shelf

			if(mu.x[k]>0 & k==1) {
				left<- tmp[[k]][1,1]
				right<- tmp[[k]][2,1]

				if(dx[j] >= right) sumfx.j[k]<- 1
					else
				if(dx[j] > left & dx[j] < right) sumfx.j[k]<- 1-(right - dx[j])/(right-left)
					else sumfx.j[k]<- 0
					
					sumfx.j[k]<- min(sumfx.j[k],mu.x[k])	
				}
			else if(mu.x[k]>0 & k==n)
				{
					# last j is a shelf
	
					left<- tmp[[k]][2,1]
					right<- tmp[[k]][3,1]

					if(dx[j] <= left) sumfx.j[k]<- 1
						else
					if(dx[j] > left & dx[j] < right) sumfx.j[k]<- (right - dx[j])/(right-left)
						else sumfx.j[k]<- 0
				
					sumfx.j[k]<- min(sumfx.j[k],mu.x[k])
				}

			else if(mu.x[k]>0){
					if(nrow(tmp[[k]])==3){# triangular
					left<- tmp[[k]][1,1]
					centre<- tmp[[k]][2,1]
					right<- tmp[[k]][3,1]

					if(dx[j] > left & dx[j] <= centre) sumfx.j[k]<- 1-(centre - dx[j])/(centre-left)
						else
					if(dx[j] > centre & dx[j] < right) sumfx.j[k]<- (right - dx[j])/(right-centre)
						else 	sumfx.j[k]<- 0
	
						sumfx.j[k]<- min(sumfx.j[k],mu.x[k])
					}#end triangular
					else if(nrow(tmp[[k]])==4)
							{# trapezoid
							left<- tmp[[k]][1,1]
							centre1<- tmp[[k]][2,1]
							centre2<- tmp[[k]][3,1]
							right<- tmp[[k]][4,1]

					if(dx[j] > left & dx[j] < centre1) sumfx.j[k]<- 1-(centre1 - dx[j])/(centre1-left)
						else
					if(dx[j] > centre2 & dx[j] < right) sumfx.j[k]<- (right - dx[j])/(right-centre2)
						else
					if(dx[j] >=centre1 & dx[j] <= centre2) sumfx.j[k]<- 1
						else 	sumfx.j[k]<- 0
					}# end trapezoid
					
						sumfx.j[k]<- min(sumfx.j[k],mu.x[k])
				}		
					
			
			} # k loop
		
		sumfx[j]<- max(sumfx.j)
		
		} #j loop
		
out[i]<- sum(sumfx*dx)/sum(sumfx)

	} # i loop

	return(cbind(fs,out))
}

#========================================================

# rescaling function

rescale.mf<- function(mf,mean,min,max){
#
# rescale the membership function mf on the unit interval
# using 2 piecewise linear transformations
# so that the mean value = 0.5	

n<- length(mf)
	for(i in 1:n){
		k<- nrow(mf[[i]])
		for(j in 1:k){
		if(mf[[i]][j,1]<=mean) mf[[i]][j,1]<- (1-((mean - mf[[i]][j,1])/(mean-min)))*0.5
			else
				mf[[i]][j,1]<- 0.5 + ((1 -((max-mf[[i]][j,1])/(max-mean)))*0.5)
		}
	}

return(mf)
}


#----------------------

# rescale to unit interval

rescale<- function(dat){
	
	min.dat<- min(dat)
	max.dat<- max(dat)
	newdat<- (dat-min.dat)/(max.dat-min.dat)
	newdat
}
