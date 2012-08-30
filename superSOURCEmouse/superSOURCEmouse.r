#working on genotype using additive allele model

#increase selection pressure
#bigger population
#rnorm for traits

simulateEvolution<-function(n=1000,generations=1000,genotype=T,genotypeSize=200,maxTimeSteps=1000,loaded=NULL) {
  
  #setup output
  output<-list()
  
  if(!is.null(loaded)){
    output<-loaded
  }
  length(output)=generations
  
  ifelse(is.null(loaded),1,length(loaded)+1)->start
  
  for(generation in start:generations) {
    currentGenerationSurvive<-list()
    currentGenerationDie<-list()
    length(currentGenerationSurvive)=n
    length(currentGenerationDie)=n
    
    k=1 #alive index
    j=1 #dead index
    dRange=NULL
    random=NULL
    genotypedRange=NULL
    genotyperandom=NULL
    sameParents=F #initialized as false for each generation for clutch
    currentNumOfOffspring=0 # counter for number of individuals
    clutchSize=0
    
    #for(i in 1:n){  
    while(k<=n){
      #this n is equal to the population size of each generation.
      #For fitness affecting clutch size, this loop needs to be transformed into a while loop
      
      if(generation==1){  
        if(genotype==F){
          dRange<-runif(1,max=2)
          random<-runif(1)
        }
        else
        { #using genotype
          sample(4,size=genotypeSize,replace=T)-1->genotypedRange
          sample(4,size=genotypeSize,replace=T)-1->genotyperandom
          
          dRange=2*sum(genotypedRange)/(genotypeSize*3)
          random=sum(genotyperandom)/(genotypeSize*3)
          
          #determine sex
          sex=sample(c("male","female"),1)
        }
      }
      else {#after generation 1, for each individual
        
        #select which mice will breed
          #give weights for each mouse based  on time to cheese
        
          ###the time to cheese will also indicate clutch size
        if(sameParents==F){
          getBreedingMice(generation,output,maxTimeSteps)->breedingMice #dataframe of male and female pairs(indices only)
          print(breedingMice)
          round(mean((maxTimeSteps-getParameters(output,generation,8)[breedingMice])/maxTimeSteps)*maxClutchSize)->clutchSize
        }

        currentNumOfOffspring=currentNumOfOffspring+1
        if(currentNumOfOffspring==clutchSize){
          sameParents=F
        }
        #simulate meiosis and sexual reproduction
        ###include clutch size here and put "reproduce" in for loop.
        reproduce(output,breedingMice,generation)->newGenotypes
        print(newGenotypes)
        #introduce mutation
        
        #continue with mouse trials
        dRange=2*sum(newGenotypes[1,])/(genotypeSize*3)
        random=sum(newGenotypes[2,])/(genotypeSize*3)
        genotypedRange=newGenotypes[1,]
        genotyperandom=newGenotypes[2,]
        
        #determine sex
        sex=sample(c("male","female"),1)
      }
      
    
      mouseTrial(detectionRange=dRange,randomness=random,plotResult=F,genotype=genotype,genotypedRange=genotypedRange,genotyperandom=genotyperandom,sex=sex)->tempResult
      if(tempResult[[1]]==T){
        currentGenerationSurvive[[k]]=tempResult
        k=k+1
      }
      else {
        currentGenerationDie[[j]]=tempResult
        j=j+1
      }
    }
    ifelse(is.null(currentGenerationSurvive[[k]]),length(currentGenerationSurvive)<-(k-1),length(currentGenerationSurvive)<-k)
    ifelse(is.null(currentGenerationDie[[j]]),length(currentGenerationDie)<-(j-1),length(currentGenerationDie)<-j)
  
    output[[generation]]$survive=currentGenerationSurvive
    output[[generation]]$die=currentGenerationDie
  }
  output
}


simulateFirstGeneration<-function(n=1000) {
  firstGenerationSurvive<-list()
  firstGenerationDie<-list()
  length(firstGenerationSurvive)=n
  length(firstGenerationDie)=n
  
  k=1
  j=1
  for(i in 1:n){
    dRange<-runif(1,max=2)
    random<-runif(1)
    
    mouseTrial(detectionRange=dRange,randomness=random,plotResult=F,maxTimeSteps=maxTimeSteps)->tempResult
    if(tempResult[[1]]==T){
      firstGenerationSurvive[[k]]=tempResult
      k=k+1
    }
    else {
      firstGenerationDie[[j]]=tempResult
      j=j+1
    }
  }
  ifelse(is.null(firstGenerationSurvive[[k]]),length(firstGenerationSurvive)<-(k-1),length(firstGenerationSurvive)<-k)
  ifelse(is.null(firstGenerationDie[[j]]),length(firstGenerationDie)<-(j-1),length(firstGenerationDie)<-j)
  
  output<-list()
  output$survive=firstGenerationSurvive
  output$die=firstGenerationDie
  output
  }

  
mouseTrial<-function(arenaSize=5,detectionRange=2,randomness=0.8,foundFood=0.05,maxTimeSteps=1000,plotResult=T,genotype=F,genotypedRange=NULL,genotyperandom=NULL,sex=NULL) {
  
  #setup simulation
  mouseStart<-c(0,0)
  foodStart<-c(0,0)
  while(dist(mouseStart[1],mouseStart[2],foodStart[1],foodStart[2])<detectionRange){
    mouseStart<-runif(2,max=arenaSize)
    foodStart <-runif(2,max=arenaSize)
  }
  
  ####faster if the matrix is defined ahead of time
  ####could be intiialized to 0s and then replaced in place
  
  
  mousePosition=data.frame(x=mouseStart[1],y=mouseStart[2])
  mousePosition<-rbind(mousePosition,mousePosition[1,]+runif(2,min=-.1,max=.1))
  j=2
  
  while(T) {
    currentMousePosition<-mousePosition[j,]
    lastMousePosition<-mousePosition[j-1,]
    
    #if mouse sees food
    if(dist(currentMousePosition[1],currentMousePosition[2],foodStart[1],foodStart[2])<detectionRange) {
      cat("Mouse sees food!\t step: ",j,"\n")
      moveTowardsFood(currentMousePosition,foodStart)->nextLocation
    }
    #if no, determine chance of next random position
    else if(runif(1)<randomness) {
      runif(2,min=-0.5,max=0.5)+currentMousePosition->nextLocation
    }
    #if no, move forward
    else {
      moveForward(currentMousePosition,lastMousePosition)->nextLocation
    }
    
    #adjust nextLocation to check for walls
    if(nextLocation[1]>arenaSize)
      nextLocation[1]=arenaSize
    if(nextLocation[1]<0)
      nextLocation[1]=0
    if(nextLocation[2]>arenaSize)
      nextLocation[2]=arenaSize
    if(nextLocation[2]<0)
      nextLocation[2]=0
    
    #write out result
    mousePosition=rbind(mousePosition,c(nextLocation[[1]],nextLocation[[2]]))
    j=j+1
    
    #breaking conditions
    #find food
    if(dist(nextLocation[[1]],nextLocation[[2]],foodStart[1],foodStart[2])<foundFood){
      cat("\nYour mouse has found the cheese ball!\t step: ",j,"\n")
      procreate=T
      break
    }
    if(j==maxTimeSteps){
      procreate=F
      cat("\nYour mouse died of starvation!\n")
      break
    }
  }
  if(plotResult){
    plot(x=mousePosition[,1],y=mousePosition[,2],type="l",xlim=c(0,arenaSize),ylim=c(0,arenaSize))
    points(x=foodStart[1],y=foodStart[2],col="red",cex=3)
    points(x=mouseStart[1],y=mouseStart[2],col="green",cex=3)
  }
  
  #return(list(procreate,detectionRange,randomness,j))
  output=list()
  output$foundCheese=procreate
  output$mousePositionx=mouseStart[1]
  output$mousePositiony=mouseStart[2]
  output$foodPositionx=foodStart[1]
  output$foodPositiony=foodStart[2]
  output$detectionRange=detectionRange
  output$randomness=randomness
  output$timeToCheese=j
  if(genotype==T){
    output$detectionRangeGenotype=genotypedRange
    output$randomnessGenotype=genotyperandom
    output$sex=sex
  }
  
  output
}




#other functions
dist<-function(x1,y1,x2,y2) {
  sqrt((x2-x1)^2+(y2-y1)^2)
}

moveTowardsFood<-function(mousePosition,foodPosition){ #->nextLocation
  #xCoordinate
  if(mousePosition[1]>foodPosition[1]){
    nextx<--runif(1)*(mousePosition[1]-foodPosition[1])
  }
  else {
    nextx<-runif(1)*(foodPosition[1]-mousePosition[1])
  }
  if(mousePosition[2]>foodPosition[2]){
    nexty<--runif(1)*(mousePosition[2]-foodPosition[2])
  }
  else {
    nexty<-runif(1)*(foodPosition[2]-mousePosition[2])
  }
  c(nextx+mousePosition[1],nexty+mousePosition[2])   
}

moveForward<-function(mousePosition,lastPosition) {#->nextLocation
  if(mousePosition[1]>lastPosition[1]){
    nextx<-runif(1)*(mousePosition[1]+lastPosition[1])
  }
  else {
    nextx<--runif(1)*(lastPosition[1]-mousePosition[1])
  }
  if(mousePosition[2]>lastPosition[2]){
    nexty<-runif(1)*(mousePosition[2]+lastPosition[2])
  }
  else {
    nexty<--runif(1)*(lastPosition[2]-mousePosition[2])
  }
  c(nextx+mousePosition[1],nexty+mousePosition[2])

}

getParameters<-function(input,prevGeneration,parameterNum){
  output<-vector(,length(input[[prevGeneration]]$survive))
  for(i in 1:length(output)){
    output[i]=input[[prevGeneration]]$survive[[i]][[parameterNum]]
  }
  output
}

distancePlot<-function(input,prevGen){
  mousex=getParameters(input=input,prevGeneration=prevGen,parameterNum=2)
  mousey=getParameters(input=input,prevGeneration=prevGen,parameterNum=3)
  foodx=getParameters(input=input,prevGeneration=prevGen,parameterNum=4)
  foody=getParameters(input=input,prevGeneration=prevGen,parameterNum=5)
  plot(0:5,0:5,type="n")
  
  for(i in 1:length(input[[prevGen]]$survive))
    lines(x=c(mousex[i],foodx[i]),y=c(mousey[i],foody[i]))
  }

getBreedingMice<-function(generation,output,maxTimeSteps=maxTimeSteps){ #->breedingMice
  which.males<-getParameters(output,prevGeneration=generation-1,"sex")=="male"
  which.females<-!which.males
  
  weights<-(maxTimeSteps-(getParameters(output,prevGeneration=generation-1,"timeToCheese")))/maxTimeSteps
  weights.males=weights
  weights.females=weights
  weights.males[!which.males]=0
  weights.females[!which.females]=0
  
  maleIndices<-sample(1:length(weights.males),1,replace=T,prob=weights.males)
  femaleIndices<-sample(1:length(weights.females),1,replace=T,prob=weights.females)
  
  data.frame(maleIndices=maleIndices,femaleIndices=femaleIndices)
}

reproduce<-function(output,breedingMice,generation){#newGenotypes
  maleDetectionRangeGenotype<<-output[[generation-1]][[1]][[breedingMice[1,1]]]["detectionRangeGenotype"][[1]]
  maleRandomnessGenotype=output[[generation-1]][[1]][[breedingMice[1,1]]]["randomnessGenotype"][[1]]
  
  #print(maleDetectionRangeGenotype)
  #print(maleRandomnessGenotype)
  
  femaleDetectionRangeGenotype=output[[generation-1]][[1]][[breedingMice[1,2]]]["detectionRangeGenotype"][[1]]
  femaleRandomnessGenotype=output[[generation-1]][[1]][[breedingMice[1,2]]]["randomnessGenotype"][[1]]
  
  matrix(data=0,nrow=2,ncol=length(output[[generation-1]][[1]][[1]]["randomnessGenotype"][[1]]))->genotypeOutput
  
  runif(2)->randomSelectionOfGenes #random numbers to be tested against 0.5 to determine whether child gets 1st half from mother or father
  halfway=ncol(genotypeOutput)/2
  
  if(randomSelectionOfGenes[1]>=0.5){
    print(maleDetectionRangeGenotype)
    print(femaleDetectionRangeGenotype)
    print(c(maleDetectionRangeGenotype[1:halfway],femaleDetectionRangeGenotype[(halfway+1):ncol(genotypeOutput)]))
    genotypeOutput[1,]=c(maleDetectionRangeGenotype[1:halfway],femaleDetectionRangeGenotype[(halfway+1):ncol(genotypeOutput)])
  }
  else {
    genotypeOutput[1,]=c(femaleDetectionRangeGenotype[1:halfway],maleDetectionRangeGenotype[(halfway+1):ncol(genotypeOutput)])
  }
  
  if(randomSelectionOfGenes[2]>=0.5){
    genotypeOutput[2,]=c(maleRandomnessGenotype[1:halfway],femaleRandomnessGenotype[(halfway+1):ncol(genotypeOutput)])
  }
  else {
    genotypeOutput[2,]=c(femaleRandomnessGenotype[1:halfway],maleRandomnessGenotype[(halfway+1):ncol(genotypeOutput)])
  }
  genotypeOutput
}

plot2parameters<-function(input){
  print(names(input[[1]][[1]][[1]]))
  as.numeric(readline("What is the index of the first parameter?\n"))->param1
  as.numeric(readline("What is the index of the second parameter?\n"))->param2
  
  maxGeneration<-length(input)
  
  start1<-getParameters(input=input,prevGeneration=1,parameterNum=param1)
  start2<-getParameters(input=input,prevGeneration=1,parameterNum=param2)
  
  end1<-getParameters(input=input,prevGeneration=maxGeneration,parameterNum=param1)
  end2<-getParameters(input=input,prevGeneration=maxGeneration,parameterNum=param2)
  
  plot(x=c(start1,end1),y=c(start2,end2),type="n",main=paste(maxGeneration,"generations",sep=" "),
       xlab=names(input[[1]][[1]][[1]])[param1],
       ylab=names(input[[1]][[1]][[1]])[param2])
  
  points(x=start1,y=start2,pch=20)
  points(x=end1,y=end2,pch=20,col="blue")
  
}
