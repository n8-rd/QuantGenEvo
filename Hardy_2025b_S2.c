//A relatively simple model of how GRN architecture affects
//evolvability. This version with a Boolean network.
//This version works a little harder for haploidy.
//And adds measures of phenotypic diversity

//function to measure genotypic sensitivity to mutation
//for now this is only messing with W, should add messing with z
function (float)mutSens(o ind){
  senses = float(length=0);
  //get the individuals GRN and phenotypic states
  pheno = ind.getValue("pheno");
  muts1 = ind.genomes.mutationsOfType(m1);
  muts2 = ind.genomes.mutationsOfType(m2);
  muts3 = ind.genomes.mutationsOfType(m3);
  muts4 = ind.genomes.mutationsOfType(m4);
  muts5 = ind.genomes.mutationsOfType(m5);
  muts6 = ind.genomes.mutationsOfType(m6);
  lu0 = size(muts1) ? muts1[0].getValue("LUout") else rep(0, nLU);
  lu1 = size(muts2) ? muts2[0].getValue("LUout") else rep(0, nLU);
  lu2 = size(muts3) ? muts3[0].getValue("LUout") else rep(0, nLU);
  luM = matrix(c(lu0,lu1,lu2),ncol=3,byrow=F);
  w0 = size(muts4) ? ind.sumOfMutationsOfType(m4) else 0.0;
  w1 = size(muts5) ? ind.sumOfMutationsOfType(m5) else 0.0;
  w2 = size(muts6) ? ind.sumOfMutationsOfType(m6) else 0.0;
  ws = c(w0,w1,w2);
  //do ten iterations of sensitivity testing
  for (i in 1:10){
    //mutation something
    coin = rbinom(1,1,0.5);
    if (coin == 0){//mutate a reg gene
      newLUout = rbinom(nLU, 1, 0.5);
      gene = sample(seq(0,2),1);
      luM[,gene] = newLUout;
    } else {//mutate a z weight
      newEff = rnorm(1, 0.0, sigma);
      weight = sample(seq(0,2),1);
      ws[weight] = ws[weight] + newEff;
    }
    states = rep(0,N);
    for (j in 1:10){
      for (s in 0:2){
        LUTstring = asString(states[0]) + asString(states[1]) + asString(states[2]);
        lui = LUi.getValue(LUTstring);
        state = luM[lui,s];
        states[s] = state;
      }
    }
    //now use the weight terms to calculate the phenotype
    AltPheno = 0.0 + states[0]*w0 + states[1]*w1 + states[2]*w2;
    delta = abs(pheno-AltPheno);
    senses = c(senses, delta);
  }
 msen = quantile(senses, 0.5);
 return(msen);
}

initialize() {
  defineConstant("K", 200); // carrying capacity 
  defineConstant("O0", 2.0); //local optimum 0 [3.0]
  defineConstant("O1", 0.0); // local optimum 1 [0.0]
  defineConstant("OG", 4.0); // global optimum displacement [5.0]
  defineConstant("tau", tu); //probability of environment shift before tx [0.1][tu]
  defineConstant("tx", 500); //generation at which global optimum goes into effect
  defineConstant("omega", 1.0); //weakness of selection [1.0, 0.3]
  defineConstant("sigma", 0.2); // sd of z weight allele effect distro [0.15]
  defineConstant("Rep", rep); //[rep]
  defineConstant("N", 3); //number of genes in the GRN layer
  LUis = Dictionary("000", 0, "001", 1, "010", 2, "011", 3, "100", 4, "101", 5, "110", 6, "111", 7);
  defineConstant("nLU", 8);
  defineConstant("LUi", LUis);
  MutRate = 0.02; //[0.01]
  initializeMutationRate(MutRate);
  initializeMutationType("m1", 0.5, "n", 0.0, 0.0); // regulatory genes
  initializeMutationType("m2", 0.5, "n", 0.0, 0.0); 
  initializeMutationType("m3", 0.5, "n", 0.0, 0.0);
  initializeMutationType("m4", 0.5, "n", 0.0, sigma); // z weights
  initializeMutationType("m5", 0.5, "n", 0.0, sigma);
  initializeMutationType("m6", 0.5, "n", 0.0, sigma);
  c(m1,m2,m3,m4,m5,m6).convertToSubstitution = F;
  c(m1,m2,m3).mutationStackPolicy = 'l';
  c(m4,m5,m6).mutationStackPolicy = 's';
  initializeGenomicElementType("g1", m1, 1.0);
  initializeGenomicElementType("g2", m2, 1.0);
  initializeGenomicElementType("g3", m3, 1.0);
  initializeGenomicElementType("g4", m4, 1.0);
  initializeGenomicElementType("g5", m5, 1.0);
  initializeGenomicElementType("g6", m6, 1.0);
  initializeGenomicElement(g1, 1, 1);
  initializeGenomicElement(g2, 2, 2);
  initializeGenomicElement(g3, 3, 3);
  initializeGenomicElement(g4, 4, 4);
  initializeGenomicElement(g5, 5, 5);
  initializeGenomicElement(g6, 6, 6);
  initializeRecombinationRate(0.0);
}

//make sure to turn off any built-in allele fitness effects
mutationEffect(m1) { return 1.0; }
mutationEffect(m2) { return 1.0; }
mutationEffect(m3) { return 1.0; }
mutationEffect(m4) { return 1.0; }
mutationEffect(m5) { return 1.0; }
mutationEffect(m6) { return 1.0; }

//regulatory mutation effect generator (for regulatory genes)
mutation(m1){
    newLUout = rbinom(nLU, 1, 0.5);
    mut.setValue("LUout", newLUout);
    return T;
}
mutation(m2){
    newLUout = rbinom(nLU, 1, 0.5);
    mut.setValue("LUout", newLUout);
    return T;
}
mutation(m3){
    newLUout = rbinom(nLU, 1, 0.5);
    mut.setValue("LUout", newLUout);
    return T;
}

//set up unstructured population
1 first() {
  sim.addSubpop("p1", K);
  inds = p1.individuals;
  //set initial environmental state
  sim.setValue("Es", 0);
  //cloning
  p1.setCloningRate(1.0);
  //with corresponding z opt
  //sim.setValue("zo", O0);
}

//let's determine the regulatory and phenotype states of each individual at birth
modifyChild(){
  //deterome the LUT outputs for each reg gene and the phenotype weights
  muts1 = child.genomes.mutationsOfType(m1);
  muts2 = child.genomes.mutationsOfType(m2);
  muts3 = child.genomes.mutationsOfType(m3);
  muts4 = child.genomes.mutationsOfType(m4);
  muts5 = child.genomes.mutationsOfType(m5);
  muts6 = child.genomes.mutationsOfType(m6);
  lu0 = size(muts1) ? muts1[0].getValue("LUout") else rep(0, nLU);
  lu1 = size(muts2) ? muts2[0].getValue("LUout") else rep(0, nLU);
  lu2 = size(muts3) ? muts3[0].getValue("LUout") else rep(0, nLU);
  //catn(nLU);
  luM = matrix(c(lu0,lu1,lu2),ncol=3,byrow=F);
  w0 = size(muts4) ? child.sumOfMutationsOfType(m4) else 0.0;
  w1 = size(muts5) ? child.sumOfMutationsOfType(m5) else 0.0;
  w2 = size(muts6) ? child.sumOfMutationsOfType(m6) else 0.0;
  //determine the activation state of each reg gene
  //the tricky part here is that all states are interdependant.
  //to get around that, we'll assign each gene a starting value, 
  //and then undergo ten cycles of interaction
  states = rep(0,N);
  for (i in 1:10){
    for (s in 0:2){
      LUTstring = asString(states[0]) + asString(states[1]) + asString(states[2]);
      lui = LUi.getValue(LUTstring);
      state = luM[lui,s];
      states[s] = state;
    }
  }
  //now use the weight terms to calculate the phenotype
  pheno = 0.0 + states[0]*w0 + states[1]*w1 + states[2]*w2;
  //and make an index of the state of the reg network
  excitation = sum(c(lu0, lu1, lu2));
  child.setValue("exi", excitation);
  child.setValue("pheno", pheno);
  return T;
}

//every generation, set the environment and apply selection
1:20000 late() {
  // remove any new mutations added to the disabled diploid genomes
  sim.subpopulations.individuals.genome2.removeMutations();
  inds = p1.individuals;
  ecoin = rbinom(1,1,tau);
  if (ecoin == 1){//switch the environmental state
      ES = sim.getValue("Es");
      newES = abs(ES-1);
      sim.setValue("Es", newES);
  }
  ES = sim.getValue("Es");
  OO = O0;
  if (ES == 1){
      OO = O1;
  }

  inds = sim.subpopulations.individuals;
  WarmUpScale = dnorm(0.0, 0.0, 10.0);
  scale = dnorm(0.0, 0.0, omega);
  for (ind in inds){
    y = dnorm(OO, ind.getValue("pheno"), omega)/(1.5*scale);
    ywm = dnorm(OO, ind.getValue("pheno"), 10.0)/WarmUpScale;//[1.2]
    z = dnorm(OO+OG, ind.getValue("pheno"), omega)/scale;
    if (sim.cycle < 100){
	    ind.fitnessScaling = ywm;
    } else if (sim.cycle >= 100 & sim.cycle < tx){
	    ind.fitnessScaling = y;
    } else {
	    ind.fitnessScaling = y + z;
    }
  }
  inds1 = p1.individuals;
}

late(){
  if (sim.cycle % 10 == 0 ){
      inds = p1.individuals;
      mv = size(inds) ? mean(inds.fitnessScaling) else 0.0;
      mp = size(inds) ? mean(inds.getValue("pheno")) else 0.0;
      vp = size(inds) ? sd(inds.getValue("pheno")) else 0.0;
      me = size(inds) ? mean(inds.getValue("exi")) else 0.0;
      //assay mutational sensitivity
      subjects = sample(inds, 10);
      ms = float(length=0);
      for (j in subjects){
        msj = mutSens(j);
        ms = c(ms, msj);
      }
      //measure the genetic diversity
      H = calcHeterozygosity(sim.subpopulations.individuals.genome1);
      mms = mean(ms);
      catn(sim.cycle + ": mv = " + mv + " mp = " + mp + " ms = " + mms + " H = " + H + " exi = " + me);
      writeFile("outs.csv", paste(c(Rep,sim.cycle,tau,mv,mp,mms,H,vp),sep=","), append=T);
      if (sim.cycle > tx){
	if (mv > 0.85){
	    sim.simulationFinished();
	    catn("Back on Track!");
	}
      }
  }
}
