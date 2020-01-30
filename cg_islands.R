library(HMM)
library(markovchain)
library(data.table)

create_simulated_data<- function(){
   
   statesNames= c("A", "C", "G", "T")
   
####### Transition Probabilities for CpG Islands
##                A	      C	      G	      T
##          A	0.180	   0.274	   0.426	   0.120
##          C	0.171	   0.368	   0.274	   0.188
##          G	0.161	   0.339	   0.375	   0.125
##          T	0.079	   0.355	   0.384	   0.182
##     station 0.155	   0.341	   0.350	   0.154
   
   cg_islands_mc = new(
      "markovchain", 
      transitionMatrix=matrix(
         c(0.180, 0.274, 0.426, 0.120, 
           0.171, 0.368, 0.274, 0.187,
           0.161, 0.339, 0.375, 0.125,
           0.079, 0.355, 0.384, 0.182
         ),
         byrow=T,
         nrow=4, 
         dimnames=list(statesNames,statesNames)
      )
   )
   
######  Transition Probabilities for non CpG Islands
   
##               A	      C	      G	      T
##          A	0.300	   0.205	   0.285	   0.210
##          C	0.322	   0.298	   0.078	   0.302
##          G	0.248	   0.246	   0.298	   0.208
##          T	0.177	   0.239	   0.292	   0.292
##    station	0.262	   0.246	   0.239	   0.253
   
   non_cg_islands_mc = new(
      "markovchain", 
      transitionMatrix=matrix(
         c(0.300, 0.205, 0.285, 0.210, 
           0.322, 0.298, 0.078, 0.302,
           0.248, 0.246, 0.298, 0.208,
           0.177, 0.239, 0.292, 0.292
         ),
         byrow=T,
         nrow=4, 
         dimnames=list(statesNames,statesNames)
      )
   )
   
   
   
   ### Generating training set and test set.
   training_set= list(seq= c(), tag=c())
   test_set= list(seq= c(), tag=c())
   
   non_cg_len= 300:1000 ## Length of non-cg islands is usually larger than CG
   cg_len= 200:500
   N <- 50
   
   
   for(i in 1:N){
      
      if(runif(1)<0.6){ ## Probablity of CG island in the sequence is 0.4
         seq= markovchainSequence(
            n=sample(non_cg_len, 1), 
            markovchain=non_cg_islands_mc, 
            include=TRUE
         )
         tag= "Non_CG"
      }
      else{
         seq= markovchainSequence(
            n=sample(cg_len, 1), 
            markovchain=cg_islands_mc, 
            include=TRUE
         )
         tag= "CG"
      }
      
      ## Creating Training Set
      training_set$seq=  c(training_set$seq, seq) 
      training_set$tag=  c(training_set$tag, rep(tag, length(seq)))
      
      if(runif(1.02)<0.6){ ## Probablity of CG island in the sequence is 0.4
         test_seq= markovchainSequence(
            n=sample(non_cg_len, 1), 
            markovchain=non_cg_islands_mc, 
            include=TRUE
         )
         tag= "Non_CG"
      }
      else{
         test_seq=  markovchainSequence(
            n=sample(cg_len, 1), 
            markovchain=cg_islands_mc, 
            include=TRUE
         )
         tag= "CG"
      }
      
      ## Creating Test Set
      test_set$seq=  c(test_set$seq, test_seq)
      test_set$tag=  c(test_set$tag, rep(tag, length(test_seq)))
   }
   return(list(training_set=training_set, test_set=test_set))
}

plot_traing_test_set_nuc<- function(training_set, test_set){
   
   dt= rbind(
      data.table(
         x= 1:length(training_set$seq),
         y= training_set$seq,
         tag="Training Set"
      ), 
      data.table(
         x= 1:length(test_set$seq),
         y= test_set$seq,
         tag="Test Set"
      )
   )
   
   ggplot(dt[order(x)][1:800], aes(x=x, y= y))+geom_tile(aes(fill= tag), alpha=0.8) + 
      xlab("Location")+ ylab("Nucleotide")
}

plot_training_test_set<- function(training_set, test_set){
   dt= rbind(
      data.table(
         x= 1:length(training_set$seq),
         y= ifelse(training_set$tag=="CG", 1, 0),
         tag="Training Set"
      ), 
      data.table(
         x= 1:length(test_set$seq),
         y= ifelse(test_set$tag=="CG", 0.99,0.01),
         tag="Test Set"
      )
   )
      
   ggplot(dt, aes(x=x, y= y)) + geom_line(aes(color= tag)) +
      xlab("Location")+ ylab("Is CG Island")
}

train_1st_order_hmm<- function(training_set){
   
   statesNames= c("A+","C+", "G+", "T+", "A-", "C-", "G-", "T-")
   symbols= c("A", "C", "G", "T")
   
   
   ### Parameters Estimation using Maximum Likelihood.
   trans_count= matrix(
      rep(1, 64), 
      nrow= 8,
      dimnames = list(statesNames, statesNames) 
   )
   
   for(i in 1: length(training_set$tag)-1){
      row_id=paste0(
         training_set$seq[i], 
         ifelse(training_set$tag[i]=="CG", "+", "-")
      )
      col_id=paste0(
         training_set$seq[i+1], 
         ifelse(training_set$tag[i+1]=="CG", "+", "-")
      )
      
      trans_count[row_id, col_id]= 
         trans_count[row_id, col_id] + 1
   }
   
   trans_prob= apply(
      trans_count, 
      1, 
      function(x) x/sum(x)
   )
   
   hmm= initHMM(
      States= statesNames, 
      Symbols= symbols,
      transProbs= trans_prob,
      emissionProbs = matrix(
         c(1, 0, 0, 0,
           0, 1, 0, 0, 
           0, 0, 1, 0,
           0, 0, 0, 1,
           1, 0, 0, 0,
           0, 1, 0, 0, 
           0, 0, 1, 0,
           0, 0, 0, 1),
         nrow = 8,
         byrow = T
      ),
      startProbs = c(0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.1, 0.1)
   )
   hmm
}

test_hmm<- function(hmm, test_set, to=NULL){
   if(is.null(to)){
      to= length(test_set$seq)
   }
   test= HMM::viterbi(hmm = hmm, observation = test_set$seq[1:to])
   dt= data.table(
      x=1:length(test[1:to]), 
      y= test[1:to], 
      is_cg= ifelse(test[1:to] %in% c("A+", "C+", "G+", "T+"), 1, 0),
      tag="prediction"
   )
   dt= rbind(
      dt, 
      data.table(
         x=1:length(test[1:to]), 
         y=test_set$tag[1:to], 
         is_cg= ifelse(test_set$tag[1:to]=="CG", 0.99, 0.01),
         tag="truth"
      )
   )
   ggplot(dt, aes(x= x, y=is_cg)) +geom_line(aes(color= tag)) + 
      xlab("Location")+ ylab("Is CG Island")
}


run_experiment<- function(){
   l_datasets= create_simulated_data()
   plot_traing_test_set_nuc(l_datasets$training_set, l_datasets$test_set)
   plot_training_test_set(l_datasets$training_set, l_datasets$test_set)
   
   hmm= train_1st_order_hmm(l_datasets$training_set)
   
   test_hmm(hmm, l_datasets$test_set)
}
