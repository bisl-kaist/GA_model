####################################################################################################
# GA_model
# R script for random network generation model incorporating the grouped attachment (GA) to resemble
# real-world topological properties of network motifs
# 
# Reference:
# Jaejoon Choi, and Doheon Lee, "Topological motifs populate complex networks through grouped
# attachment", Under review at Scientific Reports
# 
# Contact
# - Jaejoon Choi : jjchoi@biosoft.kaist.ac.kr
# - Doheon Lee : dhlee@kaist.ac.kr
####################################################################################################

# Initialize & Load libraries
rm(list=ls())
library(igraph)

graph.extension <- function(p, maxN=0){
####################################################################################################
# Function - Random network generation based on extension
#
##### Input
# q - Double. Probability for edge generation. (0<q<1)
# maxN - Integer. Number of nodes. (maxN>0)
#
##### Output
# Double. Co-neighborness of the input graph.
####################################################################################################
  g <- graph.empty(1, directed=FALSE)
  
  while(TRUE){
    temp.g <- g
    
    if(runif(1,0,1)<p){
      temp.edge <- sample(1:(length(V(temp.g))+1), 2, replace=F)
      if(temp.edge[1]<=length(V(temp.g)) && temp.edge[2]<=length(V(temp.g))){
        while(temp.g[temp.edge[1],temp.edge[2]]==1){
          temp.edge <- sample(1:(length(V(temp.g))+1), 2, replace=F)
          if(temp.edge[1]>length(V(temp.g)) || temp.edge[2]>length(V(temp.g))){
            break;
          }
        }
      }
      
      if(temp.edge[1]>length(V(temp.g)) || temp.edge[2]>length(V(temp.g))){
        temp.g <- add.vertices(temp.g, 1)
      }
      temp.g <- temp.g + edge(V(temp.g)[temp.edge[1]],V(temp.g)[temp.edge[2]])
      
      if(maxN>0){
        if(length(V(temp.g))>maxN){
          break
        }
      }
      
      g <- temp.g
    }else{
      break
    }
  }
  
  return (g)
}

ga.game <- function(n,p,q,directed=FALSE,revised=FALSE,pref=FALSE,alpha=1,a=1){
####################################################################################################
# Function - GA random network generation
#
##### Input
# n - Integer. Number of nodes. (n>0)
# p - Double. Edge Probability (density of graph). (0<p<1)
# q - Double. Groupness probability. (p<q<1)
# directed - Logical. Whether to create a directed graph. (Default : FALSE)
# revised - Logical. Whether to use revised GA model. (Default : FALSE)
# pref - Logical. Whether to use preferential GA model. (Default : FALSE)
# alpha - Initial attractiveness of the nodes. It only works if pref=TRUE. (Default : 1)
# a - Power of the preferential attachment. It only works if pref=TRUE. (Default : 1)
#
##### Output
# An graph (igraph) object
####################################################################################################
  try(if(typeof(n)!="double" || n<=0) stop("Argument n is not an integer or is not in the range n>0"))
  try(if(typeof(p)!="double" || p<=0 || p>=1) stop("Argument p is not an integer or is not in the range 0<p<1"))
  try(if(typeof(q)!="double" || q<=p || q>=1) stop("Argument q is not an integer or is not in the range p<q<1"))
  try(if(typeof(directed)!="logical") stop("Argument directed is not logical"))
  try(if(typeof(revised)!="logical") stop("Argument revised is not logical"))
  try(if(typeof(pref)!="logical") stop("Argument pref is not logical"))
  
  p_mod <- p
  g <- graph.empty(1, directed = directed)
  
  while(length(V(g))<n){
    g1 <- graph.extension(q,maxN=(n-length(V(g))))
    
    vcount.g <- length(V(g))
    vcount.g1 <- length(V(g1))
    full_length = vcount.g + vcount.g1;
    
    g_union <- graph.disjoint.union(g,g1)
    
    if(revised==TRUE){
      if(directed==FALSE){
        max.e.g1 <- vcount.g1*(vcount.g1-1)/2
        max.e.gtog1 <- vcount.g1*vcount.g
      }else{
        max.e.g1 <- vcount.g1*(vcount.g1-1)
        max.e.gtog1 <- vcount.g1*vcount.g*2
      }
      
      p_mod <- (p*(max.e.g1+max.e.gtog1)-ecount(g1))/max.e.gtog1
      if(p_mod<0){
        p_mod<-0
      }
    }
    
    if(pref==FALSE){
      for (i in 1:vcount.g){
        if(runif(1,0,1)<(p_mod/q)){
          for(j in (vcount.g+1):full_length){
            if(runif(1,0,1)<q){
              g_union <- g_union + edge(V(g_union)[i],V(g_union)[j])
            }
            if(directed==TRUE){
              if(runif(1,0,1)<q){
                g_union <- g_union + edge(V(g_union)[j],V(g_union)[i])
              }
            }
          }
        }
      }
    }else{
      sel<-round(vcount.g*p_mod/q)
      if(sel==0){
        sel <- 1
      }
      
      selected <- sample(1:vcount.g, sel, replace=FALSE, prob=(degree(g)^alpha+a))
      for(i in selected){
        for(j in (vcount.g+1):full_length){
          if(runif(1,0,1)<q){
            g_union <- g_union + edge(V(g_union)[i],V(g_union)[j])
          }
          if(directed==TRUE){
            if(runif(1,0,1)<q){
              g_union <- g_union + edge(V(g_union)[j],V(g_union)[i])
            }
          }
        }
      }
    }
    
    g<-g_union
  }
  
  return (g)
}

co_neighbor <- function(g){
####################################################################################################
# Function - Co-neighborness
#
##### Input
# g - An graph (igraph) object.
#
##### Output
# Double. Co-neighborness of the input graph.
####################################################################################################
  jindex<-similarity.jaccard(g)
  Emat<-as_adjacency_matrix(g,type="upper")
  mult<-jindex*Emat
  
  ret <- sum( mult[mult > 0] ) / ecount(g)
}
