#rm(list = ls())

TransferToMutationProfile <- function(SM, CM, SMRef = matrix(0, dim(SM)[1], dim(SM)[2]), CMRef = matrix(2, dim(CM)[1], dim(CM)[2]))
{
  rownames(SM) = paste("SNV", 1:dim(SM)[1], sep = "_");
  rownames(CM) = paste("CN", 1:dim(CM)[1], sep = "_");
  mep = vector("list", length = dim(SM)[2]);
  names(mep) = paste("subC", 1:length(mep), sep = "_");
  
  for (j in 1:dim(SM)[2])
  {
    for (i in 1:dim(SM)[1])
    {
      baseij = SMRef[i, j];
      obsij = SM[i, j];
      if (baseij != obsij)
      {
        vij = baseij:obsij;
        mep[[j]] = c(mep[[j]], paste(rownames(SM)[i], ":", vij[1:(length(vij)-1)], "->", vij[2:length(vij)], sep = ""));
      }
    }
  }
  
  for (j in 1:dim(CM)[2])
  {
    for (i in 1:dim(CM)[1])
    {
      baseij = CMRef[i, j];
      obsij = CM[i, j];
      if (baseij != obsij)
      {
        vij = baseij:obsij;
        mep[[j]] = c(mep[[j]], paste(rownames(CM)[i], ":", vij[1:(length(vij)-1)], "->", vij[2:length(vij)], sep = ""));
      }      
    }
  }
  
  Result = list(flag = TRUE, mutEvePro = mep);
  return(Result);
}


testEqual <- function(x, y)
{
  if ((length(x) == 0) && (length(y) == 0))
  {
    return(1);
  }
  return(length(intersect(x, y))/max(length(x), length(y)));
}

# This function constructs the initial tree based on two genomes
# ID1 is the index of the first genome.
# ID2 is the index of the second genome.
# mep is the mutation profile of all genome
generateInitialTree <- function(mep, ID1 = 1, ID2 = 2)
{
  tree = list(Nodes = c(list(subC_0 = c()), mep), Edges = list(), SubCinTree = c()); 
  # Tree is a list encoding the tree structure. Nodes is list of mutation sets including all nodes in the tree. 
  # Edges is a list encoding the branches in the tree. Parent is the first node in the edge and child is the second node in the edge.

  treeOri = tree
  treeV = c()
  errV = c()
  
  # Case 1
  tree$Edges = list(mep[[ID1]], mep[[ID2]]);
  names(tree$Edges)[1:2] = c(paste("subC_0", names(mep)[ID1], sep = "->"), 
                             paste("subC_0", names(mep)[ID2], sep = "->"));
  tree$negEdges = list(c(), c())
  names(tree$negEdges)[1:2] = c(paste(names(mep)[ID1], "subC_0", sep = "->"), 
                             paste(names(mep)[ID2], "subC_0", sep = "->"));
  tree$SubCinTree = names(mep)[c(ID1, ID2)];
  treeV[[1]] = tree
  errV[1] = errOfTree(tree)
  
  # Case 2
  tree = treeOri
  tree$Edges = list(mep[[ID1]], setdiff(mep[[ID2]], mep[[ID1]]));
  names(tree$Edges)[1:2] = c(paste("subC_0", names(mep)[ID1], sep = "->"), 
                             paste(names(mep)[ID1], names(mep)[ID2], sep = "->"));
  tree$negEdges = list(c(), setdiff(mep[[ID1]], mep[[ID2]]))
  names(tree$negEdges)[1:2] = c(paste(names(mep)[ID1], "subC_0", sep = "->"), 
                                paste(names(mep)[ID2], names(mep)[ID1], sep = "->"));
  tree$SubCinTree = names(mep)[c(ID1, ID2)];
  treeV[[2]] = tree
  errV[2] = errOfTree(tree)

  # Case 3
  tree = treeOri
  tree$Edges = list(mep[[ID2]], setdiff(mep[[ID1]], mep[[ID2]]));
  names(tree$Edges)[1:2] = c(paste("subC_0", names(mep)[ID2], sep = "->"), 
                             paste(names(mep)[ID2], names(mep)[ID1], sep = "->")); 
  tree$negEdges = list(c(), setdiff(mep[[ID2]], mep[[ID1]]))
  names(tree$negEdges)[1:2] = c(paste(names(mep)[ID2], "subC_0", sep = "->"), 
                                paste(names(mep)[ID1], names(mep)[ID2], sep = "->"));
  tree$SubCinTree = names(mep)[c(ID2, ID1)];  
  treeV[[3]] = tree
  errV[3] = errOfTree(tree)

  # Case 4
  if (length(intersect(mep[[ID1]], mep[[ID2]])) == 0)
  {
  } else{
    tree = treeOri
    nodeNotInTree = setdiff(names(tree$Nodes), c("subC_0", names(mep)[c(ID1, ID2)]));
    if (length(nodeNotInTree) == 0)
    {
      tree = treeOri
      tree$Nodes = c(tree$Nodes, list(intersect(mep[[ID1]], mep[[ID2]])));  
      IDk = length(tree$Nodes);
      names(tree$Nodes)[IDk] = paste("subC_", IDk-1, sep = "");
      
      tree$Edges = list(tree$Nodes[[IDk]], setdiff(mep[[ID1]], tree$Nodes[[IDk]]), 
                        setdiff(mep[[ID2]], tree$Nodes[[IDk]]));
      names(tree$Edges)[1:3] = c(paste("subC_0", names(tree$Nodes)[IDk], sep = "->"), 
                                 paste(names(tree$Nodes)[IDk], names(mep)[ID1], sep = "->"), 
                                 paste(names(tree$Nodes)[IDk], names(mep)[ID2], sep = "->"));
      tree$negEdges = list(c(), setdiff(tree$Nodes[[IDk]], mep[[ID1]]), 
                           setdiff(tree$Nodes[[IDk]], mep[[ID2]]));
      names(tree$negEdges)[1:3] = c(paste(names(tree$Nodes)[IDk], "subC_0", sep = "->"), 
                                    paste(names(mep)[ID1], names(tree$Nodes)[IDk], sep = "->"), 
                                    paste(names(mep)[ID2], names(tree$Nodes)[IDk], sep = "->"));
      tree$SubCinTree = c(names(tree$Nodes)[IDk], names(mep)[c(ID1, ID2)]);
      treeV[[4]] = tree
      errV[4] = errOfTree(tree)      
    } else {

      tree = treeOri
      equalScore = sapply(X = tree$Nodes[nodeNotInTree], FUN = testEqual, y = intersect(mep[[ID1]], mep[[ID2]]));
      IDk = which.max(equalScore);
      IDk = which(names(tree$Nodes) == nodeNotInTree[IDk]);
      tree$Edges = list(tree$Nodes[[IDk]], setdiff(mep[[ID1]], tree$Nodes[[IDk]]), 
                        setdiff(mep[[ID2]], tree$Nodes[[IDk]]));
      names(tree$Edges)[1:3] = c(paste("subC_0", names(tree$Nodes)[IDk], sep = "->"), 
                                 paste(names(tree$Nodes)[IDk], names(mep)[ID1], sep = "->"), 
                                 paste(names(tree$Nodes)[IDk], names(mep)[ID2], sep = "->"));
      tree$negEdges = list(c(), setdiff(tree$Nodes[[IDk]], mep[[ID1]]), 
                           setdiff(tree$Nodes[[IDk]], mep[[ID2]]));
      names(tree$negEdges)[1:3] = c(paste(names(tree$Nodes)[IDk], "subC_0", sep = "->"), 
                                    paste(names(mep)[ID1], names(tree$Nodes)[IDk], sep = "->"), 
                                    paste(names(mep)[ID2], names(tree$Nodes)[IDk], sep = "->"));
      tree$SubCinTree = c(names(tree$Nodes)[IDk], names(mep)[c(ID1, ID2)]);
      treeV[[4]] = tree
      errV[4] = errOfTree(tree)      
      
      tree = treeOri
      tree$Nodes = c(tree$Nodes, list(intersect(mep[[ID1]], mep[[ID2]])));  
      IDk = length(tree$Nodes);
      names(tree$Nodes)[IDk] = paste("subC_", IDk-1, sep = "");
      tree$Edges = list(tree$Nodes[[IDk]], setdiff(mep[[ID1]], tree$Nodes[[IDk]]), 
                        setdiff(mep[[ID2]], tree$Nodes[[IDk]]));
      names(tree$Edges)[1:3] = c(paste("subC_0", names(tree$Nodes)[IDk], sep = "->"), 
                                 paste(names(tree$Nodes)[IDk], names(mep)[ID1], sep = "->"), 
                                 paste(names(tree$Nodes)[IDk], names(mep)[ID2], sep = "->"));
      tree$negEdges = list(c(), setdiff(tree$Nodes[[IDk]], mep[[ID1]]), 
                           setdiff(tree$Nodes[[IDk]], mep[[ID2]]));
      names(tree$negEdges)[1:3] = c(paste(names(tree$Nodes)[IDk], "subC_0", sep = "->"), 
                                    paste(names(mep)[ID1], names(tree$Nodes)[IDk], sep = "->"), 
                                    paste(names(mep)[ID2], names(tree$Nodes)[IDk], sep = "->"));
      tree$SubCinTree = c(names(tree$Nodes)[IDk], names(mep)[c(ID1, ID2)]);
      treeV[[5]] = tree
      errV[5] = errOfTree(tree)      
    }
  }

  ID = which.min(errV)
  return(list(tree = treeV[[ID]], err = errV[ID]))
}


# This function needs to be modified.
errOfTree <- function(tree)
{
  edges = c()
  for (i in 1:length(tree$Edges))
  {
    edges = c(edges, tree$Edges[[i]])
  }
  drop = 0
  for (i in 1:length(tree$negEdges))
  {
    drop = drop + length(tree$negEdges[[i]])
  }
  err = drop + sum(duplicated(edges))
  return(err)
}  
  

# This function adds a genome to an existing tree. 
# tree is the exsiting tree.
# addG is the name of the genome that needs to be added to the tree. The genome is already in tree$Nodes.
# tree$Nodes already includes the genomes from recovery result 
# and the hidden genomes generated in the process of constructing tree.
AddOneGenomeToExistingTree <- function(tree, addG)
{
  # If when adding a genome a hidden genome is added, we must check whether this hidden genome is 
  # already included in tree$Nodes. If yes, we don't need to add this hidden genome.
  #  writeLines(paste("Add ", addG, sep = ""));

  # Check genome addID is not already in the tree
  if (is.element(addG, tree$SubCinTree))
  {
    writeLines("Error: genome already in the tree");
    Result = list(tree = NULL, th = NULL);
    return(Result);
  }
  
  # Check whether every genome is a child in one and only one edge
  allChild = sapply(strsplit(names(tree$Edges), split = "->"), function(x)x[2]);
  allParent = sapply(strsplit(names(tree$Edges), split = "->"), function(x)x[1]);
  if ((!setequal(allChild, tree$SubCinTree)) || (length(allChild) != length(unique(allChild))))
  {
    writeLines("Error: not all genomes in the tree are child in one and only one edge");
    Result = list(tree = NULL, th = NULL);
    return(Result);    
  }
  
  treeOri = tree
  index = 0
  treeV = c()
  errV = c()
  
  # Case I: add one genome to root
  tree$SubCinTree = c(tree$SubCinTree, addG);
  tree$Edges = c(tree$Edges, tree$Nodes[addG]);
  names(tree$Edges)[length(tree$Edges)] = paste("subC_0", addG, sep = "->");   
  tree$negEdges = c(tree$negEdges, list(c()));
  names(tree$negEdges)[length(tree$negEdges)] = paste(addG, "subC_0", sep = "->");
  index = index + 1
  treeV[[index]] = tree
  errV[index] = errOfTree(tree)
  
  # Case II: add as a leaf node to an existing node
  for (i in 1:length(allChild))
  {
    tree = treeOri
    tree$SubCinTree = c(tree$SubCinTree, addG);
    tree$Edges = c(tree$Edges, list(setdiff(tree$Nodes[[addG]], tree$Nodes[[allChild[i]]])));
    names(tree$Edges)[length(tree$Edges)] = paste(allChild[i], addG, sep = "->");    
    tree$negEdges = c(tree$negEdges, list(setdiff(tree$Nodes[[allChild[i]]], tree$Nodes[[addG]])));
    names(tree$negEdges)[length(tree$negEdges)] = paste(addG, allChild[i], sep = "->");
    index = index + 1
    treeV[[index]] = tree
    errV[index] = errOfTree(tree)    
  }

  # Case III: add as an itermediate node to an existing edge
  for (i in 1:length(allChild))
  {
    tree = treeOri
    tree$SubCinTree = c(tree$SubCinTree, addG);
    tree$Edges = c(tree$Edges, list(setdiff(tree$Nodes[[addG]], tree$Nodes[[allParent[i]]])), 
                   list(setdiff(tree$Nodes[[allChild[i]]], tree$Nodes[[addG]])));
    names(tree$Edges)[(length(tree$Edges)-1):length(tree$Edges)] = 
      c(paste(allParent[i], addG, sep = "->"), paste(addG, allChild[i], sep = "->"));
    tree$Edges[i] = c();
    tree$negEdges = c(tree$negEdges, list(setdiff(tree$Nodes[[allParent[i]]], tree$Nodes[[addG]])), 
                   list(setdiff(tree$Nodes[[addG]], tree$Nodes[[allChild[i]]])));
    names(tree$negEdges)[(length(tree$negEdges)-1):length(tree$negEdges)] = 
      c(paste(addG, allParent[i], sep = "->"), paste(allChild[i], addG, sep = "->"));
    tree$negEdges[i] = c();
    index = index + 1
    treeV[[index]] = tree
    errV[index] = errOfTree(tree)    
  }

  # Case IV: add as a leafe node with a hidden genome
  for (i in 1:length(allChild))
  {
    tree = treeOri
    
    # Need to use set equal function with tolerance.
    newNode = intersect(tree$Nodes[[addG]], tree$Nodes[[allChild[i]]]);
    
    if ((length(newNode) > 0) & (length(setdiff(tree$Nodes[[allChild[i]]], newNode)) > 0)
        & (length(setdiff(tree$Nodes[[addG]], newNode)) > 0))
    {
      nodeNotInTree = setdiff(names(tree$Nodes), c("subC_0", tree$SubCinTree, addG));
      if (length(nodeNotInTree) == 0)
      {
        tree = treeOri
        newNodeName = paste("subC_", length(tree$Nodes), sep = "");
        tree$Nodes = c(tree$Nodes, list(newNode));
        names(tree$Nodes)[length(tree$Nodes)] = newNodeName; 
        
        tree$SubCinTree = c(tree$SubCinTree, newNodeName, addG);
        tree$Edges = c(tree$Edges, list(setdiff(tree$Nodes[[newNodeName]], tree$Nodes[[allParent[i]]])), 
                       list(setdiff(tree$Nodes[[allChild[i]]], tree$Nodes[[newNodeName]])), 
                       list(setdiff(tree$Nodes[[addG]], tree$Nodes[[newNodeName]])));
        names(tree$Edges)[(length(tree$Edges)-2):length(tree$Edges)] = c(
          paste(allParent[i], "->", newNodeName, sep = ""), 
          paste(newNodeName, "->", allChild[i], sep = ""), 
          paste(newNodeName, "->", addG, sep = ""));
        tree$Edges[i] = c();
        tree$negEdges = c(tree$negEdges, list(setdiff(tree$Nodes[[allParent[i]]], tree$Nodes[[newNodeName]])), 
                          list(setdiff(tree$Nodes[[newNodeName]], tree$Nodes[[allChild[i]]])), 
                          list(setdiff(tree$Nodes[[newNodeName]], tree$Nodes[[addG]])));
        names(tree$negEdges)[(length(tree$negEdges)-2):length(tree$negEdges)] = c(
          paste(newNodeName, "->", allParent[i], sep = ""), 
          paste(allChild[i], "->", newNodeName, sep = ""), 
          paste(addG, "->", newNodeName, sep = ""));
        tree$negEdges[i] = c();
        index = index + 1
        treeV[[index]] = tree
        errV[index] = errOfTree(tree)    
        
      } else {
        tree = treeOri
        equalScore = sapply(X = tree$Nodes[nodeNotInTree], FUN = testEqual, y = newNode);
        IDk = which.max(equalScore);
        newNodeName = nodeNotInTree[IDk];        
        tree$SubCinTree = c(tree$SubCinTree, newNodeName, addG);
        tree$Edges = c(tree$Edges, list(setdiff(tree$Nodes[[newNodeName]], tree$Nodes[[allParent[i]]])), 
                       list(setdiff(tree$Nodes[[allChild[i]]], tree$Nodes[[newNodeName]])), 
                       list(setdiff(tree$Nodes[[addG]], tree$Nodes[[newNodeName]])));
        names(tree$Edges)[(length(tree$Edges)-2):length(tree$Edges)] = c(
          paste(allParent[i], "->", newNodeName, sep = ""), 
          paste(newNodeName, "->", allChild[i], sep = ""), 
          paste(newNodeName, "->", addG, sep = ""));
        tree$Edges[i] = c();
        tree$negEdges = c(tree$negEdges, list(setdiff(tree$Nodes[[allParent[i]]], tree$Nodes[[newNodeName]])), 
                          list(setdiff(tree$Nodes[[newNodeName]], tree$Nodes[[allChild[i]]])), 
                          list(setdiff(tree$Nodes[[newNodeName]], tree$Nodes[[addG]])));
        names(tree$negEdges)[(length(tree$negEdges)-2):length(tree$negEdges)] = c(
          paste(newNodeName, "->", allParent[i], sep = ""), 
          paste(allChild[i], "->", newNodeName, sep = ""), 
          paste(addG, "->", newNodeName, sep = ""));
        tree$negEdges[i] = c();
        index = index + 1
        treeV[[index]] = tree
        errV[index] = errOfTree(tree)  

        tree = treeOri
        newNodeName = paste("subC_", length(tree$Nodes), sep = "");
        tree$Nodes = c(tree$Nodes, list(newNode));
        names(tree$Nodes)[length(tree$Nodes)] = newNodeName;           
        tree$SubCinTree = c(tree$SubCinTree, newNodeName, addG);
        tree$Edges = c(tree$Edges, list(setdiff(tree$Nodes[[newNodeName]], tree$Nodes[[allParent[i]]])), 
                       list(setdiff(tree$Nodes[[allChild[i]]], tree$Nodes[[newNodeName]])), 
                       list(setdiff(tree$Nodes[[addG]], tree$Nodes[[newNodeName]])));
        names(tree$Edges)[(length(tree$Edges)-2):length(tree$Edges)] = c(
          paste(allParent[i], "->", newNodeName, sep = ""), 
          paste(newNodeName, "->", allChild[i], sep = ""), 
          paste(newNodeName, "->", addG, sep = ""));
        tree$Edges[i] = c();
        tree$negEdges = c(tree$negEdges, list(setdiff(tree$Nodes[[allParent[i]]], tree$Nodes[[newNodeName]])), 
                          list(setdiff(tree$Nodes[[newNodeName]], tree$Nodes[[allChild[i]]])), 
                          list(setdiff(tree$Nodes[[newNodeName]], tree$Nodes[[addG]])));
        names(tree$negEdges)[(length(tree$negEdges)-2):length(tree$negEdges)] = c(
          paste(newNodeName, "->", allParent[i], sep = ""), 
          paste(allChild[i], "->", newNodeName, sep = ""), 
          paste(addG, "->", newNodeName, sep = ""));
        tree$negEdges[i] = c();
        index = index + 1
        treeV[[index]] = tree
        errV[index] = errOfTree(tree)  
      }
    }
  }

  ID = which.min(errV)
  return(list(tree = treeV[[ID]], err = errV[ID]))
}


# Remove all edges whose length is not larger than th, unless the edge is between two original subclones.
# If th is missing, remove short edges until only the original subclone left.
ShrinkTree <- function(tree, numOriSubC, th = NULL)
{
  flag = FALSE
  edgeL = sapply(X = tree$Edges, FUN = length)
  edgeName = names(tree$Edges)
  allChild = sapply(strsplit(edgeName, split = "->"), function(x)x[2]);
  allParent = sapply(strsplit(edgeName, split = "->"), function(x)x[1]);
  oriSubC = paste("subC_", 0:numOriSubC, sep = "")
  orderID = order(edgeL, decreasing = FALSE)
  edgeNum = length(orderID)

  if (!is.null(th))
  {
    th = th * mean(edgeL)    
  }
  
  if (is.null(th))
  {
    if (length(tree$Edges) > numOriSubC)
    {
      flag = TRUE      
    }
    while (flag)
    {
      flagFor = FALSE
      for (i in 1:edgeNum)
      {
        if ((allChild[orderID[i]] %in% oriSubC) & (allParent[orderID[i]] %in% oriSubC))
        {
          next()
        }
        if ((allChild[orderID[i]] %in% oriSubC) & !(allParent[orderID[i]] %in% oriSubC))
        {
          tree = removeEdge(tree, rmEdge = edgeName[orderID[i]], rmNode = allParent[orderID[i]])
        }
        if (!(allChild[orderID[i]] %in% oriSubC) & (allParent[orderID[i]] %in% oriSubC))
        {
          tree = removeEdge(tree, rmEdge = edgeName[orderID[i]], rmNode = allChild[orderID[i]])
        }
        if (!(allChild[orderID[i]] %in% oriSubC) & !(allParent[orderID[i]] %in% oriSubC))
        {
          tree = removeEdge(tree, rmEdge = edgeName[orderID[i]], rmNode = allParent[orderID[i]])
        }
        flagFor = TRUE
        edgeL = sapply(X = tree$Edges, FUN = length)
        edgeName = names(tree$Edges)
        allChild = sapply(strsplit(edgeName, split = "->"), function(x)x[2]);
        allParent = sapply(strsplit(edgeName, split = "->"), function(x)x[1]);
        orderID = order(edgeL, decreasing = FALSE)
        edgeNum = length(orderID)
        if (edgeNum == numOriSubC)
        {
          flag = FALSE
        }
        break()
      }
      if ((i >= edgeNum) & (flag == TRUE) & (flagFor == FALSE))
      {
        flag = FALSE
      }
    }
  } else {
    if (edgeL[orderID[1]] <= th)
    {
      flag = TRUE      
    }
    while(flag)
    {
      flagFor = FALSE
      for (i in 1:edgeNum)
      {
        if (edgeL[orderID[i]] <= th)
        {
          if ((allChild[orderID[i]] %in% oriSubC) & (allParent[orderID[i]] %in% oriSubC))
          {
            next()
          }
          if ((allChild[orderID[i]] %in% oriSubC) & !(allParent[orderID[i]] %in% oriSubC))
          {
            tree = removeEdge(tree, rmEdge = edgeName[orderID[i]], rmNode = allParent[orderID[i]])
          }
          if (!(allChild[orderID[i]] %in% oriSubC) & (allParent[orderID[i]] %in% oriSubC))
          {
            tree = removeEdge(tree, rmEdge = edgeName[orderID[i]], rmNode = allChild[orderID[i]])
          }
          if (!(allChild[orderID[i]] %in% oriSubC) & !(allParent[orderID[i]] %in% oriSubC))
          {
            tree = removeEdge(tree, rmEdge = edgeName[orderID[i]], rmNode = allParent[orderID[i]])
          }
          flagFor = TRUE
          edgeL = sapply(X = tree$Edges, FUN = length)
          edgeName = names(tree$Edges)
          allChild = sapply(strsplit(edgeName, split = "->"), function(x)x[2]);
          allParent = sapply(strsplit(edgeName, split = "->"), function(x)x[1]);
          orderID = order(edgeL, decreasing = FALSE)
          edgeNum = length(orderID)
          if (edgeL[orderID[1]] > th)
          {
            flag = FALSE
          }
          break()
        }
      }
      if ((i >= edgeNum) & (flag == TRUE) & (flagFor == FALSE))
      {
        flag = FALSE
      }
    }
  }
  
  addSubC = setdiff(names(tree$Nodes), oriSubC)
  if (length(addSubC) > 0)
  {
    nodeMap = cbind(existingNode = c(oriSubC, addSubC[order(as.numeric(sapply(strsplit(addSubC, split = "_"), 
              function(x)x[2])), decreasing = FALSE)]), targetNode = c(oriSubC, 
              paste("subC", (numOriSubC+1):(length(tree$Nodes)-1), sep = "_")))
  }
  else {
    nodeMap = cbind(existingNode = oriSubC, targetNode = oriSubC)
  }
  rownames(nodeMap) = nodeMap[, "existingNode"]
  names(tree$Nodes) = nodeMap[names(tree$Nodes), "targetNode"]

  names(tree$Edges) = paste(nodeMap[sapply(strsplit(names(tree$Edges), split = "->"), function(x)x[1]), "targetNode"], 
                            nodeMap[sapply(strsplit(names(tree$Edges), split = "->"), function(x)x[2]), "targetNode"], sep = "->")
  names(tree$negEdges) = paste(nodeMap[sapply(strsplit(names(tree$negEdges), split = "->"), function(x)x[1]), "targetNode"], 
                               nodeMap[sapply(strsplit(names(tree$negEdges), split = "->"), function(x)x[2]), "targetNode"], sep = "->")
  tree$SubCinTree = nodeMap[tree$SubCinTree, "targetNode"]
  
  return(tree)
}

removeEdge <- function(tree, rmEdge, rmNode)
{

  kpNode = setdiff(strsplit(rmEdge, split = "->")[[1]], rmNode)
  tree$Edges[rmEdge] = c()
  tree$Nodes[rmNode] = c()
  ID = which(tree$SubCinTree == rmNode)
  tree$SubCinTree = tree$SubCinTree[-ID]
  tree$negEdges[paste(strsplit(rmEdge, split = "->")[[1]][2], strsplit(rmEdge, split = "->")[[1]][1], sep = "->")] = c()
  for (i in 1:length(tree$Edges))
  {
    parent = strsplit(names(tree$Edges)[i], split = "->")[[1]][1]
    child = strsplit(names(tree$Edges)[i], split = "->")[[1]][2]
    if (parent == rmNode)
    {
      tree$Edges[i] = list(setdiff(tree$Nodes[[child]], tree$Nodes[[kpNode]]))
      names(tree$Edges)[i] = paste(kpNode, child, sep = "->")
      tree$negEdges[i] = list(setdiff(tree$Nodes[[kpNode]], tree$Nodes[[child]]))
      names(tree$negEdges)[i] = paste(child, kpNode, sep = "->")
    }
    if (child == rmNode)
    {
      tree$Edges[i] = list(setdiff(tree$Nodes[[kpNode]], tree$Nodes[[parent]]))
      names(tree$Edges)[i] = paste(parent, kpNode, sep = "->")
      tree$negEdges[i] = list(setdiff(tree$Nodes[[parent]], tree$Nodes[[kpNode]]))
      names(tree$negEdges)[i] = paste(kpNode, parent, sep = "->")      
    }
  }

  return(tree)
}



plotTree = function(tree, filePath)
{
  #### useage: plotTree(Result) whrere Result is treeList[[i]]
  #### output: pdf containing tree structure
  
  N = length(tree$Nodes);
  
  adj.matrix = matrix(0, ncol=N, nrow=N); 
  rownames(adj.matrix) = names(tree$Nodes);
  colnames(adj.matrix) = names(tree$Nodes);
  
  for (i in 1:N)
  {
    a = unlist(strsplit(names(tree$Edges[i]), "->"));
    r = match(a[1],rownames(adj.matrix));
    c = match(a[2],colnames(adj.matrix));
    adj.matrix[r,c] = 1;
  }    
  
  
  if (!is.null(tree$w))
  {
    nodeNameMap = cbind(oriName = names(tree$Nodes), newName = paste(names(tree$Nodes), round(tree$w, digits = 2), sep = ":"))
    rownames(nodeNameMap) = nodeNameMap[, "oriName"]
    rownames(adj.matrix) = nodeNameMap[rownames(adj.matrix), "newName"]
    colnames(adj.matrix) = nodeNameMap[colnames(adj.matrix), "newName"]
  }
  
  tree.phylo = new("graphAM", adjMat = adj.matrix, edgemode = "directed");
  
  png(filename = paste(filePath, ".png", sep = ""), width = 1000+50*N, height = 500+40*N, units = "px");

  edgeLabel = as.character(sapply(X = tree$Edges, FUN = length))
  if (!is.null(tree$w))
  {
    names(edgeLabel) = paste(nodeNameMap[sapply(X = strsplit(names(tree$Edges), split = "->"), function(x)x[1]), "newName"],
                             nodeNameMap[sapply(X = strsplit(names(tree$Edges), split = "->"), function(x)x[2]), "newName"], sep = "~")    
  } else {
    names(edgeLabel) = paste(sapply(X = strsplit(names(tree$Edges), split = "->"), function(x)x[1]),
                             sapply(X = strsplit(names(tree$Edges), split = "->"), function(x)x[2]), sep = "~")    
  }
  tree.phylo = layoutGraph(tree.phylo, attrs = list(node = list(fillcolor = "lightblue", fontsize = 30), 
                                                    edge = list(arrowsize = 1, fontsize = 30)), edgeAttrs = list(label = edgeLabel))
  renderGraph(tree.phylo)
  dev.off();
}    


# Major function to construct tree.
# SM is the Z matrix from subclone algorithm
# CM is the L matrix from subclone algorithm
# w is the cell proportion estimate from cubclone algorithm
# prunTh is the threshold to remove edge, range is 0~1, it is the proportion of average edge length
ConstructTree <- function(SM, CM, w, prunTh = 0.1)
{
  # remove the first column from SNV and copy number matrices, because it is background
  # also change non-zero elements in SNV matrix to 1
  SM = SM[, -1, drop = FALSE]
  CM = CM[, -1, drop = FALSE]
  idNon0 = which(SM != 0)
  SM[idNon0] = 1
  
  # remove the first element from w vector, because it is the cell proportion of background
  w = w[-1]
  
  # if a column of SM matrix contains all zeros, then we just remove it and corresponding w
  rmID = intersect(which(colSums(SM) == 0), which(colSums(CM != 2) == 0))
  if (length(rmID) > 0)
  {
    SM = SM[, -rmID, drop = FALSE]
    CM = CM[, -rmID, drop = FALSE]
    w = w[-rmID]
  }

  if (dim(SM)[2] > 1)
  {
    # if there are duplicated columns in SM, we need to combine them
    M = rbind(SM, CM)
    rmFlag = rep(FALSE, dim(M)[2]);
    for (i in 1:(dim(M)[2]-1))
    {
      for (j in (i+1):dim(M)[2])
      {
        if (sum(M[, i] == M[, j]) == dim(M)[1])
        {
          rmFlag[j] = TRUE;
          w[i] = w[i] + w[j]
        }
      }
    }
    rmID = which(rmFlag);
    if (length(rmID) > 0)
    {
      SM = SM[, -rmID, drop = FALSE]
      CM = CM[, -rmID, drop = FALSE]
      w = w[-rmID]
    }
  }
  
  # generate mutation profiles
  result = TransferToMutationProfile(SM = SM, CM = CM);

  if (dim(SM)[2] == 1)
  {
    Nodes = vector("list", 2)
    names(Nodes) = c("subC_0", "subC_1")
    Nodes[[2]] = result$mutEvePro[[1]]
    Edges = vector("list", 1)
    names(Edges) = "subC_0->subC_1"
    Edges[[1]] = result$mutEvePro[[1]]
    SubCinTree = c("subC_1")
    names(SubCinTree) = "subC_1"
    negEdges = vector("list", 1)
    names(negEdges) = "subC_1->subC_0"
    wAll = c(0, w)
    names(wAll) = names(Nodes)
    tree = list(Nodes = Nodes, Edges = Edges, SubCinTree = SubCinTree, negEdges = negEdges, w = wAll)
    err = errOfTree(tree)
    allTrees = vector("list", 1)
    allTrees[[1]] = tree
    allErr = err
    return(list(allTrees = allTrees, allErr = allErr, tree = tree, err = err))    
  }

  errReturn = c()
  treeReturn = c()
  
  # Generate the best initial tree
  iniTreeV = vector("list", length(result$mutEvePro)*(length(result$mutEvePro)-1)/2)
  errV = rep(NA, length(result$mutEvePro)*(length(result$mutEvePro)-1)/2)
  index = 0
  for (i in 1:(length(result$mutEvePro)-1))
  {
    for (j in (i + 1):length(result$mutEvePro))
    {
      index = index + 1
      treeResult = generateInitialTree(mep = result$mutEvePro, ID1 = i, ID2 = j);
      errV[index] = treeResult$err
      iniTreeV[[index]] = treeResult$tree
    }
  }

  idMin = which.min(errV)
  errReturn[1] = errV[idMin]
  tree = iniTreeV[[idMin]]
  treeReturn[[1]] = tree
  
  treesToAdd = setdiff(names(tree$Nodes), c("subC_0", tree$SubCinTree));
  while(length(treesToAdd) > 0)
  {
    treeV = vector("list", length(treesToAdd))
    errV = rep(NA, length(treesToAdd))
    for (k in 1:length(treesToAdd))
    {
      treeResult = AddOneGenomeToExistingTree(tree, treesToAdd[k]);
      treeV[[k]] = treeResult$tree
      errV[k] = treeResult$err
    }

    idMin = which.min(errV)
    errReturn = c(errReturn, errV[idMin])
    tree = treeV[[idMin]]
    treeReturn[[length(errReturn)]] = tree
    treesToAdd = setdiff(names(tree$Nodes), c("subC_0", tree$SubCinTree));
  }
  
  # Shrink tree according to th
  tree = ShrinkTree(treeReturn[[length(treeReturn)]], dim(SM)[2], th = prunTh)
  tree$w = c(0, w, rep(0, length(tree$Nodes)-1-dim(SM)[2]))
  names(tree$w) = names(tree$Nodes)
  err = errOfTree(tree)
  
  return(list(allTrees = treeReturn, allErr = errReturn, tree = tree, err = err))
}


calibrateWandIndex <- function(tree)
{
  if (length(tree$Nodes) == 2)
  {
    tree$w["subC_0"] = tree$w["subC_1"]
    return(tree)
  }
  
  allParent = sapply(strsplit(names(tree$Edges), split = "->"), function(x)x[1]) 
  allChild = sapply(strsplit(names(tree$Edges), split = "->"), function(x)x[2])
  new_w = rep(NA, length(tree$w))
  names(new_w) = names(tree$w)
  wNode = vector("list", length(tree$w))
  names(wNode) = names(tree$w)
  for(i in 1:length(tree$w))
  {
    iNodes = c()
    iNodes[1] = list(names(tree$w)[i])
    index = 1
    indexID = which(allParent == iNodes[[1]])
    while(length(indexID) > 0)
    {
      index = index + 1
      iNodes[index] = list(allChild[indexID])
      indexID = which(allParent %in% iNodes[[index]])
    }
    for (j in 1:length(iNodes))
    {
      wNode[[i]] = c(wNode[[i]], iNodes[[j]])
    }
    wNode[[i]] = unique(wNode[[i]])
    new_w[i] = sum(tree$w[wNode[[i]]])
  }
  tree$w = new_w

  orderID = order(tree$w[2:length(tree$w)], decreasing = TRUE)
  nodeMap = cbind(oriNode = c("subC_0", names(tree$w)[2:length(tree$w)][orderID]), 
                  mappingNode = c("subC_0", paste("subC_", 1:(length(tree$w)-1), sep = "")))
  rownames(nodeMap) = nodeMap[, "oriNode"]
  
  names(tree$w) = nodeMap[names(tree$w), "mappingNode"]
  ID = order(as.numeric(sapply(strsplit(names(tree$w), split = "_"), function(x)x[2])), decreasing = FALSE)
  tree$w = tree$w[ID]

  names(tree$Nodes) = nodeMap[names(tree$Nodes), "mappingNode"] 
  ID = order(as.numeric(sapply(strsplit(names(tree$Nodes), split = "_"), function(x)x[2])), decreasing = FALSE)
  tree$Nodes = tree$Nodes[ID]
  
  names(tree$Edges) = paste(nodeMap[sapply(strsplit(names(tree$Edges), split = "->"), function(x)x[1]), "mappingNode"], 
                     nodeMap[sapply(strsplit(names(tree$Edges), split = "->"), function(x)x[2]), "mappingNode"], sep = "->")
  
  names(tree$negEdges) = paste(nodeMap[sapply(strsplit(names(tree$negEdges), split = "->"), function(x)x[1]), "mappingNode"], 
                        nodeMap[sapply(strsplit(names(tree$negEdges), split = "->"), function(x)x[2]), "mappingNode"], sep = "->")

  tree$SubCinTree = nodeMap[tree$SubCinTree, "mappingNode"]
  names(tree$SubCinTree) = tree$SubCinTree
  
  return(tree)
}


# library(Rgraphviz);
# 
# sampleID = c("S2", "S4", "S5", "S6", "S7")
# for (fileIndex in 1:length(sampleID))
# {
#   load(paste("output-", sampleID[fileIndex], "/", sampleID[fileIndex], ".txt.RData", sep = ""))
#   result = ConstructTree(SM = point.est$Z, CM = point.est$L, w = point.est$w, prunTh = 0.1)
#   
#   result$tree = calibrateWandIndex(result$tree)
# 
#   plotTree(tree = result$allTrees[[length(result$allTrees)]], 
#            filePath = paste("BeforeShrunk-", sampleID[fileIndex], "-err-", 
#            round(result$allErr[length(result$allErr)]/sum(sapply(X = result$allTrees[[length(result$allTrees)]]$Nodes, FUN = length)), digits = 2), sep = ""))
#   plotTree(tree = result$tree, filePath = paste("calibrateAfterShrunk-", sampleID[fileIndex], "-err-", 
#           round(result$err/sum(sapply(X = result$tree$Nodes, FUN = length)), digits = 2), sep = ""))
# }
