# Assume root is called Z; choice of priors for Z are lognormal or uniform.
# Assume other distributions are bi/multinomial for nodes.
# Assume branching probabilities are beta/Dirichlet.
# Can generate .mod or .txt files

makeJAGStree <- function(data, prior = "lognormal", filename = "JAGSmodel.mod"){

  ###############################################################################
  # Check overall data structure and function inputs
  # Check prior input
  if(!prior == "lognormal" && !prior == "uniform"){
    stop("Root prior must be 'lognormal' or 'uniform'.")
  }

  ###############################################################################
  # check structure of data here
  # make sure it is a data frame
  if(!is.data.frame(data)){
    stop('data must be dataframe type.')
  }

  # if no column x, print 'no column x'
  if(is.null(data$from)){
    stop('data must have \'from\' column')
  }
  if(is.null(data$to)){
    stop('data must have \'to\' column')
  }

  ###############################################################################
  # Generating general text for JAGS model text file
  preamble <- "# This JAGS model was created using 'makeJAGStree' in the 'JAGStree' package in R."
  warning <- "# The root may have lognormal or discretized uniform prior.
# Branching and leaf prior distributions are assumed Dirichlet and Multinomial, respectively. \n"

  ###############################################################################
  # Start "data" section string to add to
  data.sec <- "data { \n"

  ###############################################################################
  # Make "model" section string to add to
  model.sec <- "model { \n"

  ###############################################################################
  # Extract info from data, and write model

  # Find root(look for node in 'from' column that does not exist in 'to' column)
  root <- data$from[!data$from %in% data$to]

  # if theres more than one root, it's not a tree, so abort.
  if(!length(unique(root)) == 1){
    stop("Data is not tree structured.")
  }

  # Add root prior to model section
  if(prior == "lognormal"){
    model.sec <- paste(model.sec, root[1], ".cont ~ dlnorm(mu, tau); \n",
                       root[1], " <- round(", root[1], ".cont); \n", sep = "")
  }
  if(prior == "uniform"){
    model.sec <- paste(model.sec, root[1], ".cont ~ dunif(Lz, Uz); \n",
                       root[1], " <- round(", root[1], ".cont); \n", sep = "")
  }

  # Find and count children of root first.  Then iterate through additional parent
  # nodes until all distributions for Dirichlet branching and Multinomial
  # children have been set
  root.children <- data$to[data$from == root[1]]
  k <- length(root.children)

  # write the parameters of that multidimensional distribution to the 'data' section
  branch.name <- paste("p", root[1], sep = "")
  data.sec <- paste(data.sec, branch.name, ".params <- c(", sep = "")
  for (i in 1:k) {
    if (i < k){
      data.sec <- paste(data.sec, branch.name, i, ", ", sep = "")
    }else{
      data.sec <- paste(data.sec, branch.name, i, "); \n", sep = "")
    }
  }

  # create the first branching prior distribution
  if (k < 3){
    model.sec <- paste(model.sec, branch.name, " ~ dbeta(", branch.name,
                       ".params[1], ", branch.name, ".params[2]); \n",
                       sep = "")
  }else{
    model.sec <- paste(model.sec, branch.name, " ~ ddirch(", branch.name,
                       ".params); \n", sep = "")
  }

  # set multinomial/binomial distribution on children of root
  child.vec <- paste(root.children, collapse = "")
  if (k < 3){
    # use binomial with beta
    model.sec <- paste(model.sec, child.vec, "[1] ~ dbinom(", branch.name,
                       ", ", root[1], "); \n", child.vec, "[2] <- ",
                       root[1], " - ", child.vec, "[1]; \n\n", sep = "")
  }else{
    # set up parameters for multinomial (in binomial segments) with Dirichlet
    model.sec <- paste(model.sec, root[1], ".bin[1] <- ", root[1], "; \n",
                       branch.name, ".bin[1] <- ", branch.name, "[1]; \n",
                       "for (i in 2:", k, "){ \n  ", root[1], ".bin[i] <- ",
                       root[1], ".bin[i-1] - ", child.vec, "[i-1] \n  ",
                       branch.name, ".bin[i] <- ", branch.name, "[i]/(sum(",
                       branch.name, "[i:", k, "])) \n} \n", sep = "")

    # write binomial distributions
    model.sec <- paste(model.sec, "for (i in 1:", k-1, "){ \n  ", child.vec,
                       "[i] ~ dbinom(", branch.name, ".bin[i], ", root[1],
                       ".bin[i]) \n} \n", child.vec, "[", k, "] <- ",
                       root[1], ".bin[1] - sum(", child.vec, "[1:", k-1,
                       "]); \n\n", sep = "")
  }



  # one at a time, cycle through other parents (which are not the root)
  parent.list <- unique(data$from[!data$from == root[1]])
  while(length(parent.list) > 0){
    node <- parent.list[1]

    # Find and count children of current node.
    node.children <- data$to[data$from == node]
    k <- length(node.children)

    # establish siblings of current node to name it correctly in modelling code
    node.parent <- data$from[data$to == node]
    node.siblings <- data$to[data$from == node.parent]
    node.name <- paste(node.siblings, collapse = "")
    node.name <- paste(node.name, "[", which(node.siblings == node),
                       "]", sep = "")


    # write the parameters of that multidimensional distribution to the 'data' section
    branch.name <- paste("p", node, sep = "")
    data.sec <- paste(data.sec, branch.name, ".params <- c(", sep = "")
    for (i in 1:k) {
      if (i < k){
        data.sec <- paste(data.sec, branch.name, i, ", ", sep = "")
      }else{
        data.sec <- paste(data.sec, branch.name, i, "); \n", sep = "")
      }
    }

    # create the branching distribution
    if (k < 3){
      model.sec <- paste(model.sec, branch.name, " ~ dbeta(", branch.name,
                         ".params[1], ", branch.name, ".params[2]); \n",
                         sep = "")

    }else{
      model.sec <- paste(model.sec, branch.name, " ~ ddirch(", branch.name,
                         ".params); \n", sep = "")
    }

    # set multinomial distribution on children of current node
    child.vec <- paste(node.children, collapse = "")

    # set up parameters for multinomial (in binomial segments)
    if (k < 3){
      # use binomial with beta
      model.sec <- paste(model.sec, child.vec, "[1] ~ dbinom(", branch.name,
                         ", ", node.name, "); \n", child.vec, "[2] <- ",
                         node.name, " - ", child.vec, "[1]; \n \n", sep = "")
    }else{
      # set up parameters for multinomial (in binomial segments) with Dirichlet
      model.sec <- paste(model.sec, node, ".bin[1] <- ", node.name, "; \n",
                         branch.name, ".bin[1] <- ", branch.name, "[1]; \n",
                         "for (i in 2:", k, "){ \n  ", node, ".bin[i] <- ",
                         node, ".bin[i-1] - ", child.vec, "[i-1] \n  ",
                         branch.name, ".bin[i] <- ", branch.name, "[i]/(sum(",
                         branch.name, "[i:", k, "])) \n} \n", sep = "")

      # write multinomial distributions
      model.sec <- paste(model.sec, "for (i in 1:", k-1, "){ \n  ", child.vec,
                         "[i] ~ dbinom(", branch.name, ".bin[i], ", node,
                         ".bin[i]) \n} \n", child.vec, "[", k, "] <- ",
                         node, ".bin[1] - sum(", child.vec, "[1:", k-1,
                         "]); \n \n", sep = "")
    }

    # take current parent off the parent list and repeat with remaining parents
    parent.list <- parent.list[! parent.list %in% node]
  }

  # Close off data and model sections
  data.sec <- paste(data.sec, "} \n", sep = "")
  model.sec <- paste(model.sec, "} \n", sep = "")

  # finally, write the entire model to a .mod (default) or .txt file
  writeLines(c(preamble, warning, data.sec, model.sec), filename)
}
