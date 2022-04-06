##### (1) Get the context vectors for each word #####
# Inputs needed:
#   notes = vector of the set of notes as characters
#   term.vecs.mat = term embeddings in the form of "spring_...txt":
#     - the first column is the term
#     - the second column is the first component of the embedding, the third column is the second component, etc...
#   term.freq = data.frame of three columns
#     - term = the word
#     - count = the number of times the word occurs in the corpus
#     - freq = the frequency of the word in the corpus (so this column sum is equal to 1)
#   word.sample = sample of words that are to be run on.

term_embeddings_path = ''
load(term_embeddings_path)
notes = as.character(note.data$ReportText)

# maximum number of word occurences for getting contexts
max.count = 1000 

# initialize flag for notes
note.flag = rep(F, length(notes))

# context list will temporarily store the contexts for each word
context.list = list()

# finished list that stores contexts for words once they get the max# of contexts
finished.list = list()

# word sample is the unique words in "notes", found previously
# initialize the context list
for(i in 1:length(word.sample)){
  tmp = list(word.sample[i],0,'',0)
  names(tmp) = c('word','num.occurrences','context.words','context.vector')
  context.list[[i]] = tmp
}

system.time({
  # remaining words keeps track of what words still need to be processed.
  remaining.words = sapply(context.list, '[[', 1)
  
  # cycle through each note
  for(i in 1:length(notes)){
    if(i %% 5e3 == 0){
      print(sprintf('scanned %i notes',i))
    }
    
    
    # check if all words are processed (as in the max number of contexts is reached)
    if(length(context.list)==0){
      print('finished all words')
      break
    }
    
    # initialize tracker of words completed
    word.ind.complete = c()
    
    # split the note into a vector of words
    note = pre.process.note(notes[i])
    
    # cycle through each word in note and remaining list
    #for(w in 1:length(remaining.words)){ 
    for(w in which(remaining.words %in% note)){
      
      # get the context for this word
      contexts = get.contexts(note = note, term = remaining.words[w], window = 10)
      
      # if there is no context (this shouldn't happen)
      if(is.null(contexts)){note.flag[i] = T; next}
      
      # update the "num.occurrences" of the term
      context.list[[w]][[2]] = context.list[[w]][[2]] + length(contexts)
      
      # store the words in the context
      context.list[[w]][[3]] = c(context.list[[w]][[3]], unlist(contexts))
      
      # if the number of word occurrences reaches the max count of terms.
      if(context.list[[w]][[2]] >= max.count){
        
        # store the indices of words that are done
        word.ind.complete = c(word.ind.complete,w)
      }
    }
    
    # if there are words finished
    if(length(word.ind.complete)>0){
      
      # store the finished words
      finished.list = c(finished.list, context.list[word.ind.complete])
      
      # remove the word from context list
      context.list = context.list[-word.ind.complete]
      
      # update the list of remaining words
      remaining.words = sapply(context.list, '[[', 1)
    }
  }
  
  # make the final finished list of words, since some won't reach the max.count
  finished.list = c(finished.list,context.list)
  
  # get context vectors by taking the vector average of the context terms
  finished.list = lapply(finished.list, function(xx) {xx$context.vector = vector.average(xx$context.words, term.vecs.mat); return(xx)})
})

# add in feature vectors for each word (the original termvectors for that vector)
finished.list = lapply(finished.list, function(xx) {
  ind = which(term.vecs.mat$V1==xx$word)
  feature.vector = as.numeric(term.vecs.mat[ind,-1])
  xx$feature.vector = feature.vector
  return(xx)
})

##### (2) Run the multiple linear regression on all words to get A #####

# ^ loads the finished list as created above, the term.freq and term.vecs.mat matrix

# creates a formula and data.frame of the feature and context vectors.
tmp <- prep_alacarte(finished.list)
data = tmp[[1]]; formula = tmp[[2]]; rm(tmp)

# create weights that weight by the log of the count of the terms.
weights = create_log_weights(finished.list, term.freq)

# run the linear model
lm.fit = lm(formula, data=data, weights=weights)
A <- get.A.data.frame(lm.fit)
A = as.matrix(A)

# Trying a higher cutoff
weights = create_log_weights(finished.list,term.freq,min.count=1000)
lm.fit = lm(formula,data=data,weights=weights)

# get A matrix and intercept from the model
tmp = get_A_data_frame_and_intercept(lm.fit)
intercept = tmp[[1]]
A = as.matrix(tmp[[2]])

# create the new embeddings basis from A and the intercept
term.vecs.mat <- create_new_embeddings_set(finished.list, A, intercept)
