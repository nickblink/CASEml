library(bit64)
library(data.table)
library(RODBC)
library(randomForest)
library(AUC)
library(flexmix)
library(ggplot2)

##### helper functions ######

#' get acronym subset from list
#' 
#' @param acronym_list a container for all acronym data.
#' @param subset the subset of acronyms (either as numbers or characters).
#' @return an integer vector of indices of acronyms of interest
get_subset <- function(acronym_list, subset){
  # if not running on a subset of acronyms, run on all of them
  if(is.null(subset)){subset = 1:length(acronym_list)}
  
  # if subset is a set of acronym characters, get the indices
  else if(class(subset)=='character'){
    # get the acronyms
    acronyms = sapply(acronym_list, function(xx) tolower(xx$acronym))
    
    # turning character to lower
    subset = tolower(subset)
    
    # find the indices of these acronyms
    subset = which(acronyms %in% subset)
  }
  return(subset)
}

# convert patient and date_time to patient:date
create_patient_date = function(patients, date_time){
  paste0(patients, ':', sapply(date_time,
                               function(xx) gsub('-','',strsplit(xx,' ')[[1]][1])))
}

# split patient.date to patient or date.
split_patient_date <- function(patient_dates, return_type = 'patient', split_by = ':'){
  if(return_type == 'patient'){
    output = sapply(patient_dates, function(xx) strsplit(xx, split_by)[[1]][1])
  }else if(return_type == 'date'){
    output = sapply(patient_dates, function(xx) strsplit(xx, split_by)[[1]][2])
  }else{
    print('not a proper return type')
  }
  return(output)
}

# count number of acronyms in a processed note.
get_acronym_count <- function(note, acronym){
  # pre-process note if it's not already.
  if(length(note==1)){
    note_list = pre_process_note(note)
  }else{
    note_list = note
  }
  # get the number of times the acronym appears
  count = sum(note_list == acronym)
  
  # set NAs to 0
  if(is.na(count)){
    count = 0
  }
  
  # remove the name
  names(count) = NULL
  
  # return the results!
  return(count)
}

## takes in an acronym and cuis to ignore and returns the cuis for that acronym
# dict_path contains path to a dictonary with the columns for the term, the CUI grouping, and the phenotype
get_cui_list <- function(acronym, ignore_cui = NULL, dict_path ){
  # pull in the nlp dict
  nlp_dict = fread(dict_path, header = F, col.names = c('term','cui','pheno'))
  
  # get the cuis
  acr_cuis = nlp_dict$cui[nlp_dict$pheno == acronym]
  
  # unpack the cuis (some are separated by colons)
  acr_cuis = unlist(lapply(acr_cuis, function(xx) strsplit(xx,';')[[1]]))
  
  # ignore the cuis asked for
  acr_cuis = setdiff(acr_cuis, ignore_cui)
  
  return(acr_cuis)
}

# takes in the acronym list and returns the acronym list with the cui list for each acronym
get_cui_list_wrapper <- function(acronym_list){
  # just doing a for loop because this doesn't take much time anyway
  for(i in 1:length(acronym_list)){
    acronym_list[[i]]$cui_list = get_cui_list(acronym_list[[i]]$acronym, acronym_list[[i]]$longform_cui)
  }
  
  # return it!
  return(acronym_list)
}

# get the different prevalence estimates for the acronyms
estimate_prevalence <- function(acronym_list, charts, train_data, test_data, subset = NULL){
  # get the subset
  subset = get_subset(acronym_list, subset)
  
  for(i in subset){
    # pull out acronym components
    acronym = acronym_list[[i]]$acronym
    acr_cui = acronym_list[[i]]$acronym_cui
    main_cui = acronym_list[[i]]$main_cui
    acr_ICD = paste0(acronym,'_ICD')
    
    # make the data set
    D = train_data[train_data[,which(colnames(train_data)==acr_cui)]>0,
                   which(colnames(train_data) %in% c(acr_ICD,acr_cui ,main_cui))]
    D = as.data.frame(as.matrix(D))
    
    # add in patient and filter value
    D$patient = split_patient_date(rownames(D))  
    D$filter = ifelse(D$patient %in% acronym_list[[i]]$filter_pos_patients,1,0)
    
    # pull out the charts set and add ICD and cui features
    charts_sub = charts[[acronym]]
    charts_sub$ICD = test_data[match(charts_sub$patient_date, rownames(test_data)),
                               which(colnames(test_data) == acr_ICD)]
    charts_sub$main_cui = test_data[match(charts_sub$patient_date, rownames(test_data)),
                                    which(colnames(test_data) == main_cui)]
    
    # convert label to integer
    if(class(charts_sub$label) == 'factor'){
      charts_sub$label = as.integer(charts_sub$label)
      # if shifted up by one, subtract
      if(max(charts_sub$label) == 2){
        charts_sub$label = charts_sub$label - 1
      }
    }
    
    # make filter positive set
    D_pos = D[D$filter == 1,]
    
    # make the charts positive set
    charts_pos = charts_sub[charts_sub$filter == 'YES',]
    
    # get the true_prevalence  
    filterpos_true = mean(charts_pos$label)
    
    # ICD-based prevalence
    filterpos_ICD = sum(D_pos[,acr_ICD]*D_pos[,acr_cui])/sum(D_pos[,acr_cui])
    
    # ICD and cui-based prevalence
    filterpos_ICD_CUI = sum(as.integer(D_pos[,acr_ICD] + D_pos[,main_cui] > 0)*D_pos[,acr_cui])/sum(D_pos[,acr_cui])
    
    # add inverse sampling weights to charts
    charts_sub$weight = ifelse(charts_sub$filter == 'YES', 1, acronym_list[[i]]$non_filter_ratio)
    
    # label (true distribution)
    nofilter_true = sum(charts_sub$label*charts_sub$weight)/sum(charts_sub$weight)
    
    # ICD-based prevalence
    nofilter_ICD = sum(D[,acr_ICD]*D[,acr_cui])/sum(D[,acr_cui])
    
    # ICD and cui-based prevalence
    nofilter_ICD_CUI = sum(as.integer(D[,acr_ICD] + D[,main_cui] > 0)*D[,acr_cui])/sum(D[,acr_cui])
    
    # print out results
    print(sprintf('%s:',acronym))
    print('----FILTER POSITIVE SET-----')
    print(sprintf('true prevalence: %s', filterpos_true))
    print(sprintf('ICD estimated prevalence: %s', filterpos_ICD))
    print(sprintf('ICD and CUI estimated prevalence: %s', filterpos_ICD_CUI))
    print('----FULL SET (NO FILTER)-----')
    print(sprintf('true prevalence: %s', nofilter_true))
    print(sprintf('ICD estimated prevalence: %s', nofilter_ICD))
    print(sprintf('ICD and CUI estimated prevalence: %s', nofilter_ICD_CUI))
    
    # add in prevalences to acronym list
    tmp_list = list(filterpos_true, filterpos_ICD, filterpos_ICD_CUI,
                    nofilter_true, nofilter_ICD, nofilter_ICD_CUI)
    names(tmp_list) = c('filterpos_true', 'filterpos_ICD', 'filterpos_ICD_CUI',
                        'nofilter_true', 'nofilter_ICD', 'nofilter_ICD_CUI')
    acronym_list[[i]]$prevalences = tmp_list
  }
  
  # return the acronym list
  return(acronym_list)
}

# get logit score from a probability ranging 0 to 1
logit_0_to_1 <-function(p){
  if(p == 1){print('b');p = 1 - 1e-5}
  if(p == 0){print('a');p = 1e-5}
  logit = log(p/(1-p))
  return(logit)
}

# get logit score from a probability rangine -1 to 1
logit_neg1_to_1 <-function(xx){
  p = (xx + 1)/2
  if(p == 1){p = 1 - 1e-5}
  logit = log(p/(1-p))
  return(logit)
}

# get the sigmoid-distance from center
sigmoid = function(d, alpha = 0.3, r = 9){
  val = 1/(1 + exp(alpha*(d-r)))
  return(val)
}

# convert object roc to a data.frame
roc_to_df <- function(roc){
  df = data.frame(cutoffs = roc$cutoffs,
                  fpr = roc$fpr,
                  tpr = roc$tpr)
  return(df)
}

# average vectors
vector_average <- function(terms, term_vecs_mat, weights.mat = NULL, A = NULL, warnings.on=F){
  #ll = setdiff(terms,term_vecs_mat$V1)
  #if(length(ll)>0){print(paste0('missing ',length(ll),' terms from term list'))}
  matid = match(terms,term_vecs_mat$V1)
  matid = matid[!is.na(matid)]
  tmp = term_vecs_mat[matid,]
  if(nrow(tmp)==0){
    if(warnings.on){warning('no terms in term_vecs_mat')}
    return()
  }
  if(is.null(weights.mat)){
    tmp_mat = tmp[,-1]
  }else{
    matid = match(tmp$V1,weights.mat[,1])
    if(any(is.na(matid))){weights[is.na(weights)] = max(weights$weight)}
    weights = weights.mat[matid,2]
    tmp_mat = tmp[,-1]*weights
  }
  vector_average = apply(tmp_mat, 2, mean)
  if(!is.null(A)){
    names.vector = names(vector_average)
    vector_average = as.numeric(as.matrix(A)%*%vector_average)
    names(vector_average) = names.vector
  }
  return(vector_average)
}

# get cosine similarity of two vectors
cosine_similarity <- function(a,b,n.remove=0){
  a = a[(n.remove+1):length(a)]
  b = b[(n.remove+1):length(b)]
  return(sum(a*b)/sqrt(sum(a^2)*sum(b^2)))
}

#
##### word-vector/embeddings functions #####

# pre-process the note into a list of terms
pre_process_note <- function(note, list=TRUE){
  note = gsub('[[:punct:]]|[[:space:]]|[[:blank:]]',' ',note)
  note <- tolower(note)
  if(list==TRUE){
    note_list = strsplit(note,' ')[[1]]
    note_list = note_list[!(note_list %in% c("","s"))]
    return(note_list)
  }else{
    return(note)
  }
}

# used in "get_contexts" to handle multiple word terms, like "rheumatoid arthritis"
subset_index <- function(list.1,list.2){
  # list.1 = subset to check.
  # list.2 = full list
  ind = which(list.2 == list.1[1])
  index.list = NULL
  if(length(ind)==0){
    return(index.list)
  }else{
    for(ii in ind){
      tmp = TRUE
      for(j in 2:length(list.1)){
        if((length(list.2) < ii+j-1) | (list.1[j] != list.2[ii+j-1])){
          tmp = FALSE
        }
      }
      if(tmp){
        index.list = c(index.list,ii)
      }
    }
  }
  return(index.list)
}

# get the contexts for a term in a note.
get_contexts <- function(note, left_right_split = F, term='ra', window = 10){
  # if not pre-processed, then do so!
  if(length(note)==1){note = pre_process_note(note)}
  context.list = NULL
  iter = 0
  # If the term is multiple words
  if(grepl(' ',term)){
    terms = strsplit(term,' ')[[1]] 
    indices = subset_index(terms,note)
    for(loc in indices){
      iter = iter + 1
      if(left_right_split){
        left = note[c(max(0,loc-window)):(loc-1)]
        if(loc==(length(note)+1-length(terms))){
          right = c()
        }else{
          right = note[(min(length(note),loc+length(terms))):(min(length(note),loc+length(terms)-1+window))]
        }
        context.list[[iter]] = list(left = left, right = right)
      }else{
        if(loc==(length(note)+1-length(terms))){
          context = note[(max(0,loc-window)):(loc-1)]
        }else{
          context = note[c((max(0,loc-window)):(loc-1),
                           (min(length(note),loc+length(terms))):(min(length(note),loc+length(terms)-1+window)))]
        }
        context.list[[iter]] = context
      }
    }
    # if the term is a single word
  }else{
    note.short = note # in this one I delete up to the previous term
    loc_prev=0
    indices = which(note==term)
    if(length(indices)==0){
      return(context.list)
    }else{
      for(loc in indices){
        iter = iter + 1
        if(left_right_split){
          left = note[c(max(0,loc-window)):(loc-1)]
          if(loc==length(note)){
            right = c()
          }else{
            right = note[(loc+1):(min(length(note),loc+window))]
          }
          context.list[[iter]] = list(left = left, right = right)
        }else{
          if(loc==length(note)){
            context = note[(max(0,loc-window)):(loc-1)]
          }else{
            context = note[c((max(0,loc-window)):(loc-1),(loc+1):(min(length(note),loc+window)))]
          }
          context.list[[iter]] = context
        }
      }
    }
  }
  return(context.list)
}

# weight the termvector matrix by the weighting matrix
get_weighted_mat <- function(term_vecs_mat, weights.mat){
  # match the terms of both
  matid = match(term_vecs_mat$V1,weights.mat[,1])
  if(any(is.na(matid))){weights[is.na(weights)] = max(weights$weight)}
  
  # get the weights corresponding to the termvecsmat terms
  weights = weights.mat[matid,2]
  
  # weight it!
  tmp_mat = term_vecs_mat[,-1]*weights
  tmp_mat = cbind(term_vecs_mat$V1, tmp_mat)
  
  # update colnames and return
  colnames(tmp_mat) = colnames(term_vecs_mat)
  return(tmp_mat)
}

# for the notes, get the indices of the terms in the term vector matrix.
get_note_termvec_indices <- function(note, term, term.vector, window = 10){
  # pre process note
  note_list = pre_process_note(note)
  
  # get contexts for the note
  contexts = get_contexts(note_list, term = term, window = window)
  
  # skip if acronym not in the note
  if(length(contexts) == 0){return(NULL)}
  
  # make a list of the contexts
  ind.list = lapply(contexts, function(xx) {tmp = match(xx, term.vector); tmp[!is.na(tmp)]})
  ind.dt = do.call(rbind, lapply(seq_along(ind.list), 
                                 function(i) {xx = ind.list[[i]]; data.table(instance = rep(i,length(xx)), matid = xx)}))
  
  return(ind.dt)
}

# get the vectors for the contexts of a note given an acronym (or term)
get_note_vectors <- function(notes, acronym, term_vecs_mat, sigmoid_weighting = F, window = 10){

  if(xor(window == 12, sigmoid_weighting)){
    warning(sprintf('are you sure you want window %s and sigmoid weighting = %s?', window, sigmoid_weighting))
  }
  
  # get the vectors for each note
  vectors_long = data.table(notes)[,get_note_termvec_indices_finley(ReportText, term = acronym, term.vector = term_vecs_mat$V1, window = window), by='DocumentSID']
  
  # creating the doc_ N feature
  vectors_long$doc_N = paste0(vectors_long$DocumentSID, '_', vectors_long$instance)
  
  # get the unique document-instances
  unique.doc_N = unique(vectors_long$doc_N)
  
  # make the matrix for multiplication, where i is the index in the term.vector matrix and j is the doc_N
  ii = vectors_long$matid
  jj = match(vectors_long$doc_N, unique.doc_N)
  
  if(sigmoid_weighting){
    xx = sigmoid(vectors_long$dist_term)  
    # create the matrix that has the indices of the term (in the term matrix) 
    # as the row and the doc_N as the column
    mult.mat = sparseMatrix(i = ii,
                            j = jj,
                            x = xx)
  }else{
    # create the matrix that has the indices of the term (in the term matrix) 
    # as the row and the doc_N as the column
    mult.mat = sparseMatrix(i = ii,
                            j = jj,
                            x = 1)
  }
  
  # need to add 0's for extra terms not used that are in dictionary mat. Need this for matrix multiplication
  mult.mat = rbind(mult.mat, matrix(0, ncol = ncol(mult.mat), nrow = (nrow(term_vecs_mat) - nrow(mult.mat))))
  
  # create the vector averages for the term vector matrix
  doc_N.vectors = as.matrix(t(mult.mat)%*%as.matrix(term_vecs_mat[,-1]))
  
  # add rownames
  rownames(doc_N.vectors) = unique.doc_N 
  
  # return the vectors!
  return(doc_N.vectors)
}

# this function creates the longform context vector based on the cui results in data and term vector basis
get_longform_context_vector <- function(acronym_list, visit_num_cuis, notes_full, term_vecs_path, term_freq_path, term_weight_col = 'weight.alpha.mean.freq', num_samples = 1000, select_ICD_filter_pos = TRUE, weight_terms = TRUE, sigmoid_weighting = F, subset = NULL, return_type = 'acronym_list'){
  
  # setting the seed for reproducibility
  set.seed(1)
  
  # get the subset
  subset = get_subset(acronym_list, subset)
  
  # load the term vectors
  load(term_vecs_path)
  
  # load the term vector frequency table to create the weights table
  load(term_freq_path)
  
  # Create weights matrix based on frequency of terms
  weights = term_freq[,c('term',term_weight_col)]; colnames(weights) = c('term','weight')
  
  # weight the matrix
  if(weight_terms){
    term_vecs_mat = get_weighted_mat(term_vecs_mat, weights)
  } 
  
  if(select_ICD_filter_pos & !('patient' %in% colnames(visit_num_cuis))){
    # get patient list
    visit_num_cuis$patient = split_patient_date(visit_num_cuis$patient_date)
  }
  
  # create longform_vector_list
  longform_vector_list = list()
  
  # cycle through each acronym
  for(i in subset){
    print(sprintf('getting longform context vector for %s',acronym_list[[i]]$acronym))
    
    # get the longform cui and the longform name
    long.cui = acronym_list[[i]]$longform_cui
    longform = acronym_list[[i]]$longform
    
    # if only keeping filter positive patients
    if(select_ICD_filter_pos){
      visit_num_sub = visit_num_cuis[visit_num_cuis$patient %in% acronym_list[[i]]$filter_pos_patients,]
    }else{
      visit_num_sub = visit_num_cuis
    }
    
    # select indices of notes with longform
    if('cui' %in% colnames(visit_num_sub)){
      ind = which(visit_num_sub$cui==long.cui)
    }else if('V1' %in% colnames(visit_num_sub)){
      ind = which(visit_num_sub$V1==long.cui)
    }else{
      stop('cannot find the cui column in visitnums. Is it there?')
    }
    
    # if not enough notes with the longform
    if(length(ind)<num_samples){
      # print warning
      print(sprintf('not enough long forms in notes (only %s). Just using this number', length(ind)))
      
      # if there are no longforms then something is wrong
      if(length(ind)==0){return()}
    }else{
      # sample the longform instances
      ind = sample(ind,num_samples)
    }
    # get the visits to pull from SQL
    visits_to_pull = as.character(visit_num_sub$visit_num[ind])
    visit_docid = sapply(visits_to_pull, function(xx) strsplit(xx, '_')[[1]][1])

    notes = notes_full[as.character(notes_full$DocumentSID) %in% visit_docid,]
    
    if(sigmoid_weighting){
      # get the indices of all the context words for the acronym
      vectors_long = data.table(notes)[,get_note_termvec_indices_finley(ReportText, term = longform, term.vector = term_vecs_mat$V1, window = 12), by='DocumentSID']
      
      # do the sigmoid transformation
      vectors_long$sigmoid_dist = sigmoid(vectors_long$dist_term)
      
      # get the matrix of the vectors for each term in vectors_long
      tmp_mat = term_vecs_mat[vectors_long$matid,-1]
      
      # get vector weighted by sigmoid distance from term
      longform_vector = apply(tmp_mat, 2, function(xx) weighted.mean(xx, vectors_long$sigmoid_dist))
    }else{
      # get the indices of all the context words for the acronym
      vectors_long = data.table(notes)[,get_note_termvec_indices(ReportText, term = longform, term.vector = term_vecs_mat$V1),
                                       by='DocumentSID']
      
      # get the matrix of the vectors vectors for each term in vectors_long
      tmp_mat = term_vecs_mat[vectors_long$matid,-1]
      
      # create the target vector by averaging all the context vecs
      longform_vector = colMeans(tmp_mat)
    }
    
    # save the longform vector for this acronym
    acronym_list[[i]]$longform_vector = longform_vector
    longform_vector_list[[i]] = longform_vector
  }
  if(return_type == 'acronym_list'){
    return(acronym_list)
  }else if(return_type == 'longform_vector'){
    return(longform_vector_list)
  }else{
    print('select a proper return_type')
    return()
  }
}

# this creates an estimate of the longform context vector based on the ICD (silver standard label) values in the data
get_acronymICD_context_vector <- function(acronym_list, data_ICD, term_vecs_path, term_freq_path, term_weight_col = 'weight.alpha.mean.freq', patient_date_to_note_function = get_note_from_patient_dates_ORDMVP, num_samples = 200, silver_standard = 'ICD_and_NLP', subset = NULL, weight_terms = TRUE, sigmoid_weighting = F, return_type = 'acronym_list'){
  
  # set seed for reproducibility
  set.seed(1)
  
  # only keep necessary columns for data
  cols = grep('C9|C8|ICD|patient',colnames(data_ICD))
  #data_ICD = data_ICD[,cols]
  warning('not removing the extra columns like the V1 version does')
  
  # get the subset
  subset = get_subset(acronym_list, subset)
  
  # load the term vectors
  load(term_vecs_path)
  
  # load the term vector frequency table to create the weights table
  load(term_freq_path)
  
  # Create weights matrix based on frequency of terms
  weights = term_freq[,c('term',term_weight_col)]; colnames(weights) = c('term','weight')
  
  # weight the termvector matrix
  if(weight_terms){
    term_vecs_mat = get_weighted_mat(term_vecs_mat, weights)
  } 
  
  # add patient column if it doesn't exist
  if(!('patient' %in% colnames(data_ICD))){
    # get patient list
    data_ICD$patient = split_patient_date(data_ICD$patient_date)
  }
  
  # create ICD_vector_list
  ICD_vector_list = list()
  
  warning('still have a for loop around each acronym')
  # cycle through each acronym
  for(i in subset){
    print(sprintf('getting acronym context vector for %s using ICD surrogates',acronym_list[[i]]$acronym))
    
    # get the acronym name and the icd and nlp column
    acronym = acronym_list[[i]]$acronym
    acr_icd = paste0(acronym,'_ICD')
    acronym = tolower(acronym)
    acr_cui = acronym_list[[i]]$acronym_cui
    main_cui = acronym_list[[i]]$main_cui
    
    # get ICD positive and negative data sets
    if(silver_standard == 'ICD_and_NLP'){
      data_ICD_pos = data_ICD[data_ICD[,acr_icd]>0 & data_ICD[,main_cui]>0,]
      data_ICD_neg = data_ICD[data_ICD[,acr_icd]==0 & data_ICD[,main_cui]==0,]
    }else if(silver_standard == 'ICD'){
      data_ICD_pos = data_ICD[data_ICD[,acr_icd]>0,]
      data_ICD_neg = data_ICD[data_ICD[,acr_icd]==0,]
    }else if(silver_standard == 'ICD_and_filter'){
      data_ICD_pos = data_ICD[data_ICD[,acr_icd]>0,]
      data_ICD_neg = data_ICD[!(data_ICD$patient %in% acronym_list[[i]]$filter_pos_patients),]
    }else if(silver_standard == 'ICD_and_NLP_and_filter'){
      data_ICD_pos = data_ICD[data_ICD[,acr_icd]>0 & data_ICD[,main_cui]>0,]
      data_ICD_neg = data_ICD[!(data_ICD$patient %in% acronym_list[[i]]$filter_pos_patients) & data_ICD[,main_cui]==0,]
    }
    print(dim(data_ICD_pos))
    print(dim(data_ICD_neg))
    
    # create list to get the two training vectors
    tmp_list = list()
    ii = 0
    
    # cycle through ICD positive and ICD negative instances.
    for(data_name in c('data_ICD_pos','data_ICD_neg')){
      ii = ii + 1
      print(paste0('getting ',data_name))
      
      # get the data for ICD.pos or ICD.neg.
      data = get(data_name)
      
      # select indices of notes with the acronym
      ind = which(data[,acr_cui]>0)
      
      # if not enough notes with the acronym
      if(length(ind)<num_samples){
        # print warning
        print(sprintf('not enough acronyms/ICD/NLP in notes (only %s). Just using this number', length(ind)))
        
        # if there are no longforms then something is wrong
        if(length(ind)==0){next}
      }else{
        # sample the acronym instances
        ind = sample(ind,num_samples)
      }
      
      # select a sample of patient_dates with the acronym
      patient_dates = rownames(data)[ind]
      
      # pull the notes for these patient_dates
      notes = patient_date_to_note_function(patient_dates)
      #notes = get_note_from_patient_dates_ORDMVP(patient_dates)
      
      if(sigmoid_weighting){
        # get the indices of all the context words for the acronym
        vectors_long = data.table(notes)[,get_note_termvec_indices_finley(ReportText, term = acronym, term.vector = term_vecs_mat$V1, window = 12), by='DocumentSID']
        
        # do the sigmoid transformation
        vectors_long$sigmoid_dist = sigmoid(vectors_long$dist_term)
        
        # get the matrix of the vectors for each term in vectors_long
        tmp_mat = term_vecs_mat[vectors_long$matid,-1]
        
        # get vector weighted by sigmoid distance from term
        tmp_list[[ii]] = apply(tmp_mat, 2, function(xx) weighted.mean(xx, vectors_long$sigmoid_dist))
      }else{
        # get the indices of all the context words for the acronym
        vectors_long = data.table(notes)[,get_note_termvec_indices(ReportText, term = acronym, term.vector = term_vecs_mat$V1),
                                         by='DocumentSID']
        
        # get the matrix of the vectors vectors for each term in vectors_long
        tmp_mat = term_vecs_mat[vectors_long$matid,-1]
        
        # create the target vector by averaging all the context vecs
        tmp_list[[ii]] = colMeans(tmp_mat)
      }
      
    } # data_ICD_pos and data_ICD_neg
    
    # name the list for this vector and put it in the overall longform list vector.
    names(tmp_list) = c('ICD_pos_vector','ICD_neg_vector')
    
    # update the ICD vector list
    ICD_vector_list[[i]] = tmp_list
    
    # update the acronym list
    vector_names = paste0(silver_standard,'.vectors')
    
    acronym_list[[i]][[paste0(silver_standard,'_pos_vector')]] = tmp_list[[1]]
    acronym_list[[i]][[paste0(silver_standard,'_neg_vector')]] = tmp_list[[2]]
    
    #acronym_list[[i]][[vector_names]] = tmp_list
    
  } # for i in subset
  names(ICD_vector_list)[subset] = sapply(acronym_list, function(xx) xx$acronym)[subset]
  
  # return all the things
  if(return_type == 'acronym_list'){
    return(acronym_list)
  }else if(return_type == 'ICD_vector_list'){
    return(ICD_vector_list)
  }else{
    print('please select a proper return type')
    return()
  }
}

# this function finds the cosine similarity of the acronym contexts with the longforms from a given longform list (not from acronyms)
get_acronym_context_longform_similarity <- function(acronym_list, data_patient_cuis = NULL, notes = NULL, docIDs = NULL, longform_list = NULL, vector_names = NULL, term_vecs_path, term_freq_path, term_weight_col = 'weight.alpha.mean.freq', subset = NULL, sigmoid_weighting = FALSE, square_root_vec = FALSE, weight_terms = TRUE){
  # added this part on 10/11 for reproducibility
  set.seed(1)
  
  # if no longform list given (these are the target wordev vectors)
  if(is.null(longform_list) & !is.null(vector_names)){
    longform_list = lapply(acronym_list, function(xx) xx[vector_names])
  }
  
  # get the subset
  subset = get_subset(acronym_list, subset)
  
  # pull out the cuis for each acronym
  acronym_cuis <- sapply(acronym_list, function(xx) xx$acronym_cui)
  
  # return if not only one source of data
  if(is.null(data_patient_cuis) + is.null(docIDs) + is.null(notes) != 2){
    print('only one of data_patient_cuis, docIDs, or notes should be not null') 
    return()
  }
  # if data_patient_cuis are the input
  if(!is.null(data_patient_cuis)){
    # check if data_patient_cuis is too big for me right now
    if(nrow(data_patient_cuis)>10000){
      print('too many data_patient_cuis for this function to handle. It would be too slow')
      return()
    }
    
    # only keep the patient-dates that have the acronym cuis in them
    data_acronym = data_patient_cuis[data_patient_cuis$V1 %in% acronym_cuis,]
    
    # print number of unique patient-dates (# of unique notes will be ~9x higher than this)
    print(paste0('number of unique patient-dates = ',nrow(data_acronym)))
    
    # pulling all notes in the data
    notes = get_note_from_patient_dates_ORDMVP(unique(data_acronym$patient_date))
    
    # if document IDs are the input
  }else if(!is.null(docIDs)){
    notes = get_note_from_documentIDs_ORDMVP(docIDs)
  } # otherwise, notes are already created.
  
  # initial time counter
  begin0 = Sys.time()
  
  # pull in the term vector matrix
  load(term_vecs_path)
  
  # load the term vector frequency table to create the weights table
  load(term_freq_path)
  
  # Create weights matrix based on frequency of terms
  weights = term_freq[,c('term', term_weight_col)]; colnames(weights) = c('term','weight')
  
  # weight the termvector matrix
  if(weight_terms){
    term_vecs_mat = get_weighted_mat(term_vecs_mat, weights)
  }
  
  # predictions keeps track of the notes and the predictions
  predictions = NULL

  # cycle through each acronym
  for(i in subset){
    # get the acronym name
    acronym = tolower(acronym_list[[i]]$acronym)
    print(acronym)
    
    # get time for time tracking
    begin = Sys.time()
    
    # get the acronym cui used
    acronym.cui = acronym_list[[i]]$acronym_cui
    
    # if using patient_dates, just take the ones associated with the acronym of choice
    if(!is.null(data_patient_cuis)){
      # pulling out the patient dates with that acronym
      acronym_patient_dates = data_acronym$patient_date[data_acronym$V1==acronym.cui]
      
      # getting notes of all patient_dates with the acronym of interest
      notes.acronym = notes[notes$patient_date %in% acronym_patient_dates,]
    }else{
      notes.acronym = notes
    }
    
    if(sigmoid_weighting){
      # get the note vectors
      note_vectors = get_note_vectors(notes, acronym = tolower(acronym), term_vecs_mat = term_vecs_mat, sigmoid_weighting = T, window = 12)
      
      # get the cosine similarities
      if(square_root_vec){
        # get sqrt while maintaining sign
        sqrt_compress = function(xx){
          yy = sign(xx)*sqrt(abs(xx))
        }
        
        # get cosine similarities
        tmp = data.frame(sapply(longform_list[[i]], function(yy) {
          apply(note_vectors, 1, function(xx) cosine_similarity(sqrt_compress(xx), sqrt_compress(yy)))
        }))
      }else{
        tmp = data.frame(sapply(longform_list[[i]], function(yy) {
          apply(note_vectors, 1, function(xx) cosine_similarity(xx, yy))
        }))
      }
    }else{
      # get the note vectors
      note_vectors = get_note_vectors(notes.acronym, acronym = acronym, term_vecs_mat = term_vecs_mat)
      
      # get the cosine similarities
      tmp = data.frame(sapply(longform_list[[i]], function(yy) {
        apply(note_vectors, 1, function(xx) cosine_similarity(xx, yy))
      }))
    }
    
    # add in all the extra columns for analysis
    tmp$doc_N = rownames(tmp)
    tmp$DocumentSID = sapply(tmp$doc_N, function(xx) strsplit(xx,'_')[[1]][1])
    tmp$instance = sapply(tmp$doc_N, function(xx) strsplit(xx,'_')[[1]][2])
    tmp$acronym = acronym
    tmp = merge(tmp, notes[,c('DocumentSID','patient_date')], by = 'DocumentSID')
    
    # update the predictions data frame
    predictions = rbind(predictions, tmp)
    
    # print time taken for this acronym
    print(Sys.time() - begin)
    
  } # for i in 1:subset
  
  # print overall time
  print(paste0('full time taken is ',Sys.time() - begin0))
  
  # return the vector similarities
  return(predictions)
}


##### a-la-carte functions #####

# Create the a-lalcarte regression matrix and formula for regression
# finished.list
prep_alacarte <- function(finished.list){
  # gather feature data
  feature.list = lapply(finished.list, function(xx) xx$feature.vector)
  feature.data = do.call(rbind.data.frame,feature.list)
  feature.names = paste0('Vfeature_',1:(ncol(feature.data)))
  colnames(feature.data) = feature.names
  
  # gather context data
  context.list = lapply(finished.list, function(xx) xx$context.vector)
  context.data = do.call(rbind.data.frame,context.list)
  context.names = paste0('Vcontext_',1:(ncol(context.data)))
  colnames(context.data) = context.names
  
  # combine both sets
  data = cbind.data.frame(feature.data,context.data)
  
  # create formula "cbind(Y1,Y2,...) ~ X1 + X2 .."
  form.str = paste0('cbind(',
                    paste(feature.names, collapse=','),
                    ') ~ ',
                    paste(context.names, collapse=' + '))
  
  return(list(data,as.formula(form.str)))
}

# create weights for alacarte regression using log(count)
create_log_weights <- function(finished.list, term_freq, min.count=NULL){
  names = sapply(finished.list, function(xx) xx$word)
  matid = match(names,term_freq$term)
  count = term_freq$count[matid]
  weights = log(count)
  if(!is.null(min.count)){
    weights = ifelse(count>=min.count,weights,0)
  }
  return(weights)
}

# create weights for alacarte regression using a minimum count
create_cutoff_weights <- function(finished.list, term_freq, min.count){
  names = sapply(finished.list, function(xx) xx$word)
  matid = match(names,term_freq$term)
  count = term_freq$count[matid]
  return(as.integer(count>=min.count))
}

# from the linear regression output, get the A matrix and intercept
get_A_data_frame_and_intercept <- function(lm.fit){
  intercept = lm.fit$coefficients[1,]
  mat = lm.fit$coefficients[-1,] # remove (intercept)
  mat = t(mat) # transpose so we context as columns and feature vecs
  return(list(intercept,as.data.frame(mat)))
}

# create new set of embeddings given the list of contexts vectors and a-la-carte transformation
create_new_embeddings_set <- function(finished.list, A, intercept=NULL){
  names = sapply(finished.list, function(xx) xx$word)
  # get context data
  context.list = lapply(finished.list, function(xx) xx$context.vector)
  context.data = do.call(rbind,context.list)
  # Get the embeddings values
  if(is.null(intercept)){
    values = as.data.frame(context.data%*%t(as.matrix(A)))
  }else{
    int.mat = t(replicate(nrow(context.data), intercept))
    values = as.data.frame(context.data%*%t(as.matrix(A)) + int.mat)
  }
  term_vecs_mat = cbind(names,values)
  names(term_vecs_mat) = paste0('V',as.character(1:501))
  return(term_vecs_mat)
}

##### Main functions #####

# Runs the random forest CUI-ICd model. 
run_CUI_ICD_model <- function(acronym_list, data, subset = NULL, subsample_max = 1e3,
                                 method = 'RF', acronym_only_notes = TRUE,
                                 sample_proportionally = TRUE, return_type = 'acronym_list',
                                 cui_list = NULL, filter_patients = TRUE){
  
  # if patient column isn't already created 
  if(!('patient' %in% colnames(data))){
    # split patient_dates into patients
    data$patient = split_patient_date(rownames(data))
  }
  
  # setting random seed for reproducibility.
  set.seed(1)
  
  # create the indices of the subset
  subset = get_subset(acronym_list, subset)
  
  # cycle through each acronym
  for(i in subset){
    
    # get the acronym, the acronym ICD, and the acronym CUI
    acronym = acronym_list[[i]]$acronym
    acr_icd = paste0(acronym,'_ICD')
    acr_cui = acronym_list[[i]]$acronym_cui
    
    # print the acronym
    print(acronym)
    
    # only keep patients who pass the acronym filter
    if(filter_patients){
      data.sub = data[data$patient %in% acronym_list[[i]]$filter_pos_patients,]
    }else{
      data.sub = data
    }
    
    if(acronym_only_notes){
      # select only acronym notes
      data.sub = data.sub[data.sub[,acr_cui]>0,]
    }else{
      # select all notes
      data.sub = data.sub
      warning('using all notes (not just acronym ones). So the AUC comparison won\'t be fair')
    }
    
    # take a smaller subsample for running RF
    if(nrow(data.sub)>subsample_max){
      data.sub = data.sub[sample(nrow(data.sub), subsample_max),]
    }
    
    # just to see the proportion of ICDs in the notes
    print(table(data.sub[,acr_icd]))
    print(table(data.sub[,acr_icd])/nrow(data.sub))
    print(dim(data.sub))
    
    # if train and test data are not data frames, convert them to data frames.
    if(sum(class(data.sub)=='data.frame')==0){
      data.sub = as.data.frame(as.matrix(data.sub))
    }
    
    # if the cui list is not given
    if(is.null(cui_list)){
      # get the cui_list for this acronym from the acronym list
      phecode_cuis = acronym_list[[i]]$cui_list
    }else{
      # if the cui list is given, keep it
      phecode_cuis = cui_list
    }
    
    # print how many cuis are being removed because they are not in the NLP data
    print(sprintf('removing %s out of %s cuis', length(setdiff(phecode_cuis, colnames(data))), length(phecode_cuis)))
    
    # only keep cuis in the NLP data
    phecode_cuis = phecode_cuis[phecode_cuis %in% colnames(data)]
    
    # create the formula
    formula = as.formula(paste0(acr_icd,' ~ ',paste(phecode_cuis,collapse=' + ')))
    
    data.sub[,acr_icd] = as.factor(data.sub[,acr_icd])
    
    if(sample_proportionally){
      print('also sampling the final training proportionally')
      data.sub <- splitstackshape::expandRows(data.sub, acr_cui, drop = F)
    }
    
    # run random forest on the training set
    rf.model = randomForest(formula, data=data.sub)
    
    # get the auc
    rf.pred = predict(rf.model, newdata = data.sub, type='prob')[,2] 
    sprintf('training AUC = %s',AUC::auc(roc(rf.pred, data.sub[,acr_icd])))
    
    if(return_type == 'acronym_list'){
      # update the acronym_list model and AUC
      acronym_list[[i]]$ICD_model = rf.model
      acronym_list[[i]]$ICD_model_ICD_AUC = rf.AUC
      return(acronym_list)
    }else if(return_type == 'model'){
      return(rf.model)
    }else{
      print('please select a proper return type')
      return()
    }
  }
}

# clusters two models using EM-based clustering algorithm with flexmix
EM_cluster_pred_2models <- function(data, col_names, normalize = c(T,T), return_type = 'predictions', max_iter = 5){
  # normalize the columns asked for
  for(i in 1:2){
    if(normalize[i]){
      col = col_names[i]
      data[,col] <- (data[,col] - mean(data[,col]))/sd(data[,col])
      #data[,col] <- scale(data[,col])
    }
  }

  # create formula for analysis
  tmpfm <- as.formula(sprintf('%s + %s ~ 1', col_names[1], col_names[2]))
  
  # run flexmix up to five times
  iter = 1
  set.seed(1)
  tmpfit = flexmix(formula = tmpfm, data = data, k = 2,
                   model = list(FLXMRglm(family = 'gaussian'),
                                FLXMRglm(family = 'gaussian')))
  
  while((tmpfit@converged == F | length(unique(tmpfit@cluster))<2) & iter < max_iter){
    iter = iter + 1
    set.seed(iter)
    tmpfit = flexmix(formula = tmpfm, data = data, k = 2,
                     model = list(FLXMRglm(family = 'gaussian'),
                                  FLXMRglm(family = 'gaussian')))
  }
  
  # if flexmix only gives us one cluster, return
  if(tmpfit@converged == F | (iter == max_iter & length(unique(tmpfit@cluster))<2)){print('couldnt converge!');return(tmpfit)}
  print(paste0('iter = ',iter))

  #print('warning: we may not have figured out the double feature clustering issue. If the AUC is low, look at that')
  
  # get the class associated with the label
  class.pos.1 = 1 + 1*(cor(tmpfit@cluster==2, as.integer(data[,col_names[1]] > mean(data[,col_names[1]])))>0)
  class.pos.2 = 1 + 1*(cor(tmpfit@cluster==2, as.integer(data[,col_names[2]] > mean(data[,col_names[2]])))>0)
  
  if(class.pos.1 == class.pos.2){
    class.pos = class.pos.1
  }else{
    print('We have a problem. Classes dont make sense. Calling backup clustering')
           
    return()
  }
  
  # get the prediction from the flexmix class
  flexmix.prediction = tmpfit@posterior$scaled[,class.pos] 
  
  # return the prediction
  if(return_type == 'predictions'){
    return(flexmix.prediction)
    
  }else if(return_type == 'class'){
    return(as.integer(tmpfit@cluster==class.pos))
    
  }else if(return_type == 'predictions_class'){
    df = data.frame(EM_prediction = flexmix.prediction,
                    EM_class = as.integer(tmpfit@cluster==class.pos))
    return(df)
    
  }else if(return_type == 'accuracy'){
    return(mean(as.integer(tmpfit@cluster==class.pos)==predictions$label))
    
  }else if(return_type == 'object'){
    return(tmpfit)
    
    # otherwise return AUC
  }else{
    auc = AUC::auc(roc(flexmix.prediction, data$label))
    return(auc)
  }
}

# does 1D EM-clustering
EM_cluster_pred_1model <- function(data, col, normalize = TRUE, return_type = 'predictions', try_init = 'RF_prediction'){
  # normalize the columns
  if(normalize){
    data[,col] <- (data[,col] - mean(data[,col]))/sd(data[,col])
  }
  
  # create formula for analysis
  tmpfm <- as.formula(sprintf('%s ~ 1', col))
  
  ### running with random seed
  # run flexmix up to five times
  iter = 1
  tmpfit = flexmix(formula = tmpfm, data = data, k = 2,
                   model = FLXMRglm(family = 'gaussian'))
  while(length(unique(tmpfit@cluster))<2 & iter < 5){
    iter = iter + 1
    tmpfit = flexmix(formula = tmpfm, data = data, k = 2,
                     model = FLXMRglm(family = 'gaussian'))
  }
  
  # printing clustering results
  print(sprintf('random seed: variance = %s, logLik = %s', var(tmpfit@posterior$scaled[,1]), tmpfit@logLik))
  
  ### setting the initial clusters based on initial values of target column
  cluster_init = as.integer(data[,col] > 0) + 1
  tmpfit_2 = flexmix(formula = tmpfm, data = data, cluster = cluster_init,
                     model = FLXMRglm(family = 'gaussian'))
  
  # only care about it if the variance is reasonable
  if(var(tmpfit_2@posterior$scaled[,1]) > .01){
    if(var(tmpfit@posterior$scaled[,1]) < .01){
      print('replacing random flexmix with initial clusters because of variance too small in the former')
      tmpfit = tmpfit_2
    }else if(tmpfit_2@logLik > tmpfit@logLik){
      print('replacing random flexmix with initial clusters because log-likelihood is higher in the latter')
      tmpfit = tmpfit_2
    }
  }
  
  # printing clustering results
  print(sprintf('cluster-column seed: variance = %s, logLik = %s', var(tmpfit_2@posterior$scaled[,1]), tmpfit_2@logLik))
  
  ### setting initial clusters based on initial values of try_init
  if(!is.null(try_init)){
    cluster_init = as.integer(data[,try_init] > median(data[,try_init])) + 1
    tmpfit_3 = flexmix(formula = tmpfm, data = data, cluster = cluster_init,
                       model = FLXMRglm(family = 'gaussian'))
    
    # only care about it if the variance is reasonable
    if(var(tmpfit_3@posterior$scaled[,1]) > .01){
      if(var(tmpfit@posterior$scaled[,1]) < .01){
        print('replacing random flexmix with try_init clusters because of variance too small in the former')
        tmpfit = tmpfit_3
      }else if(tmpfit_3@logLik > tmpfit@logLik){
        print('replacing random flexmix with try_init clusters because log-likelihood is higher in the latter')
        tmpfit = tmpfit_3
      }
    }
    
    # printing clustering results
    print(sprintf('try_init seed: variance = %s, logLik = %s', var(tmpfit_3@posterior$scaled[,1]), tmpfit_3@logLik))
  }
  
  # if flexmix only gives us one cluster, return
  if(length(unique(tmpfit@cluster))<2){print('couldnt converge!');return()}
  
  # get the class associated with the label
  class.pos = 1 + 1*(cor(tmpfit@cluster==2, data[,col])>0)
  
  # get the prediction from the flexmix class
  flexmix.prediction = tmpfit@posterior$scaled[,class.pos]
  
  var_pred = var(flexmix.prediction)
  if(var_pred < .01){
    warning(sprintf('1D clustering variance is %s. Did not get a good spread of probabilities', var.pred))
  }
  
  # return the prediction
  if(return_type == 'predictions'){
    return(flexmix.prediction)
    
  }else if(return_type == 'class'){
    return(as.integer(tmpfit@cluster==class.pos))
    
  }else if(return_type == 'object'){
    return(tmpfit)
    
    # otherwise return AUC
  }else{
    auc = AUC::auc(roc(flexmix.prediction, data$label))
    return(auc)
  }
}

# function that runs all the components of CASEml
# to include EM averaging and 1D clustering - now using v2 of 1D clustering
CASEml_wrapper_v7 <- function(acronym_list, train_data, test_data, notes, subset = NULL, filter_pos = T, include_full_test = F, max_EM_iter = 5, cosine_transform = 'logit', RF_transform = 'none', sigmoid_weighting = F, square_root_vec = F,term_vecs_path, term_freq_path, term_weight_col = 'weight.alpha.mean.freq', vector_names = c("ICD_pos_vector", "ICD_neg_vector", "longform_vector"), add_prev_classification = T, return_all_models = F){
  
  # get the subset numerically
  subset = get_subset(acronym_list, subset)
  
  # initialize results and iter counter
  res = list()
  iter = 1
  
  for(i in subset){
    # get the acronym
    acronym = acronym_list[[i]]$acronym
    
    # for training and test, only select relevant columns
    train_sub = train_data[,which(colnames(train_data) %in% c(unlist(acronym_list[[i]][c('acronym_cui','longform_cui','cui_list')]),
                                                              paste0(acronym, '_ICD')))]
    test_sub = test_data[,which(colnames(test_data) %in% c(unlist(acronym_list[[i]][c('acronym_cui','longform_cui','cui_list')]),
                                                           paste0(acronym, '_ICD')))]
    
    # only keep data with the cui in it
    train_sub = train_sub[train_sub[,acronym_list[[i]]$acronym_cui] > 0,]
    test_sub = test_sub[test_sub[,acronym_list[[i]]$acronym_cui] > 0,]
    
    # if only filter positive patients, only keep those
    if(filter_pos){
      # get the training patient list
      train_patients = split_patient_date(rownames(train_sub))
      
      # subset by filter patients
      train_sub = train_sub[which(train_patients %in% acronym_list[[i]]$filter_pos_patients),]
      
      # if excluding test set too
      if(!include_full_test){
        test_patients = split_patient_date(rownames(test_sub))
        test_sub = test_sub[which(test_patients %in% acronym_list[[i]]$filter_pos_patients),]
      }
    }
    
    # create a data frame for train and test data
    train_sub = as.data.frame(as.matrix(train_sub))
    test_sub = as.data.frame(as.matrix(test_sub))
    
    # run the RF models
    RF_model = run_CUI_ICD_model(acronym_list, train_sub, cui_list = NULL, subset = i, return_type = 'model', filter_patients = filter_pos, subsample_max = 1e4)
    
    # predict on the the test data
    RF_pred = predict(RF_model, newdata = test_sub, type='prob')#[,2] 
    
    # convert to the format I want
    RF_pred = data.frame(patient_date = rownames(RF_pred), RF_prediction = RF_pred[,2], stringsAsFactors = F)
    
    # grab them from the test data
    if(is.null(notes)){
      # give an error if it would be caused
      if(filter_pos){
        stop('This function doesnt deal with null notes and filter positive')
      }
      
      context_pred = NULL
      # cycle through batches of 10000 because otherwise this would crash R
      for(i in 1:ceiling(nrow(test_sub)/10000)){
        print(sprintf('reading in batch %s of 10000 notes', i))
        # get the indices of notes to pull
        ind = ((i - 1)*10000 + 1):(min(nrow(test_sub), i*10000))
        
        # get the notes for this batch
        notes_subset = get_note_from_patient_dates_ORDCho(rownames(test_sub)[ind])
        
        # get the context predictions for this note set
        tmp_pred = get_acronym_context_longform_similarity(acronym_list, notes = notes_subset, vector_names = vector_names, term_vecs_path = term_vecs_path, term_freq_path = term_freq_path, term_weight_col = term_weight_col, sigmoid_weighting = sigmoid_weighting, square_root_vec = square_root_vec, subset = acronym, weight_terms = TRUE)
        
        # update the context prediction
        context_pred = rbind(context_pred, tmp_pred)
      }
    }else{
      # subset the notes to either be filter positive only or not
      if(filter_pos & (!include_full_test)){
        notes_subset <- notes[notes$VINCI_ID %in% acronym_list[[i]]$filter_pos_patients,]
        # mismatches happen when NILE acronym counter doesn't match my acronym scanner
      }else{
        notes_subset <- notes
      }
      
      # only keep notes that are in the test set
      notes_subset = notes_subset[notes_subset$patient_date %in% RF_pred$patient_date,]
      
      # run the alacarte-ICD model!
      system.time({
        context_pred = get_acronym_context_longform_similarity(acronym_list, notes = notes_subset, vector_names = vector_names, term_vecs_path = term_vecs_path, term_freq_path = term_freq_path, term_weight_col = term_weight_col, sigmoid_weighting = sigmoid_weighting, square_root_vec = square_root_vec, subset = acronym, weight_terms = TRUE)
      }) # 20s for 4k notes
    }
    
    # create the alacarte_pred
    if(cosine_transform == 'logit'){
      print('doing logit transform to cosine similarities')
      context_pred$ICD_pos_minus_neg = sapply(context_pred$ICD_pos_vector, logit_neg1_to_1) - 
        sapply(context_pred$ICD_neg_vector, logit_neg1_to_1)
      
      context_pred$LF_minus_neg = sapply(context_pred$longform_vector, logit_neg1_to_1) - 
        sapply(context_pred$ICD_neg_vector, logit_neg1_to_1)
    }else if(cosine_transform == 'exp'){
      print('doing exp transform to cosine similarities')
      context_pred$ICD_pos_minus_neg = exp(context_pred$ICD_pos_vector) - exp(context_pred$ICD_neg_vector)
      context_pred$LF_minus_neg = exp(context_pred$longform_vector) - exp(context_pred$ICD_neg_vector)
    }else if(cosine_transform == 'logit_test'){
      print('doing logit TEST transform to cosine similarities')
      context_pred$ICD_pos_minus_neg = sapply(context_pred$ICD_pos_vector, logit_neg1_to_1) - 
        sapply(context_pred$ICD_neg_vector, logit_0_to_1)
      
      context_pred$LF_minus_neg = sapply(context_pred$longform_vector, logit_neg1_to_1) - 
        sapply(context_pred$ICD_neg_vector, logit_0_to_1)
    }else{
      # create the alacarte_pred
      context_pred$ICD_pos_minus_neg = context_pred$ICD_pos_vector - context_pred$ICD_neg_vector
      context_pred$LF_minus_neg = context_pred$longform_vector - context_pred$ICD_neg_vector
    }
    
    if(RF_transform == 'logit'){
      print('doing logit transform to RF pred')
      RF_pred$RF_prediction = sapply(RF_pred$RF_prediction, logit_0_to_1)
    }else if(RF_transform == 'exp'){
      print('doing exp transform to RF pred')
      RF_pred$RF_prediction = exp(RF_pred$RF_prediction)
    }
    
    context_pred = context_pred[context_pred$patient_date %in% RF_pred$patient_date,]
    
    # merge the different predictions
    predictions = merge(context_pred, RF_pred, by = 'patient_date')
    
    # cluster using flexmix EM algorithm
    flexmix_pred_1 = EM_cluster_pred_2models(predictions, col_names = c('RF_prediction','ICD_pos_minus_neg'),
                                                return_type = 'predictions_class', normalize = c(T,T), max_iter = max_EM_iter)
    flexmix_pred_2 = EM_cluster_pred_2models(predictions, col_names = c('RF_prediction','LF_minus_neg'),
                                                return_type = 'predictions_class', normalize = c(T,T), max_iter = max_EM_iter)
    
    # if this doesn't work, try not normalizing
    if(class(flexmix_pred_1) == 'flexmix'){
      flexmix_pred_1 = EM_cluster_pred_2models(predictions, col_names = c('RF_prediction','ICD_pos_minus_neg'),
                                                  return_type = 'predictions_class', normalize = c(F,F), max_iter = max_EM_iter)
      if(class(flexmix_pred_1) == 'flexmix'){
        print('error in flexmix 1')
      }else{
        print('ran flexmix 1 without normalizing')
      }
    }
    if(class(flexmix_pred_2) == 'flexmix'){
      flexmix_pred_2 = EM_cluster_pred_2models(predictions, col_names = c('RF_prediction','LF_minus_neg'),
                                                  return_type = 'predictions_class', normalize = c(T,F), max_iter = max_EM_iter)
      if(class(flexmix_pred_2) == 'flexmix'){
        print('error in flexmix 2')
      }else{
        print('ran flexmix 2 without normalizing LF')
      }
    }

    # add cluster results to predictions
    predictions$EM_prediction_ICD = flexmix_pred_1$EM_prediction
    predictions$EM_class_ICD = flexmix_pred_1$EM_class
    
    # add cluster results to predictions
    predictions$EM_prediction_LF = flexmix_pred_2$EM_prediction
    predictions$EM_class_LF = flexmix_pred_2$EM_class
    
    # do the 1D clustering
    predictions$LF_cluster1 = EM_cluster_pred_1model(predictions, 'LF_minus_neg')
    predictions$RF_cluster1 = EM_cluster_pred_1model(predictions, 'RF_prediction')
    
    # bring together the 1D clustering
    predictions$CASEml = (predictions$LF_cluster1 + predictions$RF_prediction)/2
    predictions$EM_average = (predictions$LF_cluster1 + predictions$RF_cluster1)/2

    # do the 1D clustering
    predictions$LF_cluster1 = EM_cluster_pred_1model(predictions, 'LF_minus_neg')
    predictions$CASEml = (predictions$LF_cluster1 + predictions$RF_prediction)/2
    
    if(return_all_models == F){
      cols_to_include = c('LF_cluster1','RF_prediction','CASEml')
    }else{
      cols_to_include = c('EM_prediction_ICD','EM_prediction_LF','ICD_pos_minus_neg','LF_minus_neg','RF_prediction','EM_average','CASEml')
    }
    predictions = predictions[,c('patient_date','DocumentSID','instance','doc_N','acronym',cols_to_include)]


    if(add_prev_classification){
      # get the classifications for all predictions
      tmp <- classify_by_prev_wrapper(acronym_list[[i]]$prevalences, predictions, filter_pos = filter_pos, cols = cols_to_include)
      
      # add classifications to predictions
      predictions = cbind(predictions, tmp)
      
    }
    # save!
    res[[iter]] = predictions
    names(res)[iter] = acronym
    iter = iter + 1
  }
  
  return(res)
}


# classifies predictions to match a specified prevalence
classify_by_prev <- function(prediction, prevalence, verbose = T){
  # sort the prediction
  pred_sort = sort(prediction)
  
  # create the cutoff by finding the (1-prevalence)% entry of the predictions
  cutoff = pred_sort[ceiling(length(pred_sort)*(1-prevalence))]
  
  # make the classifications!
  class = as.integer(prediction >= cutoff)
  
  # print the cutoff
  if(verbose){
    print(sprintf('for prevalence %s: cutoff = %s', prevalence, cutoff))
  }
  
  # return it!
  return(class)
}

# a wrapper that classifies a set of columns by a list of prevalences
classify_by_prev_wrapper <- function(prevalences, predictions, filter_pos,
                                     cols = c('EM_prediction_ICD','EM_prediction_LF','ICD_pos_minus_neg','LF_minus_neg','RF_prediction')){
  
  # initialize data frame to store the results
  df = NULL
  
  # cycle through columns
  for(col in cols){
    print(col)
    # use different prevalences for the filter positive and nofilter sets
    if(filter_pos){
      # get the three classifications
      col_class_1 <- classify_by_prev(prediction = predictions[,col], 
                                      prevalence = prevalences$filterpos_true)
      col_class_2 <- classify_by_prev(prediction = predictions[,col], 
                                      prevalence = prevalences$filterpos_ICD)
      col_class_3 <- classify_by_prev(prediction = predictions[,col], 
                                      prevalence = prevalences$filterpos_ICD_CUI)
    }else{
      # get the three classifications
      col_class_1 <- classify_by_prev(prediction = predictions[,col], 
                                      prevalence = prevalences$nofilter_true)
      col_class_2 <- classify_by_prev(prediction = predictions[,col], 
                                      prevalence = prevalences$nofilter_ICD)
      col_class_3 <- classify_by_prev(prediction = predictions[,col], 
                                      prevalence = prevalences$nofilter_ICD_CUI)
    }
    
    # create the column names for the returned data frame
    col_names = c(paste0(col,'_class_prevTrue'), paste0(col,'_class_prevICD'), paste0(col,'_class_prevICDNLP'))
    
    # create the data frame of classifications
    tmp = data.frame(col_class_1, col_class_2, col_class_3)
    colnames(tmp) = col_names
    
    if(is.null(df)){
      # set the data frame equal to this run
      df = tmp
    }else{
      # update the full data frame
      df = cbind(df, tmp)
    }
  }
  # return the data frame!
  return(df)
}  

##### Performance metric functions  #####

# gets the performance of the predictions
get_prediction_performance <- function(acronym_list, predictions, charts, filter_pos = T, prediction_cols = c('ICD_pos_minus_neg','LF_minus_neg','RF_prediction','EM_prediction_ICD','EM_prediction_LF', 'EM_average','method_1'), class_cols = NULL, boot_n = 10000){
  # get the set of acronyms
  acronyms = intersect(names(predictions), names(charts))
  
  # get the class columns if not provided
  if(is.null(class_cols)){
    class_cols = grep('class_prevICDNLP', colnames(predictions[[1]]), value = T)
  }
  
  # check that the acronyms exist in each set
  if(length(predictions) != length(charts) | length(acronyms) != length(predictions)){
    print('name mismatch between predictions and charts. Continuing anyway')
  }
  
  # initialize the return list
  return_list = NULL
  
  for(acronym in acronyms){
    # initialize return data frame
    fscore_df = NULL
    accuracy_df = NULL
    auc_df = NULL
    
    # get chart review data
    chart_sub = charts[[acronym]]
    
    # create the doc_N merging column
    chart_sub$doc_N = paste0(chart_sub$DocumentSID, '_', chart_sub$instance)
    
    # only keep filter_pos patients
    if(filter_pos){
      chart_sub = chart_sub[chart_sub$filter == 'YES',]
    }
    
    # only select the columns we care about
    chart_sub = chart_sub[,c('DocumentSID', 'instance' ,'doc_N','filter','label')]
    
    # get prediction data
    pred_sub = predictions[[acronym]]
    
    # only keep predictions on the same documents
    pred_sub = pred_sub[pred_sub$DocumentSID %in% chart_sub$DocumentSID,]
    
    if(length(intersect(pred_sub$doc_N, chart_sub$doc_N)) != length(pred_sub$doc_N)){
      print('doc_N mismatch with charts')
      browser()
    }
    
    # only keep relevant columns
    pred_sub = pred_sub[,c('doc_N', prediction_cols, class_cols)]
    
    # merge the predictions and labels
    res = merge(pred_sub, chart_sub, by = 'doc_N')
    
    # get the filter/non-filter weights. This is only used if applying to non-filter positive patients
    res$weight = ifelse(res$filter == 'YES', 1,
                        acronym_list[[acronym]]$non_filter_ratio)
    
    # print acronym characteristics
    print(sprintf('acronym: %s', acronym))
    print(sprintf('N documents = %s', length(unique(res$DocumentSID))))
    print(sprintf('N instances of acronym = %s', length(res$doc_N)))
    
    roc_list = list()
    iter = 1
    
    # print all the results
    if(filter_pos){
      
      for(col in prediction_cols){
        # get the auc scores and print the point estimate
        auc_scores = bootstrap_auc(res[,col], res$label, weight = rep(1, nrow(res)), boots = boot_n)
        
        auc_df = rbind(auc_df, data.frame(method = col, AUC = auc_scores[2], AUC_L = auc_scores[1], AUC_U = auc_scores[3]))
        print(sprintf('AUC %s: %s', col, auc_scores[2]))
      }
      
      for(col in class_cols){
        # accuracy
        accuracy = mean(res[,col] == res$label)
        accuracy_df = rbind(accuracy_df, data.frame(method = col, accuracy = accuracy))
        
        # f score
        tp = sum(res[,col] == 1 & res$label == 1)
        fp = sum(res[,col] == 1 & res$label == 0)
        fn = sum(res[,col] == 0 & res$label == 1)
        
        f1 = tp/(tp + 0.5*(fp + fn))
        
        fscore_df = rbind(fscore_df, data.frame(method = col, fscore = f1))
        
        # printing results
        print(sprintf('accuracy of %s = %s; f-score = %s', col, accuracy, f1))
        
      }
      MFS = as.integer(names(which.max(table(chart_sub$label))))
      # accuracy
      accuracy_MFS = mean(res$label==MFS)
      accuracy_df = rbind(accuracy_df, data.frame(method = 'MFS', accuracy = accuracy_MFS))
      
      # f score
      fscore_MFS = get_weighted_fscore(MFS, res$label, weight = rep(1, nrow(res)))
      fscore_df = rbind(fscore_df, data.frame(method = 'MFS', fscore = fscore_MFS))
      
      # printing output
      print(sprintf('accuracy of MFS (%s) = %s; f-score is %s', MFS, accuracy_MFS, fscore_MFS))
      print('-----')
    }else{
      # Need to use the get_weighted_auc function. Test it to be normal.
      
      for(col in prediction_cols){
        auc_scores = bootstrap_auc(res[,col], res$label, weight = res$weight, boots = boot_n)
        
        auc_df = rbind(auc_df, data.frame(method = col, AUC = auc_scores[2], AUC_L = auc_scores[1], AUC_U = auc_scores[3]))
        print(sprintf('weighted AUC median %s: %s', col, auc_scores[2]))
        
      }
      for(col in class_cols){
        # accuracy
        accuracy = get_weighted_accuracy(res[,col], res$label, res$weight)
        accuracy_df = rbind(accuracy_df, data.frame(method = col, accuracy = accuracy))
        
        # f score
        f1 = get_weighted_fscore(res[,col], res$label, res$weight)
        fscore_df = rbind(fscore_df, data.frame(method = col, fscore = f1))
        
        print(sprintf('accuracy of %s = %s; f-score = %s', col, accuracy, f1))
      }
      MFS = as.integer(names(which.max(table(chart_sub$label))))
      # accuracy
      accuracy_MFS = get_weighted_accuracy(MFS, res$label, res$weight)
      accuracy_df = rbind(accuracy_df, data.frame(method = 'MFS', accuracy = accuracy_MFS))
      
      # f score
      fscore_MFS = get_weighted_fscore(MFS, res$label, res$weight)
      fscore_df = rbind(fscore_df, data.frame(method = 'MFS', fscore = fscore_MFS))
      
      # printing output
      print(sprintf('weighted accuracy of MFS (%s) = %s; weighted f-score = %s', MFS, accuracy_MFS, fscore_MFS))
      print('-----')
    }
    tmp_list = list(accuracy_df, auc_df, fscore_df)
    names(tmp_list) = c('accuracy', 'auc', 'fscore')
    return_list[[acronym]] = tmp_list
  }
  return(return_list)
}

get_weighted_auc <- function(prediction, label, weight, return_type = 'auc'){
  # convert label to a 0,1 integer
  if(class(label) == 'factor'){
    label = as.integer(label) - 1
  }
  
  if(length(unique(label)) != 2 | length(setdiff(label, c(0,1))) >0){
    stop('non-unique labels. Please check them')
  }

  junk=ROC.Est.FUN(label,prediction,yy0=0.5,fpr0=seq(0,1,0.01),
                   wgti=weight,yes.smooth=F)
  
  if(return_type == 'auc'){
    return(junk[1])
  }else if(return_type == 'roc'){
    roc=as.data.frame(matrix(junk[-1],ncol=6))
    colnames(roc) = c("cutoff","pos.rate","FPR","TPR","PPV","NPV")
    return(roc)
  }
  
}

get_weighted_accuracy <- function(class, label, weight){
  return(sum(as.integer(class == label) * weight)/sum(weight))
}

get_weighted_fscore <- function(class, label, weight){
    # f score
    tp = sum(as.integer(class == 1 & label == 1)*weight)
    fp = sum(as.integer(class == 1 & label == 0)*weight)
    fn = sum(as.integer(class == 0 & label == 1)*weight)
    
    f1 = tp/(tp + 0.5*(fp + fn))
    # weighted.ppv = sum(as.integer(class == 1 & label == 1)*weight)/sum(as.integer(class == 1)*weight)
    # weighted.sens = sum(as.integer(class == 1 & label == 1)*weight)/sum(as.integer(label == 1)*weight)
    # f1 = 2*weighted.ppv*weighted.sens/(weighted.ppv + weighted.sens)

  return(f1)
}

ROC.Est.FUN <- function(Di,yyi,yy0,fpr0=NULL,wgti=NULL,yes.smooth=F)
{
  out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- NULL
  if(is.null(wgti)){wgti=rep(1,length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));  
  mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0  
  for(k in 1:pp)
  {
    yy = yy0; 
    if(!is.null(fpr0)){
      tpr.all = S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
      fpr.all = S.FUN(yyi[,k],Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth);
      TPR = approx(c(0,fpr.all,1),c(0,tpr.all,1),fpr0,method="linear",rule=2)$y; 
      TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth), TPR); 
      yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))           
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }else{
      TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }
    out.yy = cbind(out.yy, yy)
    out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth))
    out.TPR = cbind(out.TPR,  TPR);  out.FPR  <- cbind(out.FPR,  FPR)
    PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
    out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
    #AUC <- sum((sum.I(yyi[,k],"<=",Yi=yyi[,k],Vi=Di*wgti)+sum.I(yyi[,k],"<",Yi=yyi[,k],Vi=Di*wgti))*(1-Di)*wgti/2
    #             )/(sum((1-Di)*wgti)*sum(Di*wgti))
    AUC = sum(S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
    out.AUC <- c(out.AUC, AUC)
  }
  out = c(out.AUC,out.yy,out.pp,out.FPR,out.TPR,out.PPV,out.NPV)
  out
}

S.FUN <- function(yy,Yi,Di,yes.smooth=F)
{
  if(yes.smooth){
    Y1i = Yi[Di==1]; n1 = sum(Di); bw = bw.nrd(Y1i)/n1^0.6
    c(t(rep(1/n1,n1))%*%pnorm((Y1i-VTM(yy,n1))/bw))
  }else{
    return((sum.I(yy,"<",Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  }
  ##sum.I(yy,"<=",Yi,Vi=Di)/sum(Di)
}

sum.I <- function(yy,FUN,Yi,Vi=NULL)
  ## sum_i I(yy FUN Yi)Vi
  # Vi weight
{
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  # for each distinct ordered failure time t[j], number of Xi < t[j]
  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')    
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
  if (!is.null(Vi)) {
    ## if FUN contains '=', tmpind is the order of decending
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}

Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F)
{
  yy0<-unique(sort(Yi,decreasing=T)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth) 
  return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
}

bootstrap_auc <- function(prediction, label, weight, boots = 2000, return_type = 'quantiles'){
  set.seed(1)
  raw_auc = get_weighted_auc(prediction, label, weight)
  
  # length of data
  n = length(label)
  auc_vec = c()
  
  # cycle through each bootstrap and calculate AUC
  for(i in 1:boots){
    ind = sample(n, n, replace = T)
    auc_vec = c(auc_vec, get_weighted_auc(prediction[ind], label[ind], weight[ind]))
  }
  
  # # takes the same time
  # auc_vec_2 = sapply(1:boots, function(xx){
  #   ind = sample(n, n, replace = T)
  #   get_weighted_auc(prediction[ind], label[ind], weight[ind])
  # })
  
  if(return_type == 'all'){
    return(auc_vec)
  }else{
    quants = quantile(auc_vec, probs = c(0.025, 0.5, 0.975))
    return(quants)
  }
}
