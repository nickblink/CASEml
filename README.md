# CASEml
Here you will find the code for "Binary Acronym Disambiguation in Clinical Notes from Electronic Health Records with an Application in Computational Phenotyping" by Link et. al. This code runs the CASEml algorithm, which uses visit-level information as well as word context to classify the meaning of acronyms in clinical text. This code could easily be adapted to work with different acronyms or different text.

Due to privacy concerns the data used in this study cannot be shared. Also, due to the specific nature of the EHR data in the VA and it being different from other facilities, this github does not contain the the data wrangling and cleaning code used in this project. However, some of the functions used in the data prep stage are still included, such as the functions to perform a-la-carte and create the longform vectors for the acronyms

The two files containing code are 

1) Applying_alacarte: Contains the code to run alacarte, a process to de-noise word vectors (refer to https://arxiv.org/abs/1805.05388)
2) CASEml_functions: Contains all the functions used for running CASEml. There are a lot of helper functions in here, but the main workhorse functions are:
  a) "CASEml_wrapper": which runs the CASEml algorithm
  b) "run_CUI_ICD_model": which runs the random forest ICD model using visit-level information
  c) "get_acronym_context_longform_similarity": which runs the wordvec-score model using contextual information of the acronym to predict its meaning.

Data inputs to these functions:
1) term_vecs_path: path to an RData structure "term_vecs_mat" that is a data.frame n x (d + 1), containing the n word embeddings of length d. The first column is the term and the remaining columns are the embeddings.
2) term_freq_path: path to an RData structure "term_freq" that contains the frequency of terms, such as inverse document frequency. The first column is the term and the remaining columns are different frequency measures.
3) An important data structure used throughout this code is the acronym_list, which is a list with an entry for each acronym being disambiguated (RA, MS, and MI in the referenced article). Within each, acronym list item is a sub-list with the following values

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp  3a) acronym: the acronym, e.g. "RA"
  3b) longform: the target longform, e.g. "rheumatoid arthritis"
  3c) main cui: the main CUI associated with the longform, i.e. C0003873
  3d) phecode: the phecode (from the Josh Denny mapping) associated with the main CUI, i.e. 714_1
  3e) phecode_list: the list of phecodes associated with the acronym used in the RF-CUI-ICD algorithm
  3f) acronym_cui: the cui code used to indicate acronym counts in the dataset
  3g) longform_cui: the cui code used to indicate longform counts in the dataset
  3h) icd_codes: the list of ICD codes associated with the longform
  3i) filter_pos_patients: the list of patient identifiers that pass the filter for the study (this is specific to the Link et. al paper)
  3j) ICD_pos_vector: the embedding for the target sense created by using positive ICD silver-standard labels
  3k) ICD_neg_vector: the embedding for NOT the  target sense created by using negative ICD silver-standard labels
  3l) longform_vector: the embedding for the  target sense using the longform contexts
  3m) non_filter_ratio: the ratio of patients not passing the filter to ones who do pass (used for weighting results, again this is likely specific to the Link et. al paper)
  3n) prevalences: list of the estimated prevalences of the acronym target sense class



