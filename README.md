# CASEml
Here you will find the code for "Binary Acronym Disambiguation in Clinical Notes from Electronic Health Records with an Application in Computational Phenotyping" by Link et. al. This code runs the CASEml algorithm, which uses visit-level information as well as word context to classify the meaning of acronyms in clinical text. This code could easily be adapted to work with different acronyms or different text.

Due to privacy concerns the data used in this study cannot be shared. Also, due to the specific nature of the EHR data in the VA and it being different from other facilities, this github does not contain the the data wrangling and cleaning code used in this project. However, some of the functions used in the data prep stage are still included, such as the functions to perform a-la-carte and create the longform vectors for the acronyms

To main workhorse function for CASEml_wrapper is "XXX". The inputs to the function are as follows. If the data is prepared properly for these inputs then the model should run smoothly:

Stopped at "Main Functions"

