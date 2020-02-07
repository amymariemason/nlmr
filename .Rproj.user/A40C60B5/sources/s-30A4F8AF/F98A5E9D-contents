

###################################################
# Turn table of snps into qctool ready form 

qctoolify <- function(table, chr="chr", pos="pos", filename){
  #table is the table of chromosomes and positions you want to extract
  #chr is the name of the column containing chromosone information (integrer)
  #pos is the name of the column containing position information
  #filename is where the file should be saved
  ###################
  # this function creates a space seperated list of snps ready to use in qctool
  table_temp<-table
  names(table)[names(table) == chr] <- "chr"
  names(table)[names(table) == pos] <- "pos"
  table$output<-ifelse(as.numeric(as.character(table$chr))<10, paste0("0",table$chr,":", table$pos), paste0(table$chr,":", table$pos))
  table2<-t(table$output)
  write.table(table2,filename,sep=" ",row.names=FALSE, quote =FALSE, col.names=FALSE)
  table <-table_temp
  rm(table_temp)
}

# use this to rename columns if misbehaving
library(tidyverse)
alcohol_snps<-rename(alcohol_snps, pos=Position )
alcohol_snps<-rename(alcohol_snps, chr=Chr)

##########################################
# Functions
#########################################

# create function to loop round variables


create_blank_file<- function(inputfile){
  
  #create empty data frame for variant data; includes the column headings from inputfile and an additional column 
  #called 'duplicate' 
  output<-fread(file=inputfile,nrows=1)
  output$duplicate<-0
  output<- output[0,]
  
  return(output)
}

#This function creates a data frame with the variable 'varname' from varlist - this will be 'pos' from snplist  
#in our example. It also adds an additional column called 'missing' where each variant has a value of zero.  
#Hence, will provide information on variants that are in snplist but do not appear in the statsfile.   

create_missing_report<-function(varlist, varname=names(varlist)[[1]]){
  missinglist<-as.data.frame(varlist[,varname])
  names(missinglist)<-varname
  missinglist$missing<-rep(0, nrow(missinglist)) 
  return(missinglist)
}


# this changes the proposed column types by fread by overriding them
# varlist = empty file created for the variant data
# varname = vector of variables for which you wish to change class; selects all if not specifies
# coltype = vector of new classes; by default this overrides into character

change_col_character<-function(varlist, varchange=names(varlist), coltype=rep("character", length(varchange))) {
  current<- sapply(varlist, class)
  whichcols<-which(colnames(varlist)%in%varchange)
  current[whichcols]<-coltype
  return(current)
}

# IMPROVEMENT NOTE: this in many cases would be better implemented with grep in bash. 
# IMPROVEMENT NOTE: look into iterator package https://www3.nd.edu/~steve/computing_with_data/23_data_import/data_import.html
# IMPROVEMENT NOTE: currently ADDING SPACES AROUND THE OUTSIDE SOMETIMES PREVENTS MATCHES -EPIC FILES
        #ADD OPTION TO MATCH " "+j+" " or j+"_" (currently matches latter)
# imports data into blank file, reports missing data in list
import_data<-function(inputfile, inputname, varlist, varname=names(varlist)[[1]], maxdup=5, col_change=FALSE, varchange=NULL, coltype=NULL){
  #inputfile = file to subset
  #inputname = variable in input file to search for variable
  #varlist = file containing variables to search for
  # (varname) = variable to search; defaults to first column of varlist
  # col_change = TRUE -> changes variable types to character
  #varchange = vector of variable names to restrict col_change to; default is NULL
  #coltype = vector of classes to change the column to; default is NULL
  
  # adding coltype will impose those classes as long as they are not more restrictive than the default
  # see fread for more details on colClasses
  # (maxdup) = maximum number of duplications to search for; default = 5; if some allele identifying as logical instead of character
  # try increase max dup or use col_change to hard set  
  # Duplications of the position number may  
  #occur when there has been more than one mutation i.e. G to T and G to C etc. 
  

  require(data.table)
  # create blank report files
  missingreport<-create_missing_report(varlist,varname)
  output<-create_blank_file(inputfile)  
  num_er=0 # number of errors reported
  
  # changes character types if needed
  char_read<-fread(file=inputfile,nrows=1)
  class_change <- sapply(char_read,class)
  if(col_change==TRUE){
    class_change<-change_col_character(char_read,varchange, coltype)}
  names(class_change)<-NULL
  
  # loop to read in data from outside file. The fread() function is used to check whether the variant in snplist is  
  #contained in the statsfile. To allow the code to run when the variant is missing, the tryCatch() function is used.  
  #One of two things will happen: a) it doesn't find "j", in which case it will report an error and the variant will be  
  #reported as missing; or b) if the variant is found then the fread() function will be executed and the data subseted 
  #and saved as test. 
  
  num_er =0 # error counter
  for (j in varlist[,varname]){
    # to watch progress
    print(j)
    #add spaces to ensure exact match: Note this may still pull wrong line IF value matches another var
    j1<- paste0(" ",j," ")
    j2<-paste0(j,"_")
    # clear test dummy var
    #  if (exists("test")){rm(test)}
    # test if input exists for this value; reports error if missing
    test <- tryCatch(
      # this is what I want it to do, i.e. search for position "j" and then extract maxdup rows (e.g. 5). fread() will 
      #only extract data for the first instance of this, hence we are assuming the position is in some numerical order 
      #as otherwise we may miss duplicate entires of the same variant.  
      fread(file=inputfile,nrows=maxdup, skip=as.character(j2), colClasses=class_change)
      ,
      # if error occurs
      error=function(error_message) {
        message(paste0("Error AT VALUE: ", j))
        message(error_message)
        return("ERROR")
      }
    )
    # if error, report variant as not found
    #Hence, if the variant is not found in the inputfile then report as missing by replacing with 1  
    
    if(is.character(test)){
      #Look for j in the varname and replacing missing column with 1 if not there 
      
      missingreport[grep(j, missingreport[, varname]),"missing"]<-1;
      #If there the variant is not found increase the error counting variable num_er
      num_er <- num_er +1;
    }
    # if no error, add collected lines to file, reporting duplicate matches
    if(!(is.character(test))) {
      test$duplicate<-rep(0, nrow(test));
      #Takes the col names from output (i.e. the input file) and puts the names onto test. This will only work if the 
      #variant has been found and the relevant row(s) has been extracted by the fread() function. 
      names(test)<-names(output);
      test<-as.data.frame(test)
      testSubset<-test[which(test[,inputname] == j),]
      #If there is more than one row for the position then replace the duplicate entry with a 1       
      if(nrow(testSubset)!=1) testSubset$duplicate<-rep(1, nrow(testSubset));
      output <- rbind(output, testSubset)
    }
  }  
  #Return the output, number of variants not found, and information on missing variants. 
  outlist<-list(output, num_er, missingreport)
  return(outlist)
}


#############################

# adds binary outcome column to match ordering of the .sample file
sample_all<- #ordered samplefile here
wantedoutcomes<- # list of outcomes
outcomefile<- #outcomefile here  
output<-rep(NA, nrow(sample_all))
output<-as.data.frame(output)

extract_outcome<-function(sample_list, sample_id="ID_1", outcomefile, outcome, outcome_id="eid"){
  stopifnot(!is.null(outcomefile[,outcome]))
  stopifnot(!is.null(outcomefile[,outcome_id]))
  stopifnot(!is.null(sample_list[,sample_id]))
  events<-which(outcomefile[,outcome]=="1")
  cases = unique(outcomefile[events,][[outcome_id]])
  output_temp = ifelse(as.numeric(sample_list[2:nrow(sample_list),sample_id])%in%cases, "1", "0")
  NAevents<-which(is.na(outcomefile[,outcome]))
  missing = unique(outcomefile[NAevents,][[outcome_id]])
  output_temp = ifelse(as.numeric(sample_list[2:nrow(sample_list),sample_id])%in%missing, "NA", output_temp)
  output_temp<-as.data.frame(output_temp)
  output_temp<-rbind("B", output_temp)
  names(output_temp)<-outcome
  return(output_temp)
}

# apply in loop
allframes = lapply(wanted_outcomes,
                   function(x)extract_outcome(sample_list=sample_all, 
                                              sample_id="ID_1", 
                                              outcomefile=outcomes, 
                                              outcome=x, 
                                              outcome_id="n_eid"))
