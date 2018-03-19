# !/usr/bin/Rscript

## Run analysis for MaBoSS output file probtraj.csv
## Find the trajectories for the specified conditions
## Store the trajectories in a file opt$outprefix.csv
## Plot the trajectories in opt$outprefix.pdf

## Input:
## data: MaBoss data in *_probtraj.csv file saved with read.csv()
## conditions: List of conditions which probability we want to know when they have certain value set in listOfConForNodes.
## Multiple conditions should be separated with &. Example: list <- c('Node1', 'Node2&Node3')
## listOfConForNodes: Associated value to the conditions to analyze
## Example: values <- c(1,0,0) Value 1 for Node1, and Value 0 for Node2 and Node3

install.packages(optparse)

## import arguments
library("optparse")


option_list = list(
	make_option(c("-i", "--input"), type="character", metavar="character", default=NULL,  
		help="Input file: .probTraj file from a MaBoSS simulation.
                MANDATORY"),
	make_option(c("-o", "--outprefix"), type="character", metavar="character", default="StateProbabilities.txt",
		help="Output file prefix
                [default= %default]"),
	make_option(c("-n", "--conditions"), type="character", metavar="character", default=NULL,  
		help="List of conditions (or combination of conditions, separated by &) for which probability trajectories should be extracted.
                List between double quotes, separated by comma, no space.
                e.g. \"node1,node2,node3&node4\"
                MANDATORY"),
	make_option(c("-s", "--states"), type="character", metavar="character", default=NULL,  
		help="State (0 or 1) of each node for which the probability is calculated.
                List between double quotes, separated by comma. For node combination the states can be separated by comma or by &
                If no value is given, then each state is consider to be 1. 
                e.g. \"1,0,1,0\" or \"1,0,1&0\".
                [default= \"1,1,1,...\"]")
)
opt_parser = OptionParser(option_list=option_list, usage = "%prog -i <input.probTraj> -o <output_prefix> -n <\"node1,node2,node1&node3\"> [options]")
opt = parse_args(opt_parser)



## Parse required parameters
if (!("input" %in% names(opt))){
	stop("No input file given! \n---> Specify an input file (.probtraj file from MaBoSS simulation) using '-i input.probtraj' or '--input input.protraj'")
}
if (!("conditions" %in% names(opt))){
	stop("No conditions given to analyse! Specify at least one node with '-n \"conditions1,conditions,...\"' or '--conditions \"conditions1,conditions,...\"'")
}


## format parameters
conditions <- unlist(strsplit(opt$conditions,split=","))

if (!("states" %in% names(opt))){
	states <- as.character(rep(1, length(unlist(strsplit(opt$conditions,",")))))
} else {
	states <- unlist(strsplit(opt$states,","))
}

# colN=max(count.fields(opt$input, sep = "\t"))
# colNames=c("Time","TH","ErrorTH","H","HD=0",rep(c("State","Proba","Error"), (colN-5)/3))
# data <- read.csv(opt$input,sep="\t", fill=T, col.names=colNames, header=F, skip=1)
data <- read.csv(opt$input,sep="\t", header=T)

## Function to grep values for a specific node
##############################################################################
averageProb <- function(data, conditions, states){
	#Verify input data
	if (length(conditions) != length(states)){
		stop("Length of states should be the same as the total number of conditions in conditions, considering combined conditions separately.")
	}
	## get columns that define states probabilities
	data_probCol <- !(grepl('Err', colnames(data))) & grepl('Prob', colnames(data))
	data_probCol_colnames <- colnames(data)[data_probCol]

	## get columns of states fullfiling the conditions
	conditions_proba_list <- lapply(1:length(conditions), FUN=function(conditionN){
		node <- unlist(strsplit(conditions[conditionN], '&'))
		value <- as.numeric(unlist(strsplit(states[conditionN], '&')))
		## for each node
		data_probCol_colnamesFullfillCondition_list <- lapply(1:length(node), FUN=function(nodeN){
			## for each colname
			sapply(data_probCol_colnames,FUN=function(state){
				if (value[nodeN] == 1){
					return(node[nodeN] %in% unlist(strsplit(gsub("Prob.","",gsub(".$","",state)),split="\\.\\.")))
				}
				if (value[nodeN] == 0) {
					return(!(node[nodeN] %in% unlist(strsplit(gsub("Prob.","",gsub(".$","",state)),split="\\.\\."))))
				}
			})
		})
		data_probCol_colnamesFullfillCondition <- data_probCol_colnames[Reduce('*', data_probCol_colnamesFullfillCondition_list)>0]
		if (length(data_probCol_colnamesFullfillCondition)==0){
			condition_proba <- rep(0, length(data[,"Time"]))
		} 
		if (length(data_probCol_colnamesFullfillCondition)==1){
			condition_proba <- data[,data_probCol_colnamesFullfillCondition]
		}
		if (length(data_probCol_colnamesFullfillCondition)>1){
			condition_proba <- apply(data[,data_probCol_colnamesFullfillCondition],1,sum)
		}
		return(condition_proba)
	})
	## change names
	names(conditions_proba_list) <- unlist(lapply(1:length(conditions), FUN=function(conditionN){
		node=unlist(strsplit(conditions[conditionN], '&'))
		value=as.numeric(unlist(strsplit(states[conditionN], '&')))
		return(paste0(node,".",value, collapse="_"))
	}))
	#conditions_proba_df <- data.frame(Time=data[,"Time"],do.call(cbind, conditions_proba_list))
	conditions_proba_df <- data.frame(Time=data[,"Time"],do.call(cbind, conditions_proba_list),TH=data[,"TH"])
	return(conditions_proba_df)

}




## extract trajectories

cat(paste("--- Extracting probability trajectories in... ",opt$outprefix,'.csv\n',sep=""))
conditions_proba_df <- averageProb(data, conditions, states)
#cat("---------- \n Trajectories Extracted. Saving file... \n")

## Create a file
write.table(conditions_proba_df, paste(opt$outprefix,'.csv',sep=""), quote=F,sep="\t",row.names=F)




