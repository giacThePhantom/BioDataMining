#!/usr/bin/Rscript

library(dplyr)

args <- commandArgs(trailingOnly=TRUE)

metadata <- args[1]

data <- read.csv(metadata)

#print(head(data))


disease_filter <- c('',
			'primary colorectal adenocarcinomaadenoma',
			'colorectal cancer',
			'tumor',
			'colon cancercancer',
			'tissue: colon tumor',
			'tissue: adenocarcinoma',
			'tissue: adenoma',
			'colon tumour',
			'colon cancer tissue',
			'Primary tumor resection',
			'Colon adenocarcinoma',
			'tissue: colorectal cancer biopsy',
			'Adenocarcinoma, NOS',
			'Papillary adenocarcinoma, NOS',
			'Adenocarcinoma with mixed subtypes'
			)

print("before filter for disease")
print(nrow(data))


data <- data[data$diseas %in% disease_filter, ]




print("after filter for disease")
print(nrow(data))


vitamin_d_samples <- data[complete.cases(data$vit_d_level), ]


print("Samples with vitamin D level")
print(nrow(vitamin_d_samples))

sample_name <- vitamin_d_samples$sample_name
sample_index <- c()
for(i in sample_name){
	sample_index <- append(sample_index, strsplit(i, '_')[[1]][1])
}
sample_index <- unique(sample_index)

vitamin_d_low <- vitamin_d_samples[FALSE, ]
vitamin_d_high <- vitamin_d_samples[FALSE, ]



for(i in sample_index){
	individual <- dplyr::filter(vitamin_d_samples, grepl(i, sample_name))
	vitamin_d_low <- rbind(vitamin_d_low, individual[which.min(individual$vit_d_level), ])
	vitamin_d_high <- rbind(vitamin_d_high, individual[which.max(individual$vit_d_level), ])
}


print("Samples with low vitamin D level")
print(nrow(vitamin_d_low))

write.csv(vitamin_d_low, '../../data/splitted_samples/vitamin_d_low.csv')

print("Samples with high vitamin D level")
print(nrow(vitamin_d_high))

write.csv(vitamin_d_high, '../../data/splitted_samples/vitamin_d_high.csv')


data <- data[!data$sample_name %in% vitamin_d_samples$sample_name,]
print(nrow(data))
print(head(data))

stage_samples <- data[complete.cases(data$stage),]

print(nrow(stage_samples))
print(unique(stage_samples$stage))
low_stage <- data[(stage_samples$stage == "0" |
									 stage_samples$stage == "1" |
									 stage_samples$stage == "1A" |
									 stage_samples$stage == "2" |
									 stage_samples$stage == "2A" |
									 stage_samples$stage == "2B" |
									 stage_samples$stage == "2C"
									 ), ]

high_stage <- data[(stage_samples$stage == "3" |
									 stage_samples$stage == "3A" |
									 stage_samples$stage == "3B" |
									 stage_samples$stage == "3C" |
									 stage_samples$stage == "4" |
									 stage_samples$stage == "4A" |
									 stage_samples$stage == "4B"
									 ), ]

print("Samples with low stage")
print(nrow(low_stage))
write.csv(low_stage, '../../data/splitted_samples/stage_low.csv')
print("Samples with high stage")
print(nrow(high_stage))
write.csv(high_stage, '../../data/splitted_samples/stage_high.csv')

os_data <- data[(complete.cases(data$overall_surv_months) & data$overall_surv_months != ""),]
print(os_data$overall_surv_months)

print("Sample with overall survival")
print(nrow(os_data))
write.csv(os_data, '../../data/splitted_samples/os_data.csv')


print("Intersection between overall survival and stage dataset")
print(nrow(os_data[!(os_data$sample_name %in% stage_samples$sample_name), ]))

print(nrow(stage_samples[sample(nrow(stage_samples), size=2),]))
