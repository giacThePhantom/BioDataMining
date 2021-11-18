#!/usr/bin/Rscript
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
	individual <- vitamin_d_samples[grepl(vitamin_d_samples$sample_name, i), ]
	vitamin_d_low <- rbind(vitamin_d_low, individual[min(individual$vit_d_level), ])
	vitamin_d_high <- rbind(vitamin_d_high, individual[max(individual$vit_d_level), ])
}


print("Samples with low vitamin D level")
print(nrow(vitamin_d_low))
print("Samples with high vitamin D level")
print(nrow(vitamin_d_high))



data <- data[!data$sample_name %in% vitamin_d_samples$sample_name,]
print(nrow(data))

stage_samples <- data[!(is.na(data$stage) | data$stage == ""),]

print(nrow(stage_samples))

low_stage <- data[stage_samples$stage %in% "[12]", ]
high_stage <- data[grepl(stage_samples$stage, "3"), ]

print("Samples with low stage")
print(nrow(low_stage))
print("Samples with high stage")
print(nrow(high_stage))
