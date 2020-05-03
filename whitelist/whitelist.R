#### Library_1 analysis ####

lib_1 <- read.table('/projectnb/bf528/users/group7/project4/Data/fastq_counts/SRR3879604_1.countbc.txt', header = FALSE)
#Reorder lib_1 so it's in decreasing order
lib_1 <- lib_1[order(lib_1$V2, decreasing = T),]

count_1 <- lib_1$V2
mean_count_1 <- mean(count_1)
sd_count_1 <- sd(count_1)

#Setting a threshold for the data 
threshold_1 <- mean_count_1 + 2*sd(count_1)


#Numbers from 1 to however many barcodes exist in the data
num <- seq(1, nrow(lib_1))

#Color for the plot: Red if above threshold (Kept for further analysis), blue if below threshold (should be removed)
color.plot <- sapply(count_1, function(x){
    if(x > threshold_1){
        return("Red")
    }else{
        return("Blue")
    }
})


#Plot: x- axis is barcodes numbered from 1 to 1293792 (in this case)
#y-axis is the number of UMI counts

plot_1 <- plot(num, count_1, xlab = "Number of Barcodes", ylab = "Count", col = color.plot, pch = 20, cex = 0.75)

#obtaining whitelist 
whitelist_list_1 <- unlist(sapply(count_1, function(x){
  if(x > threshold_1){
    return(x)
  }
}))


whitelist1_data<- as.data.frame(whitelist_list_1)
whitelist_lib_1 <- as.data.frame(lib_1[1:2454,1])
rownames(whitelist_lib_1) <- 1: nrow(whitelist_lib_1)
colnames(whitelist_lib_1) <- c("Barcodes")
check <- cbind(whitelist_lib_1, whitelist1_data)
write.table(whitelist_lib_1, file = "library_1_whitelist.txt", row.names = FALSE, col.names = FALSE,  quote = FALSE)
#### end ####



#### Library_2 analysis ####
lib_2 <- read.table('/projectnb/bf528/users/group7/project4/Data/fastq_counts/SRR3879605_1.countbc.txt', header = FALSE)

#Reorder lib_2 so it's in decreasing order
lib_2 <- lib_2[order(lib_2$V2, decreasing = T),]
count_2 <- lib_2$V2
mean_count_2 <- mean(count_2)
sd_count_2 <- sd(count_2)

#Setting a threshold for the data 
threshold_2 <- mean_count_2 + 2*sd(count_2)

#Numbers from 1 to however many barcodes exist in the data
num2 <- seq(1, nrow(lib_2))

#Color for the plot: Red if above threshold (Kept for further analysis), blue if below threshold (should be removed)
color.plot_2 <- sapply(count_2, function(x){
  if(x > threshold_2){
    return("Red")
  }else{
    return("Blue")
  }
})


#Plot: x- axis is barcodes numbered from 1 to 1333842 (in this case)
#y-axis is the number of UMI counts

plot_2 <- plot(num2, count_2, xlab = "Number of Barcodes", ylab = "Count", col = color.plot_2, pch = 20, cex = 0.75)

#obtaining whitelist 
whitelist_list_2 <- unlist(sapply(count_2, function(x){
  if(x > threshold_2){
    return(x)
  }
}))


whitelist2_data<- as.data.frame(whitelist_list_2)
whitelist_lib_2 <- as.data.frame(lib_2[1:2955,1])
rownames(whitelist_lib_2) <- 1: nrow(whitelist_lib_2)
colnames(whitelist_lib_2) <- c("Barcodes")
check2 <- cbind(whitelist_lib_2, whitelist2_data)
write.table(whitelist_lib_2, file = "library_2_whitelist.txt", row.names = FALSE, col.names = FALSE,  quote = FALSE)
#### end ####


#### Library_3 analysis ####
lib_3 <- read.table('/projectnb/bf528/users/group7/project4/Data/fastq_counts/SRR3879606_1.countbc.txt', header = FALSE)

#Reorder lib_3 so it's in decreasing order
lib_3 <- lib_3[order(lib_3$V2, decreasing = T),]
count_3 <- lib_3$V2
mean_count_3 <- mean(count_3)
sd_count_3 <- sd(count_3)

#Setting a threshold for the data 
threshold_3 <- mean_count_3 + 2*sd(count_3)

#Numbers from 1 to however many barcodes exist in the data
num3 <- seq(1, nrow(lib_3))

#Color for the plot: Red if above threshold (Kept for further analysis), blue if below threshold (should be removed)
color.plot_3 <- sapply(count_3, function(x){
  if(x > threshold_3){
    return("Red")
  }else{
    return("Blue")
  }
})


#Plot: x- axis is barcodes numbered from 1 to 1227152 (in this case)
#y-axis is the number of UMI counts

plot_3 <- plot(num3, count_3, xlab = "Number of Barcodes", ylab = "Count",  col = color.plot_3, pch = 20, cex = 0.75)

#obtaining whitelist 
whitelist_list_3 <- unlist(sapply(count_3, function(x){
  if(x > threshold_3){
    return(x)
  }
}))


whitelist3_data<- as.data.frame(whitelist_list_3)
whitelist_lib_3 <- as.data.frame(lib_3[1:2864,1])
rownames(whitelist_lib_3) <- 1: nrow(whitelist_lib_3)
colnames(whitelist_lib_3) <- c("Barcodes")
check3 <- cbind(whitelist_lib_3, whitelist3_data)
write.table(whitelist_lib_3, file = "library_3_whitelist.txt", row.names = FALSE, col.names = FALSE,  quote = FALSE)
#### end ####

