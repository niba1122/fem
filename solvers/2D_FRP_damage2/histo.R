data <- read.csv('out.csv', header=F, sep=",")

pdf("out.pdf")

hist(data[,1])
hist(data[,2])
