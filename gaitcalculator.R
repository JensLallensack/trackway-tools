
#settings
GAD2correction <- TRUE
slopecorrection <- TRUE

require("retistruct")
require("DescTools")
require("rlist")
require("zoo")
require("emmeans")
source("tracklib.R")

gaitcalculator <- function(a){


#find intersection
inters <- function(x){

#function: Normalize range from 0 to max(x) to 0 to 1
normalizedZero <- function(x){x/(max(x))}

#normalize GAD measurements to the range between 0 and 50
one <- x[1,1]; two <- x[1,2]; three <- x[2,1]; four <- x[2,2]
norm <- (normalizedZero(c(one,two,three,four))*50)

P1 <- c(50,norm[1])
P2 <- c(100,norm[2])
P3 <- c(50,norm[3])
P4 <- c(100,norm[4])

incl <- rad2deg(atan2L(P1,P2,P3,P4))
incl <- abs(incl)

if(incl > 1.5){
#filename <- paste("int",P1[2],P2[2],P3[2],P4[2],".svg")
#svg(filename)
#plot(rbind(P1,P2,P3,P4),pch=19)
#lines(rbind(P1,P2))
#lines(rbind(P3,P4))

if(slopecorrection == TRUE){
#check if both lines rise or fall, and correct for scale
if(P3[2] > P1[2] && P4[2] > P2[2]){
diff <- P3[2]-P1[2]
diff2 <- P4[2]-P2[2]
if(diff2 > diff){
P4[2] <- P4[2]-(diff*0.5)
}
if(diff > diff2){
P3[2] <- P3[2]-(diff2*0.5)
}
}

if(P1[2] > P3[2] && P2[2] > P4[2]){
diff <- P1[2]-P3[2]
diff2 <- P2[2]-P4[2]
if(diff2 > diff){
P2[2] <- P2[2]-(diff*0.5)
}
if(diff > diff2){
P1[2] <- P1[2]-(diff2*0.5)
}
}
}
#lines(rbind(P1,P2),col="red")
#lines(rbind(P3,P4),col="red")

#dev.off()

int <- line.line.intersection(P1,P2,P3,P4)
#print(incl)
}else{
int <- NA
}
return(c(int[[1]],incl))
}

#calculate limb phases
lphase <- function(datlist){ #single trackway
#Intersection approach
limbphases <- vector()
slopes <- vector()
if(nrow(datlist) > 1){
for(t in 1:(nrow(datlist)-1)){
for(z in t:(nrow(datlist)-1)){
int <- inters(rbind(datlist[t,],datlist[z+1,]))
limbphases <- c(limbphases,int[1])
slopes <- c(slopes,int[2])
}
}}

limbphases_raw <<- limbphases
slopes <<- slopes

#slope analysis (which slopes are most faithful)
# tab <- cbind(limbphases_raw,slopes)
# tab <- tab[tab[,1] > 0 & tab[,1] < 150 & !is.na(tab[,1]),]

#########
                                                             
limbphases <- limbphases[is.finite(limbphases)]

limbphases <- sapply(limbphases,function(x){if(x < 0){NA}else{if(x > 150){NA}else{x}}})
limbphases <- limbphases[!is.na(limbphases)]

if(length(limbphases) > 0){

#convert limbphases into standard Hildebrandt scheme (25% = lateral sequence singlefoot)
limbphases <- sapply(limbphases,function(x){50-(x-50)})

limbphases <<- limbphases
median_lp <- median(limbphases)
median_lp <<- median

mad_lp <- mad(limbphases,na.rm=TRUE)
mad_lp <<- mad_lp

number_lp <- length(limbphases)
number_lp <<- number_lp

return(c(median_lp,mad_lp,number_lp))}
else{
return(NA)
}
}


require("smatr")
if(length(a) > 1){

#ppmd vs stride evaluation; remove all non-significant ones
ppmdAp <- vector()
twname <- vector()
for(i in 1:length(a)){
Table <- corrtableSingle(a[[i]])
Table <<- Table
if(is.vector(Table) == FALSE){
if(nrow(Table) > 2){
if(any(Table[,2] < 0) == TRUE){
minneg <- min(Table[,2])
Table[,2] <- Table[,2] + minneg
}

#fit <- lmodel2(Table[,2] ~ Table[,1],range.y = "relative",range.x = "relative", nperm = 99)
fit <- NA
try(fit <- sma(Table[,2] ~ Table[,1],robust=T))

fit <<- fit

if(is.na(fit) == FALSE){
#slope <- coef(fit)[2]
#slope <- fit$regression.results$Slope[3]
slope <- fit$coef[[1]][2,1]
#pvalue <- summary(fit)$coefficients[2,4]
#pvalue <- fit$regression.results[1,5]
pvalue <- fit$pval[[1]]
}else{
slope <- 0
pvalue <- 1
}
twname <- c(twname,comment(a[[i]])[1])
if(pvalue < 0.05 & slope < 0){
significant <- TRUE
}else{
significant <- FALSE
}
#Activate if trackways should be filtered
if(significant == FALSE){
#a[[i]] <- NA
}
ppmdAp <- rbind(ppmdAp,c(slope,pvalue,significant))
}}}
ppmdAp <<- ppmdAp
rownames(ppmdAp) <- twname
colnames(ppmdAp) <- c("slope","p-value","significance")
print(ppmdAp)
}


a <- a[!is.na(a)]

datlist <- a

## Correct with GAD2

datlist <- lapply(datlist,function(x){
cbind(unlist(GAD1_B(x)),unlist(GAD3(x)),GAD2_F(x))
})

#naming
tracknames <- vector()
for(i in 1:length(a)){
tracknames <- c(tracknames,comment(a[[i]])[1])
}


for(i in 1:length(datlist)){
rownames(datlist[[i]]) <- 1:nrow(datlist[[i]])
}

datlist_raw <<- datlist
datlist_raw <- datlist

#estimate missing GADs
datlist <- lapply(datlist,function(x){
apply(x,1,function(o){
#GAD2 (third column) missing?
if(is.na(o[3]) == TRUE & anyNA(o[1:2]) == FALSE){
o[3] <- mean(c(o[1],o[2]))
}
return(o)
})})

datlist <- lapply(datlist,t)


for(i in 1:length(datlist)){
datlist[[i]] <- datlist[[i]][complete.cases(datlist[[i]][,1:3]),]
if(is.matrix(datlist[[i]]) == FALSE){
datlist[[i]] <- NA
}
}

datlist <- datlist[!is.na(datlist)]

if(length(datlist) > 0){

for(i in 1:length(datlist)){
if(is.vector(datlist[[i]]) == TRUE){
datlist[[i]] <- NA
}
if(nrow(datlist[[i]]) == 0){
datlist[[i]] <- NA
}
}
datlist <- datlist[!is.na(datlist)]

if(length(datlist) > 0){

newsel <- which(sapply(datlist,length) > 1)
datlist <- datlist[lapply(datlist,length)>1]
tracknames <- tracknames[newsel]

if(GAD2correction == TRUE){
datlist <- lapply(datlist,GAD2acc)
}

###


datlist <<- datlist


### Remove Outliers (duplicate code present)

r <- 1
for(i in 1:length(datlist)){
if(nrow(datlist[[r]]) < 2){
datlist <- datlist[-r]
tracknames <- tracknames[-r]
r <- r-1
}
r <- r+1
}

#GAD1
datlist <- lapply(datlist,function(x){
Var1 <- sdN(x[,1])*2
Mean1 <- meanN(x[,1])
apply(x,1,function(y){
if(y[1] < Mean1-Var1 | y[1] > Mean1+Var1){
#print("outlier detected GAD1")
y <- c(NA,NA,NA)
y
}else{
y
}
})})
datlist <- lapply(datlist,t)

for(i in 1:length(datlist)){
datlist[[i]] <- datlist[[i]][complete.cases(datlist[[i]]),]
}

r <- 1
for(i in 1:length(datlist)){
if(nrow(datlist[[r]]) < 2){
datlist <- datlist[-r]
tracknames <- tracknames[-r]
r <- r-1
}
r <- r+1
}

#GAD2
datlist <- lapply(datlist,function(x){
Var2 <- sdN(x[,2])*2
Mean2 <- meanN(x[,2])
apply(x,1,function(y){
if(y[2] < Mean2-Var2 | y[2] > Mean2+Var2){
#print("outlier detected GAD3")
y <- c(NA,NA,NA)
y
}else{
y
}
})})
datlist <- lapply(datlist,t)

for(i in 1:length(datlist)){
datlist[[i]] <- datlist[[i]][complete.cases(datlist[[i]]),]
}

r <- 1
for(i in 1:length(datlist)){
if(nrow(datlist[[r]]) < 2){
datlist <- datlist[-r]
tracknames <- tracknames[-r]
r <- r-1
}
r <- r+1
}

### END remove outliers


#Permutation test

ultlist <- list()
sumlist <- list()
validlist <- vector()
iis <- vector()

for(i in 1:length(datlist)){
wo <- datlist[[i]]
if(nrow(datlist[[i]]) > 3){
singtrackw <- list()
sumlistT <- vector()

validvec <- vector()

for(o in 1:nrow(datlist[[i]])){
lph <- lphase(wo)
if(all(is.na(lph)) == FALSE){
sumlistT <- rbind(sumlistT,lph)
singtrackw[[o]] <- limbphases
ultlist[[i]] <- singtrackw}
wo <- datlist[[i]]
wo <- wo[-o,]
}

if(length(sumlistT) > 0){

iis <- c(iis,i)
valid <- TRUE
me <- sumlistT[,1]
mev <- sapply(me,function(r){dif(me[1],r)})
### Change threshold here
if(any(mev > 10)){
valid <- FALSE
}

#if MAD of this trackway is too large, remove
### Change threshold here
if(sumlistT[1,2] > 30){
valid <- FALSE
}

validvec <- c(validvec,valid)
colnames(sumlistT) <- c("Median","MAD","n")
sumlistT <<- sumlistT
#print(sumlistT)
sumlist <- list.append(sumlist,sumlistT)


if(any(validvec) == FALSE){
validlist <- c(validlist,FALSE)
}else{
validlist <- c(validlist,TRUE)
}}}}

iis <<- iis
if(is.logical(iis) == FALSE){
ultlistCleaned <- list()
for(h in 1:length(iis)){
ultlistCleaned <- list.append(ultlistCleaned,ultlist[[iis[h]]])
}}else{
ultlistCleaned <- NA
}

validlist <<- validlist
ultlistCleaned <<- ultlistCleaned
sumlist <<- sumlist


if(length(sumlist) == 1){
print("==== Intersection approach (main approach 2) ====")
colnames(sumlist[[1]]) <- NULL 
print(c("Median",sumlist[[1]][1,1]))
print(c("MAD",sumlist[[1]][1,2]))
print(c("n",sumlist[[1]][1,3]))
limbphases <- ultlistCleaned[[1]][[1]]
hist(limbphases,breaks=15)

svg("histogram.svg")
hist(limbphases,breaks=15)
dev.off()

svg("trackplot.svg")
trackplot(a[[1]])
dev.off()

strideplot(a[[1]])

stridevar(a[[1]])

gadvar(a[[1]])

GADVAR <- exact_GADvar(gads,GAD2=GAD2correction)

#borrowed from the gaitcalculator-pmd.R script
print("ppmd vs stride")
Table <- corrtableSingle(a[[1]])
if(any(Table[,2] < 0) == TRUE){
minneg <- min(Table[,2])
Table[,2] <- Table[,2] - minneg
}


#fit <- lmodel2(Table[,2] ~ Table[,1],range.y = "relative",range.x = "relative", nperm = 99)
fit <- NA
try(fit <- sma(Table[,2] ~ Table[,1],robust=T))
if(is.na(fit) == FALSE){
svg("ppmd-vs-pes-str.svg")
#slope <- coef(fit)[2]
#slope <- fit$regression.results$Slope[3]
slope <- fit$coef[[1]][2,1]
#pvalue <- summary(fit)$coefficients[2,4]
#pvalue <- fit$regression.results[1,5]
pvalue <- fit$pval[[1]]
plot(Table,pch=19,asp=1)
abline(fit,col="red")
dev.off()
}else{
slope <- NA
pvalue <- NA
}

fit <<- fit
Table <<- Table
slope <<- slope
pvalue <<- pvalue
print(c("slope",slope,"p-value",pvalue))

strides <- stride(a[[1]])

#GAD vs Stride (pes stride, but manus strides when needed)
gadsr <- datlist_raw[[1]]

stri <- strides[[1]]
STRIall <- vector()
STRIDES <- vector()
for(i in 1:nrow(gadsr)){
if(is.even(i) == FALSE){
STRI <- NA
if(i > 4){
STRI <- stri[i-4]
}

M1 <- NA
M2 <- NA
if(i > 1){
M1 <- stri[i-1]
}
if(i < nrow(gadsr)){
M2 <- stri[i+1]
}
if(is.na(STRI) == TRUE){
STRI <- meanN(c(M1,M2))
}
}else{
STRI <- NA
}
STRIDES <- c(STRIDES,STRI)
STRIall <- c(STRIall,meanN(c(STRI,M1,M2)))
}

relspeed <- meanN(gadsr[,2]-gadsr[,1])/meanN(c(gadsr[,2],gadsr[,1]))

gaddd1 <- cbind(STRIDES,gadsr[,1])
gaddd1 <- gaddd1[complete.cases(gaddd1),]
gaddd2 <- cbind(STRIDES,gadsr[,3])
gaddd2 <- gaddd2[complete.cases(gaddd2),]
gaddd3 <- cbind(STRIDES,gadsr[,2])
gaddd3 <- gaddd3[complete.cases(gaddd3),]


print(c("relative speed estimator:",relspeed))

gaddd1 <<- gaddd1
gaddd2 <<- gaddd2
gaddd3 <<- gaddd3

gaddds <- rbind(gaddd1,gaddd2,gaddd3)
svg("strideVSgad.svg")
plot(gaddds,col="white",xlab="Stride",ylab="GAD")
points(gaddd1,pch=19,col="lightcoral")
abline(lm(gaddd1[,2] ~ gaddd1[,1]))
points(gaddd2,pch=19,col="royalblue1")
abline(lm(gaddd2[,2] ~ gaddd2[,1]))
points(gaddd3,pch=19,col="olivedrab3")
abline(lm(gaddd3[,2] ~ gaddd3[,1]))
dev.off()

#check for significance of slope differences

GAAD1 <- cbind(stride=gaddd1[,1],gad=gaddd1[,2],fac="GAD1")
GAAD2 <- cbind(stride=gaddd2[,1],gad=gaddd2[,2],fac="GAD2")
GAAD3 <- cbind(stride=gaddd3[,1],gad=gaddd3[,2],fac="GAD3")

allgads <- rbind(GAAD1,GAAD2,GAAD3)
allgads <- as.data.frame(allgads)
allgads$gad <- as.numeric(allgads$gad)
allgads$stride <- as.numeric(allgads$stride)
allgads$fac <- as.factor(allgads$fac)

interac <- lm(gad ~ stride*fac, data = allgads)
slopes <- lstrends(interac,"fac",var="stride")
pp <- pairs(slopes)
slopediff <- summary(pp)$p.value[2]
print(c("GAD vs stride slope (evaluation approach 2, p-value)",slopediff))
print(pp)

tw.summary <- c(n=nrow(datlist[[1]]),relspeed=relspeed,slope_ppmdVSstride=slope,pvalue=pvalue,int_appr=sumlist[[1]][1,1],gadvar_appr=GADVAR,sc.res=sc.res,scttest=scttest,pvalue_gadVSstride=slopediff)
names(tw.summary) <- c("n half-step cycles","relative speed","slope (ppmd vs stride)","pvalue slope (ppmd vs stride)","limb phase int.app.","limb phase GAD.var.app","stride comparison approach","p-value stride comparison approach","pvalue gadVSstride regression")
tw.summary <- t(tw.summary)
tw.summary <<- tw.summary

write.csv(tw.summary,file="summary.csv")

svg("gadsinseq.svg")
split.screen(c(3,1))
screen(1)
plot(1:nrow(gads),gads[,1],xlab="n",ylab="GAD50")
screen(2)
plot(1:nrow(gads),gads[,2],xlab="n",ylab="GAD25")
screen(3)
plot(1:nrow(gads),gads[,3],xlab="n",ylab="GAD0")
dev.off()

return("tw.summary")
}

#if multiple trackways …
if(length(sumlist) > 1){
#sumlist <- sumlist[lapply(sumlist,length)>1]

sumlistt <<- sumlist

tracklist <- vector()
for(i in 1:length(validlist)){
if(validlist[i] == TRUE){
tracklist <- c(tracklist,sumlist[[i]][1])
tracklist <<- tracklist
}}

tracknames <- tracknames[validlist]
tracknames <<- tracknames

print("==== Intersection approach (main approach 2) ====")
#
med <- median(tracklist)
mad <- mad(tracklist)
print(c("median:",med,"MAD:",mad,"n trackways:",length(tracklist),"n removed trackways:",sum(validlist == FALSE)))

print("all limbphases together:")
#all limbphases together
Limbphases <- vector()
for(i in 1:length(ultlistCleaned)){
if(validlist[i] == TRUE){
Limbphases <- c(Limbphases,ultlistCleaned[[i]][[1]])
}
}
Limbphases <<- Limbphases	#all limbphases together
medA <- median(Limbphases)
madA <- mad(Limbphases)
svg("hisogram_int_app.svg")
hist(Limbphases,breaks=15)
dev.off()
print(c("median all:",medA,"MAD all:",madA,"n limbphases:",length(Limbphases)))

print("==== Stride comparison approach (evaluation approach 1)")

stridediffs <- vector()
for(i in 1:length(a)){
if(nrow(a[[i]]) > 9){
stvr <- stridevar(a[[i]],quiet=TRUE)
stridediffs <- rbind(stridediffs,stvr)
}}
strot <- sum(stridediffs[,1],na.rm=TRUE)
space <- sum(stridediffs[,2],na.rm=TRUE)

if(is.numeric(strot) & is.numeric(space)){
if(nrow(stridediffs) > 3){
if(strot < space){
print("Pace can be excluded")
ttest <- t.test(stridediffs[,1],stridediffs[,2],alternative="less")
print(c("T-test: ",round(ttest$p.value,2)))
}else{if(space < strot){
print("Trot can be excluded")
ttest <- t.test(stridediffs[,2],stridediffs[,1],alternative="less")
print(c("T-test: ",round(ttest$p.value,2)))
}}}}

print("==== GAD variation approach (main approach 1) ====")

percs <- vector()

GADS <- list()
gads <- vector()
significances <- vector()
gadvvvr <- list()

for(i in 1:length(a)){
GADS[[i]] <- gadvar(a[[i]])
}


#remove empty, NA, and NULL list elements
GADS <- Filter(Negate(anyNA),GADS)
GADS <- GADS[!sapply(GADS,is.null)]
empty <- function(x){
len <- length(x)
if(len == 0){
return(TRUE)
}else{
return(FALSE)
}
}
GADS <- GADS[!sapply(GADS,empty)]

#exact GAD variation approach multiple sample
lplist <- vector()
summarylist <- vector()
Weighed_limbphase <- vector()
for(i in 1:length(GADS)){
trackwaylength <- nrow(GADS[[i]])

limbphase <- exact_GADvar(GADS[[i]],quiet=TRUE,GAD2=GAD2correction)

Weighed_limbphase <- c(Weighed_limbphase,rep(limbphase,trackwaylength))
newrow <- matrix(c(trackwaylength,limbphase),nrow=1)
rownames(newrow) <- rownames(GADS[[i]])[1]
summarylist <- rbind(summarylist,newrow)
}
Weighed_limbphase <<- Weighed_limbphase
colnames(summarylist) <- c("n step cycles","LP GAD var")
summarylist <<- summarylist
print("precise GAD variation approach:")
print(c("weighed mean:",mean(Weighed_limbphase)))
print(summarylist)
hist(Weighed_limbphase)
}}}
}
