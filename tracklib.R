require("circular")
require("BBmisc") #for data normalization, function "stridevar"

#BASE FUNCTIONS
#Calculate euclidian distance between two points (from Momocs)
ed <- function (pt1, pt2) {return(sqrt((pt1[1] - pt2[1])^2 + (pt1[2] - pt2[2])^2))}
#difference (always positive) between two points
dif <- function(pt1,pt2) {return(abs(diff(as.numeric(c(pt1,pt2)))))}
#convert degree to radians
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

#check if number is even
is.even <- function(x) {x %% 2 == 0}

#calculate coefficient of variation
varCof <- function(x){sd(x,na.rm=TRUE)/meanN(x)}

#Calculate means and standard deviations by not taking NAs into account
sdN <- function(x){sd(x,na.rm=TRUE)}

meanN <- function(x){mean(x,na.rm=TRUE)}
medianN <- function(x){median(x,na.rm=TRUE)}
meanNcirc <- function(x){circ <- suppressWarnings(circular(x)); return(suppressWarnings(meanN(circ)))}
maxN <- function(x){max(x,na.rm=TRUE)}
minN <- function(x){min(x,na.rm=TRUE)}
varN <- function(x){var(x,na.rm=TRUE)}
madN <- function(x){var(x,na.rm=TRUE)}

#only select pes measures
pes <- function(x){
if(is.list(x) == TRUE){
x <- x[[1]]
}
PES <- x[c(TRUE, FALSE)]
return(PES)
}

#only select manus measures
manus <- function(x){
if(is.list(x) == TRUE){
x <- x[[1]]
}
MANUS <- x[c(FALSE, TRUE)]
return(MANUS)
}

levtest <- function(x, y) {#
require("car")
  leveneTest(dv~gr, data=rbind(data.frame(dv=x, gr='gr1'),
  data.frame(dv=y, gr='gr2')), center='median')
}

fligner <- function(x, y) {#
  fligner.test(dv~gr, data=rbind(data.frame(dv=x, gr='gr1'),
  data.frame(dv=y, gr='gr2')))
}

normalized <- function(x){(x-min(x))/(max(x)-min(x))}

#find the difference/distance between two angles (x and y). Does not take direction into account (angles are interchangeable). Direction is taken into accound if the "abs" function is removed.
diffang <- function(x,y){abs(atan2(sin(x-y), cos(x-y)))}

#atan2 function: find angle between a line defined by two points to the vertical (y-axis). First center points on coordinate origin, second compute angle
atan2P <- function(p1,p2){
xcomp <- as.numeric(p2[1]-p1[1])
ycomp <- as.numeric(p2[2]-p1[2])
ang <- atan2(xcomp,ycomp)
return(ang)
}

#atan2 function: find angle between to lines.
#p1, p2: start and end point of first line; p3, p4: start and end point of second line.
atan2L <- function(p1,p2,p3,p4){
firstang <- atan2P(p1,p2)
secondang <- atan2P(p3,p4)
angle <- secondang-firstang
if(is.na(angle) == FALSE){
#correct for -180° to 180°:
angD <- rad2deg(angle)
if(angD > 180){angD <- (360-angD)*-1; angle <- deg2rad(angD)}
if(angD < -180){angD <- (360-angD)*-1; angle <- deg2rad(angD)}
return(angle)
}else{
return(NA)
}
}

#find average angle against the vertical for set of lines
avlineang <- function(x){	#where x is a matrix or data frame of xy-coordinates
n <- nrow(x)
nn <- nrow(x)/2		#number of lines
ang <- vector()
for(i in 1:nn*2-1){	#for every second row (the first point of each line)
ang <- c(ang,atan2P(x[i,],x[i+1,]))	#calculate the angle against the vertical using the atan2P function
}
av_ang <- mean.circular(ang,na.rm=TRUE)
return(av_ang)
}

#rotate one point relative to another
#p1: origin point, p2: point to be rotated around origin point, angle: rotation angle in rad.
#rotates counter-clockwise
rot <- function(p1,p2,angle){
#set p1 as coordinate system center, translate p2 accordingly
new_p2 <- p2
new_p2[1] <- as.numeric(p2[1]-p1[1])	#x value
new_p2[2] <- as.numeric(p2[2]-p1[2])	#y value
#apply rotation matrix
new_x <- new_p2[1]*cos(angle)-new_p2[2]*sin(angle)
new_y <- new_p2[1]*sin(angle)+new_p2[2]*cos(angle)
#translate back
new_x <- new_x+p1[1]
new_y <- new_y+p1[2]
p2_rot <- unlist(c(new_x,new_y))
return(p2_rot)
}

#Rotate two points vertically, and a third or more points accordingly
#p1, the origin around which the other two points are rotated
#p2, the point that is going to form a vertical line with the origin
#p3, the third point which stays in a constant relative geometric relationship with the other two points
#p3 MAY be a matrix of multiple points
rotvert <- function(p1,p2,p3){
angle <- atan2P(p1,p2)	#find angle to vertical
#perform rotation
p2_rot <- rot(p1,p2,angle)
#if there is only a single third point
#check whether we have a matrix/data frame of multiple points or a single point
a <- vector()
if(is.data.frame(p3) == TRUE | is.matrix(p3) == TRUE){
if(nrow(p3) > 1) {a <- "multiple"}else{a <- "single"}}else{a <- "single"}

#if multiple points
if(a == "multiple")
{
p3_rot <- matrix(rep(NA,nrow(p3)*2),ncol=2)	#initiate empty matrix of known length
for(i in 1:nrow(p3)){	#for each point, apply the rotation
p3_rot[i,] <- rot(p1,p3[i,],angle)
}}else{
p3_rot <- rot(p1,p3,angle)
}
matrix <- rbind(p1,p2_rot,p3_rot)
return(matrix)
}

#plot trackway
#accepts matrix of xy coordinates (LP1, LM1, …)
trackplot <- function(x){
# detect if trackway is upside down; mirror if necessary so that direction of travel is always up the screen
tx <- x[complete.cases(x),]
if(tx[1,2] > tx[nrow(tx),2])
{x[,2] <- x[,2]*-1}	#mirror
plot(x,asp=1,pch=19,col="white")	
points(x[seq(1,nrow(x),2),],col="blue",pch=15,cex=1) #hind, blue
points(x[seq(2,nrow(x),2),],col="red",pch=16,cex=1) #front, red
lines(rbind(x[1,],x[5,]),col="blue")
lines(rbind(x[3,],x[7,]),col="blue")
lines(rbind(x[4,],x[8,]),col="red")
lines(rbind(x[6,],x[10,]),col="red")
}


#return coordinates of surrounding footprints
#accepts "i" value (row numer of current footprint)
#only works relative to pes footprints!
#x: list of coordinates returned by the imp_trackway_calculator script
rettrack <- function(x,i){

RP0 <- rep(NA,2)
RM0 <- rep(NA,2)
LP1 <- rep(NA,2)
LM1 <- rep(NA,2)
RP1 <- rep(NA,2)
RM1 <- rep(NA,2)

LM2 <- rep(NA,2)
RP2 <- rep(NA,2)
RM2 <- rep(NA,2)
LP3 <- rep(NA,2)
LM3 <- rep(NA,2)
RP3 <- rep(NA,2)
RM3 <- rep(NA,2)

LP2 <- x[i,]

if(i > 2){
RP1 <- x[i-2,]
RM1 <- x[i-1,]
if(i > 4){
LP1 <- x[i-4,]
LM1 <- x[i-3,]
if(i > 6){
RP0 <- x[i-6,]
RM0 <- x[i-5,]
}}}

n <- nrow(x)

if(i < n){
LM2 <- x[i+1,]
if(i+1 < n){
RP2 <- x[i+2,]
if(i+2 < n){
RM2 <- x[i+3,]
if(i+3 < n){
LP3 <- x[i+4,]
if(i+4 < n){
LM3 <- x[i+5,]
if(i+5 < n){
RP3 <- x[i+6,]
if(i+6 < n){
RM3 <- x[i+7,]
}}}}}}}

ext_coords <- rbind(RP0=RP0,RMO=RM0,LP1=LP1,LM1=LM1,RP1=RP1,RM1=RM1,LP2=LP2,LM2=LM2,RP2=RP2,RM2=RM2,LP3=LP3,LM3=LM3,RP3=RP3,RM3=RM3)

if(i == 5){
ex <<- ext_coords
}

return(ext_coords)

}


#Variance comparison (compare variance between two variables)
#reads directly from bdattable object
varcom <- function(bdattable,par1,par2){
par1 <- bdattable[par1]
par2 <- bdattable[par2]
print(c("sd: ",sdN(par1[,1])," vs. ",sdN(par2[,1])))
print(c("varCof: ",varCof(par1[,1])," vs. ",varCof(par2[,1])))
vcom <- var.test(par1[,1],par2[,1])
return(vcom)
}



#Base functions
#pesonly/manusonly: (select every second element from vector); takes list of trackways
pesonly <- function(x){
for(i in 1:length(x)){
x[[i]] <- x[[i]][c(TRUE,FALSE)]
}
return(x)
}
manusonly <- function(x){
for(i in 1:length(x)){
x[[i]] <- x[[i]][c(FALSE,TRUE)]
}
return(x)
}



#CALCULATE TRACKWAY PARAMETERS

#Stride length
#takes single trackway or list of trackways
stride <- function(x){
st <- function(x){
pl <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-4)){
pl[i] <- ed(x[i,],x[i+4,])
}
pl <- unlist(pl); pl <- list(pl)
names(pl) <- comment(x)[1]
return(pl)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- st(x)}else{
fl <- sapply(x,st)}
Stride <<- fl
return(fl)
}

#Pace length
#takes single trackway or list of trackways
pace <- function(x){
pa <- function(x){
pl <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-1)){
pl[i] <- ed(x[i,],x[i+2,])
}
pl <- unlist(pl); pl <- list(pl)
names(pl) <- comment(x)[1]
return(pl)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- pa(x)}else{
fl <- sapply(x,pa)}
return(fl)
}

#footprint length
#takes single trackway or list of trackways
fl <- function(x){
flength <- function(x){
u <- 1	#counter
pl <- rep(NA,floor(nrow(x)/5)+1)
for(i in 1:((nrow(x)/5)+1)){
pl[i] <- ed(x[u,],x[u+1,])
u <- u+5
}
pl <- unlist(pl); pl <- list(pl)
names(pl) <- comment(x)[1]
return(pl)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- flength(x)}else{
fl <- sapply(x,flength)}
return(fl)
}

#footprint width
#takes single trackway or list of trackways
fw <- function(x){
fwidth <- function(x){
u <- 1	#counter
pw <- rep(NA,floor(nrow(x)/5)+1)
for(i in 1:((nrow(x)/5)+1)){
pw[i] <- ed(x[u+2,],x[u+3,])
u <- u+5
}
pw <- unlist(pw); pw <- list(pw)
names(pw) <- comment(x)[1]
return(pw)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- fwidth(x)}else{
fl <- sapply(x,fwidth)}
return(fl)
}

#pace angulation
#takes single trackway or list of trackways
pace_ang <- function(x){
pa <- function(x){
pang <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-4)){
pang[i] <- abs(rad2deg(atan2L(x[i+2,],x[i,],x[i+2,],x[i+4,])))
}
pang <- unlist(pang); pang <- list(pang)
names(pang) <- comment(x)[1]
return(pang)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- pa(x)}else{
fl <- sapply(x,pa)}
return(fl)
}

#trackway width (relative to opposite stride)
#takes single trackway or list of trackways
tw_width <- function(x){
tww <- function(x){
pang <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-4)){
vert <- rotvert(x[i,],x[i+4,],x[i+2,])
pang[i] <- dif(vert[3,1],vert[2,1])
}
pang <- unlist(pang); pang <- list(pang)
names(pang) <- comment(x)[1]
return(pang)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- tww(x)}else{
fl <- sapply(x,tww)}
return(fl)
}

#footprint rotation
#takes single trackway or list of trackways (listfull)
footrot <- function(x){
fr <- function(x){
u <- 1	#counter
pl <- rep(NA,floor(nrow(x)/5)+1)
nat <- as.data.frame(rbind(c(NA,NA),c(NA,NA)))
for(i in 1:((nrow(x)/5)+1)){
u <<- u
x <<- x
#previous stride
if(u+3 < nrow(x) & u-16 > 0){
prevstride <- rbind(x[u-16,],x[u+4,])}else{
prevstride <- nat}
#subsequent stride
if(u+23 < nrow(x)){
substride <- rbind(x[u+4,],x[u+24,])}else{
substride <- nat}
#opposite stride
if(u+13 < nrow(x) & u-6 > 0){
oppstride <- rbind(x[u-6,],x[u+14,])}else{
oppstride <- nat}
#enclosing stride (manus enclosed by pes, or pes enclosed by manus)
if(u+8 < nrow(x) & u-11 > 0){
encstride <- rbind(x[u-11,],x[u+9,])}else{
encstride <- nat}
#same names
colnames(prevstride) <- colnames(substride)
colnames(oppstride) <- colnames(substride)
colnames(encstride) <- colnames(substride)
stridetab <- rbind(prevstride,substride,oppstride,encstride)
#find mean angle of three strides
stridetab <<- stridetab
avang <- avlineang(stridetab)
#print(c("avang:",rad2deg(avang)))
#long axis of pes orientation
 if(is.na(x[u,]) == FALSE && is.na(x[u+1,]) == FALSE){
 peslong <- atan2P(x[u,],x[u+1,])}
 #check which coords are available, to be robust in case L1 or L2 are missing.
# else{
# if(is.na(x[u,]) == FALSE && is.na(x[u+4,]) == FALSE){
# peslong <- atan2P(x[u,],x[u+4,])}else{
# if(is.na(x[u+1,]) == FALSE && is.na(x[u+4,]) == FALSE){
# peslong <- atan2P(x[u+4,],x[u+1,])}
else{peslong <- NA}
# }}
r <- rad2deg(peslong-avang)
if(is.na(r) == FALSE){
if(r > 180){r <- (360-r)*-1}else{
if(r < -180){r <- (-360-r)*-1}}}
u <- u+5
pl[i] <- r
}
if(all(is.na(pl)) == TRUE){return(NA);break}
#change sign for right footprints
if(is.even(length(pl)) == FALSE){pl <- c(pl,NA)}	#even numbers of elements to create 2-column matrix
mat <- matrix(pl,ncol=2,byrow=TRUE)
mat[c(FALSE,TRUE),] <- mat[c(FALSE,TRUE),]*-1 	#select every second row and change sign
pl <- c(t(mat))	#change back to original order, as vector
pl <- na.trim(pl,sides="right")
pl <- list(pl)
names(pl) <- comment(x)[1]
return(pl)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- fr(x)}else{
fl <- sapply(x,fr)}
return(fl)
}

#Direction of travel (in degrees) at one point (footprint) in the trackway based on previous, preceeding, and opposite strides
#takes single trackway or list of trackways. Directions for pes are based on pes strides, and directions for manus on manus strides only.
#pesandmanus: if TRUE, cardinal directions for each footprint will be based on both surrounding pes and manus strides.
dot <- function(x,pesandmanus=FALSE){
DOT <- function(x){
pl <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x))){
#following stride
if(i < nrow(x)-3){
followstride <- rbind(x[i,],x[i+4,])
}else{
followstride <- rbind(c(NA,NA),c(NA,NA))
}
#preceeding stride
if(i > 4){
prevstride <- rbind(x[i-4,],x[i,])
}else{
prevstride <- rbind(c(NA,NA),c(NA,NA))
}
#opposite stride
if(i > 2 & i < nrow(x)-2){
oppstride <- rbind(x[i-2,],x[i+2,])
}else{
oppstride <- rbind(c(NA,NA),c(NA,NA))
}
### pesandmanus==TRUE
oppmanus <- rbind(c(NA,NA),c(NA,NA))
follomanus <- rbind(c(NA,NA),c(NA,NA))
surrmanus <- rbind(c(NA,NA),c(NA,NA))
if(pesandmanus == TRUE){
#surrounding manus stride (for manus: opposite pes stride)
if(i > 3 & i < nrow(x)){
surrmanus <- rbind(x[i-3,],x[i+1,])
}
#following manus stride (for manus: opposite pes stride)
if(i < nrow(x)-4){
follomanus <- rbind(x[i+1,],x[i+5,])
}
#opposite manus stride (for pes: enclosing pes stride)
if(i > 1 & i < nrow(x)-2){
oppmanus <- rbind(x[i-1,],x[i+3,])
}}
colnames(followstride) <- c("x","y")
colnames(prevstride) <- c("x","y")
colnames(oppstride) <- c("x","y")
colnames(oppmanus) <- c("x","y")
colnames(follomanus) <- c("x","y")
colnames(surrmanus) <- c("x","y")
stridetab <- rbind(followstride,prevstride,oppstride,oppmanus,follomanus,surrmanus)
stridetab <<- stridetab
pl[i] <- rad2deg(avlineang(stridetab))
}
pl <- unlist(pl); pl <- list(pl)
names(pl) <- comment(x)[1]
return(pl)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- DOT(x)}else{
fl <- sapply(x,DOT)}
return(fl)
}


#Pes-manus distance
#takes single trackway or list of trackways
pmd <- function(x){
PMD <- function(x){
pl <- rep(NA,nrow(x))	#initialise empty vector of set length
strideangT <- unlist(deg2rad(dot(x,pesandmanus=TRUE)[[1]]))
for(i in 1:(nrow(x)-1)){
if(is.even(i) == FALSE){
if(anyNA(x[i,]) == FALSE){ 
pii <- ed(x[i,],x[i+1,])
piiang <- atan2P(x[i,],x[i+1,])
strideang <- strideangT[i]
### check if negative or positive (not activated, as it makes more sense for parallel pes-manus distance)
# if(anyNA(c(piiang,strideang)) == FALSE){
# difff <- diffang(piiang,strideang)
# if(difff > deg2rad(90)){
# pii <- pii*-1
# piiang <<- piiang
# strideang <<- strideang
# }}
pl[i] <- pii
}}}
pl <- unlist(pl); pl <- list(pl)
names(pl) <- comment(x)[1]
return(pl)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- PMD(x)}else{
fl <- sapply(x,PMD)}
return(fl)
}

#parallel pes-manus distance
#takes single trackway or list of trackways
ppmd <- function(x){
PPMD <- function(x){
pl <- rep(NA,nrow(x))	#initialise empty vector of set length
strideangT <- unlist(deg2rad(dot(x,pesandmanus=TRUE)[[1]]))
strideangF <- unlist(deg2rad(dot(x,pesandmanus=FALSE)[[1]]))
for(i in 1:(nrow(x)-1)){
if(is.even(i) == FALSE){
if(anyNA(c(x[i,],x[i+1,])) == FALSE){
if(is.na(strideangF[i]) == FALSE){ 	#prefer defining the trackway midline based on pes prints only
strideang <- strideangF[i]
}else{
strideang <- strideangT[i]	#but refer to surrounding manus prints if pes prints are not available 
}
pii <- ed(x[i,],x[i+1,])
M <- rot(x[i,],x[i+1,],strideang)
ppmdd <- dif(M[2],x[i,2])
#now check if pmd is negative (behind the pes) or positive (in front of the pes), and change sign accordingly
piiang <- atan2P(x[i,],x[i+1,])
if(anyNA(c(piiang,strideang)) == FALSE){
difff <- diffang(piiang,strideang)

if(difff > deg2rad(90)){
ppmdd <- ppmdd*-1
}
pl[i] <- ppmdd
}}}}
pl <- unlist(pl); pl <- list(pl)
names(pl) <- comment(x)[1]
return(pl)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- PPMD(x)}else{
fl <- sapply(x,PPMD)}
return(fl)
}

#GAD1
#takes single trackway or list of trackways
GAD1 <- function(x){
gad1 <- function(x){
gad <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-5)){
even <- is.even(i)
if(even == FALSE){
#find center between x and y values of both the manus and pes pairs
co_pes <- c(mean(c(x[i,1],x[i+2,1])),mean(c(x[i,2],x[i+2,2])))
co_manus <- c(mean(c(x[i+3,1],x[i+5,1])),mean(c(x[i+3,2],x[i+5,2])))
gad[i] <- ed(co_pes,co_manus)
}}
gad <- unlist(gad); gad <- list(gad)
names(gad) <- comment(x)[1]
return(gad)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- gad1(x)}else{
fl <- sapply(x,gad1)}
return(fl)
}

#GAD1 one stepcycle in-between
#takes single trackway or list of trackways
GAD1_coup1 <- function(x){
gad1_coup1 <- function(x){
gad <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-9)){
even <- is.even(i)
if(even == FALSE){
#find center between x and y values of both the manus and pes pairs
co_pes <- c(mean(c(x[i,1],x[i+2,1])),mean(c(x[i,2],x[i+2,2])))
co_manus <- c(mean(c(x[i+7,1],x[i+9,1])),mean(c(x[i+7,2],x[i+9,2])))
gad[i] <- ed(co_pes,co_manus)
}}
gad <- unlist(gad); gad <- list(gad)
names(gad) <- comment(x)[1]
return(gad)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- gad1_coup1(x)}else{
fl <- sapply(x,gad1_coup1)}
return(fl)
}

#GAD2
#takes single trackway or list of trackways
GAD2 <- function(x,method="both"){
gad2 <- function(x){
#method options:
#both: Projecting mid-stride of RM1 based on both half the stride length RM1-RM2 and the position of LM2.
#opposite: Only taking into account the position of the opposite manus (LM2); this is the traditional approach following Leonardi.
#same: Only taking into account the actual stride (/2), RM1-RM2.
gad <- rep(NA,nrow(x))	#initialise empty vector of set length

for(i in 1:(nrow(x)-5)){
even <- is.even(i)
if(even == FALSE){
ex <- rettrack(x,i)	#writes coordinates for LP1, LM1 etc., with LP2 as first track
#find orientation of trackway midline
RP0 <- ex[1,]; RM0 <- ex[2,]; LP1 <- ex[3,]; LM1 <- ex[4,]; RP1 <- ex[5,]; RM1 <- ex[6,]; LP2 <- ex[7,]; LM2 <- ex[8,]; RP2 <- ex[9,]; RM2 <- ex[10,]; LP3 <- ex[11,]; LM3 <- ex[12,]; RP3 <- ex[13,]; RM3 <- ex[14,]
refstrides <- rbind(RP1,RP2,LP2,LP3,RM2,RM3)
ang <- avlineang(refstrides)	#average angle of selected strides to vertical
# Rotate into vertical position
RP2 <- rot(LP2,RP2,ang)
RM2 <- rot(LP2,RM2,ang)
LM3 <- rot(LP2,LM3,ang)
RM3 <- rot(LP2,RM3,ang)
# find points between which the GAD is measured
coo_pes <- mean(c(LP2[[2]],RP2[[2]]))
coo_same <- mean(c(RM2[[2]],RM3[[2]]))
coo_opposite <- LM3[[2]]
# Measure according to specified method
if(method == "both"){
coo_manus <- meanN(c(coo_same,coo_opposite))
}
if(method == "same"){
coo_manus <- coo_same
}
if(method == "opposite"){
coo_manus <- coo_opposite
}

gad[i] <- dif(coo_manus,coo_pes)
if(is.nan(gad[i])){
gad[i] <- NA
}
}}
gad <- unlist(gad); gad <- list(gad)
names(gad) <- comment(x)[1]
return(gad)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- gad2(x)}else{
fl <- sapply(x,gad2)}

return(fl)
}

#GAD3
#takes single trackway or list of trackways
GAD3 <- function(x){
gad3 <- function(x){
gad <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-5)){
even <- is.even(i)
if(even == FALSE){	#only uneven rows = only the pes prints, not the manus prints
ex <- rettrack(x,i)	#writes coordinates for LP1, LM1 etc., with LP2 as first track
RP0 <- ex[1,]; RM0 <- ex[2,]; LP1 <- ex[3,]; LM1 <- ex[4,]; RP1 <- ex[5,]; RM1 <- ex[6,]; LP2 <- ex[7,]; LM2 <- ex[8,]; RP2 <- ex[9,]; RM2 <- ex[10,]; LP3 <- ex[11,]; LM3 <- ex[12,]; RP3 <- ex[13,]; RM3 <- ex[14,]

refstrides1 <- rbind(LP1,LP2,LM2,LM3)
rownames(refstrides1) <- c("LP1","LP2","LM2","LM3")	#prefered set of strides
refstrides2 <- rbind(RP1,RP2,RM2,RM3)
rownames(refstrides2) <- c("RP1","RP2","RM2","RM3")	#backup set of strides (second choice)
if(anyNA(LP1,LP2) == FALSE | anyNA(LM2,LM3) == FALSE){
refstrides <- refstrides1
}else{refstrides <- refstrides2}
ang <- avlineang(refstrides)
new_LM3 <- rot(LP2,LM3,ang)		#rotate LM3 around LP2
gad[i] <- dif(new_LM3[2],LP2[2])	#difference between y-values of rotated points is the stride-parallel LR_p measurement
}}
gad <- unlist(gad); gad <- list(gad)
names(gad) <- comment(x)[1]
return(gad)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- gad3(x)}else{
fl <- sapply(x,gad3)}
return(fl)
}


#GAD4 diagonal sequence
#takes single trackway or list of trackways
GAD4 <- function(x,method="both"){
gad4 <- function(x){
#method options:
#both: Projecting mid-stride of RM1 based on both half the stride length RM1-RM2 and the position of LM2.
#opposite: Only taking into account the position of the opposite manus (LM2); this is the traditional approach following Leonardi.
#same: Only taking into account the actual stride (/2), RM1-RM2.
gad <- rep(NA,nrow(x))	#initialise empty vector of set length
negposC <- vector()
for(i in 1:(nrow(x)-3)){
even <- is.even(i)
if(even == FALSE){
ex <- rettrack(x,i)	#writes coordinates for LP1, LM1 etc., with LP2 as first track
#find orientation of trackway midline
RP0 <- ex[1,]; RM0 <- ex[2,]; LP1 <- ex[3,]; LM1 <- ex[4,]; RP1 <- ex[5,]; RM1 <- ex[6,]; LP2 <- ex[7,]; LM2 <- ex[8,]; RP2 <- ex[9,]; RM2 <- ex[10,]; LP3 <- ex[11,]; LM3 <- ex[12,]; RP3 <- ex[13,]; RM3 <- ex[14,]
refstride1 <- rbind(LM2,LM3)
refstride2 <- rbind(RP1,RP2,LP2,LP3)
if(anyNA(refstride1) == FALSE){
refstrides <- refstride1		#refstride 1 is preferred
}else{
refstrides <- refstride2
}
ang <- avlineang(refstrides)	#average angle of selected strides to vertical
# Rotate into vertical position
LP1 <- rot(LP2,LP1,ang)
LM1 <- rot(LP2,LM1,ang)
RP1 <- rot(LP2,RP1,ang)
RM1 <- rot(LP2,RM1,ang)
LM2 <- rot(LP2,LM2,ang)
RP2 <- rot(LP2,RP2,ang)
RM2 <- rot(LP2,RM2,ang)
LP3 <- rot(LP2,LP3,ang)
LM3 <- rot(LP2,LM3,ang)
# find points between which the GAD is measured
coo_pes <- mean(c(LP2[[2]],RP2[[2]]))
coo_same <- mean(c(LM2[[2]],LM3[[2]]))
coo_opposite <- RM2[[2]]
# Measure according to specified method
if(method == "both"){
coo_manus <- meanN(c(coo_same,coo_opposite))
}
if(method == "same"){
coo_manus <- coo_same
}
if(method == "opposite"){
coo_manus <- coo_opposite
}
gad[i] <- dif(coo_manus,coo_pes)
LP2 <<- LP2
RP2 <<- RP2
#midpoint between LP1 and RP1
LPmid <- c(mean(c(LP2[[1]],RP2[[1]])),mean(c(LP2[[2]],RP2[[2]])))
#orientation:
#1: trackway leading to the top
#-1: trackway leading to the bottom
neeg <- 0
try(if(LP2[[2]] > LP1[[2]]){neeg <- neeg+1}else{neeg <- neeg-1},silent=TRUE)
try(if(LP3[[2]] > LP2[[2]]){neeg <- neeg+1}else{neeg <- neeg-1},silent=TRUE)
try(if(RP2[[2]] > RP1[[2]]){neeg <- neeg+1}else{neeg <- neeg-1},silent=TRUE)
if(neeg > 0){
orientation <- 1
}else{
orientation <- -1
}
if(anyNA(c(RM2,LPmid)) == FALSE){
if(RM2[[2]] > LPmid[[2]]){
negpos <- 1
}else{
negpos <- -1
}
}else{
negpos <- NA
}
gad[i] <- gad[i]*negpos*orientation
if(is.nan(gad[i])){
gad[i] <- NA
}
}}
gad <- unlist(gad); gad <- list(gad)
names(gad) <- comment(x)[1]
return(gad)
}
if(is.data.frame(x) == TRUE | is.matrix(x) == TRUE){
fl <- gad4(x)}else{
fl <- sapply(x,gad4)}
return(fl)
}


#!!!!!!!Special functions for gait analysis!!!!!!!

# correct for GAD2 function
GAD2acc <- function(x){
if(is.null(nrow(x))){
x <- t(matrix(x))
}
for(i in 1:nrow(x)){
P1 <- c(50,x[i,1])
P2 <- c(100,x[i,2])
P3 <- c(75,x[i,3])
PM <- rbind(P1,P2,P3)
mod1 <- lm(PM[,2] ~ PM[,1])
res <- mod1$residuals
for(o in 1:ncol(x)){
x[i,o] <- x[i,o]-res[o]
}}
return(x)
}

#sort strides according to step cycles (LP1-LP2, RM1-RM2, etc.)
#accepts single trackway
sortstride <- function(x){
x <- unlist(stride(x))
#expand vector
xx <- c(NA,NA,x,NA,NA)
Stride <- rep(NA,length(xx))
for(i in 1:length(x)){
if(is.na(x[i]) == FALSE){
if(is.even(i) == TRUE){  #if manus
Stride[i] <- xx[i+2]
}else{	#if pes
Stride[i+2] <- xx[i+2]
}
}}
Stride <- na.trim(Stride,sides="right")	#remove trailing NAs and leading NAs

svg("strideplot.svg")
plot(Stride,type="b",pch=19)
xes <- 1:length(Stride)
rightpes <- Stride[seq(1,length(Stride),4)]
rightpesX <- xes[seq(1,length(xes),4)]
rightmanus <- Stride[seq(1,length(Stride),4)+1]
rightmanusX <- xes[seq(1,length(xes),4)+1]
leftpes <- Stride[seq(1,length(Stride),4)+2]
leftpesX <- xes[seq(1,length(xes),4)+2]
leftmanus <- Stride[seq(1,length(Stride),4)+3]
leftmanusX <- xes[seq(1,length(xes),4)+3]

points(rightpesX,rightpes,pch=15,col="red",cex=1.5)
points(rightmanusX,rightmanus,pch=16,col="red",cex=1.5)
points(leftpesX,leftpes,pch=15,col="blue",cex=1.5)
points(leftmanusX,leftmanus,pch=16,col="blue",cex=1.5)

fit <- lm(Stride ~ poly(xes,round(length(xes)/3.7),raw=TRUE))
fit <<- fit
Stride <<- Stride
lines(xes,predict(fit,data.frame(x=xes)),col="red",lwd=2)
dev.off()
return(fit)
}

# autstepcomp <- function(){
# source("gaitcalculator.R")
# svg("hist-before.svg")
# gaitcalculator(a,stepcomp=FALSE)
# dev.off()
# svg("trackway.svg")
# trackplot(a[[1]])
# dev.off()
# svg("strides-before.svg")
# correct <- stepcomp(a[[1]],dev=3)
# dev.off()
# svg("strides-after.svg")
# sortstride(correct)
# dev.off()
# svg("hist-after.svg")
# gaitcalculator(a)
# dev.off()
# }

#step compensation corrector
stepcomp <- function(x,dev=2){
fit <- sortstride(x)
res <- fit$residuals
Stride <- unlist(stride(x))
#fill vector of residuals with NAs to match original (Stride) data
o <- 0
restmp <- res
for(i in 1:length(res)){
o <- o+1
while(names(res[i]) > o){
if(i > 1){
restmp <- c(restmp[1:(o-1)],NA,restmp[o:length(restmp)])
o <- o+1
}else{
restmp <- c(NA,restmp)
o <- o+1
}}}
res <- restmp
res <<- res
#detect compensation events
newcoords <- x
checkvec <- rep(FALSE,length(res))
z <- 0 	#markes if a footprint has already been moved and should not moved twice
for(i in 1:length(res)){
LP2N <- vector()
LP2 <- vector()
if(i > 2){
if(any(z == i) == FALSE){#was this footprint already adjusted? If yes, skip it to avoid amplification.
if(i < length(res)-3 & is.na(res[i]) == FALSE & is.na(res[i+4]) == FALSE){
mres <- meanN(c(abs(res[i-1]),abs(res[i-2]),abs(res[i+1]),abs(res[i+2])))
#mres2 <- meanN(c(abs(res[i+4-1]),abs(res[i+4-2]),abs(res[i+4+1]),abs(res[i+4+2])))
if(is.na(mres) == FALSE){
if(abs(res[i]) > mres*dev){
# if(is.na(mres) == FALSE & is.na(mres2) == FALSE){
# if(abs(res[i]) > mres*dev | abs(res[i+4]) > mres2*dev){
if((res[i] > 0 & res[i+4] < 0) | (res[i] < 0 & res[i+4] > 0)){

checkvec[i] <- TRUE
checkvec <<- checkvec
####implement shortening the second stride if negative, otherwise we get error in turning trackways
#if(res[i] > 0){		#if positiv, shorten first stride
if(is.even(i) == FALSE){	#if pes
#calculate back in original footprint order (non-stepcycle compliant order)
LP1 <- x[i-2,]	#original coordinates
LP2 <- x[i+2,]
}else{	#if manus
LP1 <- x[i,]	#original coordinates
LP2 <- x[i+4,]
}
#magnitude of shortening: the stride length that deviates less
wmax <- which(c(abs(res[i]),abs(res[i+4])) == max (c(abs(res[i]),abs(res[i+4]))))
if(wmax == 1){
mag <- res[i+4]
}else{
mag <- res[i]
}
if(res[i] > 0){			#if stride too large, shorten it
xdev <- (((LP2[1]-LP1[1])/Stride[i])*abs(mag))*-1		#new x coordinate
ydev <- (((LP2[2]-LP1[2])/Stride[i])*abs(mag))*-1		#new y coordinate
}else{					#else if too short, make it longer
xdev <- ((LP2[1]-LP1[1])/Stride[i])*abs(mag)		#new x coordinate
ydev <- ((LP2[2]-LP1[2])/Stride[i])*abs(mag)		#new y coordinate
}
#if(LP2[2] > LP1[2]){								#if the trackway goes down the screen and not up
LP2N <- cbind(LP2[1]+xdev,LP2[2]+ydev)
# }else{
# LP2N <- cbind(LP2[1]+xdev,LP2[2]+ydev)
# }
if(is.even(i) == FALSE){
newcoords[i+2,] <- NA								#write new coords
}else{
newcoords[i+4,] <- NA
}
z <- c(z,i+4)


}
}

}}}}}
#}
return(newcoords)
}

#get table with parallel pes-manus distance and associated mean stride
corrtable <- function(x){
strides <- unlist(stride(x))
ppmds <- unlist(ppmd(x))
Table <- vector()
for(i in 1:(nrow(x)-1)){
#first check if all required data are available
if(i > 4 & i < nrow(x) & anyNA(c(x[i,],x[i+1,],x[i-3,],x[i-4,])) == FALSE & is.even(i) == FALSE){
Table <- rbind(Table,c(mean(c(strides[i-3],strides[i-4])),ppmds[i]))
}else{
Table <- rbind(Table,c(NA,NA))
}
Table <- rbind(Table,c(NA,NA))
colnames(Table) <- c("mean str","ppmd")
}
Table <- Table[complete.cases(Table),]
return(Table)
}


#get table with parallel pes-manus distance and associated stride. Pes stride only
corrtableSingle <- function(x){
strides <- unlist(stride(x))
ppmds <- unlist(ppmd(x))
Table <- vector()
for(i in 1:(nrow(x)-1)){
#first check if all required data are available
if(i > 4 & i < nrow(x) & anyNA(c(x[i,],x[i+1,],x[i-3,],x[i-4,])) == FALSE & is.even(i) == FALSE){
Table <- rbind(Table,c(strides[i-4],ppmds[i]))
}else{
Table <- rbind(Table,c(NA,NA))
}
Table <- rbind(Table,c(NA,NA))
colnames(Table) <- c("pes str","ppmd")
}
Table <- Table[complete.cases(Table),]
return(Table)
}


#GAD2_B
#Calculate GAD2 from single pes to half-distance between manus pair
GAD2_B <- function(x){
#method options:
#both: Projecting mid-stride of RM1 based on both half the stride length RM1-RM2 and the position of LM2.
#opposite: Only taking into account the position of the opposite manus (LM2); this is the traditional approach following Leonardi.
#same: Only taking into account the actual stride (/2), RM1-RM2.
gad <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-5)){
even <- is.even(i)
if(even == FALSE){
ex <- rettrack(x,i)	#writes coordinates for LP1, LM1 etc., with LP2 as first track
#find orientation of trackway midline

RP0 <- ex[1,]; RM0 <- ex[2,]; LP1 <- ex[3,]; LM1 <- ex[4,]; RP1 <- ex[5,]; RM1 <- ex[6,]; LP2 <- ex[7,]; LM2 <- ex[8,]; RP2 <- ex[9,]; RM2 <- ex[10,]; LP3 <- ex[11,]; LM3 <- ex[12,]; RP3 <- ex[13,]; RM3 <- ex[14,]

refstride1 <- rbind(RM2,RM3)
refstride2 <- rbind(RP1,RP2,LP2,LP3)
if(anyNA(refstride1) == FALSE){
refstrides <- refstride1		#refstride 1 is preferred
}else{
refstrides <- refstride2
}
ang <- avlineang(refstrides)	#average angle of selected strides to vertical

# Rotate into vertical position
RP2 <- rot(LP2,RP2,ang)
RM2 <- rot(LP2,RM2,ang)
LM3 <- rot(LP2,LM3,ang)
RM3 <- rot(LP2,RM3,ang)

# find points between which the GAD is measured
coo_pes <- LP2[[2]]
coo_manus <- mean(c(RM2[[2]],LM3[[2]]))
gad[i] <- dif(coo_manus,coo_pes)
if(is.nan(gad[i])){
gad[i] <- NA
}
}}
return(gad)
}

#GAD3_B
#Calculate GAD3 based on four footprints
GAD3_B <- function(x){
gad <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-5)){
even <- is.even(i)
if(even == FALSE){
#find center between x and y values of both the manus and pes pairs
co_pes <- c(mean(c(x[i,1],x[i+2,1])),mean(c(x[i,2],x[i+2,2])))
co_manus <- c(mean(c(x[i+5,1],x[i+7,1])),mean(c(x[i+5,2],x[i+7,2])))
gad[i] <- ed(co_pes,co_manus)
}}
return(gad)
}

#GAD3_0
#GAD3 assuming that the manus is directly ahead of the pes with no step cycle in-between
GAD3_0 <- function(x){
gad <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-1)){
even <- is.even(i)
if(even == FALSE){	#only uneven rows = only the pes prints, not the manus prints
ex <- rettrack(x,i)	#writes coordinates for LP1, LM1 etc., with LP2 as first track
RP0 <- ex[1,]; RM0 <- ex[2,]; LP1 <- ex[3,]; LM1 <- ex[4,]; RP1 <- ex[5,]; RM1 <- ex[6,]; LP2 <- ex[7,]; LM2 <- ex[8,]; RP2 <- ex[9,]; RM2 <- ex[10,]; LP3 <- ex[11,]; LM3 <- ex[12,]; RP3 <- ex[13,]; RM3 <- ex[14,]
refstrides1 <- rbind(LP1,LP2,LM1,LM2)
rownames(refstrides1) <- c("LP1","LP2","LM1","LM2")	#prefered set of strides
refstrides2 <- rbind(RP1,RP2,RM2,RM3)
rownames(refstrides2) <- c("RP1","RP2","RM2","RM3")	#backup set of strides (second choice)
if(anyNA(LP1,LP2) == FALSE | anyNA(LM1,LM2) == FALSE){
refstrides <- refstrides1
}else{refstrides <- refstrides2}
ang <- avlineang(refstrides)
new_LM2 <- rot(LP2,LM2,ang)		#rotate LM3 around LP2
gad[i] <- dif(new_LM2[2],LP2[2])	#difference between y-values of rotated points is the stride-parallel LR_p measurement
}}
return(gad)
}

#GAD3_0_B
#Calculate GAD3_0 based on four footprints
GAD3_0_B <- function(x){
gad <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-5)){
even <- is.even(i)
if(even == FALSE){
#find center between x and y values of both the manus and pes pairs
co_pes <- c(mean(c(x[i,1],x[i+2,1])),mean(c(x[i,2],x[i+2,2])))
co_manus <- c(mean(c(x[i+1,1],x[i+3,1])),mean(c(x[i+1,2],x[i+3,2])))
gad[i] <- ed(co_pes,co_manus)
}}
return(gad)
}

#LR (stride-parallel distance between LP1 and RM1, the two feet that hit the ground simultaneously in a trot)
GAD1_B <- function(x){
lrp <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-3)){

even <- is.even(i)
if(even == FALSE){	#only uneven rows = only the pes prints, not the manus prints
#check which prints are available for reference
#the first two are required in any case:
ex <- rettrack(x,i)	#writes coordinates for LP1, LM1 etc., with LP2 as first track
RP0 <- ex[1,]; RM0 <- ex[2,]; LP1 <- ex[3,]; LM1 <- ex[4,]; RP1 <- ex[5,]; RM1 <- ex[6,]; LP2 <- ex[7,]; LM2 <- ex[8,]; RP2 <- ex[9,]; RM2 <- ex[10,]; LP3 <- ex[11,]; LM3 <- ex[12,]; RP3 <- ex[13,]; RM3 <- ex[14,]
refstrides1 <- rbind(RM1,RM2,LP1,LP2)	#prefered set of strides
rownames(refstrides1) <- c("RM1","RM2","LP1","LP2")
refstrides2 <- rbind(RP1,RP2,LM2,LM3)	#backup set of strides (second choice)
rownames(refstrides2) <- c("RP1","RP2","LM2","LM3")
if(anyNA(RM1,RM2) == FALSE | anyNA(LP1,LP2) == FALSE){
refstrides <- refstrides1
}else{refstrides <- refstrides2}
ang <- avlineang(refstrides)

if(is.na(ang) == TRUE){
allrefs <- rbind(RM1=RM1,RM2=RM2,LP1=LP1,LP2=LP2,RP1=RP1,RP2=RP2,LM2=LM2,LM3=LM3,LP2=LP2,LP3=LP3,RM2=RM2,RM3=RM3)

refstrides <- vector()
for(o in 1:(nrow(allrefs)-1)){
if(is.even(o) == FALSE){
if(anyNA(c(allrefs[o,],allrefs[o+1,])) == FALSE){
refstrides <- rbind(refstrides,allrefs[o,],allrefs[o+1,])
}}}
if(is.logical(refstrides) == FALSE){
ang <- avlineang(refstrides)	#find average angle against the vertical of selected strides
}else{
ang <- NA
}}

refstrides <<- refstrides

new_RM2 <- rot(LP2,RM2,ang)		#rotate RM2 around LP2
lrp[i] <- dif(new_RM2[2],LP2[2])	#difference between y-values of rotated points is the stride-parallel LR_p measurement.

#check if positive or negative
negpos <- 0

if(anyNA(c(LP3,LP2,RM2)) == FALSE){
LP2_LP3 <- ed(LP2,LP3)
RM2_LP3 <- ed(RM2,LP3)
if(LP2_LP3 > RM2_LP3){
negpos <- negpos+1
}else{negpos <- negpos-1}
}

if(anyNA(c(LM3,LP2,RM2)) == FALSE){
LP2_LM3 <- ed(LP2,LM3)
RM2_LM3 <- ed(RM2,LM3)
if(LP2_LM3 > RM2_LM3){
negpos <- negpos+2
}else{negpos <- negpos-2}
}

if(anyNA(c(RM2,LP2,RP3)) == FALSE){
LP2_RP3 <- ed(LP2,RP3)
RM2_RP3 <- ed(RM2,RP3)
if(LP2_RP3 > RM2_RP3){
negpos <- negpos+2
}else{negpos <- negpos-2}
}

if(anyNA(c(LP1,RM2,LP2)) == FALSE){
LP2_LP1 <- ed(LP2,LP1)
RM2_LP1 <- ed(RM2,LP1)
if(LP2_LP1 < RM2_LP1){
negpos <- negpos+2
}else{negpos <- negpos-2}
}

if(anyNA(c(RP1,RM2,LP2)) == FALSE){
LP2_RP1 <- ed(LP2,RP1)
RM2_RP1 <- ed(RM2,RP1)
if(LP2_RP1 < RM2_RP1){
negpos <- negpos+1
}else{negpos <- negpos-1}
}

if(anyNA(c(RM2,LP2,LM1)) == FALSE){
LP2_LM1 <- ed(LP2,LM1)
RM2_LM1 <- ed(RM2,LM1)
if(LP2_LM1 < RM2_LM1){
negpos <- negpos+1
}else{negpos <- negpos-1}
}

if(negpos < 0){		#if index is negative, need to switch sign
lrp[i] <- lrp[i]*(-1)}}
}
return(lrp)
}


#GAD50 based on midstride points
GAD50_midstride <- function(x){
lrp <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-7)){
even <- is.even(i)
if(even == FALSE){	#only uneven rows = only the pes prints, not the manus prints
LP1 <- x[i,]; LP2 <- x[i+4,]; RM1 <- x[i+3,]; RM2 <- x[i+7,]

if(anyNA(c(LP1,LP2,RM1,RM2)) == FALSE){
LP1 <<- LP1
LP2 <<- LP2
RM1 <<- RM1
RM2 <<- RM2

refstrides <- rbind(LP1,LP2,RM1,RM2)
ang <- avlineang(refstrides)	#find average angle against the vertical of selected strides
midpes <- c(mean(c(LP1[[1]],LP2[[1]])),mean(c(LP1[[2]],LP2[[2]])))
midmanus_raw <- c(mean(c(RM1[[1]],RM2[[1]])),mean(c(RM1[[2]],RM2[[2]])))
midmanus <- rot(midpes,midmanus_raw,ang)
lrp[i] <- dif(midmanus[[2]],midpes[2])
}
}
}
return(lrp)
}


#GAD50 based on midstride points
GAD0_midstride <- function(x){
lrp <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 1:(nrow(x)-9)){
even <- is.even(i)
if(even == FALSE){	#only uneven rows = only the pes prints, not the manus prints
LP1 <- x[i,]; LP2 <- x[i+4,]; LM2 <- x[i+5,]; LM3 <- x[i+9,]

if(anyNA(c(LP1,LP2,LM2,LM3)) == FALSE){
refstrides <- rbind(LP1,LP2,LM2,LM3)
ang <- avlineang(refstrides)	#find average angle against the vertical of selected strides
midpes <- c(mean(c(LP1[[1]],LP2[[1]])),mean(c(LP1[[2]],LP2[[2]])))
midmanus_raw <- c(mean(c(LM2[[1]],LM3[[1]])),mean(c(LM2[[2]],LM3[[2]])))
midmanus <- rot(midpes,midmanus_raw,ang)
lrp[i] <- dif(midmanus[[2]],midpes[2])
}
}
}
return(lrp)
}




#LR with no step cycle in-between
GAD1_B_0 <- function(x){
lrp <- rep(NA,nrow(x))	#initialise empty vector of set length
for(i in 2:(nrow(x))){
even <- is.even(i)
if(even == FALSE){	#only uneven rows = only the pes prints, not the manus prints
#check which prints are available for reference
#the first two are required in any case:
ex <- rettrack(x,i)	#writes coordinates for LP1, LM1 etc., with LP2 as first track
RP0 <- ex[1,]; RM0 <- ex[2,]; LP1 <- ex[3,]; LM1 <- ex[4,]; RP1 <- ex[5,]; RM1 <- ex[6,]; LP2 <- ex[7,]; LM2 <- ex[8,]; RP2 <- ex[9,]; RM2 <- ex[10,]; LP3 <- ex[11,]; LM3 <- ex[12,]; RP3 <- ex[13,]; RM3 <- ex[14,]
refstrides1 <- rbind(RM0,RM1,LP1,LP2)	#prefered set of strides
rownames(refstrides1) <- c("RM0","RM1","LP1","LP2")
refstrides2 <- rbind(RP1,RP2,LM1,LM2)	#backup set of strides (second choice)
rownames(refstrides2) <- c("RP1","RP2","LM1","LM2")
if(anyNA(RM0,RM1) == FALSE | anyNA(LP1,LP2) == FALSE){
refstrides <- refstrides1
}else{refstrides <- refstrides2}
ang <- avlineang(refstrides)	#find average angle against the vertical of selected strides
new_RM1 <- rot(LP2,RM1,ang)		#rotate RM2 around LP2
lrp[i] <- dif(new_RM1[2],LP2[2])	#difference between y-values of rotated points is the stride-parallel LR_p measurement.
}}
return(lrp)
}

GAD1_M <- function(x){
A <- GAD1(x)
B <- GAD1_B(x)
GAD <- cbind(A,B)
GAD <- apply(GAD,1,meanN)
for(z in 1:length(GAD)){
if(is.nan(GAD[z]) == TRUE){
GAD[z] <- NA
}}
return(GAD)
}

GAD2_M <- function(x){
A <- GAD2(x)
B <- GAD2_B(x)
GAD <- cbind(A,B)
GAD <- apply(GAD,1,meanN)
for(z in 1:length(GAD)){
if(is.nan(GAD[z]) == TRUE){
GAD[z] <- NA
}}
return(GAD)
}

GAD3_M <- function(x){
A <- GAD3(x)
B <- GAD3_B(x)
GAD <- cbind(A,B)
GAD <- apply(GAD,1,meanN)
for(z in 1:length(GAD)){
if(is.nan(GAD[z]) == TRUE){
GAD[z] <- NA
}}
return(GAD)
}

GAD2_F <- function(x){
#calculate GAD2 based on five points

gad <- rep(NA,nrow(x))	#initialise empty vector of set length
if(nrow(x) > 7){
for(i in 1:(nrow(x)-7)){
even <- is.even(i)
if(even == FALSE){
#find center between x and y values of both the manus and pes pairs
co_pes <- c(mean(c(x[i,1],x[i+2,1])),mean(c(x[i,2],x[i+2,2])))

#midpoint RM1-RM2, taking LM2 into account
midRM <- c(mean(c(x[i+3,1],x[i+7,1])),mean(c(x[i+3,2],x[i+7,2])))
co_manus <- c(mean(c(midRM[1],x[i+5,1])),mean(c(midRM[2],x[i+5,2])))
gad[i] <- ed(co_pes,co_manus)
}}}
return(gad)
}


#produce matrix of all measurements
trackmatrix <- function(x){
tm <- as.data.frame(cbind(stride(x),pace(x),GAD1(x),GAD1_B(x),GAD2(x),GAD2_B(x),GAD3(x),GAD3_B(x)))
colnames(tm) <- c("stride","pace","GAD1","GAD1_B","GAD2","GAD2_B","GAD3","GAD3_B")

GAD1_M <- rep(NA,nrow(tm))
GAD2_M <- GAD1_M
GAD3_M <- GAD1_M

tm <<- tm

for(i in 1:nrow(tm)){
GAD1_M[i] <- meanN(c(tm[[i,"GAD1"]],tm[[i,"GAD1_B"]]))	#GAD1, mean of both approaches (GAD1 and GAD1_B)
GAD2_M[i] <- meanN(c(tm[[i,"GAD2"]],tm[[i,"GAD2_B"]]))
GAD3_M[i] <- meanN(c(tm[[i,"GAD3"]],tm[[i,"GAD3_B"]]))
}

tm <- as.data.frame(cbind(tm[,"stride"],tm[,"pace"],tm[,"GAD1"],tm[,"GAD1_B"],GAD1_M,tm[,"GAD2"],tm[,"GAD2_B"],GAD2_M,tm[,"GAD3"],tm[,"GAD3_B"],GAD3_M))
colnames(tm) <- c("stride","pace","GAD1","GAD1_B","GAD1_M","GAD2","GAD2_B","GAD2_M","GAD3","GAD3_B","GAD3_M")

return(tm)
}

#get simple GAD1 (GAD1_B) associated with median stride
corrtable2 <- function(x){
strides <- unlist(stride(x))
GADs <- unlist(GAD1(x))
Table <- vector()
for(i in 1:(nrow(x)-3)){
i <<- i
#first check if all required data is available
if(is.even(i) == FALSE & i > 4){	#if pes …
#calculate median stride based on 6 strides (if all are available)
Table <- rbind(Table,c(mean(c(strides[i-4],strides[i-1])),GADs[i]))
}else{
Table <- rbind(Table,c(NA,NA))
}
colnames(Table) <- c("median str","GAD1")
}
Table <- Table[complete.cases(Table),]
return(Table)
}

#Stride comparison approach (Hypothesis 1)

stridevar <- function(x,gads,quiet=FALSE){
strides <- unlist(stride(x))
diffpace <- vector()
difftrot <- vector()
for(i in 1:(nrow(x)-8)){
if(is.even(i) == FALSE){
#pace
ddpace <- dif(strides[i],strides[i+5])
diffpace <- c(diffpace,ddpace)
#trot
ddtrot <- dif(strides[i],strides[i+3])
difftrot <- c(difftrot,ddtrot)
}}
bindd <- cbind(difftrot,diffpace)
bindd <- bindd[complete.cases(bindd),]
bindd <<- bindd
if(is.vector(bindd == TRUE)){
bindd <- matrix(bindd,ncol=2)}
strot <- sum(bindd[,1],na.rm=TRUE)
space <- sum(bindd[,2],na.rm=TRUE)

if(quiet == FALSE){
if(nrow(bindd) > 3){
if(strot > 0 & space > 0){
print("==== Stride comparison approach (Evaluation approach 1) ====")
Perc <- 100-((100/max(c(strot,space)))*min(c(strot,space)))
#print(c("Trot: ",strot,"n:",length(difftrot)))
#print(c("Pace: ",space,"n:",length(diffpace)))
print(c("diff percent: ",round(Perc,2)))
sc.res <- NA
if(strot < space){
print("Trot more likely than Pace")
ttest <- t.test(bindd[,1],bindd[,2],alternative="less")
scttest <<- ttest$p.value
sc.res <<- "trot"
print(c("T-test: ",round(ttest$p.value,2)))
}else{if(space < strot){
print("Pace more likely than Trot")
ttest <- t.test(bindd[,1],bindd[,2],alternative="greater")
scttest <<- ttest$p.value
sc.res <<- "pace gait"
print(c("T-test: ",round(ttest$p.value,2)))
}
}
}}else{
print("Data insufficient")
}
}

bindd <- bindd/meanN(strides)

#normalize differences to range 1:0
# n <- nrow(bindd)
# norm <- normalize(c(bindd[,1],bindd[,2]),method="range")
# bindd <- cbind(norm[1:n],norm[(n+1):length(norm)])
return(bindd)
}

#Distance variation approach (Hypothesis 2).
#Based on MAD (deviation around median), non-parametric
distvar <- function(x){
#LP1-RM1
trot <- vector()
pace <- vector()
for(i in 1:(nrow(x)-5)){
if(is.even(i) == FALSE){
distT <- unlist(ed(x[i,],x[i+3,]))
trot <- c(trot,distT)
distP <- unlist(ed(x[i,],x[i+5,]))
pace <- c(pace,distP)
}}
trotpace <- cbind(trot,pace)
trotpace <<- trotpace
trotpace <- trotpace[complete.cases(trotpace),]
trotvar <- mad(unlist(trotpace[,1]))
pacevar <- mad(unlist(trotpace[,2]))
trotvarV <- var(unlist(trotpace[,1]))
pacevarV <- var(unlist(trotpace[,2]))
print(c("Trot","mad",trotvar,"var",trotvarV))
print(c("Pace:","mad",pacevar,"var",pacevarV))
if(trotvar < pacevar){print("Trot!")}else{if(pacevar < trotvar){print("Pace!")}}
#Levene's test
levene <- levtest(trotpace[,1],trotpace[,2])
levene <<- levene
print(c("Levene test",levene$p.value))
}


#Distance variation approach (Hypothesis 2)
#As "distvar", but based on stride-parallel measurements
#Based on MAD (deviation around median), non-parametric

distvar2 <- function(x){
gads <- cbind(unlist(GAD1_B(x)),unlist(GAD3(x)))
gads <- gads[complete.cases(gads),]
gad1var <- mad(gads[,1])
gad3var <- mad(gads[,2])
print(c("GAD1:",gad1var))
print(c("GAD3:",gad3var))
if(gad1var < gad3var){print("Trot!")}
if(gad3var < gad1var){print("Pace!")}
}







##GAD var

gadvar <- function(x,GAD2=TRUE){
if(GAD2==TRUE){
gads <- cbind(unlist(GAD1_B(x)),unlist(GAD3(x)),unlist(GAD2(x)))
}
if(GAD2==FALSE){
gads <- cbind(unlist(GAD1_B(x)),unlist(GAD3(x)))
}
output <- list()

if(GAD2==TRUE){
#GAD2 (third column) missing?
for(m in 1:nrow(gads)){
if(is.na(gads[m,3]) == TRUE & anyNA(gads[m,1:2]) == FALSE){
gads[m,3] <- mean(c(gads[m,1],gads[m,2]))
}}
}

#remove NA
gads <- gads[complete.cases(gads),]
if(is.vector(gads) == FALSE){
if(nrow(gads) > 3){
if(GAD2==TRUE){
#correct using GAD2
gads <- GAD2acc(gads)
#resort into logical order
gads <- cbind(gads[,1],gads[,3],gads[,2])
}

gads <<- gads

}
}
}




#GAD distribution vs stride
gaddistri<- function(x){
require("MASS")

gads1 <- vector()
gads3 <- vector()
for(i in 1:length(x)){
gads <- cbind(unlist(GAD1_B(x[[i]])),unlist(GAD3(x[[i]])),unlist(stride(x[[i]])))
#gads <<- gads

gads <- apply(gads,2,scale)

#data frame

one <- data.frame(cbind(gads[,1],gads[,3],"GAD1"))
two <- data.frame(cbind(gads[,2],gads[,3],"GAD2"))
gadS <- rbind(one,two)
gadS[,3] <- as.factor(gadS[,3])
colnames(gadS) <- c("GAD","Stride","Factor")
gadS <- gadS[complete.cases(gadS),]
gadS[,1] <- as.numeric(gadS[,1])
gadS[,2] <- as.numeric(gadS[,2])
mod1 <- aov(GAD ~ Stride*Factor,data=gadS)
summ <- summary(mod1)
print(summary(mod1))
print(c(comment(x[[i]])[1],"p:",summ[[1]]$"Pr(>F)"[3]))

Med <- medianN(gads[,3])
gads[,1] <- gads[,1]/Med
gads[,2] <- gads[,2]/Med
#GAD1
gadsB1 <- vector()
if(nrow(gads) > 4){
for(i in 1:nrow(gads)){
if(is.even(i) == FALSE){
if(i > 4){im4 <- gads[i-4,3]}else{im4 <- NA}
stridesT <- meanN(c(im4,gads[i-1,3]))
gadsB1 <- rbind(gadsB1,c(gads[i,1],stridesT))
}}}
gads1 <- rbind(gads1,gadsB1[complete.cases(gadsB1),])
#GAD3
gadsB3 <- vector()
if(nrow(gads) > 4){
for(i in 1:nrow(gads)){
if(is.even(i) == FALSE){
if(i > 4){im4 <- gads[i-4,3]}else{im4 <- NA}
if(i < nrow(gads)){im5 <- gads[i+1,3]}else{im5 <- NA}
stridesT <- meanN(c(im4,im5))
gadsB3 <- rbind(gadsB3,c(gads[i,2],stridesT))
}}}
gads3 <- rbind(gads3,gadsB3[complete.cases(gadsB3),])
}
lm1 <- rlm(gads1[,1] ~ gads1[,2])
lm3 <- rlm(gads3[,1] ~ gads3[,2])
svg("gaddistri.svg")
plot(c(gads1[,2],gads3[,2]),c(gads1[,1],gads3[,1]),col="white",xlab="Stride",ylab="GAD / Stride")
points(gads1[,2],gads1[,1],pch=19,col="darkgreen")
points(gads3[,2],gads3[,1],pch=19,col="red")
abline(lm1)
abline(lm3)
dev.off()
#print(summary(lm3))
lmGAD3 <<- lm1

}


#PCA analysis
PCAgait <- function(x){
GAD1 <- unlist(GAD1_B(x[[1]]))
GAD1 <<- GAD1
GAD3 <- unlist(GAD3(x))
Stride <- unlist(stride(x))
Stride <<- Stride
PPMD <- unlist(ppmd(x))

dattable <- vector()
if(length(GAD1) > 4){
for(i in 2:(length(GAD1)-1)){
if(is.even(i) == FALSE){
strides <- mean(c(Stride[i-1],Stride[i+1]))
dattable <- rbind(dattable,c(GAD1[i],GAD3[i],strides))
}}}
dattable <- dattable[complete.cases(dattable),]
colnames(dattable) <- c("GAD1 / Median","GAD3 / Median","Stride")
dattable <<- dattable
pca <- prcomp(dattable,scale=TRUE)
print(pca)
print(summary(pca))
pca <<- pca
}

#PMD-Stride comparison
PMDStride <- function(x){
PMD <- pmd(a)
STRIDE <- stride(a)
datt <- vector()
for(i in 1:length(PMD)){datt <- rbind(datt,cbind(STRIDE[[i]],PMD[[i]]))}
DD <- datt[complete.cases(datt),]
DDD <- apply(DD,1,function(x){x[2]/x[1]})
DD[,2] <- DDD
plot(DD)
}

#exact GADvar

exact_GADvar <- function(GADS,quiet=FALSE,GAD2=TRUE){
x <- GADS
#excepts GADS list of corrected gads
#if GAD2=FALSE: Only two columns should be present (GAD1 and GAD3)
#check if all gads are normally distributed:
shaplist <- apply(gads,2,shapiro.test)
if(GAD2==TRUE){
plist <- c(shaplist[[1]]$p.value,shaplist[[2]]$p.value,shaplist[[3]]$p.value)
}
if(GAD2==FALSE){
plist <- c(shaplist[[1]]$p.value,shaplist[[2]]$p.value)
}
if(any(plist < 0.05) == TRUE){
distribution <- "non-normal"}else{
distribution <- "normal"}

#extended GAD list, with 50 columns, one for each limb phase:
if(GAD2==TRUE){
gadtable <- vector()
for(i in 1:nrow(x)){
interval <- (x[i,2]-x[i,1])/24
firsthalf <- vector(length=24)
firsthalf[1] <- x[i,1]
for(o in 2:length(firsthalf)){
firsthalf[o] <- firsthalf[o-1]+interval
}
firsthalf <- c(firsthalf,x[i,2])
interval2 <- (x[i,3]-x[i,2])/25
secondhalf <- vector(length=25)
secondhalf[1] <- x[i,2]+interval2
for(p in 2:length(secondhalf)){
secondhalf[p] <- secondhalf[p-1]+interval2
}
gadtable <- rbind(gadtable,c(firsthalf,secondhalf))
}
}
if(GAD2==FALSE){
gadtable <- vector()
for(i in 1:nrow(x)){
interval <- (x[i,2]-x[i,1])/49
bothhalfs <- vector(length=49)
bothhalfs[1] <- x[i,1]
for(o in 2:length(bothhalfs)){
bothhalfs[o] <- bothhalfs[o-1]+interval
}
bothhalfs <- c(bothhalfs,x[i,2])
gadtable <- rbind(gadtable,bothhalfs)
}
}
tttt <<- gadtable
require("DescTools")
gadtable <- Rev(gadtable,2)     #reverse order of columns, so that we have 0=pace and 50=trot
colnames(gadtable) <- 1:50
mads <- apply(gadtable,2,madN)
madsmin <- which(mads == min(mads))
mads <<- mads
var <- apply(gadtable,2,varN)
varmin <- which(var == min(var))
VAR <<- var
if(quiet == FALSE){
print("exact GADvar (main approach 1):")
print(c("var",varmin,"mad:",madsmin))
print(c("Distribution:",distribution))
}

if(distribution == "normal"){
return(varmin)
}else{
return(madsmin)
}

#Plot:
# xa <- 1:50
# plot(xa,mads,pch=19,col="black")
# points(xa,sqrt(VAR),pch=19,col="red")
# points(xa,sqrt(IQR),pch=19,col="blue")
}


#Stride plot
strideplot <- function(x){
x <- unlist(stride(x))
#expand vector
xx <- c(x,rep(NA,10))
medstride <- vector()
for(i in 1:length(x)){
if(is.even(i) == FALSE){
ps <- xx[i]
ms1 <- xx[i+3]
ms2 <- xx[i+5]
footprints <- c(ps,ms1,ms2)
if(sum(is.na(footprints)) > 1){
medstride <- c(medstride,NA)
}else{
medstride <- c(medstride,meanN(footprints))
}}
else{
medstride <- c(medstride,NA)
}
}
Stride <<- x
medstride <<- medstride
#expand vector
xxx <- c(NA,NA,x,NA,NA)
Stride <- rep(NA,length(xxx))
for(i in 1:length(x)){
if(is.na(x[i]) == FALSE){
if(is.even(i) == TRUE){  #if manus
Stride[i] <- xxx[i+2]
}else{	#if pes
Stride[i+2] <- xxx[i+2]
}
}}
Stride <- Stride[3:((length(Stride)-2))]
medstride <- c(rep(NA,1),medstride[1:(length(medstride)-1)])
xes <- 1:length(medstride)
xes <<- xes
forplot <- cbind(xes,medstride)
forplot <- forplot[complete.cases(forplot),]
Stride <<- Stride
svg("strideplot.svg")
plot(Stride,col="white",xlab="footprint number (along trackway)",ylab="stride length (cm)")
points(forplot,type="l",lwd="2",col="black")
points(Stride,type="b",pch=19,col="forestgreen")
xes2 <- 1:length(Stride)
xes2 <<- xes2
rightpes <- Stride[seq(1,length(Stride),4)]
rightpesX <- xes[seq(1,length(xes),4)]
rightmanus <- Stride[seq(1,length(Stride),4)+1]
rightmanusX <- xes[seq(1,length(xes),4)+1]
leftpes <- Stride[seq(1,length(Stride),4)+2]
leftpesX <- xes[seq(1,length(xes),4)+2]
leftmanus <- Stride[seq(1,length(Stride),4)+3]
leftmanusX <- xes[seq(1,length(xes),4)+3]
points(rightpesX,rightpes,pch=15,col="lightcoral",cex=1.5)
points(rightmanusX,rightmanus,pch=16,col="lightcoral",cex=1.5)
points(leftpesX,leftpes,pch=15,col="royalblue1",cex=1.5)
points(leftmanusX,leftmanus,pch=16,col="royalblue1",cex=1.5)
dev.off()
return(medstride)
}


#Sensitivity analysis. Add random normally distributed noise to coordinates
sensitivitytrack <- function(x,variant="one"){
stri <- medianN(unlist(stride(x)))
sp <- stri/100
S <- list()

if(variant == "one"){
noises <- c(1,sp,sp*2,sp*3,sp*5,sp*10,sp*15,sp*20,sp*25,sp*30,sp*35,sp*40,sp*45,sp*50)
S[[1]] <- x
for(i in 2:length(noises)){
S[[i]] <- apply(x,c(1,2),function(f){
rnorm(1,mean=f,sd=noises[i])
})
}
}
if(variant == "two"){
spp <- sp*20
noises <- rep(spp,50)
for(i in 1:length(noises)){
S[[i]] <- apply(x,c(1,2),function(f){
rnorm(1,mean=f,sd=noises[i])
})
}}


return(S)
}

#quick function to perform f-test on gads
fga <- function(x,y){
var.test(gads[,x],gads[,y],alternative="less")$p.value
}



















































