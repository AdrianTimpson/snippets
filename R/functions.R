#-----------------------------------------------------------------------------------------------
data.in.polygon <- function(data,kml.path,polygons=NULL,index=NULL){
	# returns the data in a particular polygon, use the full kml file path of the polygon
	# alternatively, if kml.path is NULL, polygons can be used: a list of matrixes of two columns (no names): long,lat
	# index: can be used to select particular polygons only in the kml
	require(maptools)
	require(splancs)
	require(rgdal)
	data.region <- NULL
	if(!is.null(kml.path)&!is.null(polygons))stop('Which one do you want? the kml.path or the polygons? Make one of them NULL')
	if(is.null(polygons))polygons <- getKMLcoordinates(kml.path,ignoreAltitude=T)
	if(is.null(index))index <- 1:length(polygons)
	for(n in index){
		polygon <- polygons[[n]]
		if(!'Longitude'%in%names(data))stop('Longitude is not in the data')		
		if('Longitude'%in%names(data)){
			data <- subset(data,!is.na(Longitude))
			georefs <- pip(pts=data.frame(x=data$Longitude,y=data$Latitude),poly=polygon,bound=T)
			index <- data$Longitude%in%georefs$x & data$Latitude%in%georefs$y
			}

		data.region <- rbind(data.region,data[index,])	
		}	
return(data.region)}
#-----------------------------------------------------------------------------------------------
summary.matrix.maker <- function(start,end,nbin,N=10000,prior=c(0.5,0.5),start.date,end.date,presence,absence,adaptive){
	# take any count data type (LP presence/absence, or milk/non.milk, or domestic/wild)
	# take date ranges for each sample
	# generate summary matrices that account for chronological uncertainty, and small sample sizes

	# a few data consistency checks
	if(length(start.date)!=length(end.date))stop('E1')
	if(length(start.date)!=length(presence))stop('E2')
	if(length(start.date)!=length(absence))stop('E3')
	if(sum(end.date>start.date)!=0)stop('E4')

	# matrices to store various statistics
	sample.mat <- perc.mat <- matrix(, N, nbin)

	if(!adaptive)bins <- seq(start, end, length.out=nbin+1)
	if(adaptive){
		x <- runif(nrow(d)*100,end.date,start.date)
		x <- x[x>end & x<start]
		bins <- rev(round(as.numeric(quantile(x, probs = seq(0,1,length.out=nbin+1)))))
		bins[1] <- start
		bins[length(bins)] <- end
		}

	for(n in 1:N){
		# random date from each phase or sample
		dates <- runif(length(end.date), end.date, start.date)

		# loop for each timeslice
		for(b in 1:nbin){

			# subset (index) of data in each timeslice
			i <- dates<bins[b] & dates>bins[b+1]

			# presence and absence
			P <- presence[i]
			A <- absence[i]

			# if data is not phased (each row is a sample, such as LP)
			if(!is.numeric(P)){
				P <- sum(P)
				A <- sum(A)
				}

			# presence percentage including the prior and small sample size uncertainty
			raw.perc <- rbeta(length(P),P+prior[1],A+prior[2])
	
			# weight proportions in each phase according to sample size
			weight <- 1 - 1/(sqrt(P+A+1))

			# weighted average percentage
			perc <- sum(raw.perc*weight)/sum(weight)

			# store results
			if(sum(i)!=0)perc.mat[n,b] <- perc

			if(!is.numeric(P)){
				# if data is not phased (each row is a sample, such as LP), we also want to include the small sample size uncertainty
				P <- sum(P)
				A <- sum(A)

				perc.mat[n,b] <- rbeta(1,P+prior[1],A+prior[1])
				}

			# if no samples fall within a bin
			if(length(P)==0)P <- 0
			if(length(A)==0)A <- 0
			sample.mat[n,b] <- sum(P)+sum(A)
			}
		}

		NS <- colMeans(sample.mat)

	# HPD and MAP
	require(LaplacesDemon)
	res <- data.frame(matrix(,ncol(perc.mat),7));names(res) <- c('startBP','endBP','MAP','lower95','upper95','lower50','upper50')
	for(n in 1:ncol(perc.mat)){
		if(NS[n]>0){
			x <- perc.mat[,n]
			x <- x[!is.na(x)]
			if(length(x)>10){
				int50 <- p.interval(x,prob=0.50,MM=F, HPD=F)
				int95 <- p.interval(x,prob=0.95,MM=F, HPD=F)
				res$MAP[n] <- mean(p.interval(x,prob=0.1,MM=F, HPD=F))
				res$lower95[n] <- int95[1,1]
				res$upper95[n] <- int95[1,2]
				res$lower50[n] <- int50[1,1]
				res$upper50[n] <- int50[1,2]
				}
			}
		}

	res$startBP <- bins[1:nbin]
	res$endBP <- bins[2:(nbin+1)]

return(res)}
#-----------------------------------------------------------------------------------------
write.csv.utf8.BOM <- function(df, filename){
	 con <- file(filename, "w")
	tryCatch({
	for (i in 1:ncol(df))
 	df[,i] = iconv(df[,i], to = "UTF-8") 
	writeChar(iconv("\ufeff", to = "UTF-8"), con, eos = NULL, nchar=1)
	write.csv(df, file = con, na='NULL', row.names=FALSE)
	},finally = {close(con)})
	}
#-----------------------------------------------------------------------------------------------
getKMLnames <- function(file){	
	names <- readLines(file)
	names <- names[grep('\t\t\t<name>',names)]
	names <- gsub('\t\t\t<name>','',names)
	names <- gsub('</name>','',names)
	names <- gsub('&amp;','&',names,fixed=T)
return(names)}
#-----------------------------------------------------------------------------------------------
generate.colours.for.a.numeric <- function(x){
	if(max(x,na.rm=T)<1.1)res <- generate.colours.for.zero.to.one(x)
	if(max(x,na.rm=T)>5)res <- generate.colours.for.integers(x)
return(res)}
#--------------------------------------------------------------------------------------------------
generate.colours.for.zero.to.one <- function(x){
	posts <- round(quantile(x,na.rm=T),3)
	N <- length(posts)-1
	key <- code <- col <- c()
	for(n in 1:N){
		lower <- posts[n]
		upper <- posts[n+1]
		key[n] <- paste(lower,upper,sep='=>')
		i <- x>=lower & x<upper
		code[i] <- n
		}
	cols <- colorRampPalette(c("red", "blue"))(N)
	for(n in 1:(N))col[code==n] <- cols[n]
	legend <- data.frame(key=key,col=cols)
return(list(col=col,legend=legend))}
#--------------------------------------------------------------------------------------------------
generate.colours.for.integers <- function(x){
	posts <- floor(quantile(x[!x%in%c(0,1,2)],na.rm=T))	
	N <- length(posts)-1
	posts[N+1] <- posts[N+1]+1
	key <- code <- col <- c()
	for(n in 1:N){
		lower <- posts[n]
		upper <- posts[n+1]
		key[n] <- paste(lower,upper,sep='=>')
		i <- x>=lower & x<upper
		code[i] <- n+2
		}
	code[x==1] <- 1
	code[x==2] <- 2
	cols <- colorRampPalette(c("red", "blue"))(N+2)
	for(n in 1:(N+2))col[code==n] <- cols[n]
	legend <- data.frame(key=c(1,2,key),col=cols)
return(list(col=col,legend=legend))}
#--------------------------------------------------------------------------------------------------
summary.maker <- function(d){

	x <- as.data.frame(table(d$SiteID)); names(x) <- c('SiteID','count')
	x <- merge(x,unique(d[,1:3]),by='SiteID')
	x$code[x$count==1] <- 1
	x$code[x$count==2] <- 2
	posts <- floor(unique(quantile(x$count[!x$count%in%c(1,2)])))
	N <- length(posts)-1
	posts[N+1] <- posts[N+1]+1
	key <- c()
	for(n in 1:N){
		lower <- posts[n]
		upper <- posts[n+1]
		key[n] <- paste(lower,upper,sep='=>')
		i <- x$count>=lower & x$count<upper
		x$code[i] <- n+2
		}
	cols <- colorRampPalette(c("red", "blue"))(N+2)
	for(n in 1:(N+2))x$col[x$code==n] <- cols[n]
	legend <- c(1,2,key)
return(list(summary=x,cols=cols,legend=legend))}
#--------------------------------------------------------------------------------------------------
slc <- function(x,y,ax,ay,input='rad'){
	# inputs required in rad or deg
	# calculate shortest distance between a single point (x,y) and all points(xa,ya) using spherical law of cosines
	# returns distances in radians, therefore only needs multiplying by radius of earth 6378.1 to convert to km
	if(length(ax)!=length(ay))stop('ax must be same length as ay')
	if(length(ax)==0)return(0); if(length(ax)>0){
	if(input=='deg'){x <- x*pi/180; y <- y*pi/180; ax <- ax*pi/180; ay <- ay*pi/180;}
	step.1 = sin(y) * sin (ay) + cos(y) * cos(ay) * cos(ax - x)
	
	# floating point bullshit, as sometimes 1 is greater than 1 (Rinferno!) 
	step.1[step.1>1]=1
	step.1[step.1<(-1)]=-1
	dist <- acos(step.1)	
return(dist)}}
#-----------------------------------------------------------------------------------------------
plot.frequency <- function(res,xlab='x',ylab='y',ylim=c(0,100),xaxt='n',...){
	plot(NULL,xlim=c(max(res$startBP),min(res$endBP)),ylim=ylim,xlab=xlab,ylab=ylab,xaxt=xaxt,...)
	for(n in 1:nrow(res)){
		polygon(x=c(res$startBP[n],res$endBP[n],res$endBP[n],res$startBP[n]), y=100 * c(res$lower95[n],res$lower95[n],res$upper95[n],res$upper95[n]), col='lightgrey', border=NA)
		polygon(x=c(res$startBP[n],res$endBP[n],res$endBP[n],res$startBP[n]), y=100 * c(res$lower50[n],res$lower50[n],res$upper50[n],res$upper50[n]), col='grey', border=NA)
		lines(x=c(res$startBP[n],res$endBP[n]),y=100 * c(res$MAP[n],res$MAP[n]),col='black',lwd=1)
		}
	}
#-----------------------------------------------------------------------------------------------
plot.data.density <- function(res,xlab='x',ylab='y'){
	plot(NULL,xlim=c(max(res$bins),min(res$bins)),ylim=c(0,max(res$ns)),xlab=xlab,ylab=ylab)
	for(b in 1:(length(res$bins)-1)){
		polygon(x=c(res$bins[b],res$bins[b+1],res$bins[b+1],res$bins[b]), y=c(0,0,res$ns[b],res$ns[b]), col='lightgrey')
		#text(x=mean(c(res$bins[b],res$bins[b+1])),y=max(res$ns)*0.8,labels=round(res$P[b],1),col='blue',cex=0.5)
		#text(x=mean(c(res$bins[b],res$bins[b+1])),y=max(res$ns)*0.7,labels=round(res$A[b],1),col='red',cex=0.5)
		}
	}
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------