microarray2soft<-function(samplenames,sampleinfo,seriesnames,seriesinfo,datadir=NULL,infodir=NULL,writedir=NULL,softname=NULL,expressionmatrix=NULL,verbose=TRUE) {

	##parameters

	nbsamples<-length(samplenames)
	nbseries<-length(seriesnames)
	

	##defaults

	if (is.null(datadir)) datadir<-getwd()
	if (is.null(infodir)) infodir<-datadir
	if (is.null(writedir)) writedir<-datadir
	if (is.null(softname)) softname<-paste(paste(unlist(sapply(seriesinfo,function(x){strsplit(x,'\\.')[[1]][1]})),collapse='_'),'.soft',sep='')

	##full filenames

	fullsampleinfo<-file.path(infodir,sampleinfo)
	fullseriesinfo<-file.path(infodir,seriesinfo)
	fullsoftname<-file.path(writedir,softname)
	if (softname=="") fullsoftname<-"" #write to standard output

	if (!is.null(expressionmatrix)) {
		if (dirname(expressionmatrix)!='.') fullexpressionmatrix<-expressionmatrix
		else fullexpressionmatrix<-file.path(datadir,expressionmatrix)
	}


	##check if output file already exists

	answer<-''
	if (file.exists(fullsoftname)) {
		while (!(answer %in% c('y','n'))) {
		answer<-readline(paste('File ',fullsoftname,' already exists.\nDo you want to overwrite it? (y/n) ',sep=''))
		if (answer=='n') stop('Terminated. No output written.')
		if (answer=='y') file.remove(fullsoftname)
		}
	}
	

	##read-in sampleinfo and perform consistency checks

	if (verbose) message('Loading sampleinfo...\n')
	sampletable<-read.delim(fullsampleinfo,colClasses='character') # quote='' would keep double-quotes

	if (sum(samplenames %in% sampletable$SAMPLE)!=nbsamples) stop('One or more sample names absent from sampleinfo')

	celnames<-sampletable$Sample_supplementary_file[match(samplenames,sampletable$SAMPLE)]
	fullcelnames<-file.path(datadir,celnames)
	if (sum(file.exists(fullcelnames))!=nbsamples) stop('Non-existing raw data files:\n',paste(fullcelnames[!file.exists(fullcelnames)],collapse='\n'))	


	##read-in seriesinfo and perform consistency checks (ONLY 1 SERIES IMPLEMENTED SO FAR)

	if (verbose) message('Loading seriesinfo...\n')
	seriestable<-read.delim(fullseriesinfo,colClasses='character')

	if (sum(seriesnames %in% seriestable$SERIES)!=nbseries) stop('One or more series names absent from seriesinfo')

	if (dim(seriestable)[1]!=1) stop('More than one series defined in ',fullseriesinfo)
	if (seriesnames!=seriestable$SERIES) stop('Mismatch between series names and Series_names in seriesinfo')
	seriessampleid<-strsplit(seriestable$Series_sample_id,';')[[1]]
	if (!setequal(seriessampleid,samplenames)) stop('Mismatch between sample names and Series_sample_id in seriesinfo')


	##read-in CEL files and normalize

	readexpressionmatrix<-NA
	writeexpressionmatrix<-NA

	if (is.null(expressionmatrix)) {readexpressionmatrix<-0; writeexpressionmatrix<-0}
	else if (!file.exists(fullexpressionmatrix)) {readexpressionmatrix<-0; writeexpressionmatrix<-1}
	else {readexpressionmatrix<-1; writeexpressionmatrix<-0}

	#if (is.null(expressionmatrix) | (if (!is.null(expressionmatrix)) {!file.exists(fullexpressionmatrix)})) {
	
	if (!readexpressionmatrix) {
		if (verbose) message('Loading CEL files and calculating normalized expression values...\n')
		ab<-ReadAffy(filenames=fullcelnames)
		dset<-rma(ab)
		drma<-exprs(dset)
		#if (!is.null(expressionmatrix) & !file.exists(fullexpressionmatrix)) write.table(drma,file=fullexpressionmatrix,quote=FALSE,sep='\t',row.names=TRUE,col.names=TRUE)
		if (writeexpressionmatrix) {
			if (verbose) message(paste('Writing expression matrix to ',fullexpressionmatrix,'...\n',sep=''))
			write.table(drma,file=fullexpressionmatrix,quote=FALSE,sep='\t',row.names=TRUE,col.names=TRUE)
		}
	}
	else {
		if (verbose) message(paste('Loading expression matrix in ',fullexpressionmatrix,'...\n',sep=''))
		drma<-read.delim(fullexpressionmatrix,header=TRUE,quote='')
		drmacolnames<-sapply(celnames,function(x){if (!is.na(charmatch(substr(x,1,1),0:9))) paste('X',x,sep='') else x}) #NB: CEL file name '12..CEL' and 'X12..CEL' are thus considered to match the same sample!
		if (!setequal(drmacolnames,colnames(drma))) stop('Column names in expression matrix do not match the names of raw microarray data files given in sampleinfo file.\n') 
		if (sum(drmacolnames==colnames(drma))!=nbsamples) {
			drma<-drma[match(drmacolnames,colnames(drma)),]
			if (verbose) message('Column order in expression matrix was different than sample order specified in input. Expression matrix was reordered.\n')
		}	
	}

	##define utility functions

	write.soft.append<-function(...) write(...,file=fullsoftname,append=TRUE)
	paste.soft<-function(...) paste(...,sep=' = ')


	##write sample info and expression values to soft file

	sep=' = '
	label_sample<-c('Sample_title','Sample_supplementary_file','Sample_source_name','Sample_organism','Sample_characteristics','Sample_molecule','Sample_extract_protocol','Sample_label','Sample_label_protocol','Sample_hyb_protocol','Sample_scan_protocol','Sample_data_processing','Sample_description','Sample_platform_id')

	samplesind<-match(samplenames,sampletable$SAMPLE)
	
	if (verbose) message(paste('Starting to write to ',fullsoftname,'...\n',sep=''))
	for (i in 1:nbsamples) {

		if (verbose) message(paste('Writing sampleinfo for ',samplenames[i],' (',sampletable$Sample_supplementary_file[samplesind][i],')...\n',sep=''))

		write.soft.append(paste.soft('^SAMPLE',sampletable$SAMPLE[samplesind][i]))
		
		labelvalue_sample<-paste('!',sapply(label_sample,function(x){paste.soft(x,sampletable[[x]][samplesind][i])}),sep='') #sub('\"$','',sub('^\"','', )) is an unelegant way of removing possible starting and ending double-quote
		write.soft.append(labelvalue_sample)
		
		write.soft.append('#ID_REF =')
		write.soft.append(paste.soft('#VALUE','RMA-calculated signal intensity')) 
		
		if (verbose) message(paste('Writing expression data for ',colnames(drma)[i],'...\n',sep=''))
		write.soft.append('!Sample_table_begin')

		write.soft.append(paste('ID_REF','VALUE',sep='\t'))
		write.table(as.matrix(drma)[,i],file=fullsoftname,append=TRUE,quote=FALSE,row.names=TRUE,col.names=FALSE,sep='\t')  ##as.matrix necessary so that row names are properly written for data.frame

		write.soft.append('!Sample_table_end')

	}
	

	##write series to soft file (SO FAR IMPLEMENTED FOR 1 SERIES ONLY)

	label_series<-c('Series_title','Series_summary','Series_type','Series_overall_design')

	write.soft.append(paste.soft('^SERIES',seriestable$SERIES))

	labelvalue_series<-paste('!',sapply(label_series,function(x){paste.soft(x,seriestable[[x]])}),sep='')
	write.soft.append(labelvalue_series)

	##write.soft.append(paste.soft('!Series_summary',seriestable$Series_summary))

	seriescontributor<-strsplit(seriestable$Series_contributor,';')[[1]]
	write.soft.append(paste.soft('!Series_contributor',seriescontributor))

	write.soft.append(paste.soft('!Series_sample_id',seriessampleid))




}