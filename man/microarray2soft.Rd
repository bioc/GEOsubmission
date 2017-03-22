\name{microarray2soft}
\alias{microarray2soft}
\title{Prepares microarray data for submission to GEO}
\description{
Converts a microarray dataset and the corresponding sample information
into a SOFT file to be used for GEO submission.} 
\usage{
microarray2soft(samplenames, sampleinfo, seriesnames, seriesinfo,
                datadir=NULL, infodir=NULL, writedir=NULL, softname=NULL,
                expressionmatrix=NULL, verbose=TRUE)}
\arguments{
  \item{samplenames}{character vector. Sample names.}
  \item{sampleinfo}{character string. Name of a tab-delimited file with
    sample information. Each line in the file corresponds to a
    sample. Sample names are stored under column header SAMPLE. Other
    column headers will be looked up for a match against the sample
    labels required by the SOFT format.} 
  \item{seriesnames}{charcter vector. Series name. Works only with 1
    series at present time.} 
  \item{seriesinfo}{character string. Name of a tab-delimited file with
    series information. Only 1 line (ie. 1 series) allowed at present
    time (in addition to the header). Series names are stored under
    column header SERIES. Other column headers will be looked up for a
    match against the series labels required by the SOFT format.} 
  \item{datadir}{character string. Name of the directory (with
    associated path) containing the microarray data files. Defaults to
    the working directory.} 
  \item{infodir}{character string. (Full path to) Directory containing
    the text file with sample info (sampleinfo) and series info
    (seriesinfo). Defaults to datadir.} 
  \item{writedir}{character string. (Full path to) Directory to write
    the SOFT file to. Defaults to datadir.} 
  \item{softname}{character string. SOFT file name. Defaults to the name
    of the series info file, with suffix .soft. If "", print to the
    standard output connection.} 
  \item{expressionmatrix}{character string. Name of a tab-delimited text
    file with normalized expression values, by default in datadir (but
    if the string contains "/", it is considered a full path and it is
    read or written at the corresponding location). If expressionmatrix
    is non-null and a corresponding file exists, it is used as the
    source of gene expression values to be written in the SOFT
    file. Column names in this file should match the names of raw data
    files given in sampleinfo. If expressionmatrix is non-null but no
    corresponding file is found, normalized expression values calculated
    from the microarray data files (only Affymetrix CEL files at the
    moment) are written to this file.} 
  \item{verbose}{logical. If FALSE suppress all messages (useful to
    write the SOFT file to standard output). Defaults to TRUE.} 
}
\details{
Converts a microarray dataset and the corresponding sample information
into a SOFT file that is ready to be used for direct deposit into GEO. 

Sample information and series information are taken from two separate
plain text files. The column names in these two files should match the
labels defined by the SOFT format (see
http://www.ncbi.nlm.nih.gov/geo/info/soft2.html). This allows for the
extraction of the information needed to write the SOFT format file. 

Currently only a subset of SOFT keywords are accepted but they are
expected to suffice for standard submissions (in particular for 1-color
microarrays). 

Normalized expression values are read from a tab-delimited
file. Alternatively, in the case of an Affymetrix dataset, raw
microarray data files (CEL files) can be processed directly (they are
RMA-normalized) and written to the SOFT file. 

See the example below or the vignette for the format of sampleinfo,
seriesinfo and expressionmatrix. 
}
\author{Alexandre Kuhn, alexandre.m.kuhn@gmail.com}
\note{The SOFT file format is defined at
  http://www.ncbi.nlm.nih.gov/geo/info/soft2.html}

\examples{
##location of example files for sampleinfo and seriesinfo (i.e. provided
##with this vignette)

dataDirectory <- system.file(package='GEOsubmission', 'extdata')

## Define a directory and a file to write the generated example SOFT
## file out to (for the sake of this demo) 
soft_example_fullpath <- tempfile(pattern='soft_example')
soft_example_name <- basename(soft_example_fullpath)
soft_example_dir <- dirname(soft_example_fullpath)

## Example with a dataset (series) named "neuronalCultures", comprised
## of 2 samples (named "1" and "2")

## This example is not run since the CEL files corresponding to sample
## "1" and "2" (and indicated in "sampleInfo.txt") do not exist (which
## would be purposedly detected by microarray2soft and cause it to
## stop). 
\dontrun{
microarray2soft(c('1','2'), 'sampleInfo.txt', 'neuronalCultures',
                'seriesInfo.txt', datadir=dataDirectory,
                writedir=soft_example_dir, softname=soft_example_name)

}

## The files "sampleInfo.txt" (as sampleInfo) and "seriesInfo.txt" (as
## seriesInfo) are provided as example files in the "extdata" directory of
## the package installation directory. Their content can be printed to the
## R console with 
\dontrun{
dataDirectory<-system.file(package='GEOsubmission','extdata')
read.delim(file.path(dataDirectory,'sampleInfo.txt'))
read.delim(file.path(dataDirectory,'seriesInfo.txt'))
}

## Example specifying a separate file containing normalized expression
## values (here the example file "expressionNormalized.txt") 
microarray2soft(c('1','2'), 'sampleInfo.txt', 'neuronalCultures',
                'seriesInfo.txt', datadir=dataDirectory,
                writedir=soft_example_dir, softname=soft_example_name, 
                expressionmatrix='expressionNormalized.txt')


## The file "expressionNormalized.txt" (as expressionmatrix) is provided
## as example file in the "extdata" directory of the package installation
## directory. It can be printed to the R console with 
\dontrun{
dataDirectory<-system.file(package='GEOsubmission','extdata')
read.delim(file.path(dataDirectory,'expressionNormalized.txt')
}

## Same as the preceding example but write the example SOFT file to the
##console (note the use of verbose=FALSE to suppress messages) 
microarray2soft(c('1','2'), 'sampleInfo.txt', 'neuronalCultures',
                'seriesInfo.txt', datadir=dataDirectory,
                softname='', expressionmatrix='expressionNormalized.txt',  
                verbose=FALSE)


\dontshow{
## Clean-up
unlink(soft_example_fullpath)
}
}
