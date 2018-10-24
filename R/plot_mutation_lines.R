#' Function to plot the lineplots with all clinical and genomic information
#' @param blasts_to_plot output from `import_blast_for_lineplots()`
#' @param indels_to_plot output from `import_indels_for_lineplots()`
#' @param snvs_to_plot output from `import_snvs_for_lineplots()`
#' @param cnas_to_plot output from `import_cnas_for_lineplots()`
#' @param clones_to_plot output from `import_clones_for_lineplots()`
#' @param spread spread


plot_mutation_lines = function(blasts_to_plot, indels_to_plot, snvs_to_plot, cnas_to_plot, clones_to_plot, spread=0.1) {

  samples = colnames(blasts_to_plot$y_matrix)
  time = gsub('\\..+$', '', gsub('^D[0-9]*\\.', '', samples))
  samples = samples[order(time!='Screen', time)]
  if ( class(indels_to_plot) != 'NULL' ) samples = intersect(samples, colnames(indels_to_plot$y_matrix))
  if ( class(snvs_to_plot) != 'NULL' ) samples = intersect(samples, colnames(snvs_to_plot$y_matrix))
  if ( class(cnas_to_plot) != 'NULL' ) samples = intersect(samples, colnames(cnas_to_plot$y_matrix))
  #order the samples here, the order of samples is used in the plot.
  N = length(samples)
  totalYmx = NULL
  labels = data.frame(y=numeric(0), label=character(0), col=character(0))
  plot(0, type='n', xlim=c(0.5, N+2), ylim=c(-0.15,1.1), frame=F, xlab='', xaxt='n', ylab='clonality', main=patient)
  if ( class(indels_to_plot) != 'NULL' ) {
    indelCol = superFreq::mcri('red')
    yMx = indels_to_plot$y_matrix[,samples,drop=F]
    for ( row in 1:nrow(yMx) ) lines(1:N+spread*(0.5-row/nrow(yMx)), yMx[row,], lwd=2, col=indelCol)
    if ( ncol(yMx) == 1 ) points(1+spread*(0.5-row(yMx)/nrow(yMx)), yMx, cex=1.5, pch=16, col=indelCol)
    indels_labels = data.frame(y=yMx[,N], label=indels_to_plot$mutations, col=indelCol)
    labels = rbind(labels, indels_labels)
    totalYmx = yMx
  }
  if ( class(clones_to_plot) != 'NULL' ) {
    clonesCol = superFreq::mcri('grey')
    yMx = clones_to_plot$y_matrix[,samples,drop=F]
    for ( row in 1:nrow(yMx) ) lines(1:N+spread*(0.5-row/nrow(yMx)), yMx[row,], lwd=2, col=clonesCol)
    if ( ncol(yMx) == 1 ) points(1+spread*(0.5-row(yMx)/nrow(yMx)), yMx, cex=1.5, pch=16, col=clonesCol)
    clones_labels = data.frame(y=yMx[,N], label=clones_to_plot$mutations, col=clonesCol)
    labels = rbind(labels, clones_labels)
    totalYmx = rbind(totalYmx, yMx)
  }
  if ( class(blasts_to_plot) != 'NULL' ) {
    blastCol = 'black'
    yMx = blasts_to_plot$y_matrix[,samples,drop=F]
    for ( row in 1:nrow(yMx) ) lines(1:N+spread*(0.5-row/nrow(yMx)), yMx[row,], lwd=2, col=blastCol)
    if ( ncol(yMx) == 1 ) points(1+spread*(0.5-row(yMx)/nrow(yMx)), yMx, cex=1.5, pch=16, col=blastCol)
    blasts_labels = data.frame(y=yMx[,N], label=blasts_to_plot$mutations, col=blastCol)
    labels = rbind(labels, blasts_labels)
    totalYmx = rbind(totalYmx, yMx)
  }
  if ( class(snvs_to_plot) != 'NULL' ) {
    snvCol = superFreq::mcri('blue')
    yMx = snvs_to_plot$y_matrix[,samples,drop=F]
    for ( row in 1:nrow(yMx) ) lines(1:N+spread*(0.5-row/nrow(yMx)), yMx[row,], lwd=2, col=snvCol)
    if ( ncol(yMx) == 1 ) points(1+spread*(0.5-row(yMx)/nrow(yMx)), yMx, cex=1.5, pch=16, col=snvCol)
    snvs_labels = data.frame(y=yMx[,N], label=snvs_to_plot$mutations, col=snvCol)
    labels = rbind(labels, snvs_labels)
    totalYmx = rbind(totalYmx, yMx)
  }
  if ( class(cnas_to_plot) != 'NULL' ) {
    cnaCol = superFreq::mcri('green')
    yMx = cnas_to_plot$y_matrix[,samples,drop=F]
    for ( row in 1:nrow(yMx) ) lines(1:N+spread*(0.5-row/nrow(yMx)), yMx[row,], lwd=2, col=cnaCol)
    if ( ncol(yMx) == 1 ) points(1+spread*(0.5-row(yMx)/nrow(yMx)), yMx, cex=1.5, pch=16, col=cnaCol)
    cnas_labels = data.frame(y=yMx[,N], label=cnas_to_plot$mutations, col=cnaCol)
    labels = rbind(labels, cnas_labels)
    totalYmx = rbind(totalYmx, yMx)
  }
  labels$y = superFreq:::spreadPositions(labels$y, 0.05)
  text(N+1, labels$y, labels$label, col=as.character(labels$col))
  text(1:N, -0.1, samples, srt=15)

  rownames(totalYmx) = labels$label
  colnames(totalYmx) = samples
  invisible(as.data.frame(totalYmx))
}

