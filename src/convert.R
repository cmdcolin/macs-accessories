
# Convert two digit numbers to roman numerals
convertFileSacCer3<-function(filename) {
  
  filetext<-readLines(filename)
  writeLines(filetext,sprintf("%s.bak",filename))
  filetext <- paste(filetext,collapse="\n")
  
  for(i in 1:16) {
    str<-sprintf("%02d",i)
    filetext<-str_replace_all(filetext, 
      sprintf("chr%s.fs",str), sprintf("chr%s",romanNum(str)))
    printf("Finished chr%02d\n", i)
  }
  filetext<-str_replace_all(filetext,  "chrmt.fs", "chrM")
  printf("Finished chrmt\n", i)
  writeLines(filetext,filename)
  
}


# Convert to digits to roman numerals
convertFileSacCer3_mod<-function(filename) {
  
  filetext<-readLines(filename)
  writeLines(filetext,sprintf("%s.bak",filename))
  filetext <- paste(filetext,collapse="\n")
  
  # Reverse order so that chr11-chr16 gets converted before chr1
  for(i in 16:1) {
    str<-sprintf("%d",i)
    filetext<-str_replace_all(filetext, 
      sprintf("chr%s",str), sprintf("chr%s",romanNum(str)))
    printf("Finished chr%d\n", i)
  }
  
  writeLines(filetext,paste(filename))
}


# Convert to digits to roman numerals
convertFileSacCer3_mod2<-function(filename) {
  
  filetext<-readLines(filename)
  writeLines(filetext,sprintf("%s.bak",filename))
  filetext <- paste(filetext,collapse="\n")
  
  # Reverse order so that chr11-chr16 gets converted before chr1
  for(i in 16:1) {
    str<-sprintf("%02d",i)
    filetext<-str_replace_all(filetext, 
                              sprintf("chr%s",str), sprintf("chr%s",romanNum(str)))
    printf("Finished chr%d\n", i)
  }
  
  writeLines(filetext,paste(filename))
}


convertFileS288C<-function(filename) {
  filetext<-readLines(filename)
  writeLines(filetext,sprintf("%s.bak",filename))
  filetext <- paste(filetext,collapse="\n")
  
  for(i in 1:16) {
    str<-sprintf("%02d",i)
    filetext<-str_replace_all(filetext, 
      sprintf("chr%s.fsa",str), sprintf("Chr%s",str))
    if(debug==TRUE)
      printf("Finished chr%02d\n", i)
  }
  filetext<-str_replace_all(filetext,  "chrmt.fsa", "Chrmt")
  if(debug==TRUE)
    printf("Finished chrmt\n", i)
  
  writeLines(filetext,filename)
  
}


convert<-function() {
  for(file in c(list.files(recursive=TRUE,pattern="*.wig$"),
                list.files(recursive=TRUE,pattern="*.bed$"))) 
  {
    printf("Converting %s\n",file)
    # Uses sacCer3 chromosome
    if(saccer==TRUE)
      convertFileSacCer3(file)
    # Uses S288c chromosome
    else if(s288c==TRUE)
      convertFileS288C(file)
  }
}



romanNum<-function(str) {
  num=as.numeric(str)
  if(num==1) "I"
  else if(num==2) "II"
  else if(num==3) "III"
  else if(num==4) "IV"
  else if(num==5) "V"
  else if(num==6) "VI"
  else if(num==7) "VII"
  else if(num==8) "VIII"
  else if(num==9) "IX"
  else if(num==10) "X"
  else if(num==11) "XI"
  else if(num==12) "XII"
  else if(num==13) "XIII"
  else if(num==14) "XIV"
  else if(num==15) "XV"
  else if(num==16) "XVI"
  else if(str=="mt") "M"
  else "Error"
}


###! Unused

#extended regex
#str_replace_all(text, "(chr)([0-9a-z]{2})(.fsa)", sprintf("\\1%s","\\1",getRoman("\\2")))
convertFileOld<-function(filename) {
  con<-file(filename)
  open(con)
  lines<-readLines(con,-1)
  outlines<-sapply(lines, function(cline) {
    x<-str_match(cline,"(.*chr)([0-9a-z]{2})(.fsa)(.*)")
    if(sum(is.na(x))==0) {
      sprintf("%s%s%s",x[2],romanNum(x[3]),x[5])
    }
    else
      cline
  })
  writeLines(outlines,filename)
  close(con)
}




# Write NormDiff
# temp1<-read.table("data/HS959combined_normdiff.txt",skip=1)
# temp2<-data.frame(Chr=temp1$V2,temp1$V3,
#                   Start=temp1$V3,End=temp1$V3+10,Score=temp1$V5)
# write.table(temp2,"out_normdiff.txt",sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)
