#Convert files

s288c=TRUE
saccer=FALSE

#pwd<-getwd
#cd<-setwd



# Uses roman numerals
convertFileSacCer3<-function(filename) {
  
  filetext<-readLines(filename)
  filetext <- paste(filetext,collapse="\n")
  
  for(i in 1:16) {
    str<-sprintf("%02d",i)
    filetext<-str_replace_all(filetext, 
                              sprintf("chr%s.fsa",str), sprintf("chr%s",romanNum(str)))
    printf("Finished chr%02d\n", i)
  }
  filetext<-str_replace_all(filetext,  "chrmt.fsa", "chrM")
  printf("Finished chrmt\n", i)
  writeLines(filetext,filename)
  
}

convertFileS288C<-function(file) {
  
  filetext<-readLines(file)
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
  
  writeLines(filetext,file)
  
}


convert<-function() {}
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
  if(str=="01") "I"
  else if(str=="02") "II"
  else if(str=="03") "III"
  else if(str=="04") "IV"
  else if(str=="05") "V"
  else if(str=="06") "VI"
  else if(str=="07") "VII"
  else if(str=="08") "VIII"
  else if(str=="09") "IX"
  else if(str=="10") "X"
  else if(str=="11") "XI"
  else if(str=="12") "XII"
  else if(str=="13") "XIII"
  else if(str=="14") "XIV"
  else if(str=="15") "XV"
  else if(str=="16") "XVI"
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