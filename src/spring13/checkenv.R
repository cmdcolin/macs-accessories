# use CodeTools to analyze environment for errors

require(codetools)
checkUsageEnv(globalenv(),all=TRUE)