 #ins tall.packages("seqinr")
#library(seqinr)

#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#library(Biostrings)

sequence1='GAATTGCTGAC'
sequence2='GTATTGGTGGGAC'
myScoringMat=nucleotideSubstitutionMatrix(match=1,mismatch=-1,baseOnly=TRUE)
gapOpen=2
gapExtend=0
gapExtend=1
myAlignment <- pairwiseAlignment(sequence1, sequence2,    substitutionMatrix = myScoringMat, gapOpening = gapOpen, gapExtension = gapExtend, type="global", scoreOnly = FALSE) 
score=myAlignment@score

seq2=getSequence(sequence2)#向量

len2=length(seq2)
acount=0
tcount=0
ccount=0
gcount=0
for (x in seq2){
if(x=='A'){acount=acount+1}
if(x=='T'){tcount=tcount+1}
if(x=='C'){ccount=ccount+1}
if(x=='G'){gcount=gcount+1}
}#统计seq2中ATCG出现的数目


times=10000#实验次数
score_list=rep(0,times)#存储每次实验后得到的比对得分
count=0#统计score_list中超过score的得分的数目

for (k in 1:times){
random_num=floor(runif(len2)*len2)#长为len2的向量，存储0~len2-1的整数
new_seq2=rep('',len2)#存储shuffle后的序列（向量）

for (i in 1:length(seq2) ){
if (random_num[i]<acount){#0~acount-1
new_seq2[i]='A'
}else {
if (random_num[i]<acount+tcount){#acount~acount+tcount
new_seq2[i]='T'
}else{
if (random_num[i]<acount+tcount+ccount){
new_seq2[i]='C'
}else{
if (random_num[i]<acount+tcount+ccount+gcount){
new_seq2[i]='G'
}}}}}
new_sequence2=paste(new_seq2,collapse='')#shuffle后的序列（字符串）

myAlignment <- pairwiseAlignment(sequence1, new_sequence2,    substitutionMatrix = myScoringMat, gapOpening = gapOpen, gapExtension = gapExtend, type="global", scoreOnly = FALSE)
score_list[k]=myAlignment@score
if (score_list[k]>=score){
count=count+1
}
}


r=count/times

r



