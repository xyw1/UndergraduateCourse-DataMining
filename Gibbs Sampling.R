library(seqinr)


#######
#######
#自定义函数，将A T C G 分别映射到 1 2 3 4作为数组索引
base_index <- function(base)
{
if (base == 'A') {x=1}
else if (base == 'T') {x=2}
else if (base == 'C') {x=3}
else if (base == 'G') {x=4}
return (x)
}
########
########

seq <- function(n){
  seq1='TTGCCACAAAATAATCCGCCTTCGCAAATTGACCTACCTCAATAGCGGTAGAAAAACGCACCACTGCCTGACAG'
  seq2='GTAAGTACCTGAAAGTTACGGTCTGCGAACGCTATTCCACTGCTCCTTTATAGGTACAACAGTATAGTCTGA'
  seq3='CCACACGGCAAATAAGGAGTAACTCTTTCCGGGTATGGGTATACTTCAGCCAATAGCCGAGAATACTGCCATT'
  seq4='CCATACCCGGAAAGAGTTACTCCTTATTTGCCGTGTGGTTAGTCGCTTTACATCGGTAAGGGTAGGGATTTT'
  seq5='AAACTATTAAGATTTTTATGCAGATGGGTATTAAGGAGTATTCCCCATGGGTAACATATTAATGGCTCTTA'
  seq6='TTACAGTCTGTTATGTGGTGGCTGTTAATTATCCTAAAGGGGTATCTTAGGAATTTACTT'
  seq1=getSequence(seq1)
  seq2=getSequence(seq2)
  seq3=getSequence(seq3)
  seq4=getSequence(seq4)
  seq5=getSequence(seq5)
  seq6=getSequence(seq6)
  if (n==1){seq= seq1}
  else if (n==2){seq= seq2}
  else if (n==3){seq= seq3}
  else if (n==4){seq= seq4}
  else if (n==5){seq= seq5}
  else if (n==6){seq= seq6}
  return (seq)
}

W=16




times=100
F=rep(0,times)
N=6



c <- matrix (0.01,nrow=4,ncol=W)#计数器初始化为0.01以防最后的0影响Q[x]和P[x]的计算结果
q <- matrix (0.01,nrow=4,ncol=W)
bc =rep(0.01,4)
p =rep(0,4)
a=rep(0,N)
for (z in 1:N)
  #每次循环选择一个作为z
{
  #########1.predictive update step##########
  
  for (k in 1:N)
  {
    if (k!=z)
      #一条序列为z，其余为k
    {	
      
      lenk=length(seq(k))
      a[k]=ceiling( runif(1)*(lenk-W+1 ))#在seq[k]的1~lenk-W+1位置处选择一个作为a[k]
      
      for ( pos in 1:lenk )
      {
        # 遍历seq[k],pos在a[k]后W长度范围内(即pattern instance范围内)时更新计数器c(4行W列)
        # pos在pattern instance范围之外时更新计数器bc(四列，对应ATCG)
        baseIndex <- base_index(seq(k)[pos])
        pos_in_W=pos-a[k]+1
        if (pos>=a[k] && pos<=(a[k]+W-1)){c[baseIndex,pos_in_W] <- c[baseIndex,pos_in_W] + 1} #pattern中每个位置碱基出现次数
        else {bc[baseIndex] <- bc[baseIndex] + 1 } #pattern之外每个碱基出现的次数
      }
    }
  }
  #for (i in 1:4){cat(c[i,])
  #cat('\n')}
  #cat(c)
  #cat('\n')
  #lenc=length(c)
  #cat(lenc)
  #cat('\n')
  #cat(bc)






for (time in 1:times)
{


		##########2.Sampling step##############
	
		#计算N条seq的长度
		lens <- rep(0,N)
		for (i in 1:N){lens[i]=length(seq(i))}

		M=sum(lens)-N*W
		cat(M)
		#计算Ax
		lenz=length(seq(z))
		#cat(z)
		#cat(seq(z))
		#cat(lenz)#74
		#cat('\n')
		x_range=lenz-W+1
		#cat(x_range)#59
		Q=rep(0 , x_range)
		P=rep(0 , x_range)
		A=rep(0 , x_range)
		
		for (x in 1:x_range)
		{
			pattern=seq(z)[x:(x+W-1)]
			tmp_p=rep(0,W)
			
			for (tmp_pos in 1:W){tmp_p[tmp_pos]=c[base_index(pattern[tmp_pos]),tmp_pos]/(N-1)}
			lenp=length(tmp_p)
			#cat(lenp)
			#cat('\n')
			#cat(tmp_p)
			#cat('\n')
			Q[x]=prod(tmp_p)
			for (tmp_pos in 1:W){tmp_p[tmp_pos]=bc[base_index(pattern[tmp_pos])]/M}
			P[x]=prod(tmp_p)
			#cat(Q[x])
			#cat('\n')
			#cat(P[x])
			#cat('\n')
			A[x]=Q[x]/P[x]
			#cat(A[x])
			#cat('\n')
		}
		#归一化Ax
		sum_A=sum(A)
		#lenA=length(A)
		#cat(lenA)
		#cat('\n')
		#cat(sum_A)
		#cat('\n')
		for (tmp_pos in 1:x_range){A[tmp_pos]=A[tmp_pos]/sum_A}
		#cat(A)
		#cat('\n')
		#以Ax为权重随机生成a[z]
		tmp_A=rep(0 , x_range+1)
		tmp_A[1]=0
		for (tmp_pos in 1:x_range){tmp_A[tmp_pos+1]=sum(A[1:tmp_pos])}
		#cat(tmp_A)
		rand_num=runif(1)
		#cat(rand_num)
		#cat('\n')
		for (tmp_pos in 1:(x_range-1)){
		  #cat(tmp_pos)
		  if( (rand_num >tmp_A[tmp_pos]) & (rand_num<tmp_A[tmp_pos+1])) {a[z]=tmp_pos}  }
	

	#根据c和bc计算q和p和f
	f <- matrix (0,nrow=4,ncol=W)
	for (j in 1:W)
	{
		#p[j]=bc[j]/M
		#cat('pj')
		#cat(p[j])
		#cat('\n')
		for (i in 1:4)
		{
		  #cat('i:')
			#q[i,j]=c[i,j]/(N-1)
			#cat('cij' )
			#cat('\n')
			#cat(c[i,j])
			#cat('\n')
			#cat('qij:')
			#cat(j)
			#cat('\n')
			f[i,j]=c[i,j]*( log10((c[i,j]*M)/(bc[i]*(N-1))) )
			cat(f[i,j])
			cat('\n')
			cat((c[i,j]*M)/(bc[i]*(N-1)))
			cat('\n')
			#cat('log')
			#cat(q[i,j]/p[j])
			#cat('qij')
			#cat(q[i][j])
			#cat('\n')
			#cat('pj')
			#cat(p[j])
			#cat('\n')
			#cat('cij')
			#cat(c[i,j])
			#cat('\n')
			#cat('fij')
		  #cat(f[i,j])
		  #cat('\n')
		}
	}
#for (i in 1:4){cat(q[i,])
  #cat('\n')
  #cat(c[i,])
  #cat('\n')}
#cat (p)
	for (j in 1:W)
	{
		for (i in 1:4){F[time]=F[time]+f[i,j]}
	}

}}
cat(c)
cat('\n')
cat(bc)
cat('\n')
	cat(F)
	cat('\n')
	cat(A)
	cat(a)
