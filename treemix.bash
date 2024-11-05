# TREEMIX进行基因流的检测
# TREEMIX官方manual在https://bitbucket.org/nygcresearch/treemix/downloads/

#先安装两个依赖软件boost（v>1.42）和gsl,如果解压不成功就先下到本地再上传到服务器
#gsl
wget -c http://ftp.club.cc.cmu.edu/pub/gnu/gsl/gsl-latest.tar.gz
tar -zxf gsl-latest.tar.gz
cd gsl-2.7.1/
./configure --prefix=/data01/wangyf/software/gsl
make
make install
#添加环境变量
export PATH=$PATH:/data01/wangyf/software/gsl/bin
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/data01/wangyf/software/gsl/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data01/wangyf/software/gsl/lib
export GSL_LD=/data01/wangyf/software/gsl/lib

#boost
wget -c https://www.boost.org/users/download/#repository
#下载最新的：
#boost_1_74_01.tar.gz
tar -xzvf boost_1_74_01.tar.gz
#step 1
cd boost_1_74
#step 2
bootstrap.sh –with-libraries=all –with-toolset=gcc 
#step3
./b2 install --prefix=/data01/wangyf/software/boost_1_83_0


#安装treemix
wget https://bitbucket.org/nygcresearch/treemix/downloads/treemix-1.13.tar.gz
tar -xvf treemix-1.13.tar.gz
cd treemix-1.13
./configure --prefix=/data01/wangyf/software/treemix-1.13
# checking for cblas_dgemm in -lgslcblas... no
# configure: error: could not find GSL BLAS
./configure LDFLAGS=-L/data01/wangyf/software/gsl/lib CPPFLAGS=-I/data01/wangyf/software/gsl/include
make
make install
#Treemix是假设你的SNP是不连锁的，并且它并不喜欢SNP VCF中有缺失的数据
#所以需要先对ld和maxmissing 1进行过滤

##########################################################################################################################################

# ld pruning
./ldPruning.sh RefTdalaica.Tdalaica_17.snpgap.xindels.biallele.depth.nomiss.minq.maf.vcf.gz 0.5
# https://speciationgenomics.github.io/Treemix/


# run treemix 
#m为migration edge值，i为重复运行次数
# -i 指定基因频率输入文件
# -o 指定输出文件前缀
# -root指定外类群(指定的是居群名称)，多个用逗号分隔；最好指定，否则后面plot_tree画树没找到更换外类群的参数会很麻烦
# -m 为the number of migration edges即基因渗入的次数
# -k 500 因为SNP之间不是独立位点，为了避免连锁不平衡，用k参数指定SNP数量有连锁，比如这里指定用1000个SNP组成的blocks评估协方差矩阵
# -noss 关闭样本量校正。TreeMix计算协方差会考虑每个种群的样本量，有些情况(如果有种群的样本只有1个)会过度校正，可以关闭。


for m in {1..20}
   do
   for i in {1..10}
      do
      treemix \
         -i samples.max1.LDpruned.treemix.frq.gz \
         -o test.${i}.${m} \
         -m ${m} \
         -bootstrap -k 1000 -noss \
		 -root DN
      done 
done

# 使用optM推测最佳的treemix migration edge值
optM在服务器R4.3.2版本下安装成功了
https://cran.r-project.org/web/packages/OptM/readme/README.html
https://github.com/carolindahms/TreeMix/blob/main/Step2%264_TreeMix.R

# run using various linear modeling estimates rather than the ad hoc statistic
library(OptM)
folder <- file.path(path="E:\\Rworkspace\\treemix\\usa-every-sample")  #path to files of TreeMix replicates with different migration edges (m) to test
test.linear = optM(folder, method = "linear", tsv="linear.txt")   #test m: produces a list of which the $out dataframe suggests optimum m based on multiple linear models
plot_optM(test.linear, method = "linear")                         #shows changes in log likelihood for different m, and suggests optimum m as 'change points'




# plotting
setwd("E:/Rworkspace/treemix/1.mirrorcarp-18samples/")

library(RColorBrewer)
library(R.utils)
source("E:/Rworkspace/plotting_funcs.R")



pdf("treemix1.m1.pdf", width = 5.5, height = 5)
plot_tree("usa-sv-every.1.1")
dev.off()
