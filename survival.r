#!/usr/bin/env  Rscript
library(docopt)
"Usage: survival.r  -i <file> -o <dir> [ --model <string> --Rlib <dir> --width <int> --survivalroc_time <numeric> --survivalroc_method <string> ]
Options:
   -i, --input <file>                    输入文件，第一列样本名，第二列生存时间，第三列是否截尾（0截尾,1去世），第四列及之后的数据为要分析的表型
                                         注意：列名尽可能不要出现‘-’减号字符，会报错
   -o, --output_dir <dir>                结果输出目录,需要提前建立
   --model <string>                      分析方式, KM/log-rank/COX/My.stepwise.coxph  [default: KM]
   
                                         KM : Kaplan-Meier plots to visualize survival curves（根据生存时间分布，估计生存率以及中位生存时间，以生存曲线方式展示，从而分析生存特征，一般用Kaplan-Meier法，还有寿命法）
                                              该分析只需要数据前3列

                                         log-rank : Log-rank test to compare the survival curves of two or more groups（通过比较两组或者多组之间的的生存曲线，一般是生存率及其标准误，从而研究之间的差异，一般用log rank检验）
                                                    必须是分类表型，可以是数字、字符。该分析对第四列及之后的表型分别进行分析
                                                    注意：每一组中的样本数量不宜过少，否则会报错

                                         COX  : Cox proportional hazards regression to describe the effect of variables on survival（用Cox风险比例模型来分析变量对生存的影响，可以两个及两个以上的因素，很常用）
                                                该分析对第四列及之后的表型汇总建模

                                         My.stepwise.coxph : 选择一组与生存相关性较强的一组变量，构成COX模型。注：最终只会在桌面打印一组变量名称，能够较好的构成COX模型。

                                         survivalROC : 生存资料的ROC曲线分析。预测第四列及之后的数据与样本的生存状态的ROC曲线。第四列及之后的数据应该是连续性变量。

    --survivalroc_time <numeric>         ntime生存期的ROC曲线， 如果是survivalROC分析的话，必须填写该内容。
    --survivalroc_method <string>        survivalROC分析方法 NNE/KM [default: NNE]
                                         
   --width <int>                         pdf图像宽度 [default: 7]
   --Rlib <dir>                          R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts               <- docopt(doc, version = '生存分析软件 \n')
input              <- opts$input
output_dir         <- opts$output_dir
model              <- opts$model
survivalroc_method <- opts$survivalroc_method
survivalroc_time   <- as.numeric(opts$survivalroc_time)
width              <- as.integer(opts$width)
Rlib               <- opts$Rlib
.libPaths(Rlib)


# input = 'data2.txt'
# output_dir = "./"
# model      = 'COX'

# cat(input,'\n')
# cat(output_dir,'\n')
# cat(model,'\n')

library("survival")
library("survminer")
library("My.stepwise")
library("survivalROC")
if(model != 'KM' & model != 'log-rank' & model != 'COX')
{
    message("无法识别'model'字符，请确认输入是否正确")
    q()
}


# 读入数据
input_data = read.table(input, head = T, row.names = 1, check.names = F, sep = '\t')
# 对前两列换名字
colnames(input_data)[1:2] = c('time', 'status') 

##################
# （1）所有样本放在一起，检查生存率
##################
if(model == 'KM')
{   
    cat('Model : KM \n')
    # 取出分析的信息，并去掉缺失，以及空字符
    input_data <- input_data[, c('time', 'status')]
    input_data <- input_data[complete.cases(input_data), ]  
    input_data <- input_data[apply(input_data, 1, function(x){sum(x=='')}) == 0, ]
 
    # 拟合
    fit <- survfit(Surv(time, status)~1, data = input_data)
    # 拟合详细结果输出
    surv_info <- surv_summary(fit)
    write.table(surv_info, paste0(output_dir, '/survival.KM.surv_info.txt'), row.names = F, quote = F, sep = '\t')

    # 生存曲线绘图
    pdf(paste0(output_dir, '/survival.KM.surv.pdf'), width = width)
    P <- ggsurvplot(fit, data =input_data,
       conf.int = TRUE,
       risk.table = TRUE, # Add risk table
       linetype = "strata", # Change line type by groups
       surv.median.line = "hv", # Specify median survival
       ggtheme = theme_bw(), # Change ggplot2 theme
       palette = c("#2E9FDF")
       )
    print(P, newpage= FALSE)
    dev.off()
}

##################
# （2）计算分组间差异（对每一个表型都计算一遍）
##################
if(model == 'log-rank')
{   
    cat('Model : log-rank \n')
    time <- input_data$time
    status <- input_data$status

    # 对每一个分层表型进行计算
    pdf(paste0(output_dir, '/survival.LogRank.surv.pdf'), width = width)
    pvalue_data = matrix(0,ncol(input_data) - 2, 2)
    colnames(pvalue_data) = c('Pheno', 'Pvalue')
    
    count = 0
    for(pheno_name in colnames(input_data)[3:ncol(input_data)])
    {   
        message("[process Log-Rank] ", pheno_name)
        count = count + 1
        # 创建分析数据
        new_data <- data.frame(time = time, status = status, variable = input_data[[pheno_name]])
        # colnames(new_data)[3] = pheno_name 
        new_data <- new_data[complete.cases(new_data), ]
        new_data <- new_data[apply(new_data, 1, function(x){sum(x=='')}) == 0, ]
        # 只有一组数据，无法分析
        if(length(unique(new_data[, 3])) == 1 )
        {
            pvalue_data[count, ] = c(pheno_name, 'NA')
            next
        }

        # 拟合
        # formula <- as.formula( paste('Surv(time, status) ~ ', pheno_name, sep = '') )
        # fit <- survfit(formula, data = new_data)  # 注：不能这么做，否则后续绘图会报错
        # fit <- survfit(as.formula( paste('Surv(time, status) ~ ', pheno_name, sep = '') ), data = new_data)  
        # surv_diff <- survdiff(as.formula( paste('Surv(time, status) ~ ', pheno_name, sep = '') ), data = new_data) # 差异分析
        
        # 由于列名存在‘-’字符，不支持，所以表型统一命名为variable,ggsurvplot中添加legend.title
        fit <- survfit(Surv(time, status) ~ variable, data = new_data)
        surv_diff <- survdiff(Surv(time, status) ~ variable, data = new_data) # 差异分析
        pvalue    <- pchisq(surv_diff$chisq, length(surv_diff$n)-1, lower.tail = FALSE) # 获取P值
        pvalue_data[count, ] = c(pheno_name, pvalue)

        # 拟合详细结果输出
        surv_info <- surv_summary(fit, data = new_data)
        write.table(surv_info, paste0(output_dir, '/survival.LogRank.surv_info.', pheno_name, '.txt' ), row.names = F, quote = F, sep = '\t')

        # 生存曲线绘图
        p <- ggsurvplot(fit, data = new_data,
                    pval = TRUE, conf.int = TRUE,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "strata", # Change line type by groups
                    surv.median.line = "hv", # Specify median survival
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    legend.title=pheno_name,
                    palette = "npg" #杂志nature的配色
                    # palette = "aaas", #杂志Science的配色
                    # palette = "jco", #按jco杂志配色方案

                    )

        is_new_page = TRUE
        if(count == 1) is_new_page = FALSE  # 第一幅图要设置newpage = FALSE, 否则会多一个空白页
        result = tryCatch(print(p, newpage= is_new_page),
                            error=function(e){cat("ERROR :",conditionMessage(e),"\n")}
                            ) # 注部分图可能无法绘制，报错
        
    }
    # 输出组间差异P值
    write.table(pvalue_data, paste0(output_dir, '/survival.LogRank.pvalue..txt' ), row.names = F, quote = F, sep = '\t')

    dev.off()
}

##################
# （3）COX分析
##################
if(model == 'COX')
{   
    cat('Model : COX \n')
    new_data <- input_data[complete.cases(input_data), ]
    new_data <- new_data[apply(new_data, 1, function(x){sum(x=='')}) == 0, ]
    pheno_names <- colnames(input_data)[3:ncol(input_data)]

    # 拟合
    formula <- as.formula( paste0('Surv(time, status) ~ ', paste(pheno_names, collapse = ' + ') ))
    fit <- coxph(as.formula( paste0('Surv(time, status) ~ ', paste(pheno_names, collapse = ' + ') )), data = new_data)

    # 结果输出
    result <- summary(fit)$coefficients
    result <- cbind(row.names(result), result)
    colnames(result)[1] = 'Pheno'
    write.table(result, paste0(output_dir, '/survival.cox.txt' ) , row.names = F, quote = F, sep = '\t')
}

##################
# （4）My.stepwise.coxph分析
##################
if(model == 'My.stepwise.coxph')
{   
    cat('Model : My.stepwise.coxph \n')
    new_data <- input_data[complete.cases(input_data), ]
    new_data <- new_data[apply(new_data, 1, function(x){sum(x=='')}) == 0, ]

    time_var <- colnames(input_data)[1]
    event_var <- colnames(input_data)[2]
    pheno_names <- colnames(input_data)[3:ncol(input_data)]

    # 分析
    My.stepwise.coxph(Time = time_var, Status = event_var, variable.list = pheno_names, data = new_data)
    message("最后一行输出的变量名称就是最优组合结果")
}


##################
# （5）survivalROC分析
##################
if(model == 'survivalROC')
{   
    cat('Model : survivalROC \n')
    if(is.null(time_survivalroc))
    {
        message("survivalROC模型下，必须填写 time_survivalroc")
        q()
    }

    time <- input_data$time
    status <- input_data$status

    # 对每一个分层表型进行计算
    pdf(paste0(output_dir, '/survival.survivalROC.pdf'), width = width)
    
    count = 0
    for(pheno_name in colnames(input_data)[3:ncol(input_data)])
    {   
        message("[process survivalROC] ", pheno_name)
        count = count + 1
        # 创建分析数据
        new_data <- data.frame(time = time, status = status, variable = input_data[[pheno_name]])
        # colnames(new_data)[3] = pheno_name 
        new_data <- new_data[complete.cases(new_data), ]
        new_data <- new_data[apply(new_data, 1, function(x){sum(x=='')}) == 0, ]
        # 只有一组数据，无法分析
        if(length(unique(new_data[, 3])) == 1 )
        {
            pvalue_data[count, ] = c(pheno_name, 'NA')
            next
        }

        my_surv = survivalROC(new_data$time, new_data$status, new_data$variable, entry = NULL, predict.time = survivalroc_time, cut.values = NULL, method = survivalroc_method, lambda = NULL, span = NULL, window = "symmetric")
        ## NNE法
        plot(my_surv$FP, my_surv$TP, ## x=FP,y=TP
            type="l",col="#BC3C29FF", ##线条设置
            xlim=c(0,1), ylim=c(0,1),   
            xlab="FP",  
            ylab="TP",
            main=paste0("Time dependent ROC :", pheno_name))## \n换行符
        abline(0,1,col="gray",lty=2)##线条颜色
        info = paste0("AUC of ", survivalroc_method, " = ", round(my_surv$AUC,3))
        legend(0.6, info,
                 x.intersp=1, y.intersp=0.8,
                 lty= 1 ,lwd= 2,col=c("#BC3C29FF"),
                 bty = "n",# bty框的类型
                 seg.len=1,cex=0.8)# 
        
    }
    dev.off()
}
