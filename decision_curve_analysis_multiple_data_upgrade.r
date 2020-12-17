#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: decision_curve_analysis_multiple_data.r -p <prefix> -f <file> -n <string> [--xmin <numeric>  --xmax <numeric>  --ymax <numeric> --pdf_width <numeric> --pdf_height <numeric>   --rlib <dir> ]

Options:
    -f, --files <file>            分析的数据文件列表。包含表头，第一列样本名，第二列分组信息（只能有1、0）,第三列往后是用于数据分析的所有自变量。多个文件时，用“逗号”分隔。注：文件的列名字不能有空格、‘-’、‘+’等特殊符号
    -n, --names <string>          与files中的文件一一对应，用于图像展示的命名。
    -p, --output_prefix <prefix>  结果输出前缀，示例：./abc
    --xmin <numeric>              DCA绘图中，x轴最小值 [default: 0]
    --xmax <numeric>              DCA绘图中，x轴最大值 [default: 1]
    --ymax <numeric>              DCA绘图中，y轴最大值 [default: None]
    --pdf_width <numeric>         DCA图的PDF宽度 [default: 12]
    --pdf_height <numeric>        DCA图的PDF宽度 [default: 9]
    --rlib <dir>                  R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='李澎鹏，DCA分析\n')
files             <- opts$files
names             <- opts$names
output_prefix     <- opts$output_prefix
xmin              <- as.numeric(opts$xmin) 
xmax              <- as.numeric(opts$xmax) 
ymax              <- opts$ymax
pdf_width         <- as.numeric(opts$pdf_width) 
pdf_height        <- as.numeric(opts$pdf_height)
rlib              <- opts$rlib

files = unlist(strsplit(files, ","))  # 拆分文件列表
names = unlist(strsplit(names, ","))  # 拆分

# 加载R包
if(!is.null(rlib)) .libPaths(rlib)
library(DecisionCurve)
set.seed(91)


# 每个数据文件做DCA拟合分析
fit_decision_curve <- function(input_file){
    
    message("process file : ", input_file)
    # 读取
    input_data <- read.table(input_file,
                         header = TRUE,
                         sep = "\t",
                         quote = "",
                         check.names = FALSE,
                         )
    # 去掉缺失值
    input_data <- input_data[complete.cases(input_data), ]

    # 前两列设定固定名字，方便统一
    colnames(input_data)[1] <- 'Sample' 
    colnames(input_data)[2] <- 'Group'

    # 检查变量名称是否有特殊符号
    independent_variables <- colnames(input_data)[3:ncol(input_data)] 
    not_allow_character <- grep("-", independent_variables)
    if(length(not_allow_character) > 0){
        message("[ERROR] 输入文件表头含有字符 ‘-’， 请修改名称 : ", input_file)
        q()
    }

    # 拟合
    formula <- as.formula(paste('Group ~', paste(independent_variables , collapse = ' + ') ) )
    fit <- decision_curve(formula, 
        data = input_data, 
        family = binomial(link ='logit'),
        thresholds= seq(0,1, by = 0.01),
        confidence.intervals = 0.95,
        study.design = 'case-control'
    )
    fit
}

# （1）拟合
fit_list = lapply(files, fit_decision_curve)
names(fit_list) = names  # 设定列表名称


#（2）绘图
mycol <- c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)
mycol <- colors()[mycol]
color_list = c('red', 'blue', 'green', mycol)

dca_pdf_file = paste(output_prefix, ".pdf", sep = "", collapse = NULL)
pdf(dca_pdf_file, width = pdf_width, height = pdf_height)
message('plot: ', dca_pdf_file, "\n")

# ymax是否定义？
if(ymax == 'None')
{
    plot_decision_curve(fit_list, 
           curve.names = names(fit_list),
           cost.benefit.axis = FALSE,
           col = color_list[1:length(fit_list)],
           confidence.intervals = FALSE,
           standardize = FALSE,
           xlim = c(xmin, xmax),
           )
}else {
    ymax = as.numeric(ymax)
    plot_decision_curve(fit_list, 
           curve.names = names(fit_list),
           cost.benefit.axis = FALSE,
           col = color_list[1:length(fit_list)],
           confidence.intervals = FALSE,
           standardize = FALSE,
           xlim = c(xmin, xmax),
           ylim = c(-0.05, ymax),
  
           )
}

dev.off()



# (3) 数据输出
txt_file = paste(output_prefix, ".decision_curve.txt", sep = "", collapse = NULL)
red_pdf_file = paste(output_prefix, ".reduction.pdf", sep = "", collapse = NULL)
pdf(red_pdf_file, width = pdf_width, height = pdf_height)

output_data = list()
count = 0

message('write: ', txt_file, "\n")
message('plot: ', red_pdf_file, "\n")
for(name in names(fit_list))
{   
    count = count + 1
    fit_result <- fit_list[[name]]
    data_all <- subset(fit_result$derived.data, is.element(model, c("All", "None"))) # all 数据
    data_fit <- fit_result$derived.data # 拟合数据
    thresholds_count <- 101

    thresholds =   100 * data_all[, 'thresholds'][1:thresholds_count]
    all_sNB    =   data_all[, 'NB'][1:thresholds_count]  
    fit_sNB    =   data_fit[, 'NB'][1:thresholds_count] 
    reduction  = (fit_sNB - all_sNB) * 100 / (thresholds / (100 - thresholds)) # 公式参考地址 https://www.mskcc.org/sites/default/files/node/4511/documents/v3-worked-example-of-decision-curve-analysis-using-r.pdf


    output_data[[paste(name, 'thresholds', sep = '.')]] = thresholds
    output_data[[paste(name, 'ALL', sep = '.')]] = all_sNB
    output_data[[paste(name, 'Fit', sep = '.')]] = fit_sNB
    output_data[[paste(name, 'Reduction', sep = '.')]] = reduction

    # plot reduction
    #cut off reductions < 0
    reduction <- pmax(0, reduction) 
    if(count == 1) plot(thresholds, reduction, col = color_list[count], type="l", lwd=2, xlab="Threshold probability (%)", ylab="Net reduction in interventions per 100 patients") 
    if(count > 1) lines(thresholds, reduction, col = color_list[count], type="l", lwd=2, xlab="Threshold probability (%)", ylab="Net reduction in interventions per 100 patients") 
}
legend("bottomright", legend = names(fit_list), text.col = color_list[1:count] )
dev.off()
output_data = as.data.frame(output_data, check.names = FALSE)

write.table(output_data, 
    txt_file, 
    quote = FALSE, 
    sep = "\t", 
    row.names = FALSE,
    )




# # 系统 -> 初始化
# R_VERSION <- "3.4.0"
# R_SYS_VERSION <- R.version.string #返回R版本
# if(!grepl(R_VERSION, R_SYS_VERSION)){
#     cat(" -- [ERR] 必须使用 R", R_VERSION, "版本运行此脚本,", R_SYS_VERSION, "\n")
#     q()
# }
# .libPaths("/home/ganb/soft/R-3.4.0/lib64/R/library")

# # 检测 -> 脚本输入
# init <- function(){
#     tmp <- commandArgs();
#     tmp <- tmp[grep("--file=", tmp)]
#     tmp <- substring(tmp, 8, nchar(tmp))
#     tmp_pos <- gregexpr("[\\/]", tmp)
#     tmp_pos <- tmp_pos[[1]][length(tmp_pos[[1]])]
#     return(list(substring(tmp, 1, tmp_pos), substring(tmp, tmp_pos+1, nchar(tmp))))
# }
# SCRIPT <- init();
# SCRIPTDIR <- SCRIPT[[1]]
# SCRIPT <- SCRIPT[[2]]

# # 检测 -> 脚本输入
# ARGS <- commandArgs(T)
# if(length(ARGS) < 3){
#     cat("
# Program: Decision Curve Analysis，DCA
# Version: v1.0
# Contact: 李澎鹏

# Usage:   Rscript", SCRIPT, "OUTDIR PREFIX INPUTFile1 INPUTFile2 ...

# Options:

#          OUTDIR          结果文件输出目录
#          PREFIX          输出结果前缀
#          XMAX            x轴最大值
#          INPUTFile1      输入文件,包含表头，第一列样本名，第二列分组信息（只能有1、0）,第三列往后是用于逻辑回归的所有自变量         
#     \n");
#     q()
# }
# OUTDIR <- normalizePath(ARGS[1])
# PREFIX <- ARGS[2] 
# XMAX   <- ARGS[3] 
# INPUTFile_list <- ARGS[4:length(ARGS)]  # 数据路径

# ###################################################################### 主程序
# library(DecisionCurve)

# # (1) 数据读取并拟合

# fit_decision_curve <- function(input_file){
    
#     # 获取文件
#     input_prefix = input_file[1]
#     input_file = input_file[2]

#     # 读取
#     Data <- read.table(input_file,
#                          header = TRUE,
#                          sep = "\t",
#                          quote = "",
#                          check.names = FALSE,
#                          )
#     colnames(Data)[1] <- 'Sample' 
#     colnames(Data)[2] <- 'Group'
#     col_names <- colnames(Data)
#     independent_variables <- col_names[3:length(col_names)]

#     not_allow_character <- grep("-", independent_variables)
#     if(length(not_allow_character) > 0){
#         cat('[ERROR] independent variables has unsupport character [-] in ', input_prefix, ' : ', independent_variables[not_allow_character], "\n")
#         q()
#     }
#     formula <- as.formula(paste('Group ~', paste(independent_variables , collapse = ' + ') ) )

#     # 拟合
#     fit <- decision_curve(formula, 
#         data = Data, 
#         family = binomial(link ='logit'),
#         thresholds= seq(0,1, by = 0.01),
#         confidence.intervals = 0.95,
#         study.design = 'case-control'
#     )

#     fit
# }

# # 转list, 获取文件名前缀
# INPUTFile_list2 <- lapply(INPUTFile_list, function(input_file){
#             tmp_pos <- gregexpr("[\\/]", input_file)
#             tmp_pos <- tmp_pos[[1]][length(tmp_pos[[1]])]
#             file_name = substring(input_file, tmp_pos+1, nchar(input_file)) 
#             tmp_pos <- gregexpr("[.]", file_name)
#             tmp_pos <- tmp_pos[[1]][length(tmp_pos[[1]])]
#             file_prefix = substring(file_name, 1, tmp_pos - 1)
#             c(file_prefix, input_file) })
# # 修改列表名称
# names(INPUTFile_list2) = lapply(INPUTFile_list2, function(x){x[1]}) 

# fit_list = lapply(INPUTFile_list2, fit_decision_curve)

# #（2）绘图
# mycol <- c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)
# mycol <- colors()[mycol]
# color_list = c('red', 'blue', 'green', mycol)

# pdf_file = paste(OUTDIR, '/', PREFIX, ".pdf", sep = "", collapse = NULL)
# pdf(pdf_file, width = 12, height = 9)
# plot_decision_curve(fit_list, 
#            curve.names = names(fit_list),
#            cost.benefit.axis = FALSE,
#            col = color_list[1:length(fit_list)],
#            confidence.intervals = FALSE,
#            standardize = FALSE,
#            xlim = c(0, as.numeric(XMAX))
#            )
# dev.off()

# cat('plot: ', pdf_file, "\n")

# # (3) 数据输出
# txt_file = paste(OUTDIR, '/', PREFIX, ".decision_curve.txt", sep = "", collapse = NULL)
# output_data = list()
# name_lists = names(fit_list)
# pdf_file = paste(OUTDIR, '/', PREFIX, ".reduction.pdf", sep = "", collapse = NULL)
# pdf(pdf_file, width = 12, height = 9)
# count = 0
# for(name in name_lists)
# {   
#     count = count + 1
#     fit_result <- fit_list[[name]]
#     data_all <- subset(fit_result$derived.data, is.element(model, c("All", "None"))) # all 数据
#     data_fit <- fit_result$derived.data # 拟合数据
#     thresholds_count <- 101

#     thresholds =   100 * data_all[, 'thresholds'][1:thresholds_count]
#     all_sNB    =   data_all[, 'NB'][1:thresholds_count]  
#     fit_sNB    =   data_fit[, 'NB'][1:thresholds_count] 
#     reduction  = (fit_sNB - all_sNB) * 100 / (thresholds / (100 - thresholds)) # 公式参考地址 https://www.mskcc.org/sites/default/files/node/4511/documents/v3-worked-example-of-decision-curve-analysis-using-r.pdf


#     output_data[[paste(name, 'thresholds', sep = '.')]] = thresholds
#     output_data[[paste(name, 'ALL', sep = '.')]] = all_sNB
#     output_data[[paste(name, 'Fit', sep = '.')]] = fit_sNB
#     output_data[[paste(name, 'Reduction', sep = '.')]] = reduction

#     # plot reduction
#     #cut off reductions < 0
#     reduction <- pmax(0, reduction) 
#     if(count == 1) plot(thresholds, reduction, col = color_list[count], type="l", lwd=2, xlab="Threshold probability (%)", ylab="Net reduction in interventions per 100 patients") 
#     if(count > 1) lines(thresholds, reduction, col = color_list[count], type="l", lwd=2, xlab="Threshold probability (%)", ylab="Net reduction in interventions per 100 patients") 
# }
# legend("bottomright", legend = name_lists, text.col = color_list[1:count] )
# dev.off()
# output_data = as.data.frame(output_data, check.names = FALSE)

# write.table(output_data, 
#     txt_file, 
#     quote = FALSE, 
#     sep = "\t", 
#     row.names = FALSE,
#     )
# cat('write: ', 