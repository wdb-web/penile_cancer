library(readxl)
library(dplyr)
library(dplyr)
library(readxl)
library(stringr)
#library(memisc)
library(progressr)
library(lavaan)
library(dplyr)
library(readxl)
library(stringr)
library(memisc)
library(mediation)
setwd("/home/wangdongbin/work/word/data/2024/T007_GDMU")

path="/home/wangdongbin/work/word/data/2024/T007_GDMU/one_run/中介分析/组学数据汇总.xlsx"
"M14-粪便差异微生物"     "M14-广靶代谢组差异物质" "M14-粪便广靶代谢组"    
"广靶脂质组"    
readxl::excel_sheets (path)

gz_TG=read_excel(path,sheet="广靶脂质组")%>% data.frame

gz_TG_meta=gz_TG[,c(stringr::str_detect(colnames(gz_TG),"^A") %>% which,stringr::str_detect(colnames(gz_TG),"^B") %>% which)]%>%
    t%>%data.frame %>%  tibble::rownames_to_column("line")%>%
  mutate(#sample=stringr::str_extract(line,"[A|B].*"),
         group=stringr::str_extract(line,"A|B"))%>%{.[,colnames(.)!="line"]}

colnames(gz_TG_meta)[1:length(gz_TG$sig_Metabolites)]=gz_TG$sig_Metabolites

M2otu=read_excel(path,sheet="M14-粪便差异微生物")[-1,]%>% data.frame %>% na.omit
M2otu_data=M2otu[,c(stringr::str_detect(colnames(M2otu),"^A") %>% which,stringr::str_detect(colnames(M2otu),"^B") %>% which)]%>%
    t%>%data.frame %>%  tibble::rownames_to_column("line")%>%
  mutate(#sample=stringr::str_extract(line,"[A|B].*"),
         group=stringr::str_extract(line,"A|B"))%>%{.[,colnames(.)!="line"]}
colnames(M2otu_data)[1:length(M2otu$Sample)]=M2otu$Sample


M2_fengbian=read_excel(path,sheet="M14-粪便广靶代谢组")[-1,-1]%>% data.frame %>% na.omit
M2_fengbian_data=M2_fengbian[,c(stringr::str_detect(colnames(M2_fengbian),"A") %>% which,stringr::str_detect(colnames(M2_fengbian),"B") %>% which)]%>%
    t%>%data.frame %>%  tibble::rownames_to_column("line")%>%
  mutate(#sample=stringr::str_extract(line,"[A|B].*"),
         group=stringr::str_extract(line,"A|B"))%>%{.[,colnames(.)!="line"]}
colnames(M2_fengbian_data)[1:length(M2_fengbian$Sample)]=M2_fengbian$Sample



options(encoding ="UTF-8")


data=list(fengbian=(M2_fengbian_data),
         # xueye=(M2_xueye),
          gz_TG=(gz_TG_meta),
         otu=(M2otu_data))
data2=names(data)%>% lapply( function(.x){
  col=ncol(data[[.x]])
  a=data[[.x]]
  colnames(a)[-col] =paste0( .x,"_",colnames(data[[.x]])[-col])
  return(a)
})
names(data2)=names(data)
data=data2
x="gz_TG"
y="otu"
corP <- function(x,y,full) {
  cat(x,y)
  x_col=ncol(data[[x]])
  y_col=ncol(data[[y]])
  
  join_data=dplyr::full_join(data[[x]],data[[y]],by="group")
  lapply(colnames(data[[x]])[-x_col], function(.x){
    lapply(colnames(data[[y]])[-y_col], function(.y) {
      k=cor.test(join_data[,.x]%>%unlist() %>% as.numeric  ,join_data[,.y]%>%unlist() %>% as.numeric)
      data.frame(a=.x,b=.y,cor=k[["estimate"]][["cor"]],p=k$p.value)
    })%>%do.call(rbind,.)
  }
  )%>%do.call(rbind,.)
}
names(data)
data1=corP(x="fengbian",y="gz_TG")%>%#rbind(corP("xueye","pizhi"))%>%
  #rbind(corP("fengbian","pizhi"))%>%
  rbind(corP("gz_TG","otu"))%>%
  rbind(corP("fengbian","otu"))%>%na.omit()
data1->p
corP <- function(x, y, threshold = 0, pval_threshold = 0.05) {
  x_col <- ncol(data[[x]])
  y_col <- ncol(data[[y]])
  join_data <- dplyr::full_join(data[[x]], data[[y]], by = "group")
  
  results <- lapply(colnames(data[[x]])[-x_col], function(.x) {
    lapply(colnames(data[[y]])[-y_col], function(.y) {
      k <- cor.test(join_data[[.x]]%>% as.numeric, join_data[[.y]]%>% as.numeric, use = "complete.obs")
      if (abs(k$estimate) >= threshold && k$p.value <= pval_threshold) {
        data.frame(a = .x, b = .y, cor = k$estimate, p = k$p.value)
      }
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  return(results)
}

# 筛选显著相关的变量
data1 <- bind_rows(
  corP("fengbian", "gz_TG"),
  corP("gz_TG", "otu"),
  corP("otu", "fengbian")
) %>% na.omit()






paste3 <- function(out, input, sep = ".") {
  
  full_join(out, input,by="group")
}
datajoin=data%>%purrr::reduce(paste3)
all_combin=data%>% lapply( function(.x){
  x_col=ncol(.x)
  data.frame(sample=colnames(.x)[-x_col],group="A")
  
})%>%purrr::reduce(paste3)

all_combin=all_combin[,c(4,1,3)]
colnames(data1)[1:2]=colnames(all_combin)[1:2]
all_combin=left_join(all_combin,data1)%>% na.omit
colnames(data1)[1:2]=colnames(all_combin)[2:3]
all_combin=left_join(all_combin,data1,by = join_by(sample.x, sample.y))%>% na.omit

# 定义你的 do.data 函数
do.data = function(.x) {
  try({


          # cat(.x,.y,.z)
      #data.frame(data[["a"]][,.x],data[["b"]][,.y],data[["c"]][,.z])->datass
      # k=PROCESS(data1, y=.x, x=.z,meds =.y,
      #           #   covs="family_inc", 
      #           ci="boot", nsim=100, seed=1)
      datass = datajoin[,.x|>unlist()] |> na.omit()
    #colnames(line_data) = c("X", "M1", "M2", "Y")
      colnames(datass)=c("X","Y","M")
      model.m=lm(M~X,datass)
      model.y=lm(Y~X+M,datass)
      set.seed(1)
      summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 100))
      colnames(datass)=c("X","M","Y")
      model.m=lm(M~X,datass)
      model.y=lm(Y~X+M,datass)
      summary2=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 100))
      data.frame(group=paste0(.x,collapse= "==>"),
                 Pval_mediate=summary$d.avg.p,
                 Pval_direct=summary$z.avg.p,
                 Pval_mediate_inverse=summary2$d.avg.p,
                 Pval_direct_inverse=summary2$z.avg.p
      )
  })
}




# 创建进度条进度器
# 设置并行处理
future::plan("multisession", workers =30)  # Reduce the number of workers

options(future.globals.maxSize = 200 * 1024^3)

# 创建进度条
handlers(global = TRUE)
handlers("progress")
with_progress({
p <- progressr::progressor(steps = nrow(all_combin))
# 使用 future_apply 并传递进度条
get1 = future.apply::future_apply(all_combin[,1:3], 1, function(.x) {
    p()  # Update the progress bar
    result <- do.data(.x)
    return(result)
})
})
# 打印结果
saveRDS(get1,"get14.RDS")
get_gene_meta2=do.call(rbind,get1)




library(ggplot2)
library(ggalluvial)
##install.packages("ggalluvial")
library(gridExtra)
get_gene_meta2$Qval_mediate=p.adjust(get_gene_meta2$Pval_mediate,method = "BH")
get_gene_meta2$Qval_mediate_inverse=p.adjust(get_gene_meta2$Pval_mediate_inverse,method = "BH")
writexl::write_xlsx(get_gene_meta2,"中介M14.xlsx")

#dzs_get=dzs_get[which(dzs_get$Pval_mediate<0.05&dzs_get$Pval_mediate_inverse>0.05|dzs_get$Pval_mediate>0.05&dzs_get$pval_mediate_inverse<0.05),]
get_gene_meta2=get_gene_meta2[which((get_gene_meta2$Qval_mediate<0.05&get_gene_meta2$Pval_mediate_inverse>0.05)),]
get_gene_meta2$input= get_gene_meta2$group%>%str_extract("^.*?(?===>)")
get_gene_meta2$line= get_gene_meta2$group%>%str_extract("(?<=>).*?(?===>)")
get_gene_meta2$output= get_gene_meta2$group%>%str_remove(".*==>")
get_gene_meta2->net
net$fre=1
library(showtext)
showtext::showtext_auto()

p=ggplot(get_gene_meta2,aes(axis1 = net$input, axis2 = net$line,
                            axis3 = net$output,y= net$fre))+
  # scale_x_discrete(limits = c("Diet", "Microbiome", "Meta")) +
  geom_alluvium(aes(fill = get_gene_meta2$line),cex=0.4)+
  geom_stratum(width = 1/12, fill = "black", color = "grey")+
  #geom_stratum() +
  geom_stratum(width = 0.1) +
  
  theme_minimal()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  
  theme(legend.position="none",axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),panel.grid=element_blank())#+

ggsave("中介M14.pdf",p,width = 15,height = 12)

# 脂质
library(readxl)
library(dplyr)
library(dplyr)
library(readxl)
library(stringr)
#library(memisc)
library(progressr)
library(lavaan)
library(dplyr)
library(readxl)
library(stringr)
library(memisc)
library(mediation)
setwd("/data/project/wangdongbin/data/2024/T007_GDMU/one_run")

path="/data/project/wangdongbin/data/2024/T007_GDMU/one_run/中介分析/组学数据汇总.xlsx"
#"M14-粪便差异微生物"     "M14-广靶代谢组差异物质" "M14-粪便广靶代谢组"    "广靶脂质组"    
readxl::excel_sheets (path)

gz_TG=read_excel(path,sheet="M14-广靶代谢组差异物质")[-1,]%>% data.frame %>% na.omit

gz_TG_meta=gz_TG[,c(stringr::str_detect(colnames(gz_TG),"A") %>% which,stringr::str_detect(colnames(gz_TG),"B") %>% which)]%>%
    t%>%data.frame %>%  tibble::rownames_to_column("line")%>%
  mutate(#sample=stringr::str_extract(line,"[A|B].*"),
         group=stringr::str_extract(line,"A|B"))%>%{.[,colnames(.)!="line"]}

colnames(gz_TG_meta)[1:length(gz_TG$`...1`)]=gz_TG$`...1`

M2otu=read_excel(path,sheet="M14-粪便差异微生物")[-1,]%>% data.frame %>% na.omit
M2otu_data=M2otu[,c(stringr::str_detect(colnames(M2otu),"^A") %>% which,stringr::str_detect(colnames(M2otu),"^B") %>% which)]%>%
    t%>%data.frame %>%  tibble::rownames_to_column("line")%>%
  mutate(#sample=stringr::str_extract(line,"[A|B].*"),
         group=stringr::str_extract(line,"A|B"))%>%{.[,colnames(.)!="line"]}
colnames(M2otu_data)[1:length(M2otu$Sample)]=M2otu$Sample


M2_fengbian=read_excel(path,sheet="M14-粪便广靶代谢组")[-1,-1]%>% data.frame %>% na.omit
M2_fengbian_data=M2_fengbian[,c(stringr::str_detect(colnames(M2_fengbian),"A") %>% which,stringr::str_detect(colnames(M2_fengbian),"B") %>% which)]%>%
    t%>%data.frame %>%  tibble::rownames_to_column("line")%>%
  mutate(#sample=stringr::str_extract(line,"[A|B].*"),
         group=stringr::str_extract(line,"A|B"))%>%{.[,colnames(.)!="line"]}
colnames(M2_fengbian_data)[1:length(M2_fengbian$Sample)]=M2_fengbian$Sample



options(encoding ="UTF-8")


data=list(fengbian=(M2_fengbian_data),
         # xueye=(M2_xueye),
          gz_TG=(gz_TG_meta),
         otu=(M2otu_data))
data2=names(data)%>% lapply( function(.x){
  col=ncol(data[[.x]])
  a=data[[.x]]
  colnames(a)[-col] =paste0( .x,"_",colnames(data[[.x]])[-col])
  return(a)
})
names(data2)=names(data)
data=data2
x="gz_TG"
y="otu"
corP <- function(x,y,full) {
  cat(x,y)
  x_col=ncol(data[[x]])
  y_col=ncol(data[[y]])
  
  join_data=dplyr::full_join(data[[x]],data[[y]],by="group")
  lapply(colnames(data[[x]])[-x_col], function(.x){
    lapply(colnames(data[[y]])[-y_col], function(.y) {
      k=cor.test(join_data[,.x]%>%unlist() %>% as.numeric  ,join_data[,.y]%>%unlist() %>% as.numeric)
      data.frame(a=.x,b=.y,cor=k[["estimate"]][["cor"]],p=k$p.value)
    })%>%do.call(rbind,.)
  }
  )%>%do.call(rbind,.)
}
names(data)
data1=corP(x="fengbian",y="gz_TG")%>%#rbind(corP("xueye","pizhi"))%>%
  #rbind(corP("fengbian","pizhi"))%>%
  rbind(corP("gz_TG","otu"))%>%
  rbind(corP("fengbian","otu"))%>%na.omit()
data1->p
corP <- function(x, y, threshold = 0, pval_threshold = 0.05) {
  x_col <- ncol(data[[x]])
  y_col <- ncol(data[[y]])
  join_data <- dplyr::full_join(data[[x]], data[[y]], by = "group")
  
  results <- lapply(colnames(data[[x]])[-x_col], function(.x) {
    lapply(colnames(data[[y]])[-y_col], function(.y) {
      k <- cor.test(join_data[[.x]]%>% as.numeric, join_data[[.y]]%>% as.numeric, use = "complete.obs")
      if (abs(k$estimate) >= threshold && k$p.value <= pval_threshold) {
        data.frame(a = .x, b = .y, cor = k$estimate, p = k$p.value)
      }
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  return(results)
}

# 筛选显著相关的变量
data1 <- bind_rows(
  corP("fengbian", "gz_TG"),
  corP("gz_TG", "otu"),
  corP("otu", "fengbian")
) %>% na.omit()






paste3 <- function(out, input, sep = ".") {
  
  full_join(out, input,by="group")
}
datajoin=data%>%purrr::reduce(paste3)
all_combin=data%>% lapply( function(.x){
  x_col=ncol(.x)
  data.frame(sample=colnames(.x)[-x_col],group="A")
  
})%>%purrr::reduce(paste3)

all_combin=all_combin[,c(4,1,3)]
colnames(data1)[1:2]=colnames(all_combin)[1:2]
all_combin=left_join(all_combin,data1)%>% na.omit
colnames(data1)[1:2]=colnames(all_combin)[2:3]
all_combin=left_join(all_combin,data1,by = join_by(sample.x, sample.y))%>% na.omit

# 定义你的 do.data 函数
do.data = function(.x) {
  try({


          # cat(.x,.y,.z)
      #data.frame(data[["a"]][,.x],data[["b"]][,.y],data[["c"]][,.z])->datass
      # k=PROCESS(data1, y=.x, x=.z,meds =.y,
      #           #   covs="family_inc", 
      #           ci="boot", nsim=100, seed=1)
      datass = datajoin[,.x|>unlist()] |> na.omit() %>%apply(2,as.numeric) %>% data.frame
    #colnames(line_data) = c("X", "M1", "M2", "Y")
      colnames(datass)=c("X","Y","M")
      model.m=lm(M~X,datass)
      model.y=lm(Y~X+M,datass)
      set.seed(1)
      summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 100))
      colnames(datass)=c("X","M","Y")
      model.m=lm(M~X,datass)
      model.y=lm(Y~X+M,datass)
      summary2=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 100))
      data.frame(group=paste0(.x,collapse= "==>"),
                 Pval_mediate=summary$d.avg.p,
                 Pval_direct=summary$z.avg.p,
                 Pval_mediate_inverse=summary2$d.avg.p,
                 Pval_direct_inverse=summary2$z.avg.p
      )
  })
}




# 创建进度条进度器
# 设置并行处理
future::plan("multisession", workers =30)  # Reduce the number of workers

options(future.globals.maxSize = 200 * 1024^3)

# 创建进度条
handlers(global = TRUE)
handlers("progress")
with_progress({
p <- progressr::progressor(steps = nrow(all_combin))
# 使用 future_apply 并传递进度条
get1 = future.apply::future_apply(all_combin[,1:3], 1, function(.x) {
    p()  # Update the progress bar
    result <- do.data(.x)
    return(result)
})
})
# 打印结果
saveRDS(get1,"广靶代谢组差异物质.RDS")
get_gene_meta2=do.call(rbind,get1)




library(ggplot2)
library(ggalluvial)
##install.packages("ggalluvial")
library(gridExtra)
get_gene_meta2$Qval_mediate=p.adjust(get_gene_meta2$Pval_mediate,method = "BH")
get_gene_meta2$Qval_mediate_inverse=p.adjust(get_gene_meta2$Pval_mediate_inverse,method = "BH")
writexl::write_xlsx(get_gene_meta2,"中介M14.xlsx")

#dzs_get=dzs_get[which(dzs_get$Pval_mediate<0.05&dzs_get$Pval_mediate_inverse>0.05|dzs_get$Pval_mediate>0.05&dzs_get$pval_mediate_inverse<0.05),]
get_gene_meta2=get_gene_meta2[which((get_gene_meta2$Qval_mediate<0.05&get_gene_meta2$Pval_mediate_inverse>0.05)),]
get_gene_meta2$input= get_gene_meta2$group%>%str_extract("^.*?(?===>)")
get_gene_meta2$line= get_gene_meta2$group%>%str_extract("(?<=>).*?(?===>)")
get_gene_meta2$output= get_gene_meta2$group%>%str_remove(".*==>")
get_gene_meta2->net
net$fre=1
library(showtext)
showtext::showtext_auto()

p=ggplot(get_gene_meta2,aes(axis1 = net$input, axis2 = net$line,
                            axis3 = net$output,y= net$fre))+
  # scale_x_discrete(limits = c("Diet", "Microbiome", "Meta")) +
  geom_alluvium(aes(fill = get_gene_meta2$line),cex=0.4)+
  geom_stratum(width = 1/12, fill = "black", color = "grey")+
  #geom_stratum() +
  geom_stratum(width = 0.1) +
  
  theme_minimal()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  
  theme(legend.position="none",axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),panel.grid=element_blank())#+

ggsave("中介M14.pdf",p,width = 15,height = 12)

