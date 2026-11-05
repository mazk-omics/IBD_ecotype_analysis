# 直接读取已有官方绘图函数的颜色；本脚本只导出色表，不改动EcoTyper源码。
# ecotyper_interpretation环境运行，颜色与冻结原图的State/CellType/Ecotype语义一致。
source('software/ecotyper/pipeline/lib/misc.R')
cts = sort(unique(read.delim('06_validation/paper_figures_20261012/source_data/frozen_states.tsv')$CellType))
type_colors = data.frame(Kind='CellType', Label=cts, Color=get_colors(length(cts),column='CellType'))
state_colors = data.frame(Kind='State', Label=sprintf('S%02d',1:7), Color=get_colors(7,column='State'))
eco_colors = data.frame(Kind='Ecotype', Label=paste0('E',1:7), Color=get_colors(7,column='Ecotype'))
colors = rbind(type_colors,state_colors,eco_colors)
write.table(colors,'06_validation/paper_figures_20261012/source_data/official_annotation_colors.tsv',sep='\t',row.names=FALSE,quote=FALSE)
