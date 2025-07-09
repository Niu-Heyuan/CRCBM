data_CRC<-gistic@data
cnv.summary_CRC<-gistic@cnv.summary
cytoband.summary_CRC<-gistic@cytoband.summary
gene.summary_CRC<-gistic@gene.summary
gis.scores_CRC<-gistic@gis.scores
summary_CRC<-gistic@summary
data_list <- list(
  data_CRC= data_CRC,
  cnv_summary_CRC = cnv.summary_CRC,
  cytoband_summary_CRC = cytoband.summary_CRC,
  gene_summary_CRC= gene.summary_CRC,
  gis_scores_CRC = gis.scores_CRC,
  summary_CRC= summary_CRC)
# 绘制基因组图
gisticChromPlot(
  gistic = gistic,
  fdrCutOff = 0.05,
  markBands = "all",
  color = NULL,
  ref.build = "hg38",
  cytobandOffset = 0.03,
  txtSize = 0.8,
  cytobandTxtSize = 0.6,
  mutGenesTxtSize = 0.6
)
# 绘制GISTIC气泡图
gisticBubblePlot(
  gistic = gistic,
  color = NULL,
  markBands = NULL,
  fdrCutOff = 0.05,
  log_y = TRUE,
  txtSize = 3
)
# 绘制GISTIC条形图
gisticBubblePlot(
  gistic = gistic,
  color = NULL,
  markBands = NULL,
  fdrCutOff = 0.05,
  log_y = TRUE,
  txtSize = 3
)
