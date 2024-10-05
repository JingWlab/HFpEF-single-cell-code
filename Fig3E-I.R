##### Fig 3E
N <- read.table("/NASdata/wuwq/HFpEF/ALL/HFp_CM_reg.txt",  # 读取的数据文件名称，这里文件是放在工作目录下
                header=T, # 数据集第一行为变量名
                row.names=1, # 第一列为行名
                sep="\t") # 指定分隔符号
N<-as.matrix(N)
p <- pheatmap(N, scale="none",
              border="black", # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              legend = T, # 添加图例
              legend_breaks=c(1,3,5),
              legend_labels = c(1,3,5),
              fontsize_row = 5, # 分别设置横向和纵向字体大小
              fontsize_col = 8,
              clustering_distance_rows = "euclidean", # 设置聚类的距离类型
              clustering_method="centroid",
              fontsize_number = 5,
              color =  colorRampPalette(c("white","pink","firebrick3"))(100))
print(p)


##### Fig 3F
#main color of the circle
top_languagesQ<-c(
  'Ppargc1a',
  'Atf6',
  'E2f6',
  'Mitf')
top_language_colorsQ <- c(
  '#efb306',
  '#e8351e',
  '#cd023d',
  '#0f8096',
  '#000000'
)
# name of the 4 circles
names(top_language_colorsQ) <- c(
  'Ppargc1a',
  'Atf6',
  'E2f6',
  'Mitf',
  'All'
)
#lines from inner to outer
edges1Q<-read.table("/NASdata/wuwq/HFpEF/ALL/edges1Q.txt", header=TRUE,sep = "")
edges1Q$total_code<-as.numeric(edges1Q$total_code)
edges2Q<-read.table("/NASdata/wuwq/HFpEF/ALL/edges2Q.txt", header=TRUE,sep = "")
edges2Q$from<- ""
edges2Q$total_code<- as.numeric(edges2Q$total_code)
edgesQ <- bind_rows(edges1Q, edges2Q)
vertices1Q<-read.table("/NASdata/wuwq/HFpEF/ALL/vertices1Q.txt", header=TRUE,sep = "")
vertices1Q$total_code<-as.numeric(vertices1Q$total_code)
vertices1Q$level<-as.numeric(vertices1Q$level)
vertices2Q<-read.table("/NASdata/wuwq/HFpEF/ALL/vertices2Q.txt", header=TRUE,sep = "")
vertices2Q$total_code<-as.numeric(vertices2Q$total_code)
vertices2Q$level<-as.numeric(vertices2Q$level)
vertices3Q <- tibble(
  node = '', language = NA, total_code = 0, level = 3
)
verticesQ <- bind_rows(vertices1Q, vertices2Q, vertices3Q) %>%
  mutate(
    radius = total_code**(1.8), # scaling circles
    language = factor(language, names(top_language_colorsQ))
  ) %>%
  arrange(level, language, node)

graphQ <- graph_from_data_frame(edgesQ, vertices = verticesQ)
# create custom layout by updating existing circle layout
layoutQ <- create_layout(graphQ, layout = 'circle')
outer_circleQ <- layoutQ %>%
  filter(level == 1) %>%
  mutate(language = factor(language, names(top_language_colorsQ))) %>%
  arrange(language, desc(name)) %>%
  mutate(
    x = cos((row_number() - 1) / 295 * 2 * pi),
    y = sin((row_number() - 1) / 295 * 2 * pi)
  )
# positioning circle centers manually by specifying polar coords
anglesQ <- c(30, 120, 210, 300, 0)
radiiQ <- c(0.5, 0.5, 0.5, 0.5, 0)
centersQ <- tibble(
  x = radiiQ * cos(anglesQ / 180 * pi),
  y = radiiQ * sin(anglesQ / 180 * pi)
)
inner_circleQ <- bind_cols(centersQ, select(filter(layoutQ, level != 1), -x, -y))
layoutQ[] <- bind_rows(outer_circleQ, inner_circleQ) %>% 
  arrange(.ggraph.index)
ggraph(layoutQ) +
  geom_edge_diagonal(
    aes(edge_color = node1.language, edge_alpha = as.factor(main)),
    edge_width = 0.2, show.legend = FALSE
  ) +
  geom_node_point(
    aes(size = radius*10000, color = language),
    alpha = 0.5, show.legend = FALSE
  ) +
  geom_node_text(
    aes(
      x = 1.0175* x,
      y = 1.0175 * y,
      label = name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = !(name %in% top_languagesQ)
    ),
    size = 2, hjust = 'outward', family = 'Oswald'
  ) +
  scale_edge_color_manual(values = top_language_colorsQ) +
  scale_color_manual(values = top_language_colorsQ) +
  scale_size_area(max_size = 70) +
  scale_edge_alpha_manual(values = c(0.6, 5)) +
  coord_fixed() +
  theme_void() +
  theme(
    text = element_text(family = 'Oswald'),
    legend.position = c(0.645, 0.51),
    plot.title = element_text(
      face = 'bold', hjust = 0.5, size = 5, margin = margin(t = 45, b = 3)
    ),
    plot.subtitle = element_text(
      face = 'plain', hjust = 0.5, size = 13, margin = margin(t = 5, b = 3)),
    plot.caption = element_text(
      face = 'plain', color = '#dedede', size = 8, hjust = 1,
      margin = margin(b = 20)
    )
  )


##### Fig 3G
N <- read.table("/NASdata/wuwq/HFpEF/ALL/heatmap.txt",  
                header=T, 
                row.names=1, 
                sep="\t") 
N<-as.matrix(N)
p <- pheatmap(N, scale="none",
              border="white", 
              cluster_cols = F, 
              cluster_rows = F,
              legend = T, 
              fontsize_row = F, 
              fontsize_col = 8,
              clustering_distance_rows = "euclidean", 
              clustering_method="centroid",
              fontsize_number = 5,
              color =  colorRampPalette(c("#4169e1","#87cefa","#e0ffff","white"))(100)) 
print(p)


##### Fig 3H
#main color of the circle
Rtop_languagesQ<-c(
  'Ppargc1a',
  'Rxrg',
  'Mitf',
  'Uqcrb'
)
#Rtop_language_colorsQ <- c(
#  '#efb306',
#  '#4e54ac',
#  '#0f8096',
#  '#17a769',
#  '#000000'
#)
Rtop_language_colorsQ <- c(
  '#e8351e',
  '#4e54ac',
  '#cd023d',
  '#17a769',
  '#000000'
)
# name of the 4 circles
names(Rtop_language_colorsQ) <- c(
  'Ppargc1a',
  'Rxrg',
  'Mitf',
  'Uqcrb',
  'All'
)
#lines from inner to outer
Redges1Q<-read.table("/NASdata/wuwq/HFpEF/ALL/edges1Q_HFr.txt", header=TRUE,sep = "")
Redges1Q$total_code<-as.numeric(Redges1Q$total_code)
Redges2Q<-read.table("/NASdata/wuwq/HFpEF/ALL/edges2Q_HFr.txt", header=TRUE,sep = "")
Redges2Q$from<- ""
Redges2Q$total_code<- as.numeric(Redges2Q$total_code)
RedgesQ <- bind_rows(Redges1Q, Redges2Q)
Rvertices1Q<-read.table("/NASdata/wuwq/HFpEF/ALL/vertices1Q_HFr.txt", header=TRUE,sep = "")
Rvertices1Q$total_code<-as.numeric(Rvertices1Q$total_code)
Rvertices1Q$level<-as.numeric(Rvertices1Q$level)
Rvertices2Q<-read.table("/NASdata/wuwq/HFpEF/ALL/vertices2Q_HFr.txt", header=TRUE,sep = "")
Rvertices2Q$total_code<-as.numeric(Rvertices2Q$total_code)
Rvertices2Q$level<-as.numeric(Rvertices2Q$level)
Rvertices3Q <- tibble(
  node = '', language = NA, total_code = 0, level = 3
)
RverticesQ <- bind_rows(Rvertices1Q, Rvertices2Q, Rvertices3Q) %>%
  mutate(
    radius = total_code**(1.8), # scaling circles
    language = factor(language, names(Rtop_language_colorsQ))
  ) %>%
  arrange(level, language, node)

graphQ <- graph_from_data_frame(RedgesQ, vertices = RverticesQ)


# create custom layout by updating existing circle layout
layoutQ <- create_layout(graphQ, layout = 'circle')
###此处乘的是最终去重后的node数289
outer_circleQ <- layoutQ %>%
  filter(level == 1) %>%
  mutate(language = factor(language, names(Rtop_language_colorsQ))) %>%
  arrange(language, desc(name)) %>%
  mutate(
    x = cos((row_number() - 1) / 289 * 2 * pi),
    y = sin((row_number() - 1) / 289 * 2 * pi)
  )
# positioning circle centers manually by specifying polar coords
anglesQ <- c(10,100,190,280, 0)
radiiQ <- c(0.5, 0.5, 0.5, 0.5, 0)
centersQ <- tibble(
  x = radiiQ * cos(anglesQ / 180 * pi),
  y = radiiQ * sin(anglesQ / 180 * pi)
)
inner_circleQ <- bind_cols(centersQ,  select(filter(layoutQ, level != 1), -x, -y))
layoutQ[] <- bind_rows(outer_circleQ, inner_circleQ) %>% 
  arrange(.ggraph.index)
ggraph(layoutQ) +
  geom_edge_diagonal(
    aes(edge_color = node1.language, edge_alpha = as.factor(main)),
    edge_width = 0.2, show.legend = FALSE
  ) +
  geom_node_point(
    aes(size = radius*10000, color = language),
    alpha = 0.5, show.legend = FALSE
  ) +
  geom_node_text(
    aes(
      x = 1.0175* x,
      y = 1.0175 * y,
      label = name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = !(name %in% Rtop_languagesQ)
    ),
    size = 2, hjust = 'outward', family = 'Oswald'
  ) +
  scale_edge_color_manual(values = Rtop_language_colorsQ) +
  scale_color_manual(values = Rtop_language_colorsQ) +
  scale_size_area(max_size = 70) +
  scale_edge_alpha_manual(values = c(0.6, 5)) +
  coord_fixed() +
  theme_void() +
  theme(
    text = element_text(family = 'Oswald'),
    legend.position = c(0.645, 0.51),
    plot.title = element_text(
      face = 'bold', hjust = 0.5, size = 5, margin = margin(t = 45, b = 3)
    ),
    plot.subtitle = element_text(
      face = 'plain', hjust = 0.5, size = 13, margin = margin(t = 5, b = 3)),
    plot.caption = element_text(
      face = 'plain', color = '#dedede', size = 8, hjust = 1,
      margin = margin(b = 20)
    )
  )


##### Fig 3I
### 基于SCENIC分析得到的单个细胞binary activity, 根据CM1-8亚群的barcode，在excel中将CM3亚群中细胞binary activity降序
### 作图在excel 完成
