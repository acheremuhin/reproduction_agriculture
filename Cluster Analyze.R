library("tidyverse") # коллекция пакетов от Hadley Wickham

library("lmtest") # тесты в линейных моделях
library("sandwich") # оценки ковариационной матрицы робастные к гетероскедастичности
library("erer") # подборка пакетов для эмпирических исследований
library("AUC") # подсчёт показателя AUC
library("mfx") # для предельных эффектов в logit/probit
library("estimatr") # модели с робастными стандартными ошибками
library("e1071") # SVM
library("lubridate")  # работа с датами
library("quantmod")

library("GGally") # матрица диаграмм рассеяния
library("lattice") # конкурент ggplot2
library("ggplot2") # базовый графический пакет
library("vcd") # мозаичный график
library("hexbin") # график из шестиугольников
library("sjPlot") # визуализация результатов МНК
library("factoextra") # визуализация для метода главных компонент и не только
library("sm") # + графики
library("vioplot") # + графики
library("beeswarm")
library("lattice")

library("reshape2") # длинные <-> широкие таблицы
library("psych") # описательные статистики
library("skimr") # описательные статистики

library("glmnet") # LASSO
library("HSAUR")
library("sgof")
library("car") # для тестирования линейных гипотез, подсчёта vif
library("spikeslab")
library("quantreg") 

library("spikeslab") # байесовская регрессия пик-плато
library("quantreg") # квантильная регрессия
library("MCMCpack") # набор моделей с байесовским подходом

library("devtools") # разработка пакетов

library("caret") # подбор параметров с помощью кросс-валидации
library("AER")
library("ivpack") # инструментальные переменные

library("zoo") # нерегулярные временные ряды
library("xts") # еще ряды
library("forecast") # ARMA, экспоненциальное сглаживание
library("rugarch") # хорош для GARCH

Base_cluster <- Date
# Base_cluster_non <- Base_cluster[Base_cluster$X1>quantile(Base$X1, p=c(0.05)) & Base_cluster$X1<quantile(Base$X1, p=c(0.95)) &
#                                Base_cluster$X2>quantile(Base$X2, p=c(0.05)) & Base_cluster$X2<quantile(Base$X2, p=c(0.95)) &
#                                 Base_cluster$X3>quantile(Base$X3, p=c(0.05)) & Base_cluster$X3<quantile(Base$X3, p=c(0.95)) &
#                                  Base_cluster$X4>quantile(Base$X4, p=c(0.05)) & Base_cluster$X4<quantile(Base$X4, p=c(0.95)) &
#                                  Base_cluster$X5>quantile(Base$X5, p=c(0.05)) & Base_cluster$X5<quantile(Base$X5, p=c(0.95)) &
#                                  Base_cluster$X6>quantile(Base$X6, p=c(0.05)) & Base_cluster$X6<quantile(Base$X6, p=c(0.95)) &
#                                  Base_cluster$X7A>quantile(Base$X7A, p=c(0.05)) & Base_cluster$X7A<quantile(Base$X7A, p=c(0.95)) &
#                                  Base_cluster$X7B>quantile(Base$X7B, p=c(0.05)) & Base_cluster$X7B<quantile(Base$X7B, p=c(0.95)) &
#                                  Base_cluster$X8>quantile(Base$X8, p=c(0.05)) & Base_cluster$X8<quantile(Base$X8, p=c(0.95)) &
#                                  Base_cluster$X9A>quantile(Base$X9A, p=c(0.05)) & Base_cluster$X9A<quantile(Base$X9A, p=c(0.95)) &
#                                  Base_cluster$X9B>quantile(Base$X9B, p=c(0.05)) & Base_cluster$X9B<quantile(Base$X9B, p=c(0.95)) &
#                                  Base_cluster$X10>quantile(Base$X10, p=c(0.05)) & Base_cluster$X10<quantile(Base$X10, p=c(0.95)) &
#                                  Base_cluster$X11>quantile(Base$X11, p=c(0.05)) & Base_cluster$X11<quantile(Base$X11, p=c(0.95)) &
#                                 Base_cluster$X12>quantile(Base$X12, p=c(0.05)) & Base_cluster$X12<quantile(Base$X12, p=c(0.95)) &
#                                  Base_cluster$X13>quantile(Base$X13, p=c(0.05)) & Base_cluster$X13<quantile(Base$X13, p=c(0.95)) &
#                                  Base_cluster$X14>quantile(Base$X14, p=c(0.05)) & Base_cluster$X14<quantile(Base$X14, p=c(0.95)) &
#                                  Base_cluster$X15>quantile(Base$X15, p=c(0.05)) & Base_cluster$X15<quantile(Base$X15, p=c(0.95)),]


Base_cluster_non <- Base_cluster[Base_cluster$X1>quantile(Base_cluster$X1, p=c(0.05)) & Base_cluster$X1<quantile(Base_cluster$X1, p=c(0.95)) &
                                   Base_cluster$X2>quantile(Base_cluster$X2, p=c(0.05)) & Base_cluster$X2<quantile(Base_cluster$X2, p=c(0.95)) &
                                   Base_cluster$X3>quantile(Base_cluster$X3, p=c(0.05)) & Base_cluster$X3<quantile(Base_cluster$X3, p=c(0.95)) &
                                   Base_cluster$X4>quantile(Base_cluster$X4, p=c(0.05)) & Base_cluster$X4<quantile(Base_cluster$X4, p=c(0.95)) &
                                   Base_cluster$X5>quantile(Base_cluster$X5, p=c(0.05)) & Base_cluster$X5<quantile(Base_cluster$X5, p=c(0.95)),]

Base <- Base_cluster_non[,1:5]
# Определяем максимальное количество кластеров по дисперсии
k.max <- 10 # максимальное число кластеров
wss <- sapply(1:k.max, function(k)
{kmeans(Base, k)$tot.withinss})
plot(1:k.max, wss, type="b", pch = 19, frame = FALSE,
     xlab="Number of cluster",
     ylab="Total sum square")
# Формируем график с помощью fviz_nbclust()
library(factoextra)
library(cluster)
fviz_nbclust(Base, kmeans, method = "wss")  # Построение графика зависимости дисперсии от количества кластеров
# Находим оптимальное число кластеров по статистике разрыва (Grap)
  # k-means
set.seed(123)
gap_stat <- clusGap(Base, FUN = kmeans, nstart = 10, # k-means
                 K.max = 10, B = 50)
          # Печать и визуализация результатов
print(gap_stat, method = "firstmax")
fviz_gap_stat(gap_stat)  # Находим оптимальное число кластеров по статистике разрыва (Grap)
  # k-medoids
set.seed(123)
gap_stat <- clusGap(Base, FUN = pam, 
                    K.max = 20, B = 50)
print(gap_stat, method = "globalSEmax")
fviz_gap_stat(gap_stat)
# Оптимальное разбиение через диаграмму силуэтов 
(k.pam <- pam(Base, k=3))
fviz_silhouette(silhouette(k.pam)) # Диаграмма силуэтов (чем выше среднее, тем больше качество силуэтов)
fviz_nbclust(Base, kmeans, method = "silhouette")  # Определяем оптимальное разбиение через среднюю величину силуэтов (чем больше)
fviz_cluster(pam(Base, 3), stand=FALSE)
fviz_cluster(clara(Base, 3), stand=FALSE,
             ellipse.type = "norm", ellipse.level = 0.7) 
# Оптимальное количество кластеров
library(NbClust)
nb <- NbClust(Base, distance = "minkowski", min.nc = 2,
              max.nc = 10, method = "average", index ="all")
nb$Best.nc    # Оптимальная схема объединения в кластеры
fviz_nbclust(nb) + theme_minimal()   # График
# Смотрим вручную какие-то кластеры и их статистику
clus <- kmeans(Base, centers = 3)
#clus <- fviz_cluster(pam(Base, 3), stand=FALSE)
#hist(as.numeric(clus$data$cluster))
clus$centers
clus$size
clus$cluster
# Разбиение методом k-средних получилось нечетким, поэтому разбиваем через вот эту диаграмму
fviz_cluster(pam(Base, 3), stand=FALSE,geom = c("point"), pointsize = 0.8)
aaa <- fviz_cluster(pam(Base, 3), stand=FALSE)
Base_with_cluster <- cbind(aaa$data$cluster,Base_cluster_non[,c(1:5,18:40)])
glimpse(Base_with_cluster)
aab<-Base_with_cluster[Base_with_cluster$`aaa$data$cluster`=="3",2:6]
mx1<-c(mean(aab[,1]),mean(aab[,2]),mean(aab[,3]),mean(aab[,4]),mean(aab[,5]))
sx1<-c(sd(aab[,1]),sd(aab[,2]),sd(aab[,3]),sd(aab[,4]),sd(aab[,5]))
Base_reg_1<-Base_with_cluster[Base_with_cluster$`aaa$data$cluster`=="1",-1]
Base_reg_2<-Base_with_cluster[Base_with_cluster$`aaa$data$cluster`=="2",-1]
Base_reg_3<-Base_with_cluster[Base_with_cluster$`aaa$data$cluster`=="3",-1]
library(leaps)
library(MASS)
library(dplyr)
# Step-wise regression X1
Base_reg_1_1<-mutate(Base_reg_1[,c(1,6:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P5=as.factor(P5),P6=as.factor(P6),P7=as.factor(P7),P8=as.factor(P8),P9=as.factor(P9),P10=as.factor(P10),P11=as.factor(P11),P12=as.factor(P12))
Base_reg_1_1
full.model <- lm(X1 ~., data = Base_reg_1_1)
step.model_1 <- stepAIC(full.model, direction = "both", 
                      trace = FALSE)
summary(step.model_1)
Base_reg_1_2<-mutate(Base_reg_1[,c(2,6:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P5=as.factor(P5),P6=as.factor(P6),P7=as.factor(P7),P8=as.factor(P8),P9=as.factor(P9),P10=as.factor(P10),P11=as.factor(P11),P12=as.factor(P12))
Base_reg_1_2
full.model <- lm(X2 ~., data = Base_reg_1_2)
step.model_2 <- stepAIC(full.model, direction = "both", 
                        trace = FALSE)
step.model_2_2 <- lm(X2 ~ Z4+P8+Z11, data = Base_reg_1_2)
summary(step.model_2_2)
Base_reg_1_3<-mutate(Base_reg_1[,c(3,6:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P5=as.factor(P5),P6=as.factor(P6),P7=as.factor(P7),P8=as.factor(P8),P9=as.factor(P9),P10=as.factor(P10),P11=as.factor(P11),P12=as.factor(P12))
Base_reg_1_3
full.model <- lm(X3 ~., data = Base_reg_1_3)
step.model_3 <- stepAIC(full.model, direction = "both", 
                        trace = FALSE)
step.model_3_2 <- lm(X3 ~ Z3+P5+P8, data = Base_reg_1_3)
summary(step.model_3_2)
Base_reg_1_4<-mutate(Base_reg_1[,c(4,6:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P5=as.factor(P5),P6=as.factor(P6),P7=as.factor(P7),P8=as.factor(P8),P9=as.factor(P9),P10=as.factor(P10),P11=as.factor(P11),P12=as.factor(P12))
Base_reg_1_4
full.model <- lm(X4 ~., data = Base_reg_1_4)
step.model_4 <- stepAIC(full.model, direction = "both", 
                        trace = FALSE)
step.model_4_2 <- lm(X4 ~ P8, data = Base_reg_1_4)
summary(step.model_4_2)
Base_reg_1_5<-mutate(Base_reg_1[,c(5,6:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P5=as.factor(P5),P6=as.factor(P6),P7=as.factor(P7),P8=as.factor(P8),P9=as.factor(P9),P10=as.factor(P10),P11=as.factor(P11),P12=as.factor(P12))
Base_reg_1_5
full.model <- lm(X5 ~., data = Base_reg_1_5)
step.model_5 <- stepAIC(full.model, direction = "both", 
                        trace = FALSE)
summary(step.model_5)
# Step-wise regression X2
Base_reg_2_1<-mutate(Base_reg_2[,c(1,6:15,19:20,22,24:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P8=as.factor(P8),P9=as.factor(P9),P11=as.factor(P11))
glimpse(Base_reg_2_1)
full.model <- lm(X1 ~., data = Base_reg_2_1)
step.model_6 <- stepAIC(full.model, direction = "both", 
                        trace = FALSE)
step.model_6_2 <- lm(X1 ~ Z1+Z5, data = Base_reg_2_1)
summary(step.model_6_2)
Base_reg_2_2<-mutate(Base_reg_2[,c(2,6:15,19:20,22,24:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P8=as.factor(P8),P9=as.factor(P9),P11=as.factor(P11))
full.model <- lm(X2 ~., data = Base_reg_2_2)
step.model_7 <- stepAIC(full.model, direction = "both", 
                        trace = FALSE)
summary(step.model_7)
Base_reg_2_3<-mutate(Base_reg_2[,c(3,6:15,19:20,22,24:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P8=as.factor(P8),P9=as.factor(P9),P11=as.factor(P11))
full.model <- lm(X3 ~., data = Base_reg_2_3)
step.model_8 <- stepAIC(full.model, direction = "both", 
                        trace = FALSE)
step.model_8_2 <- lm(X3 ~ Z5+Z7, data = Base_reg_2_3)
summary(step.model_8_2)
Base_reg_2_4<-mutate(Base_reg_2[,c(4,6:15,19:20,22,24:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P8=as.factor(P8),P9=as.factor(P9),P11=as.factor(P11))
full.model <- lm(X4 ~., data = Base_reg_2_4)
step.model_9 <- stepAIC(full.model, direction = "both", 
                        trace = FALSE)
step.model_9_2 <- lm(X4 ~ Z5+Z6+P8+Z11, data = Base_reg_2_4)
summary(step.model_9_2)
Base_reg_2_5<-mutate(Base_reg_2[,c(5,6:15,19:20,22,24:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P8=as.factor(P8),P9=as.factor(P9),P11=as.factor(P11))
full.model <- lm(X5 ~., data = Base_reg_2_5)
step.model_10 <- stepAIC(full.model, direction = "both", 
                        trace = FALSE)
step.model_10_2 <- lm(X5 ~ 0+Z5+Z9, data = Base_reg_2_5)
summary(step.model_10_2)
# Step-wise regression X3
Base_reg_3_1<-mutate(Base_reg_3[,c(1,6:15,17,19:20,22,24:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P6=as.factor(P6),P8=as.factor(P8),P9=as.factor(P9),P11=as.factor(P11))
full.model <- lm(X1 ~., data = Base_reg_3_1)
step.model_11 <- stepAIC(full.model, direction = "both", 
                        trace = FALSE)
step.model_11_2 <- lm(X1 ~ Z1+Z2+Z4+Z6+Z7+Z11, data = Base_reg_3_1)
summary(step.model_11_2)
Base_reg_3_2<-mutate(Base_reg_3[,c(2,6:15,17,19:20,22,24:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P6=as.factor(P6),P8=as.factor(P8),P9=as.factor(P9),P11=as.factor(P11))
full.model <- lm(X2 ~., data = Base_reg_3_2)
step.model_12 <- stepAIC(full.model, direction = "both", 
                         trace = FALSE)
step.model_12_2 <- lm(X2 ~ Z2+Z4+P8+Z7+Z8+Z9, data = Base_reg_3_2)
summary(step.model_12_2)
Base_reg_3_3<-mutate(Base_reg_3[,c(3,6:15,17,19:20,22,24:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P6=as.factor(P6),P8=as.factor(P8),P9=as.factor(P9),P11=as.factor(P11))
full.model <- lm(X3 ~., data = Base_reg_3_3)
step.model_13 <- stepAIC(full.model, direction = "both", 
                         trace = FALSE)
step.model_13_2 <- lm(X3 ~ Z1+Z2+Z6, data = Base_reg_3_3)
summary(step.model_13_2)
Base_reg_3_4<-mutate(Base_reg_3[,c(4,6:15,17,19:20,22,24:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P6=as.factor(P6),P8=as.factor(P8),P9=as.factor(P9),P11=as.factor(P11))
full.model <- lm(X4 ~., data = Base_reg_3_4)
step.model_14 <- stepAIC(full.model, direction = "both", 
                         trace = FALSE)
step.model_14_2 <- lm(X4 ~ Z7+Z8, data = Base_reg_3_4)
summary(step.model_14_2)
Base_reg_3_5<-mutate(Base_reg_3[,c(5,6:15,17,19:20,22,24:28)],P1=as.factor(P1),P2=as.factor(P2),P3=as.factor(P3),P4=as.factor(P4),P6=as.factor(P6),P8=as.factor(P8),P9=as.factor(P9),P11=as.factor(P11))
full.model <- lm(X5 ~., data = Base_reg_3_5)
step.model_15 <- stepAIC(full.model, direction = "both", 
                         trace = FALSE)
step.model_15_2 <- lm(X5 ~ Z1+Z2+Z6, data = Base_reg_3_5)
summary(step.model_15_2)
