rm(list = ls())

data=read.csv(file.choose(),header=T,sep=";",dec=",")
head(data)
str(data)

# MATRIKS DATA Y
Y <- as.matrix(data[,11:19])
head(Y)
summary(Y)

# IDENTIFIKASI OUTLIER
# Menghitung rata-rata dan matriks kovarians
mu <- colMeans(Y)  
cov_matrix <- cov(Y)  

# Menghitung jarak mahalanobis untuk setiap observasi
library(MASS)
mahalanobis_dist <- mahalanobis(Y, center = mu, cov = cov_matrix)

# Box-Plot untuk Jarak Mahalanobis
boxplot(mahalanobis_dist, main="Intervensi Sensitif", ylab="Mahalanobis Distance")

n=22
# STANDARISASI DATA
y1 <- (Y[,1]-mean(Y[,1]))/(sd(Y[,1]))*sqrt((n-1)/n)
y2 <- (Y[,2]-mean(Y[,2]))/(sd(Y[,2]))*sqrt((n-1)/n)
y3 <- (Y[,3]-mean(Y[,3]))/(sd(Y[,3]))*sqrt((n-1)/n)
y4 <- (Y[,4]-mean(Y[,4]))/(sd(Y[,4]))*sqrt((n-1)/n)
y5 <- (Y[,5]-mean(Y[,5]))/(sd(Y[,5]))*sqrt((n-1)/n)
y6 <- (Y[,6]-mean(Y[,6]))/(sd(Y[,6]))*sqrt((n-1)/n)
y7 <- (Y[,7]-mean(Y[,7]))/(sd(Y[,7]))*sqrt((n-1)/n)
y8 <- (Y[,8]-mean(Y[,8]))/(sd(Y[,8]))*sqrt((n-1)/n)
y9 <- (Y[,9]-mean(Y[,9]))/(sd(Y[,9]))*sqrt((n-1)/n)
datastandar <- data.frame(x1=c(y1), x2=c(y2), x3=c(y3), x4=c(y4),
                          x5=c(y5),x6=c(y6),x7=c(y7),x8=c(y8),x9=c(y9))
X <- as.matrix(datastandar);X
round_X <- round(X,4);round_X

# NILAI EIGEN DAN VEKTOR EIGEN
Xt <- t(X)
XtX <- Xt%*%X
round(XtX,3)
eigen(XtX)

# MATRIKS L
Nilai_eigen <- eigen(XtX)$value
L <- diag(Nilai_eigen^(1/2)) ; L

# MATRIKS A
Vektor_eigen <- eigen(XtX)$vector
A <- Vektor_eigen ; A
round(A,3)
At <- t(A) ; At

# MATRIKS U
Linv <- solve(L)
U <- X%*%A%*%Linv ;U
round(U,3)

# IDENTIFIKASI PERSENTASE KERAGAMAN DATA
pca<- prcomp(X, scale = TRUE) # ini pake X (data yg distandar) atau Y (data asli)
sum_pca <- summary(pca)
sum_pca

# IDENTIFIKASI KORELASI VARIABEL ASLI DAN KOMPONEN UTAMA
var_diag <- diag(XtX) ; var_diag
sqrt_nilaieigen=sqrt(Nilai_eigen)
rho_matrix = t(t(Vektor_eigen) * sqrt_nilaieigen) / (var_diag) ; rho_matrix
round(rho_matrix,3)

# PROPORSI VARIANSI KOMPONEN UTAMA
q1 <- Nilai_eigen[1] / sum(Nilai_eigen);q1
q2 <- Nilai_eigen[2] / sum(Nilai_eigen);q2
q3 <- Nilai_eigen[3] / sum(Nilai_eigen);q3
q4 <- Nilai_eigen[4] / sum(Nilai_eigen);q4
q5 <- Nilai_eigen[5] / sum(Nilai_eigen);q5
q6 <- Nilai_eigen[6] / sum(Nilai_eigen);q6
q7 <- Nilai_eigen[7] / sum(Nilai_eigen);q7
q8 <- Nilai_eigen[8] / sum(Nilai_eigen);q8
q9 <- Nilai_eigen[9] / sum(Nilai_eigen);q9
q <- data.frame(PC = c("PC1","PC2",
                       "PC3","PC4","PC5","PC6", "PC7", "PC8","PC9"),q =
                  c(q1,q2,q3,q4,q5,q6,q7,q8,q9))

# VARIANS KUMULATIF
PC1 <- q1
PC2 <- PC1+q2
PC3 <- PC2+q3
PC4 <- PC3+q4
PC5 <- PC4+q5
PC6 <- PC5+q6
PC7 <- PC6+q7
PC8 <- PC7+q8
PC9 <- PC7+q9
PC <- data.frame(Cummulative =
                   c(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9))                 
cbind(q,PC)

# MATRIKS BARIS (G) DAN MATRIKS KOLOM (H)
alpha <- 0.5
G <- U%*%L^alpha; G
H <- A%*%L^alpha; H
G9 <- G[,1:9] ; G9
H9 <- H[,1:9] ; H9
round(G9,3)
round(H9,3)

#======== IDENTIFIKASI HASIL BIPLOT =========
# KORELASI ANTAR OBJEK DAN VARIABEL
# install.packages("geometry")
library(geometry)
korelasi_obvar <- function(X, Y) {
  p <- nrow(X)
  q <- nrow(Y)
  r <- matrix(, nrow = p, ncol = q)
  for (i in 1:p) {
    for (j in 1:q) {
      r[i, j] <- dot(X[i, ], Y[j, ]) / (sqrt(sum(X[i, ]^2)) *
                                          sqrt(sum(Y[j, ]^2)))
    }
  }
  print(r)
}
obvar<- korelasi_obvar(G9, H9)
obvar_3 <- round(obvar,3);obvar_3

# PEMBOBOTAN VARIABEL BERDASARKAN PCA
# Bobot formatif
sqrtvektor = sqrt(Vektor_eigen^2)
round(sqrtvektor,4)
rata_vektor <- colMeans(sqrtvektor)
bobot_formatif <- sweep(sqrtvektor, 2, rata_vektor, FUN = "/") * q1
round(bobot_formatif,4)

# Bobot terstandar
jumlah_bobot <- colSums(bobot_formatif)
round(jumlah_bobot,4)
bobot_terstandar <- sweep(bobot_formatif, 2, jumlah_bobot, FUN = "/")
round(bobot_terstandar,4)

# Bobot akhir
bobot_akhir <- sweep(bobot_terstandar, 2, q$q, FUN = "*")
print(bobot_akhir)
round(bobot_akhir,4)
round(colSums(bobot_akhir),4) #proporsi varians
bobot_variabel = round(rowSums(bobot_akhir),4) ; bobot_variabel
