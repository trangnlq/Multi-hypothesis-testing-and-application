# nhap thu vien va du lieu
library(multtest, verbose = FALSE)
library(dplyr)
data(golub)
# quan sat du lieu
golub[1:5, 1:5]
golub.gnames[1:15, ]
golub.cl
# tinh thong ke t
teststat = mt.teststat(golub, golub.cl)
teststat = round(teststat, 2)
teststat[1:100]
# tinh p-value ban dau
rawp = 2*(1-pnorm(abs(teststat)))
rawp = round(rawp, 2)
rawp[1:100]
df_rawp = data.frame(rawp)
# ham hieu chinh Bonferroni
Bonferroni = function(vector){
  m = length(vector)
  for (i in 1:m){
    vector[i] = min(vector[i]*m, 1)
  }
  return(vector)
}
adjp_B = Bonferroni(rawp)
adjp_B = round(adjp_B, 2)
df_adjp_B = data.frame(adjp_B)
# ham hieu chinh Holm
Holm = function(vector){
  m = length(vector)
  sorted = sort(vector, decreasing = FALSE)
  match = match(vector, sorted)
  for (i in 1:m){
    sorted[i] = min(sorted[i]*(m+1-i), 1)
  }
  sorted = sorted[match]
  return(sorted)
}
adjp_H = Holm(rawp)
adjp_H = round(adjp_H, 2)
df_adjp_H = data.frame(adjp_H)
# ham hieu chinh Benjamini_Hochberg
Benjamini_Hochberg = function(vector){
  m = length(vector)
  sorted = sort(vector, decreasing = FALSE)
  match = match(vector, sorted)
  for (i in 1:m){
    sorted[i] = sorted[i]*(m/i)
  }
  sorted = sorted[match]
  return(sorted)
}
adjp_BH = Benjamini_Hochberg(rawp)
adjp_BH = round(adjp_BH, 2)
df_adjp_BH = data.frame(adjp_BH)

# in ra p-value hieu chinh 
summarize = bind_cols(df_rawp,df_adjp_B,df_adjp_H,df_adjp_BH)

summarize[1:100, ]

# ham bac bo 
reject = function(vector, alpha){
  m = length(vector)
  for (i in m:1){
    if (vector[i]>alpha) { 
      vector[i] = 0
    }
    else {
      vector[i] = 1
    }
  }
  return(vector)
}
reject_raw = reject(rawp,0.05)
reject_B = reject(adjp_B,0.05)
reject_H = reject(adjp_H,0.05)
reject_BH = reject(adjp_BH,0.05)
df_reject_raw = data.frame(reject_raw)
df_reject_B = data.frame(reject_B)
df_reject_H = data.frame(reject_H)
df_reject_BH = data.frame(reject_BH)
rejected = bind_cols(df_reject_raw, df_reject_B, df_reject_H, df_reject_BH)
print("Ma tran bac bo gia thuyet:")
rejected[1:100, ]

# in ra gen co bieu hien khac
no_reject_B = 0
genes_different_B = c()
for (i in 1:m){
  if (reject_B[i] == 1) {
    no_reject_B = no_reject_B + 1
    genes_different_B[no_reject_B] = i
  }
  else next
}
genes_different_B = golub.gnames[genes_different_B, ]
print("Cac gen co bieu hien khac nhau thu duoc theo phuong phap Bonferroni:")
genes_different_B[1:100, ]
print("So luong gen co bieu hien khac nhau thu duoc theo phuong phap Bonferroni:")
no_reject_B

no_reject_H = 0
genes_different_H = c()
for (i in 1:m){
  if (reject_H[i] == 1) {
    no_reject_H = no_reject_H + 1
    genes_different_H[no_reject_H] = i
  }
  else next
}
genes_different_H = golub.gnames[genes_different_H, ]
print("Cac gen co bieu hien khac nhau thu duoc theo phuong phap Holm:")
genes_different_H[1:100, ]
print("So luong gen co bieu hien khac nhau thu duoc theo phuong phap Holm:")
no_reject_H

no_reject_BH = 0
genes_different_BH = c()
for (i in 1:m){
  if (reject_BH[i] == 1) {
    no_reject_BH = no_reject_BH + 1
    genes_different_BH[no_reject_BH] = i
  }
  else next
}
genes_different_BH = golub.gnames[genes_different_BH, ]
print("Cac gen co bieu hien khac nhau thu duoc theo phuong phap Benjamini - Hochberg:")
genes_different_BH[1:100, ]
print("So luong gen co bieu hien khac nhau thu duoc theo phuong phap Benjamini - Hochberg:")
no_reject_BH
