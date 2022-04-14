library(Seurat)
library(reshape2)
library(lme4)

#Inspired and adapated from Paulsen, Velasco, Kedaigle and Pegoni et al., 2022
#Demonstrating the script through the example of hCS30DIV dataset 
#Load the previously saved hCS30DIV Celltype file

hCS30D <- readRDS(file = "~/hCS30DIVCelltype.rds")

#Calculating cell type progression 
#OrgIdent contains Ctrl vs Mut information
#orig.ident contains the identity of the organoids C1, C2, C3, M1, M2, M3

data1 = hCS30D@data1.data[,c("celltype", "OrgIdent", "orig.ident")] 

# set Ctrl as base  
data1$OrgIdent =  factor(data1$OrgIdent, levels = c("Ctrl", "Mut")) 
data1$orig.ident = factor(paste0(data1$OrgIdent, "_", data1$orig.ident))

data1 = data1[data1$celltype != "Unknown", ] #remove unknown cells

# Change cell types name compatible with regression 
data1$celltype = factor(gsub(" ", ".", data1$celltype))
data1$celltype = factor(gsub("-", ".", data1$celltype))
data1$celltype = factor(gsub("/", ".", data1$celltype))

lmm_celltype <- function(data1, form = "(1|orig.ident)", celltype = "celltype", test = "OrgIdent") {
  dat = data1
  levs = levels(dat[,celltype])
  ret = c()
  for (i in levs){
    typ = as.character(i)
    dat[paste("celltype", typ, sep= "_")]=as.numeric(dat[,celltype]==i)# column of 1st for cell type of interest, 0s otherwise
    
    form_cur=as.formula(paste(paste("celltype",typ,sep="_"),paste(test,form,sep = " + "),sep = "~ ")) #model as fixed effect
    form_cur2 = as.formula(paste(paste("celltype",typ,sep = "_"),form,sep = "~ "))# model without treatment included
    print(form_cur)
    res=glmer(form_cur,data = dat, family = "binomial")
    print(form_cur2)
    res2 = glmer(form_cur2, data = dat, family = "binomial")
    vari = paste(test,as.character(levels(dat[,test]))[2],sep = "")
    
    #in the first model, get coefficient related to Mutant (Ctrl vs Mut) effect
    coef = fixef(res)[vari]
    OR = exp(coef)
    print("start CI")
    CI = confint(res,parm=vari,method="Wald")
    CI_OR = exp(CI)
    
    #compare two models to get p value for each cell type 
    anov = anova(res, res2)
    pv = anov$"Pr(>Chisq)"[2]
    res=as.numeric(c(coef,OR,CI,CI_OR,pv))
    ret=rbind(ret,res)
  }
  ret=data.frame(ret)
  rownames(ret) = as.character(levs)
  colnames(ret) = c("coef","OR", "CI_coef_low", "CI_coef_high", "CI_OR_low", "CI_OR_high", "pval")
  return(ret)
}

ret = lmm_celltype(data1)
ret$adj.pval = p.adjust(ret$pval, method ="BH") # get adjusted p value amlunt of regressions we ran 

save(ret, file = "~/hCS30DCellComposition.LMM.NoUnk.Robj")
write.table(ret, file ="~/hCS30DCellComposition.LMM.NoUnk.txt", sep = "\t", quote=F,row.names =T)

load("~/hCS30DCellComposition.LMM.NoUnk.Robj")
ret1m1 = ret
ret1m1$dataset = "RNA"
ret1m1$celltype = rownames(ret1m1)
ret1m1$celltype = factor(ret1m1$celltype, levels = rev(c("Radial Glia", "Immature Interneurons", "Cycling Progenitors", "Deep Layer Neurons", "Early Cajal-Retzius Neurons", "Intermediate Progenitors", "Astrocytes", "Other")))
ret1m1$sig = ret1m1$adj.pval<0.05
