object testLoad {
  val h = "barcode,doublet_probability,doublet_score,orig.ident,nCount_RNA,nFeature_RNA,percent.mt,brainregion,modality,oriBarcode,sublib,sex,rep,exp,mouse_reads_percentage,anno,isdlt,L1UMAP1,L1UMAP2,L1,L2UMAP1,L2UMAP2,L2,L1_2,L3UMAP1,L3UMAP2,L3,L1_2_3,L4UMAP1,L4UMAP2,L4,L1_2_3_4,L5UMAP1,L5UMAP2,L5,L1_2_3_4_5,cl,supertype_id_label,subclass_id_label,L5r,isNeuFromL1,class_id_label"
  val s =
    "A01:A1:B2:04,0.0150687809611147,0.0002804410995583,mousebrain,485,327,1.44329896907216,HYP,H3K27ac,4,A01,Female,FemaleB,Exp1,1,mouse,FALSE,6.2512092590332,-8.76627922058106,13,0.31991320848465,-1.3445326089859,20,1320,7.85275506973267,8.74591541290283,0,132000,,,0,13200000,,,0,1320000000,cl-14990,1175 Ependymal NN_1,323 Ependymal NN,L5-1320000000,FALSE,30 Astro-Epen"

  val hvec = h.split(",")
  val svec = s.split(",")
  Range(0, hvec.length).map(i => s"$i:${hvec(i)}:${svec(i)}")

}
