lesion_filling <- function(t1_path, flair_path, lesion_map_path, t1_brain_mask_path, output_path, npar=TRUE, print=TRUE){
t1 <- antsImageRead(t1_path, 3)
t1_brain_mask <- antsImageRead(t1_brain_mask_path, 3)
flair <- antsImageRead(flair_path, 3)
lesion_binary <- antsImageRead(lesion_map_path, 3)
t1_bias_corrected <- antsImageClone(t1)
N3BiasFieldCorrection(3, t1, t1_bias_corrected)
segs1 <- Atropos(d=3, x= t1_brain_mask, a = t1_bias_corrected, m = "[0.3,1x1x1]", p = "Socrates")
segs2 <- Atropos(d=3, x= t1_brain_mask, a = t1_bias_corrected, m = "[0.3,1x1x1]", p = "Socrates", i = segs1$probabilityimages)
wm_probability_map <- antsImageClone(segs2$probabilityimages[[3]])
wm_probability_map[wm_probability_map>0.7] <- 1
wm_binary_map <- wm_probability_map
antsRegOut <- antsRegistration(fixed=t1_bias_corrected, moving=flair, typeofTransform=c('Rigid'), outprefix=paste(output_path, "/flair2t1", sep = ""))
lesion_binary[lesion_binary>0] <- 1
lesionMapT1space <- antsApplyTransforms(fixed=t1_bias_corrected, moving=lesion_binary, transformlist = antsRegOut$fwdtransforms)
lesionMapT1space[lesionMapT1space>0] <- 1
nawm <- antsImageClone(wm_probability_map)
ImageMath(3, nawm, '-', wm_binary_map, lesionMapT1space)
nawm[nawm>0.3] <- 1
masked_t1 <- maskImage(t1_bias_corrected, nawm)
nawm_median_intenstiy <- median(as.array(masked_t1[masked_t1!=0]))
lesionMapT1spaceArray <- as.array(lesionMapT1space)
t1_bias_corrected_array <- as.array(t1_bias_corrected)
t1_bias_corrected_array[lesionMapT1spaceArray==1] <- nawm_median_intenstiy
t1_bias_corrected_image <- as.antsImage(t1_bias_corrected_array)
return(t1_bias_corrected_image)
}
