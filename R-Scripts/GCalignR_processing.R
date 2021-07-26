#GCalignR processing
library("GCalignR")



#get data source starting with all known compounds

peak_data<- paste0("Output/OutputData/GCAlign_Ready_RT.Peakarea_KnownCompounds.txt")

#if plot = T, a histogram of peaks is plotted
check_input(data = peak_data,plot = T)  


#identify parameter min_diff_peak2peak based on average distance betweenneightbopuring peaks 

peak_interspace(data = peak_data, rt_col_name = "RT",
                quantile_range = c(0, 0.8), quantiles = 0.05)



#align Chromatograms 


peak_data_aligned <- align_chromatograms(data = peak_data, # input data
                                         rt_col_name = "RT", # retention time variable name 
                                         reference = "arizon_c_3_49", # choose samples with the most peaks  
                                         max_linear_shift = 0, # max. shift for linear corrections
                                         max_diff_peak2mean = 0.03, # max. distance of a peak to the mean across samples
                                         min_diff_peak2peak = 0.06, # min. expected distance between peaks (2*peak2mean)
                                         delete_single_peak = F, # delete peaks that are present in just one sample 
                                         permute = T,
                                         write_output = NULL) # add variable names to write aligned data to text files


gc_heatmap(peak_data_aligned,threshold = 0.07) 


plot(peak_data_aligned,which_plot = "all")
print(peak_data_aligned)




write.csv(peak_data_aligned$aligned$RT,"Output/OutputData/Aligned_RT_known.csv")
write.csv(peak_data_aligned$aligned$Area,"Output/OutputData/Aligned_RT_Area_known.csv")


#get data source starting with all compounds

peak_data<- paste0("Output/OutputData/GCalign_Ready_RT.Peakarea_AllCompounds.txt")

#if plot = T, a histogram of peaks is plotted
check_input(data = peak_data,plot = T)  


#identify parameter min_diff_peak2peak based on average distance betweenneightbopuring peaks 

peak_interspace(data = peak_data, rt_col_name = "RT",
                quantile_range = c(0, 0.8), quantiles = 0.05)



#align Chromatograms 


peak_data_aligned <- align_chromatograms(data = peak_data, # input data
                                         rt_col_name = "RT", # retention time variable name 
                                         reference = "arizon_c_3_49", # choose samples with the most peaks  
                                         max_linear_shift = 0, # max. shift for linear corrections
                                         max_diff_peak2mean = 0.03, # max. distance of a peak to the mean across samples
                                         min_diff_peak2peak = 0.06, # min. expected distance between peaks (2*peak2mean)
                                         delete_single_peak = F, # delete peaks that are present in just one sample 
                                         permute = T,
                                         write_output = NULL) # add variable names to write aligned data to text files

#use plots to assess groupings
gc_heatmap(peak_data_aligned,threshold = 0.07) 



plot(peak_data_aligned,which_plot = "all")
print(peak_data_aligned)



#write to csv
write.csv(peak_data_aligned$RT,"Output/OutputData/Aligned_RT_All_compounds.csv")
write.csv(peak_data_aligned$aligned$Area,"Output/OutputData/Aligned_RT_Area_All_compounds.csv")

