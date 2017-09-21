    README

    Katie Shakman
    9/21/2017

    To use the Ch2 Version in this folder: 

 1. Decide if you want to analyze the full field of view or just
    an ROI within the images.  
    If you want an ROI, then 
    run makeROI_for_Ch2_Version.m
    It will first prompt you for a folder to look in for images, which it 
    will then average together and present in a figure window.  
    You can select an ROI on this figure (double-click to close the 
    polygon).  It will then save the ROI, which you will load for the next 
    step. 
 
 2. Run new_FluorAnalysis_batch_Sept2017_Ch2GCaMP.m
    You can run this on the entire image field of view, or you can load an 
    ROI made in Step 1 (and give a name for the ROI), and the code will 
    analyze only the part of the image within that ROI.  
    The code should be run from within an imaging folder which contains 
    subfolders for each trial, all with the same conditions.  For ex, you
    may choose to run it in a folder called "bitter", which contains five
    subfolders representing different trials when bitter was presented.  
    Each subfolder should begin with "TSeries", and should contain a series
    of images.  
    
 3. Run PlotAll_DelFoverF_ofType_zScore.m
    This should be run from within the Data_Analysis/... folder generated 
    in Step 2.  It will analyze the data from that folder and save figures
    in that same Data_Analysis/... location.  
 
 